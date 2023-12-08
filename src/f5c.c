/* @file f5c.c
**
** f5c interface implementation
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "f5c.h"
#include "f5cmisc.h"



#include <sys/wait.h>
#include <unistd.h>

/*
todo :
Error counter for consecutive failures in the skip unreadable mode
not all the memory allocations are needed for eventalign mode
*/

//read bed file
static inline char **read_bed_regions(char *bedfile, int64_t *count){

    FILE *bedfp = fopen(bedfile,"r");
    F_CHK(bedfp,bedfile);

    char* buffer = (char*)malloc(sizeof(char) * (100)); //READ+newline+nullcharacter
    MALLOC_CHK(buffer);

    int64_t reg_capcacity = 1024;
    int64_t reg_i = 0;
    char **reg_list = (char **)malloc(reg_capcacity * sizeof(char *));
    MALLOC_CHK(reg_list);


    size_t bufferSize = 100;
    ssize_t readlinebytes = 0;
    int64_t line_no = 0;


    while ((readlinebytes = getline(&buffer, &bufferSize, bedfp)) != -1) {

        char *ref = (char *)malloc(sizeof(char)*readlinebytes);
        MALLOC_CHK(ref);
        int64_t beg=-1;
        int64_t end=-1;

        //TODO can optimised though strtok etc later
        int ret=sscanf(buffer,"%s\t%ld\t%ld",ref,&beg, &end);
        if(ret!=3 || end<beg){
            ERROR("Malformed bed entry at line %ld",line_no);
            exit(EXIT_FAILURE);
        }

        if(reg_i>=reg_capcacity){
            if(reg_capcacity>1000000){
                WARNING("The region bed file has over %ld regions. To reduce memory usage, you may consider merging bed regions.",reg_i);
            }
            reg_capcacity=reg_capcacity*2;
            reg_list = (char **)realloc((void *)reg_list,reg_capcacity * sizeof(char *));
            MALLOC_CHK(reg_list);

        }

        reg_list[reg_i] = (char *)malloc(sizeof(char)*readlinebytes);
        sprintf(reg_list[reg_i],"%s:%ld-%ld",ref, beg, end);
        reg_i++;


        free(ref);

        line_no++;
    }

    fclose(bedfp);
    free(buffer);
    *count = reg_i;

    return reg_list;
}

/* initialise the core data structure */
core_t* init_core(const char* bamfilename, const char* fastafile,
                  const char* fastqfile, const char* tmpfile, opt_t opt,double realtime0, int8_t mode, char *eventalignsummary, char *slow5file) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    if(opt.num_iop > 1){
        init_iop(core,opt);
    }

    // load bam file
    core->m_bam_fh = sam_open(bamfilename, "r");
    NULL_CHK(core->m_bam_fh);

    // load bam index file
    core->m_bam_idx = sam_index_load(core->m_bam_fh, bamfilename);
    if(core->m_bam_idx==NULL){
        ERROR("could not load the .bai index file for %s", bamfilename);
        fprintf(stderr, "Please run 'samtools index %s'\n", bamfilename);
        exit(EXIT_FAILURE);
    }

    // read the bam header
    core->m_hdr = sam_hdr_read(core->m_bam_fh);
    NULL_CHK(core->m_hdr);

    // If processing a region of the genome, get clipping coordinates
    core->clip_start = -1;
    core->clip_end = -1;

    core->reg_list=NULL; //region list is NULL by default
    core->reg_i=0;
    core->reg_n=0;
    core->itr = NULL;

    if(opt.region_str == NULL){
        core->itr = sam_itr_queryi(core->m_bam_idx, HTS_IDX_START, 0, 0);
        if(core->itr==NULL){
            ERROR("%s","sam_itr_queryi failed. A problem with the BAM index?");
            exit(EXIT_FAILURE);
        }
    }
    else{
        //determine if .bed
        int region_str_len = strlen(opt.region_str);
        if(region_str_len>=4 && strcmp(&(opt.region_str[region_str_len-4]),".bed")==0 ){
            STDERR("Fetching the list of regions from file: %s", opt.region_str);
            WARNING("%s", "Loading region windows from a bed file is an experimental option and not yet throughly tested.");
            WARNING("%s", "When loading windows from a bed file, output is based on reads that are unclipped. Also, there may be repeated entries when regions overlap.");
            int64_t count=0;
            char **reg_l = read_bed_regions(opt.region_str, &count);
            core->reg_list = reg_l;
            core->reg_i = 0;
            core->reg_n = count;
        }
        else{
            STDERR("Iterating over region: %s\n", opt.region_str);
            core->itr = sam_itr_querys(core->m_bam_idx, core->m_hdr, opt.region_str);
            if(core->itr==NULL){
                ERROR("sam_itr_querys failed. Please check if the region string you entered [%s] is valid",opt.region_str);
                exit(EXIT_FAILURE);
            }
            hts_parse_reg(opt.region_str, &(core->clip_start) , &(core->clip_end));
        }
    }


    //open the bam file for writing skipped ultra long reads
    core->ultra_long_tmp=NULL; //todo :  at the moment this is used to detect if the load balance mode is enabled. A better method in the opt flags.
    if(tmpfile!=NULL){
        core->ultra_long_tmp = sam_open(tmpfile, "wb");
        NULL_CHK(core->ultra_long_tmp);

        //write the header to the temporary file
        int ret_sw=sam_hdr_write(core->ultra_long_tmp,core->m_hdr);
        NEG_CHK(ret_sw);
    }

    if(opt.flag & F5C_WR_RAW_DUMP){
        core->raw_dump = fopen("f5c.tmp.bin","wb");
        F_CHK(core->raw_dump,"f5c.tmp.bin");
    }
    if(opt.flag & F5C_RD_RAW_DUMP){
        core->raw_dump = fopen("f5c.tmp.bin","rb");
        F_CHK(core->raw_dump,"f5c.tmp.bin");
    }
    if(opt.flag & F5C_RD_SLOW5){
        if(opt.flag & F5C_RD_RAW_DUMP){
            STDERR("%s","Option --slow5 is incompatible with --read-dump");
            exit(EXIT_FAILURE);
        }
       if(opt.flag & F5C_WR_RAW_DUMP){
            STDERR("%s","Option --slow5 is incompatible with --write-dump");
            exit(EXIT_FAILURE);
       }
        core->sf = slow5_open(slow5file,"r");
        if (core->sf == NULL) {
            STDERR("Error opening SLOW5 file %s\n",slow5file);
            exit(EXIT_FAILURE);
        }
        int ret=slow5_idx_load(core->sf);
        if(ret<0){
            STDERR("Error in loading SLOW5 index for %s\n",slow5file);
            exit(EXIT_FAILURE);
        }
    }

    // reference file
    core->fai = fai_load(fastafile);
    NULL_CHK(core->fai);

    // readbb
    core->readbb = new ReadDB;
    core->readbb->load(fastqfile, (opt.flag & F5C_RD_SLOW5)? 1 : 0);

    //model
    core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //4096 is 4^6 which is hardcoded now
    MALLOC_CHK(core->model);
    core->cpgmodel = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER_METH); //15625 is 4^6 which os hardcoded now
    MALLOC_CHK(core->cpgmodel);

    //load the model from files
    uint32_t kmer_size=0;
    uint32_t kmer_size_meth=0;
    if (opt.model_file) {
        kmer_size=read_model(core->model, opt.model_file, MODEL_TYPE_NUCLEOTIDE);
    } else {
        if(opt.flag & F5C_RNA){
            INFO("%s","builtin RNA nucleotide model loaded");
            kmer_size=set_model(core->model, MODEL_ID_RNA_NUCLEOTIDE);
        }
        else{
            kmer_size=set_model(core->model, MODEL_ID_DNA_NUCLEOTIDE);
        }
    }
    if (opt.meth_model_file) {
        kmer_size_meth=read_model(core->cpgmodel, opt.meth_model_file, MODEL_TYPE_METH);
    } else {
        kmer_size_meth=set_model(core->cpgmodel, MODEL_ID_DNA_CPG);
    }
    if( mode==0 && kmer_size != kmer_size_meth){
        ERROR("The k-mer size of the nucleotide model (%d) and the methylation model (%d) should be the same.",kmer_size,kmer_size_meth);
        exit(EXIT_FAILURE);
    }
    core->kmer_size = kmer_size;

    core->opt = opt;

    //realtime0
    core->realtime0=realtime0;

    core->load_db_time=0;
    core->process_db_time=0;

    core->db_bam_time=0;
    core->db_fasta_time=0;
    core->db_fast5_time=0;
    core->db_fast5_open_time=0;
    core->db_fast5_read_time=0;

    core->event_time=0;
    core->align_time=0;
    core->est_scale_time=0;
    core->meth_time=0;

    core->output_time=0;

    //cuda stuff
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        init_cuda(core);
    }
#endif

    core->sum_bases=0;
    core->total_reads=0; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)
    core->bad_fast5_file=0; //empty fast5 path returned by readdb, could not open fast5
    core->ultra_long_skipped=0;
    core->qc_fail_reads=0;
    core->failed_calibration_reads=0;
    core->failed_alignment_reads=0;

    //eventalign related
    core->mode = mode;
    core->read_index=0;
    if(mode==1){
        if(eventalignsummary!=NULL){
            core->event_summary_fp = fopen(eventalignsummary,"w");
            F_CHK(core->event_summary_fp,eventalignsummary);
        }
        else{
            core->event_summary_fp =NULL;
        }

        if(core->opt.flag & F5C_SAM){
            core->sam_output = hts_open("-", "w");
        }
    }


    return core;
}

/* free the core data structure */
void free_core(core_t* core,opt_t opt) {
    free(core->model);
    free(core->cpgmodel);
    delete core->readbb;
    fai_destroy(core->fai);

    if(core->itr){
        sam_itr_destroy(core->itr);
    }
    if(core->reg_list){
        for(int64_t i=0;i<core->reg_n;i++){
            free(core->reg_list[i]);
        }
        free(core->reg_list);
    }

    bam_hdr_destroy(core->m_hdr);
    hts_idx_destroy(core->m_bam_idx);
    sam_close(core->m_bam_fh);

    if(core->ultra_long_tmp!=NULL){
        sam_close(core->ultra_long_tmp);
    }
    if(core->opt.flag&F5C_WR_RAW_DUMP || core->opt.flag&F5C_RD_RAW_DUMP){
        fclose(core->raw_dump);
    }
    if(core->opt.flag & F5C_RD_SLOW5){
        slow5_idx_unload(core->sf);
        slow5_close(core->sf);
    }
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        free_cuda(core);
    }
#endif
    //eventalign related
    if(core->mode==1 && core->event_summary_fp!=NULL){
        fclose(core->event_summary_fp);
    }
    if(core->mode==1 && core->opt.flag & F5C_SAM){
        hts_close(core->sam_output);
    }
    if(opt.num_iop > 1){
        free_iop(core,opt);
    }

    free(core);
}

/* initialise a data batch */
db_t* init_db(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_bam_rec = core->opt.batch_size;
    db->n_bam_rec = 0;

    db->bam_rec = (bam1_t**)(malloc(sizeof(bam1_t*) * db->capacity_bam_rec));
    MALLOC_CHK(db->bam_rec);

    int32_t i = 0;
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        db->bam_rec[i] = bam_init1();
        NULL_CHK(db->bam_rec[i]);
    }

    db->fasta_cache = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->fasta_cache);
    db->read = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->read);
    db->read_len = (int32_t*)(malloc(sizeof(int32_t) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_len);
    db->read_idx = (int64_t*)(malloc(sizeof(int64_t) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_idx);

    db->sig = (signal_t**)malloc(sizeof(signal_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->sig);

    db->et = (event_table*)malloc(sizeof(event_table) * db->capacity_bam_rec);
    MALLOC_CHK(db->et);

    db->scalings =
        (scalings_t*)malloc(sizeof(scalings_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->scalings);

    db->event_align_pairs =
        (AlignedPair**)malloc(sizeof(AlignedPair*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_align_pairs);
    db->n_event_align_pairs =
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_align_pairs);

    db->event_alignment = (event_alignment_t**)malloc(
        sizeof(event_alignment_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_alignment);
    db->n_event_alignment =
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_alignment);

    db->events_per_base =
        (double*)malloc(sizeof(double) * db->capacity_bam_rec);
    MALLOC_CHK(db->events_per_base);

    db->base_to_event_map =
        (index_pair_t**)malloc(sizeof(index_pair_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->base_to_event_map);

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->read_stat_flag);

    db->site_score_map = (std::map<int, ScoredSite> **)malloc(sizeof(std::map<int, ScoredSite> *) * db->capacity_bam_rec);
    MALLOC_CHK(db->site_score_map);

    for (i = 0; i < db->capacity_bam_rec; ++i) {
        db->site_score_map[i] = new std::map<int, ScoredSite>;
        NULL_CHK(db->site_score_map[i]);
    }

    db->total_reads=0;
    db->bad_fast5_file=0;
    db->ultra_long_skipped=0;

    //eventalign related
    if(core->mode==1){
        db->eventalign_summary = (EventalignSummary *)malloc(sizeof(EventalignSummary) * db->capacity_bam_rec);
        MALLOC_CHK(db->eventalign_summary);

        db->event_alignment_result = (std::vector<event_alignment_t> **)malloc(sizeof(std::vector<event_alignment_t> *) * db->capacity_bam_rec);
        MALLOC_CHK(db->event_alignment_result);

        db->event_alignment_result_str = (char **)malloc(sizeof(char *) * db->capacity_bam_rec);
        MALLOC_CHK(db->event_alignment_result_str);

        for (i = 0; i < db->capacity_bam_rec; ++i) {
            db->event_alignment_result[i] = new std::vector<event_alignment_t> ;
            NULL_CHK(db->event_alignment_result[i]);
            (db->eventalign_summary[i]).num_events=0; //done here in the same loop for efficiency
            db->event_alignment_result_str[i] = NULL;
        }

    }
    else{
        db->eventalign_summary = NULL;
        db->event_alignment_result = NULL;
        db->event_alignment_result_str = NULL;
    }

    return db;
}

/* load a data batch from disk */
ret_status_t load_db(core_t* core, db_t* db) {
    if(core->opt.num_iop == 1 || (core->opt.flag & F5C_RD_SLOW5) ){
        return load_db1(core,db);
    }
    else{
        if (core->opt.flag & F5C_PRINT_RAW) {
            ERROR("%s","Printing data unsupported with --iop");
            exit(EXIT_FAILURE);
        }
        if (core->opt.flag & F5C_RD_RAW_DUMP){
            ERROR("%s","Reading from raw dump is unsupported with --iop");
            assert(0);
        }
        if(core->opt.flag & F5C_WR_RAW_DUMP){
            ERROR("%s","Writing to raw dump is unsupported with --iop");
            exit(EXIT_FAILURE);
        }
        return load_db2(core,db);
    }
}


#ifdef WORK_STEAL
static inline int32_t steal_work(pthread_arg_t* all_args, int32_t n_threads)
{

	int32_t i, c_i = -1;
	int32_t k;
	for (i = 0; i < n_threads; ++i){
        pthread_arg_t args = all_args[i];
        //fprintf(stderr,"endi : %d, starti : %d\n",args.endi,args.starti);
		if (args.endi-args.starti > STEAL_THRESH) {
            //fprintf(stderr,"gap : %d\n",args.endi-args.starti);
            c_i = i;
            break;
        }
    }
    if(c_i<0){
        return -1;
    }
	k = __sync_fetch_and_add(&(all_args[c_i].starti), 1);
    //fprintf(stderr,"k : %d, end %d, start %d\n",k,all_args[c_i].endi,all_args[c_i].starti);
	return k >= all_args[c_i].endi ? -1 : k;
}
#endif

void* pthread_single(void* voidargs) {
    int32_t i;
    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;

#ifndef WORK_STEAL
    for (i = args->starti; i < args->endi; i++) {
        args->func(core,db,i);
    }
#else
    pthread_arg_t* all_args = (pthread_arg_t*)(args->all_pthread_args);
    //adapted from kthread.c in minimap2
	for (;;) {
		i = __sync_fetch_and_add(&args->starti, 1);
		if (i >= args->endi) {
            break;
        }
		args->func(core,db,i);
	}
	while ((i = steal_work(all_args,core->opt.num_thread)) >= 0){
		args->func(core,db,i);
    }
#endif

    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}


void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int)){

    if (core->opt.num_thread == 1) {
        int i;
        for (i = 0; i < db->n_bam_rec; i++) {
            func(core,db,i);
        }
    }
    else{
        //create threads
        pthread_t tids[core->opt.num_thread];
        pthread_arg_t pt_args[core->opt.num_thread];
        int32_t t, ret;
        int32_t i = 0;
        int32_t num_thread = core->opt.num_thread;
        int32_t step = (db->n_bam_rec + num_thread - 1) / num_thread;
        //todo : check for higher num of threads than the data
        //current works but many threads are created despite

        //set the data structures
        for (t = 0; t < num_thread; t++) {
            pt_args[t].core = core;
            pt_args[t].db = db;
            pt_args[t].starti = i;
            i += step;
            if (i > db->n_bam_rec) {
                pt_args[t].endi = db->n_bam_rec;
            } else {
                pt_args[t].endi = i;
            }
            pt_args[t].func=func;
        #ifdef WORK_STEAL
            pt_args[t].all_pthread_args =  (void *)pt_args;
        #endif
            //fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);

        }

        //create threads
        for(t = 0; t < core->opt.num_thread; t++){
            ret = pthread_create(&tids[t], NULL, pthread_single,
                                    (void*)(&pt_args[t]));
            NEG_CHK(ret);
        }

        //pthread joining
        for (t = 0; t < core->opt.num_thread; t++) {
            int ret = pthread_join(tids[t], NULL);
            NEG_CHK(ret);
        }
        
    }
}


void event_single(core_t* core,db_t* db, int32_t i) {

    if(db->sig[i]->nsample>0){

        float* rawptr = db->sig[i]->rawptr;
        float range = db->sig[i]->range;
        float digitisation = db->sig[i]->digitisation;
        float offset = db->sig[i]->offset;
        int32_t nsample = db->sig[i]->nsample;

        // convert to pA
        float raw_unit = range / digitisation;
        for (int32_t j = 0; j < nsample; j++) {
            rawptr[j] = (rawptr[j] + offset) * raw_unit;
        }

        int8_t rna=0;
        if (core->opt.flag & F5C_RNA){
            rna=1;
        }
        db->et[i] = getevents(db->sig[i]->nsample, rawptr, rna);

        // if(db->et[i].n/(float)db->read_len[i] > 20){
        //     fprintf(stderr,"%s\tevents_per_base\t%f\tread_len\t%d\n",bam_get_qname(db->bam_rec[i]), db->et[i].n/(float)db->read_len[i],db->read_len[i]);
        // }

        //get the scalings
        db->scalings[i] = estimate_scalings_using_mom(
            db->read[i], db->read_len[i], core->model, core->kmer_size, db->et[i]);

        //If sequencing RNA, reverse the events to be 3'->5'
        if (rna){
            event_t *events = db->et[i].event;
            size_t n_events = db->et[i].n;
            for (size_t i = 0; i < n_events/2; ++i) {
                event_t tmp_event = events[i];
                events[i]=events[n_events-1-i];
                events[n_events-1-i]=tmp_event;
            }
        }

        //allocate memory for the next alignment step
        db->event_align_pairs[i] = (AlignedPair*)malloc(
                    sizeof(AlignedPair) * db->et[i].n * 2); //todo : find a good heuristic to save memory
        MALLOC_CHK(db->event_align_pairs[i]);
     
    }
    
    else{
        db->et[i].n = 0;
        db->et[i].event = NULL;
        db->event_align_pairs[i] = NULL;
    }

}

void scaling_single(core_t* core, db_t* db, int32_t i){

    db->event_alignment[i] = NULL;
    db->n_event_alignment[i] = 0;
    db->events_per_base[i] = 0; //todo : is double needed? not just int8?

    int32_t n_kmers = db->read_len[i] - core->kmer_size + 1;

    if (db->n_event_align_pairs[i] > 0) {

        db->base_to_event_map[i]=(index_pair_t*)(malloc(sizeof(index_pair_t) * n_kmers));
        MALLOC_CHK(db->base_to_event_map[i]);

        // prepare data structures for the final calibration

        db->event_alignment[i] = (event_alignment_t*)malloc(
            sizeof(event_alignment_t) * db->n_event_align_pairs[i]);
        MALLOC_CHK(db->event_alignment[i]);

        // for (int j = 0; j < n_event_align_pairs; ++j) {
        //     fprintf(stderr, "%d-%d\n",event_align_pairs[j].ref_pos,event_align_pairs[j].read_pos);
        // }


        //todo : verify if this n is needed is needed
        db->n_event_alignment[i] = postalign(
            db->event_alignment[i],db->base_to_event_map[i], &db->events_per_base[i], db->read[i],
            n_kmers, db->event_align_pairs[i], db->n_event_align_pairs[i], core->kmer_size, db->read_len[i], db->et[i].n);

        //fprintf(stderr,"n_event_alignment %d\n",n_events);

        // run recalibration to get the best set of scaling parameters and the residual
        // between the (scaled) event levels and the model.

        // internally this function will set shift/scale/etc of the pore model
        bool calibrated = recalibrate_model(
            core->model, core->kmer_size, db->et[i], &db->scalings[i],
            db->event_alignment[i], db->n_event_alignment[i], 1, core->opt.min_num_events_to_rescale);

        // QC calibration
        if (!calibrated || db->scalings[i].var > MIN_CALIBRATION_VAR) {
            //     events[strand_idx].clear();
            free(db->event_alignment[i]);
            //free(db->event_align_pairs[i]);
            db->read_stat_flag[i] |= FAILED_CALIBRATION;
            return;
        }

        free(db->event_alignment[i]);

    } else {
        db->base_to_event_map[i] = NULL;
        // Could not align, fail this read
        // this->events[strand_idx].clear();
        // this->events_per_base[strand_idx] = 0.0f;
        //free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_ALIGNMENT;
        return;
    }

    // Filter poor quality reads that have too many "stays"

    if (db->events_per_base[i] > 5.0) {
        //     events[0].clear();
        //     events[1].clear();
        //free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_QUALITY_CHK;
        return;
    }


}

/* align a single read specified by index i (perform ABEA for a single read) */
//note that this is used in f5c.cu and thus modifications must be done with care
void align_single(core_t* core, db_t* db, int32_t i) {

    if(db->sig[i]->nsample>0){ //if a good read
        if (db->sig[i]->nsample && (db->et[i].n)/(float)(db->read_len[i]) < AVG_EVENTS_PER_KMER_MAX){
            db->n_event_align_pairs[i] = align(
                    db->event_align_pairs[i], db->read[i], db->read_len[i], db->et[i],
                    core->model, core->kmer_size, db->scalings[i], db->sig[i]->sample_rate);
                //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
        }
        else{//todo : too many avg events per base - oversegmented
            db->n_event_align_pairs[i]=0;
            if(core->opt.verbosity > 0){
                STDERR("Skipping over-segmented read %s with %f events per base",bam_get_qname(db->bam_rec[i]), (db->et[i].n)/(float)(db->read_len[i]));
            }
        }
    }
    else{ //if a bad read (corrupted/missing slow5 record)
       db->n_event_align_pairs[i]=0;
    }
  //  fprintf(stderr,"align single  read_id %d n_pairs %d\n",i, db->n_event_align_pairs[i]);
}

/* align a data batch (perform ABEA for a data batch) */
void align_db(core_t* core, db_t* db) {
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        //STDERR("%s","Performing on cuda");
        align_cuda(core, db);
    }
#endif

    if (core->opt.flag & F5C_DISABLE_CUDA) {
        //fprintf(stderr, "cpu\n");
        pthread_db(core, db, align_single);
    }
}


void align_db_FPGA_by_len(core_t * core, db_t* db)
{
  
    const int n_bam_rec = db->n_bam_rec;
    int kmer_size = core->kmer_size;

    int LENGTH_RANGE = 20000;

    int bucket = 0;
    int reads_num[40];
    int max_event_length[40];
    int max_read_length[40];
    int min_event_length[40];
    int min_read_length[40];
    int max_n_bands[40];

    for(int i=0; i< 40; i++)
    {
        reads_num[i] = 0;
        max_event_length[i] = 0;
        max_read_length[i] = 0;
        min_event_length[i] = INT_MAX;
        min_read_length[i] =  INT_MAX;
        max_n_bands[i] = 0;
    }
    int reads_id[40][BUCKET_SIZE];
    
 
    int ultra_long_reads_num = 0;
    int ultra_long_reads_id[n_bam_rec];


    int cpu_reads_id[n_bam_rec];
    int cpu_reads_num = 0;

    //ultra-long read on CPU
    for(int i = 0; i < n_bam_rec; i++)
    {

        int event_length = db->et[i].n;
        int read_length = db->read_len[i];

       if(db->sig[i]->nsample>0 && (event_length>=120000 || read_length >=60000 ) && event_length/(float)read_length < AVG_EVENTS_PER_KMER_MAX)
        {
            ultra_long_reads_id[ultra_long_reads_num] = i;
            ultra_long_reads_num++;
            core->ultra_long_reads_num++;
           //core->statis_event_length[event_length/10000]++;
        }

    }

    //pipeline kernel to process ultra-long reads
    //ultra_long_pipeline(core, db, ultra_long_reads_id, ultra_long_reads_num);


    //pthread to call ultra-long kernel
    pthread_arg_t *tmparg_ultra = NULL;
    pthread_t tid_ultra =  align_ultra_long_async(&tmparg_ultra, core, db, ultra_long_reads_id, ultra_long_reads_num);



    //fill BUCKET data
    for(int i = 0; i < n_bam_rec; i++)
    {
        int event_length = db->et[i].n;
        int read_length = db->read_len[i];
        int length = event_length + read_length;
        int n_bands = event_length + read_length - kmer_size + 1 + 2;
        int bucket_id = 999;
    
        bucket_id = length/LENGTH_RANGE;
 
        if(db->sig[i]->nsample>0 && (event_length < 120000 && read_length < 60000) && event_length/(float)read_length < AVG_EVENTS_PER_KMER_MAX)
        {
         
            reads_id[bucket_id][reads_num[bucket_id]] = i;
            if(read_length > max_read_length[bucket_id])
                max_read_length[bucket_id] = read_length;

            if(event_length > max_event_length[bucket_id])
                max_event_length[bucket_id] = event_length;

            if(read_length < min_read_length[bucket_id])
                min_read_length[bucket_id] = read_length;

            if(event_length < min_event_length[bucket_id])
                min_event_length[bucket_id] = event_length;

            if(n_bands > max_n_bands[bucket_id])
                max_n_bands[bucket_id] = n_bands; 

            reads_num[bucket_id]++;
          //  core->statis_event_length[event_length/10000]++;

            if(reads_num[bucket_id]==BUCKET_SIZE)
            {
                fprintf(stderr, "(%d, %d): bucket %d is full, reads_num %d max_bands %d max_read %d min_read %d max_event %d min_event %d \n", db->batch_idx, bucket, bucket_id,reads_num[bucket_id], max_n_bands[bucket_id], max_read_length[bucket_id], min_read_length[bucket_id], max_event_length[bucket_id], min_event_length[bucket_id] );
                core->total_num_bands+= max_n_bands[bucket_id];
                fill_bucket(core, db,bucket, kmer_size, max_read_length[bucket_id] , max_event_length[bucket_id], reads_num[bucket_id], reads_id[bucket_id]); 
                reads_num[bucket_id] = 0;
                max_read_length[bucket_id] = 0;
                max_event_length[bucket_id] = 0;
                min_read_length[bucket_id] = INT_MAX;
                min_event_length[bucket_id] = INT_MAX;
                max_n_bands[bucket_id] = 0;
                bucket++;
                //pthread
            }
        }
        
    }

   // fprintf(stderr,"\n NOT FULL:\n");
    for(int bucket_id=0; bucket_id<40; bucket_id++)
    {
        if(reads_num[bucket_id]>0)
        {
            fprintf(stderr, "(%d, %d): bucket %d, reads_num %d max_bands %d max_read %d min_read %d max_event %d min_event %d \n", db->batch_idx, bucket, bucket_id,reads_num[bucket_id], max_n_bands[bucket_id], max_read_length[bucket_id], min_read_length[bucket_id], max_event_length[bucket_id], min_event_length[bucket_id] );

            core->total_num_bands+=max_n_bands[bucket_id];
            fill_bucket(core, db, bucket, kmer_size, max_read_length[bucket_id] , max_event_length[bucket_id], reads_num[bucket_id], reads_id[bucket_id]); 
            bucket++;
            
        }

    }
 
     align_ultra_long_async_join(tmparg_ultra, tid_ultra);


}

void transfer_model_to_FPGA(core_t * core)
{
    fprintf(stderr,"open device\n");
    core->device = xrt::device(0);
    std::string binaryFile = "/home/centos/Methy/Methy_b100_unified_block/binary_container_Methy_b100_unified_block_v6.awsxclbin";

    fprintf(stderr,"load binary\n");
    auto uuid = core->device.load_xclbin(binaryFile);

    fprintf(stderr,"load kernel\n");
    core->kernel_0 = xrt::kernel(core->device, uuid, "align_top_stream:{align_top_stream_1}");
    core->kernel_1 = xrt::kernel(core->device, uuid, "align_top_stream:{align_top_stream_2}");
    core->kernel_2 = xrt::kernel(core->device, uuid, "align_top_stream_no_group:{align_top_stream_no_group_1}");

    size_t event_mean_size_bytes = sizeof(float) * BUCKET_SIZE/2 * EVENT_LEN;
    size_t read_size_bytes = sizeof(char) * BUCKET_SIZE/2 * READ_LEN;
    size_t n_align_size_bytes=sizeof(int)* BUCKET_SIZE/2;
    size_t aligned_pairs_size_bytes=sizeof(int) * BUCKET_SIZE/2 * 2 * N_BANDS;
    size_t trace_table_out_size_bytes =  32 * 8 * N_BANDS; //8 reads in pipeline
    size_t left_out_size_bytes = sizeof(int) * 8 * N_BANDS;
    size_t model_data_size_bytes = sizeof(float) * MAX_NUM_KMER * 4;
    size_t length_int_size_bytes=sizeof(int) * BUCKET_SIZE/2 * 2;
    size_t scailing_size_bytes = sizeof(float) * BUCKET_SIZE/2 * 4;
    size_t lp_size_bytes = sizeof(float) * BUCKET_SIZE/2 * 3;

    //size for ultra long
    size_t ultra_event_mean_size_bytes = sizeof(float) * ULTRA_EVENT_LEN;
    size_t ultra_read_size_bytes = sizeof(char) * ULTRA_READ_LEN;
    size_t ultra_n_align_size_bytes = sizeof(int);
    size_t ultra_aligned_pairs_size_bytes = sizeof(int) * ULTRA_N_BANDS*2;
    size_t ultra_trace_table_out_size_bytes = 32 * ULTRA_N_BANDS;
    size_t ultra_left_out_size_bytes = sizeof(int) * ULTRA_N_BANDS;
    size_t ultra_model_data_size_bytes = sizeof(float) * MAX_NUM_KMER * 4;
    size_t ultra_length_int_size_bytes = sizeof(int) * 2;
    size_t ultra_scailing_size_bytes = sizeof(float) * 4;
    size_t ultra_lp_size_bytes = sizeof(float) * 3;

    fprintf(stderr,"Allocate Buffer in Global Memory\n");
    double bo_start = realtime();

    core->event_mean_bo_0 = xrt::bo(core->device, event_mean_size_bytes, core->kernel_0.group_id(1));
    core->read_bo_0 = xrt::bo(core->device, read_size_bytes, core->kernel_0.group_id(2));
    core->n_align_bo_0 = xrt::bo(core->device, n_align_size_bytes, core->kernel_0.group_id(3));
    core->aligned_pairs_bo_0 = xrt::bo(core->device, aligned_pairs_size_bytes, core->kernel_0.group_id(4));
    core->model_data_bo_0 = xrt::bo(core->device, model_data_size_bytes, core->kernel_0.group_id(5));
    core->trace_table_bo_0 = xrt::bo(core->device,trace_table_out_size_bytes,core->kernel_0.group_id(6));
    core->left_out_bo_0 = xrt::bo(core->device, left_out_size_bytes, core->kernel_0.group_id(7));
    core->length_int_bo_0 = xrt::bo(core->device, length_int_size_bytes, core->kernel_0.group_id(8));
    core->scaling_info_bo_0 = xrt::bo(core->device, scailing_size_bytes, core->kernel_0.group_id(9));
    core->lp_info_bo_0 = xrt::bo(core->device, lp_size_bytes, core->kernel_0.group_id(10));

    core->event_mean_0 = core->event_mean_bo_0.map<float*>();
    core->reads_0 = core->read_bo_0.map<char*>();
    core->model_data_0 = core->model_data_bo_0.map<float*>();
    core->length_int_0 = core->length_int_bo_0.map<int*>();
    memset(core->length_int_0, 0, BUCKET_SIZE/2 * 2 * sizeof(int));
    core->scaling_info_0 =  core->scaling_info_bo_0.map<float*>();
    core->lp_info_0 = core->lp_info_bo_0.map<float*>();
    core->n_align_0 = core->n_align_bo_0.map<int*>();
    core->aligned_pairs_0 = core->aligned_pairs_bo_0.map<int*>();
    
        
    core->event_mean_bo_1 = xrt::bo(core->device, event_mean_size_bytes, core->kernel_1.group_id(1));
    core->read_bo_1 = xrt::bo(core->device, read_size_bytes, core->kernel_1.group_id(2));
    core->n_align_bo_1 = xrt::bo(core->device, n_align_size_bytes, core->kernel_1.group_id(3));
    core->aligned_pairs_bo_1 = xrt::bo(core->device, aligned_pairs_size_bytes, core->kernel_1.group_id(4));
    core->model_data_bo_1 = xrt::bo(core->device, model_data_size_bytes, core->kernel_1.group_id(5));
    core->trace_table_bo_1=xrt::bo(core->device,trace_table_out_size_bytes,core->kernel_1.group_id(6));
    core->left_out_bo_1 = xrt::bo(core->device, left_out_size_bytes, core->kernel_1.group_id(7));
    core->length_int_bo_1 = xrt::bo(core->device, length_int_size_bytes, core->kernel_1.group_id(8));
    core->scaling_info_bo_1 = xrt::bo(core->device, scailing_size_bytes, core->kernel_1.group_id(9));
    core->lp_info_bo_1 = xrt::bo(core->device, lp_size_bytes, core->kernel_1.group_id(10));

    core->event_mean_1 = core->event_mean_bo_1.map<float*>();
    core->reads_1 = core->read_bo_1.map<char*>();
    core->n_align_1 = core->n_align_bo_1.map<int*>();
    core->aligned_pairs_1 = core->aligned_pairs_bo_1.map<int*>();
    core->model_data_1 = core->model_data_bo_1.map<float*>();
    core->length_int_1 = core->length_int_bo_1.map<int*>();
    memset(core->length_int_1, 0, BUCKET_SIZE/2 * 2 * sizeof(int));
    core->scaling_info_1 =  core->scaling_info_bo_1.map<float*>();
    core->lp_info_1 = core->lp_info_bo_1.map<float*>();
    
    core->event_mean_bo_2 = xrt::bo(core->device, ultra_event_mean_size_bytes, core->kernel_2.group_id(1));
    core->read_bo_2 = xrt::bo(core->device, ultra_read_size_bytes, core->kernel_2.group_id(2));
    core->n_align_bo_2 = xrt::bo(core->device, ultra_n_align_size_bytes, core->kernel_2.group_id(3));
    core->aligned_pairs_bo_2 = xrt::bo(core->device, ultra_aligned_pairs_size_bytes, core->kernel_2.group_id(4));
    core->model_data_bo_2 = xrt::bo(core->device, ultra_model_data_size_bytes, core->kernel_2.group_id(5));
    core->trace_table_bo_2 = xrt::bo(core->device, ultra_trace_table_out_size_bytes, core->kernel_2.group_id(6));
    core->left_out_bo_2 = xrt::bo(core->device, ultra_left_out_size_bytes, core->kernel_2.group_id(7));
    core->length_int_bo_2 = xrt::bo(core->device, ultra_length_int_size_bytes, core->kernel_2.group_id(8));
    core->scaling_info_bo_2 = xrt::bo(core->device, ultra_scailing_size_bytes, core->kernel_2.group_id(9));
    core->lp_info_bo_2 = xrt::bo(core->device, ultra_lp_size_bytes, core->kernel_2.group_id(10));

    core->event_mean_2 = core->event_mean_bo_2.map<float*>();
    core->reads_2 = core->read_bo_2.map<char*>();
    core->n_align_2 = core->n_align_bo_2.map<int*>();
    core->aligned_pairs_2 = core->aligned_pairs_bo_2.map<int*>();
    core->model_data_2 = core->model_data_bo_2.map<float*>();
    core->length_int_2 = core->length_int_bo_2.map<int*>();
    memset(core->length_int_2, 0, 2 * sizeof(int));
    core->scaling_info_2 =  core->scaling_info_bo_2.map<float*>();
    core->lp_info_2 = core->lp_info_bo_2.map<float*>();

    double bo_end = realtime();

    core->new_bo_time += bo_end - bo_start;
    
   for(int i=0; i<MAX_NUM_KMER; i++)
   {
        core->model_data_0[i] = core->model[i].level_mean;
        core->model_data_0[i + MAX_NUM_KMER] = core->model[i].level_stdv;
        core->model_data_0[i + 2*MAX_NUM_KMER] = core->model[i].level_log_stdv;
        core->model_data_0[i + 3*MAX_NUM_KMER] = ((float)1)/core->model_data_0[i + MAX_NUM_KMER];//inverse_stdv

        core->model_data_1[i] = core->model[i].level_mean;
        core->model_data_1[i + MAX_NUM_KMER] = core->model[i].level_stdv;
        core->model_data_1[i + 2*MAX_NUM_KMER] = core->model[i].level_log_stdv;
        core->model_data_1[i + 3*MAX_NUM_KMER] = ((float)1)/core->model_data_1[i + MAX_NUM_KMER];//inverse_stdv

        core->model_data_2[i] = core->model[i].level_mean;
        core->model_data_2[i + MAX_NUM_KMER] = core->model[i].level_stdv;
        core->model_data_2[i + 2*MAX_NUM_KMER] = core->model[i].level_log_stdv;
        core->model_data_2[i + 3*MAX_NUM_KMER] = ((float)1)/core->model_data_2[i + MAX_NUM_KMER];//inverse_stdv

    }

    fprintf(stderr,"transfer model data\n");
   core->model_data_bo_0.sync(XCL_BO_SYNC_BO_TO_DEVICE);
   core->model_data_bo_1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
   core->model_data_bo_2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    auto run_0 = core->kernel_0(1 ,NULL,NULL,NULL,NULL, core->model_data_bo_0,NULL, NULL, NULL, NULL, NULL, NULL, NULL, 1);
    auto run_1 = core->kernel_1(1 ,NULL,NULL,NULL,NULL, core->model_data_bo_1,NULL, NULL, NULL, NULL, NULL, NULL, NULL,1);
    auto run_2 = core->kernel_2(1 ,NULL,NULL,NULL,NULL, core->model_data_bo_2,NULL, NULL, NULL, NULL, NULL, 100, 100, 1);
    run_0.wait();
    run_1.wait();
   run_2.wait();
   fprintf(stderr,"transfer model data successfully\n");

}

void ultra_long_pipeline(core_t* core, db_t* db, int* ultra_long_reads_id, int ultra_long_reads_num)
{
    xrt::kernel kernel = core->kernel_0;
    xrt::device device = core->device;
    xrt::bo event_mean_bo = core->event_mean_bo_0;
    xrt::bo read_bo = core->read_bo_0;
    xrt::bo n_align_bo = core->n_align_bo_0;
    xrt::bo aligned_pairs_bo = core->aligned_pairs_bo_0;
    xrt::bo trace_table_bo = core->trace_table_bo_0;
    xrt::bo left_out_bo = core->left_out_bo_0;
    xrt::bo length_int_bo = core->length_int_bo_0;
    xrt::bo scaling_info_bo = core->scaling_info_bo_0;
    xrt::bo lp_info_bo = core->lp_info_bo_0;

    int event_offset;
    int reads_offset;
    int n_bands;
    int kmer_size = core->kmer_size;
    float* event_mean = core->event_mean_0;
    char* reads = core->reads_0;
    int* length_int = core->length_int_0;
    memset(length_int, 0, BUCKET_SIZE/2 * 2 * sizeof(int));
    float* scaling_info = core->scaling_info_0;
    float* lp_info =core->lp_info_0;
    int* n_align = core->n_align_0;
    int* aligned_pairs = core->aligned_pairs_0;
    int batch_num = 1;

    for(int i=0; i<ultra_long_reads_num; i++)
    {
        double fill_start = realtime();
        int idx = ultra_long_reads_id[i];
        int read_length = db->read_len[idx];
        int event_length = db->et[idx].n;

        if(event_length%64!=0)
            event_offset = (event_length/64 + 1) * 64;
        else
            event_offset = event_length;

        if(read_length%64!=0)
            reads_offset = (read_length/64 + 1) * 64;
        else
            reads_offset = read_length;

        n_bands = event_offset + reads_offset;

        fprintf(stderr,"reads_offset %d event_offset %d n_bands %d\n", reads_offset, event_offset, n_bands);
     
        memcpy(&reads[0], db->read[idx], read_length);
        
        for(int j=0; j< event_length; j++)
            event_mean[j] = db->et[idx].event[j].mean; 

        length_int[0] =  read_length;
        length_int[BATCH_SIZE] = event_length;

        scaling_info[0] =  db->scalings[idx].log_var;
        scaling_info[BATCH_SIZE] = db->scalings[idx].scale;
        scaling_info[2 * BATCH_SIZE] = db->scalings[idx].shift;
        scaling_info[3 * BATCH_SIZE] = db->scalings[idx].var;

        int n_kmers = read_length - kmer_size + 1;
        float events_per_kmer = (float)event_length / (float)n_kmers;
        float p_stay = (float)1.0 - ((float)1.0 / (events_per_kmer + (float)1.0));
        float epsilon = 1e-10;
        lp_info[0] = log(epsilon);
        lp_info[BATCH_SIZE]= log(p_stay);
        lp_info[2 * BATCH_SIZE] = log((float)1.0 - exp(lp_info[0]) - exp(lp_info[BATCH_SIZE])); 

        double fill_end = realtime();
        core->fill_bucket_time_ultra += fill_end - fill_start;

        double sync_to_device_start = realtime();
        fprintf(stderr,"sync to device\n");
        size_t event_mean_size_bytes = sizeof(float) * event_offset;
        size_t read_size_bytes = sizeof(char) *  reads_offset;//
        size_t length_int_size_bytes=sizeof(int) * batch_num * BATCH_SIZE * 2;
        size_t scaling_size_bytes = sizeof(float) * batch_num * BATCH_SIZE * 4;
        size_t lp_size_bytes = sizeof(float) * batch_num * BATCH_SIZE * 3;//
        size_t n_align_size_bytes=sizeof(int);
        size_t aligned_pairs_size_bytes=sizeof(int) * batch_num * BATCH_SIZE * 2 * n_bands;

        event_mean_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE,event_mean_size_bytes, 0 );
        read_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, read_size_bytes, 0);
        length_int_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, length_int_size_bytes, 0);
        scaling_info_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, scaling_size_bytes, 0);
        lp_info_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, lp_size_bytes, 0);

        double sync_to_device_end = realtime();
        core->sync_to_device_ultra += sync_to_device_end - sync_to_device_start;

        double kernel_start = realtime();
        core->call_kernel_0_num++;
        fprintf(stderr,"run kernel\n");
        auto run = kernel(0, event_mean_bo,read_bo,n_align_bo,aligned_pairs_bo,NULL,trace_table_bo, left_out_bo, length_int_bo,scaling_info_bo,lp_info_bo,event_offset,reads_offset,batch_num);

        run.wait();
        double kernel_end = realtime();
        core->kernel_run_ultra += kernel_end - kernel_start;

        double sync_from_device_start = realtime(); 
        fprintf(stderr,"sync to host\n");

        aligned_pairs_bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE, aligned_pairs_size_bytes,0);
        n_align_bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE, n_align_size_bytes, 0);
        double sync_from_device_end = realtime();
        core->sync_from_device_ultra += sync_from_device_end - sync_from_device_start;
    }

}


pthread_t align_ultra_long_async(pthread_arg_t **pt_args_ptr,core_t* core, db_t* db, int* ultra_long_reads_id, int ultra_long_reads_num) {

    assert(*pt_args_ptr==NULL);
    *pt_args_ptr = (pthread_arg_t *)malloc(sizeof(pthread_arg_t));
    pthread_arg_t *pt_args=*pt_args_ptr;
    MALLOC_CHK(pt_args);
    pt_args->core = core;
    pt_args->db = db;
    pt_args->starti = 0;
    pt_args->endi = ultra_long_reads_num;
    pt_args->ultra_long_reads_id = ultra_long_reads_id;

    pthread_t tid;
    int ret = pthread_create(&tid, NULL, align_ultra_long,(void*)(pt_args));
    NEG_CHK(ret);

    return tid;
}

void* align_ultra_long(void* voidargs){

    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    int* ultra_long_reads_id   = args->ultra_long_reads_id;
    int ultra_long_reads_num = args->endi;

    xrt::kernel kernel = core->kernel_2;
    xrt::device device = core->device;
    xrt::bo event_mean_bo = core->event_mean_bo_2;
    xrt::bo read_bo = core->read_bo_2;
    xrt::bo n_align_bo = core->n_align_bo_2;
    xrt::bo aligned_pairs_bo = core->aligned_pairs_bo_2;
    xrt::bo trace_table_bo = core->trace_table_bo_2;
    xrt::bo left_out_bo = core->left_out_bo_2;
    xrt::bo length_int_bo = core->length_int_bo_2;
    xrt::bo scaling_info_bo = core->scaling_info_bo_2;
    xrt::bo lp_info_bo = core->lp_info_bo_2;

    int event_offset;
    int reads_offset;
    int n_bands;
     int kmer_size = core->kmer_size;
    float* event_mean = core->event_mean_2;
    char* reads = core->reads_2;
    int* length_int = core->length_int_2;
    float* scaling_info = core->scaling_info_2;
    float* lp_info =core->lp_info_2;
    int* n_align = core->n_align_2;
    int* aligned_pairs = core->aligned_pairs_2;
    int batch_num = 1;

    for(int i=0; i<ultra_long_reads_num; i++)
    {
        double fill_start = realtime();
        int idx = ultra_long_reads_id[i];
        int read_length = db->read_len[idx];
        int event_length = db->et[idx].n;

        if(event_length%64!=0)
            event_offset = (event_length/64 + 1) * 64;
        else
            event_offset = event_length;

        if(read_length%64!=0)
            reads_offset = (read_length/64 + 1) * 64;
        else
            reads_offset = read_length;

        n_bands = event_offset + reads_offset + 16;

     //   fprintf(stderr,"reads_offset %d event_offset %d n_bands %d\n", reads_offset, event_offset, n_bands);
     
        memcpy(&reads[0], db->read[idx], read_length);
        
        for(int j=0; j< event_length; j++)
            event_mean[j] = db->et[idx].event[j].mean; 

        length_int[0] =  read_length;
        length_int[1] = event_length;

        scaling_info[0] =  db->scalings[idx].log_var;
        scaling_info[1] = db->scalings[idx].scale;
        scaling_info[2] = db->scalings[idx].shift;
        scaling_info[3] = db->scalings[idx].var;

        int n_kmers = read_length - kmer_size + 1;
        float events_per_kmer = (float)event_length / (float)n_kmers;
        float p_stay = (float)1.0 - ((float)1.0 / (events_per_kmer + (float)1.0));
        float epsilon = 1e-10;
        lp_info[0] = log(epsilon);
        lp_info[1]= log(p_stay);
        lp_info[2] = log((float)1.0 - exp(lp_info[0]) - exp(lp_info[1])); 

        double fill_end = realtime();
        core->fill_bucket_time_ultra += fill_end - fill_start;

        double sync_to_device_start = realtime();
    //    fprintf(stderr,"sync to device\n");
        size_t event_mean_size_bytes = sizeof(float) * event_offset;
        size_t read_size_bytes = sizeof(char) *  reads_offset;
        size_t n_align_size_bytes=sizeof(int);
        size_t aligned_pairs_size_bytes=sizeof(int)* 2 * n_bands;
        size_t length_int_size_bytes=sizeof(int) * 2;
        size_t scaling_size_bytes = sizeof(float) * 4;
        size_t lp_size_bytes = sizeof(float) * 3;//

        event_mean_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE,event_mean_size_bytes, 0 );
        read_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, read_size_bytes, 0);
        length_int_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, length_int_size_bytes, 0);
        scaling_info_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, scaling_size_bytes, 0);
        lp_info_bo.sync(XCL_BO_SYNC_BO_TO_DEVICE, lp_size_bytes, 0);

        double sync_to_device_end = realtime();
        core->sync_to_device_ultra += sync_to_device_end - sync_to_device_start;

        double kernel_start = realtime();
        core->call_kernel_2_num++;
     //   fprintf(stderr,"run kernel\n");
        auto run = kernel(0, event_mean_bo,read_bo,n_align_bo,aligned_pairs_bo,NULL,trace_table_bo, left_out_bo, length_int_bo,scaling_info_bo,lp_info_bo,event_offset,reads_offset,batch_num);

        run.wait();
        double kernel_end = realtime();
        core->kernel_run_ultra += kernel_end - kernel_start;

        double sync_from_device_start = realtime(); 
      //  fprintf(stderr,"sync to host\n");
        aligned_pairs_bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE, aligned_pairs_size_bytes,0);
        n_align_bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE, n_align_size_bytes, 0);
        double sync_from_device_end = realtime();
        core->sync_from_device_ultra += sync_from_device_end - sync_from_device_start;

       /* //result from device
        if(read_length == 81616 && event_length == 171579)
        {
            fprintf(stderr,"read_len %d event len %d n_pairs %d\n", read_length, event_length, n_align[0]);
            for(int j=0;j<n_bands;j++)
               fprintf(stderr,"%d %d\n",aligned_pairs[j*2], aligned_pairs[j*2+1]);
        }*/
    }

    return NULL;

}


void align_ultra_long_async_join(pthread_arg_t *pt_args, pthread_t tid) {
    int ret = pthread_join(tid, NULL);
    NEG_CHK(ret);
    assert(pt_args);
    free(pt_args);
    return;
}

pthread_t align_cpu_reads_async(pthread_arg_t **pt_args_ptr,core_t* core, db_t* db, int* cpu_reads_id, int cpu_reads_num)
{
    assert(*pt_args_ptr==NULL);
    *pt_args_ptr = (pthread_arg_t *)malloc(sizeof(pthread_arg_t));
    pthread_arg_t *pt_args=*pt_args_ptr;
    MALLOC_CHK(pt_args);
    pt_args->core = core;
    pt_args->db = db;
    pt_args->starti = 0;
    pt_args->endi = cpu_reads_num;
    pt_args->ultra_long_reads_id = cpu_reads_id;

    pthread_t tid;
    int ret = pthread_create(&tid, NULL, align_cpu_reads,(void*)(pt_args));
    NEG_CHK(ret);

    return tid;

}

void* align_cpu_reads(void* voidargs){

    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    int* cpu_reads_id   = args->ultra_long_reads_id;
    int cpu_reads_num = args->endi;

    double cpu_start = realtime();

    //create threads
    pthread_t tids[core->opt.num_thread];
    pthread_arg_t pt_args[core->opt.num_thread];
    int32_t t, ret;
    int32_t i = 0;
    int32_t num_thread = core->opt.num_thread;
    int32_t step = (cpu_reads_num + num_thread - 1) / num_thread;
    //todo : check for higher num of threads than the data
    //current works but many threads are created despite

    //set the data structures
    for (t = 0; t < num_thread; t++) {
        pt_args[t].core = core;
        pt_args[t].db = db;
        pt_args[t].starti = i;
        i += step;
        if (i > cpu_reads_num) {
            pt_args[t].endi = cpu_reads_num;
        } else {
            pt_args[t].endi = i;
        }
        pt_args[t].func=align_single;
        pt_args[t].ultra_long_reads_id = cpu_reads_id;

    }

    //create threads
    for(t = 0; t < core->opt.num_thread; t++){
        ret = pthread_create(&tids[t], NULL, pthread_cpu_reads_single,
                                (void*)(&pt_args[t]));
        NEG_CHK(ret);
    }

    //pthread joining
    for (t = 0; t < core->opt.num_thread; t++) {
        int ret = pthread_join(tids[t], NULL);
        NEG_CHK(ret);
    }

    double cpu_end = realtime();

}

void* pthread_cpu_reads_single(void* voidargs) {

    double realtime1 = realtime();

    int32_t i,j;
    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;

    for (i = args->starti; i < args->endi; i++) {
        j=args->ultra_long_reads_id[i];
        args->func(core,db,j);
    }

    pthread_exit(0);
}

void align_cpu_reads_async_join(pthread_arg_t *pt_args, pthread_t tid) {
    int ret = pthread_join(tid, NULL);
    NEG_CHK(ret);
    assert(pt_args);
    free(pt_args);
    return;
}

//two align top functions for regular reads
void fill_bucket(core_t* core, db_t* db, int bucket, int kmer_size, int max_read_length , int max_event_length, int reads_num, int reads_id[BUCKET_SIZE])
{
    xrt::device device = core->device;
    xrt::kernel kernel_0 = core->kernel_0;
    xrt::kernel kernel_1 = core->kernel_1;

    xrt::bo event_mean_bo_0 = core->event_mean_bo_0;
    xrt::bo read_bo_0 = core->read_bo_0;
    xrt::bo n_align_bo_0 = core->n_align_bo_0;
    xrt::bo aligned_pairs_bo_0 = core->aligned_pairs_bo_0;
    xrt::bo trace_table_bo_0 = core->trace_table_bo_0;
    xrt::bo left_out_bo_0 = core->left_out_bo_0; 
    xrt::bo length_int_bo_0 = core->length_int_bo_0;
    xrt::bo scaling_info_bo_0 = core->scaling_info_bo_0;
    xrt::bo lp_info_bo_0 = core->lp_info_bo_0;

    xrt::bo event_mean_bo_1 = core->event_mean_bo_1;
    xrt::bo read_bo_1 = core->read_bo_1;
    xrt::bo n_align_bo_1 = core->n_align_bo_1;
    xrt::bo aligned_pairs_bo_1 = core->aligned_pairs_bo_1;
    xrt::bo trace_table_bo_1 = core->trace_table_bo_1;
    xrt::bo left_out_bo_1 = core->left_out_bo_1; 
    xrt::bo length_int_bo_1 = core->length_int_bo_1;
    xrt::bo scaling_info_bo_1 = core->scaling_info_bo_1;
    xrt::bo lp_info_bo_1 = core->lp_info_bo_1;

    size_t event_offset;
    size_t reads_offset;
    size_t n_bands;
    float* event_mean_0 = core->event_mean_0;
    char* reads_0 = core->reads_0;
    int* length_int_0 = core->length_int_0;
    memset(length_int_0, 0, BUCKET_SIZE/2 * 2 * sizeof(int));
    float* scaling_info_0 = core->scaling_info_0;
    float* lp_info_0 = core->lp_info_0;
    int* n_align_0 = core->n_align_0;
    int* aligned_pairs_0 = core->aligned_pairs_0;
    

    float* event_mean_1 = core->event_mean_1;
    char* reads_1 = core->reads_1;
    int* length_int_1 = core->length_int_1;
    memset(length_int_1, 0, BUCKET_SIZE/2 * 2 * sizeof(int));
    float* scaling_info_1 = core->scaling_info_1;
    float* lp_info_1 = core->lp_info_1;
    int* n_align_1 = core->n_align_1;
    int* aligned_pairs_1 = core->aligned_pairs_1;
  

   if(max_event_length%64!=0)
        event_offset = (max_event_length/64 + 1) * 64;
    else
        event_offset = max_event_length;

   if(max_read_length%64!=0)
        reads_offset = (max_read_length/64 + 1) * 64;
    else
        reads_offset = max_read_length;

    n_bands = event_offset + reads_offset + 2;

    int kernel_num = 2;
    int reads_num_0 = reads_num/2;
    int reads_num_1 = reads_num - reads_num_0;
    int batch_num_0 = reads_num_0%BATCH_SIZE==0?reads_num_0/BATCH_SIZE:reads_num_0/BATCH_SIZE + 1;
    int batch_num_1 = reads_num_1%BATCH_SIZE==0?reads_num_1/BATCH_SIZE:reads_num_1/BATCH_SIZE + 1;
    if(reads_num<=8)
    {
        kernel_num = 1;
        reads_num_0 = reads_num;
        reads_num_1 = 0;
        batch_num_0 = 1;
        batch_num_1 = 0;
    }

 //  fprintf(stderr," reads_num %d batch_num %d event_offset %ld reads_offset %ld n_bands %ld\n", reads_num, batch_num, event_offset, reads_offset, n_bands);
   double fill_start = realtime();

   for(int i=0; i < reads_num_0; i++ )
   {
      int idx = reads_id[i];
      int read_length = db->read_len[idx];
      int event_length = db->et[idx].n;
     
      memcpy(&reads_0[(int64_t)i * reads_offset], db->read[idx], read_length);

      for(int j=0; j< event_length; j++)
            event_mean_0[(int64_t)i * event_offset + j] = db->et[idx].event[j].mean; 

       int start = (i/BATCH_SIZE) * BATCH_SIZE * 2;
       int id  = i%BATCH_SIZE;
        
       length_int_0[start + id] =  read_length;
       length_int_0[start + id + BATCH_SIZE] = event_length;

       start = (i/BATCH_SIZE) * BATCH_SIZE * 4;

       scaling_info_0[start + id] =  db->scalings[idx].log_var;
       scaling_info_0[start + id + BATCH_SIZE] =db->scalings[idx].scale;
       scaling_info_0[start + id + 2 * BATCH_SIZE] = db->scalings[idx].shift;
       scaling_info_0[start + id + 3 * BATCH_SIZE] = db->scalings[idx].var;
     
       
       int n_kmers = read_length - kmer_size + 1;
       float events_per_kmer = (float)event_length / (float)n_kmers;
       float p_stay = (float)1.0 - ((float)1.0 / (events_per_kmer + (float)1.0));
       float epsilon = 1e-10;
       start  = (i/BATCH_SIZE) * BATCH_SIZE * 3;
       lp_info_0[start + id] = log(epsilon);
       lp_info_0[start + id + BATCH_SIZE]= log(p_stay);
       lp_info_0[start + id + 2 * BATCH_SIZE] = log((float)1.0 - exp(lp_info_0[start + id]) - exp(lp_info_0[start + id + BATCH_SIZE])); 
   }

    for(int i=0; i < reads_num_1; i++ )
    {
      int idx = reads_id[reads_num_0 + i];
      int read_length = db->read_len[idx];
      int event_length = db->et[idx].n;
     
      memcpy(&reads_1[(int64_t)i * reads_offset], db->read[idx], read_length);

      for(int j=0; j< event_length; j++)
            event_mean_1[(int64_t)i * event_offset + j] = db->et[idx].event[j].mean; 

       int start = (i/BATCH_SIZE) * BATCH_SIZE * 2;
       int id  = i%BATCH_SIZE;
        
       length_int_1[start + id] =  read_length;
       length_int_1[start + id + BATCH_SIZE] = event_length;

       start = (i/BATCH_SIZE) * BATCH_SIZE * 4;

       scaling_info_1[start + id] =  db->scalings[idx].log_var;
       scaling_info_1[start + id + BATCH_SIZE] =db->scalings[idx].scale;
       scaling_info_1[start + id + 2 * BATCH_SIZE] = db->scalings[idx].shift;
       scaling_info_1[start + id + 3 * BATCH_SIZE] = db->scalings[idx].var;
     
       
       int n_kmers = read_length - kmer_size + 1;
       float events_per_kmer = (float)event_length / (float)n_kmers;
       float p_stay = (float)1.0 - ((float)1.0 / (events_per_kmer + (float)1.0));
       float epsilon = 1e-10;
       start  = (i/BATCH_SIZE) * BATCH_SIZE * 3;
       lp_info_1[start + id] = log(epsilon);
       lp_info_1[start + id + BATCH_SIZE]= log(p_stay);
       lp_info_1[start + id + 2 * BATCH_SIZE] = log((float)1.0 - exp(lp_info_1[start + id]) - exp(lp_info_1[start + id + BATCH_SIZE])); 
   }

   double fill_end = realtime();

   core->fill_bucket_time += fill_end - fill_start;

   double sync_to_device_start = realtime();
   size_t event_mean_size_bytes = sizeof(float) * batch_num_0 * BATCH_SIZE * event_offset;
   size_t read_size_bytes = sizeof(char) * batch_num_0 * BATCH_SIZE * reads_offset;//
   size_t length_int_size_bytes=sizeof(int) * batch_num_0 * BATCH_SIZE * 2;
   size_t scaling_size_bytes = sizeof(float) * batch_num_0 * BATCH_SIZE * 4;
   size_t lp_size_bytes = sizeof(float) * batch_num_0 * BATCH_SIZE * 3;//

//   fprintf(stderr,"kernel 0 sync to device\n");
   
   event_mean_bo_0.sync(XCL_BO_SYNC_BO_TO_DEVICE,event_mean_size_bytes, 0 );
   read_bo_0.sync(XCL_BO_SYNC_BO_TO_DEVICE, read_size_bytes, 0);
   length_int_bo_0.sync(XCL_BO_SYNC_BO_TO_DEVICE, length_int_size_bytes, 0);
   scaling_info_bo_0.sync(XCL_BO_SYNC_BO_TO_DEVICE, scaling_size_bytes, 0);
   lp_info_bo_0.sync(XCL_BO_SYNC_BO_TO_DEVICE, lp_size_bytes, 0);
   double sync_to_device_end = realtime();
   core->sync_to_device += sync_to_device_end - sync_to_device_start;

  double kernel_start = realtime();
   core->call_kernel_0_num++;
//  fprintf(stderr,"kernel 0 run\n");
   auto run_0 = kernel_0(0, event_mean_bo_0,read_bo_0,n_align_bo_0,aligned_pairs_bo_0,NULL,trace_table_bo_0, left_out_bo_0, length_int_bo_0,scaling_info_bo_0,lp_info_bo_0,event_offset,reads_offset,batch_num_0);
   

   if(batch_num_1 > 0)
   {
   //    fprintf(stderr,"kernel 1 sync to device\n");
        double sync_to_device_start = realtime();
        event_mean_size_bytes = sizeof(float) * batch_num_1 * BATCH_SIZE * event_offset;
        read_size_bytes = sizeof(char) * batch_num_1 * BATCH_SIZE * reads_offset;//
        length_int_size_bytes=sizeof(int) * batch_num_1 * BATCH_SIZE * 2;
        scaling_size_bytes = sizeof(float) * batch_num_1 * BATCH_SIZE * 4;
        lp_size_bytes = sizeof(float) * batch_num_1 * BATCH_SIZE * 3;//

        event_mean_bo_1.sync(XCL_BO_SYNC_BO_TO_DEVICE,event_mean_size_bytes, 0 );
        read_bo_1.sync(XCL_BO_SYNC_BO_TO_DEVICE, read_size_bytes, 0);
        length_int_bo_1.sync(XCL_BO_SYNC_BO_TO_DEVICE, length_int_size_bytes, 0);
        scaling_info_bo_1.sync(XCL_BO_SYNC_BO_TO_DEVICE, scaling_size_bytes, 0);
        lp_info_bo_1.sync(XCL_BO_SYNC_BO_TO_DEVICE, lp_size_bytes, 0);
        double sync_to_device_end = realtime();
        core->sync_to_device += sync_to_device_end - sync_to_device_start;

        core->call_kernel_1_num++;

     //   fprintf(stderr,"kernel 1 run\n");
        auto run_1= kernel_1(0, event_mean_bo_1,read_bo_1,n_align_bo_1,aligned_pairs_bo_1,NULL,trace_table_bo_1, left_out_bo_1, length_int_bo_1,scaling_info_bo_1,lp_info_bo_1,event_offset,reads_offset,batch_num_1);
        run_1.wait();
   }
   
    run_0.wait();
    double kernel_end = realtime();
    core->kernel_run += kernel_end - kernel_start;

    double sync_from_device_start = realtime();
    size_t n_align_size_bytes=sizeof(int)* batch_num_0 * BATCH_SIZE;
    size_t aligned_pairs_size_bytes=sizeof(int) *batch_num_0 * BATCH_SIZE * 2 * n_bands;
  //  fprintf(stderr,"kernel0 sync to host\n");
    aligned_pairs_bo_0.sync(XCL_BO_SYNC_BO_FROM_DEVICE, aligned_pairs_size_bytes, 0);
    n_align_bo_0.sync(XCL_BO_SYNC_BO_FROM_DEVICE, n_align_size_bytes, 0);


    if(batch_num_1 > 0)
    {
    //    fprintf(stderr,"kernel1 sync to host\n");
        n_align_size_bytes=sizeof(int)* batch_num_1 * BATCH_SIZE;
        aligned_pairs_size_bytes=sizeof(int) *batch_num_1 * BATCH_SIZE * 2 * n_bands;
        aligned_pairs_bo_1.sync(XCL_BO_SYNC_BO_FROM_DEVICE, aligned_pairs_size_bytes, 0);
        n_align_bo_1.sync(XCL_BO_SYNC_BO_FROM_DEVICE, n_align_size_bytes, 0);
    }

    double sync_from_device_end = realtime();
    core->sync_from_device += sync_from_device_end - sync_from_device_start;

   /* for(int i=0; i<reads_num_0 ; i++)
    {
        int idx = reads_id[i];
        int read length = db->read len[idx];
        int event length = db->et[idx].n;
    
        fprintf(stderr,"read length %d event length %d n pairs %d\n", read_length, event_length, n_align_0[i]);
    }
*/
   

    //pthread fill in struct
     pthread_reverse_aligned_pairs(core,  db, reads_num_0, reads_num_1, n_bands,reads_id);

}

void pthread_reverse_aligned_pairs(core_t* core, db_t* db, int reads_num_0, int reads_num_1, int n_bands,int* reads_id){

        //create threads
        pthread_t tids[core->opt.num_thread];
        pthread_arg_rev pt_args[core->opt.num_thread];
        int32_t t, ret;
        int32_t i = 0;
        int reads_num = reads_num_0 + reads_num_1;
        int32_t num_thread = core->opt.num_thread;
        int32_t step = (reads_num + num_thread - 1) / num_thread;
        //todo : check for higher num of threads than the data
        //current works but many threads are created despite

        //set the data structures
        for (t = 0; t < num_thread; t++) {
            pt_args[t].core = core;
            pt_args[t].db = db;
            pt_args[t].reads_id = reads_id;
            pt_args[t].reads_num_0 = reads_num_0;
            pt_args[t].reads_num_1 = reads_num_1;
            pt_args[t].n_bands = n_bands;
            pt_args[t].starti = i;
            i += step;
            if (i > reads_num) {
                pt_args[t].endi = reads_num;
            } else {
                pt_args[t].endi = i;
            }

        }

        //create threads
        for(t = 0; t < core->opt.num_thread; t++){
            ret = pthread_create(&tids[t], NULL, pthread_single_reverse, (void*)(&pt_args[t]));
            NEG_CHK(ret);
        }

        //pthread joining
        for (t = 0; t < core->opt.num_thread; t++) {
            int ret = pthread_join(tids[t], NULL);
            NEG_CHK(ret);
        }
        
    
}

void* pthread_single_reverse(void* voidargs) {
    int32_t i;
    pthread_arg_rev* args = (pthread_arg_rev*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    int* reads_id = args->reads_id;
    int reads_num_0 = args->reads_num_0;
    int reads_num_1 = args->reads_num_1;
    int n_bands = args->n_bands;

    int* n_align_0 = core->n_align_0;
    int* aligned_pairs_0 = core->aligned_pairs_0;
    int* n_align_1 = core->n_align_1;
    int* aligned_pairs_1 = core->aligned_pairs_1;

    for (i = args->starti; i < args->endi; i++) {
        int idx = reads_id[i];
        int read_length = db->read_len[idx];
        int event_length = db->et[idx].n;

        if(i<reads_num_0)
        {

            db->n_event_align_pairs[idx] = n_align_0[i];
             int outIndex = db->n_event_align_pairs[idx]-1;
         //    
            if(db->n_event_align_pairs[idx]>0)
            {     
                if(read_length==33320&&event_length==66345)
                     fprintf(stderr,"1. read %d event %d npairs %d outIndex %d\n", read_length, event_length, db->n_event_align_pairs[idx] , outIndex);
                for(int j=read_length+event_length+2; j>=0; j--)
                {
                    int batch = i/BATCH_SIZE;
                    int read = i%BATCH_SIZE;
                    int ref_pos = aligned_pairs_0[batch * BATCH_SIZE * 2 *n_bands + j*BATCH_SIZE*2 + read* 2];
                    int read_pos = aligned_pairs_0[batch * BATCH_SIZE * 2 *n_bands + j*BATCH_SIZE*2 + read*2 + 1];
                    if(ref_pos>=0 && outIndex>=0)
                    {
                        db->event_align_pairs[idx][outIndex].ref_pos = ref_pos;
                        db->event_align_pairs[idx][outIndex].read_pos = read_pos;
                        if(read_length==33320&&event_length==66345)
                            fprintf(stderr,"%d %d\n",ref_pos, read_pos);
                        outIndex--;
                    }
                }
            }
        }else
        {

            db->n_event_align_pairs[idx] = n_align_1[i-reads_num_0];
            int outIndex = db->n_event_align_pairs[idx]-1;
          //   fprintf(stderr,"event %d npairs %d outIndex %d\n", event_length, db->n_event_align_pairs[idx] , outIndex);
            if(db->n_event_align_pairs[idx] > 0)
            {     
                if(read_length==33320&&event_length==66345)
                    fprintf(stderr,"2. read %d event %d npairs %d outIndex %d\n", read_length, event_length, db->n_event_align_pairs[idx] , outIndex);
                for(int j=read_length+event_length+2; j>=0; j--)
                {
                    int batch = i/BATCH_SIZE;
                    int read = i%BATCH_SIZE;
                    int ref_pos = aligned_pairs_1[batch * BATCH_SIZE * 2 *n_bands + j*BATCH_SIZE*2 + read* 2];
                    int read_pos = aligned_pairs_1[batch * BATCH_SIZE * 2 *n_bands + j*BATCH_SIZE*2 + read*2 + 1];
               
                    if(ref_pos>=0&&outIndex>=0)
                    {
                        db->event_align_pairs[idx][outIndex].ref_pos = ref_pos;
                        db->event_align_pairs[idx][outIndex].read_pos = read_pos;
                        if(read_length==33320&&event_length==66345)
                            fprintf(stderr,"%d %d\n",ref_pos, read_pos);
                        outIndex--;
                    }
                }
            }
        }
    }

    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}



void eventalign_single(core_t* core, db_t* db, int32_t i){
    realign_read(db->event_alignment_result[i], &(db->eventalign_summary[i]),core->event_summary_fp, db->fasta_cache[i],core->m_hdr,
                  db->bam_rec[i],db->read_len[i],
                  i,
                  core->clip_start,
                  core->clip_end,
                  &(db->et[i]), core->model,core->kmer_size, db->base_to_event_map[i],db->scalings[i],db->events_per_base[i],db->sig[i]->sample_rate);

    char* qname = bam_get_qname(db->bam_rec[i]);
    char* contig = core->m_hdr->target_name[db->bam_rec[i]->core.tid];
    std::vector<event_alignment_t> *event_alignment_result = db->event_alignment_result[i];
    int8_t print_read_names = (core->opt.flag & F5C_PRINT_RNAME) ? 1 : 0;
    int8_t scale_events = (core->opt.flag & F5C_SCALE_EVENTS) ? 1 : 0;
    int8_t collapse_events = (core->opt.flag & F5C_COLLAPSE_EVENTS) ? 1 : 0;
    int8_t write_samples = (core->opt.flag & F5C_PRINT_SAMPLES) ? 1 : 0;
    int8_t write_signal_index = (core->opt.flag & F5C_PRINT_SIGNAL_INDEX) ? 1 : 0;
    int8_t sam_output = (core->opt.flag & F5C_SAM) ? 1 : 0;

    if(sam_output==0){
        db->event_alignment_result_str[i] = emit_event_alignment_tsv(0,&(db->et[i]),core->model,core->kmer_size, db->scalings[i],*event_alignment_result, print_read_names, scale_events, write_samples, write_signal_index, collapse_events,
                   db->read_idx[i], qname, contig, db->sig[i]->sample_rate, db->sig[i]->rawptr);

    }
}

void meth_single(core_t* core, db_t* db, int32_t i){
    if(!db->read_stat_flag[i]){
        if(core->mode==0){
            calculate_methylation_for_read(db->site_score_map[i], db->fasta_cache[i], db->bam_rec[i], db->read_len[i], db->et[i].event, db->base_to_event_map[i],
            db->scalings[i], core->cpgmodel, core->kmer_size, db->events_per_base[i]);
        }
        else if (core->mode==1){
            eventalign_single(core,db,i);
        }
    }
}


void process_single(core_t* core, db_t* db,int32_t i) {
    event_single(core,db,i);
   align_single(core, db,i);
   scaling_single(core,db,i);
   meth_single(core,db,i);
}

/* completely process a data batch
   (all steps: event detection, adaptive banded event alignment, ...., HMM) */
void process_db(core_t* core, db_t* db) {

    double process_start = realtime();

    if((core->opt.flag&F5C_SEC_PROF) || (!(core->opt.flag & F5C_DISABLE_CUDA))){

        double realtime0=core->realtime0;

        double event_start = realtime();
        pthread_db(core,db,event_single);
        double event_end = realtime();
        core->event_time += (event_end-event_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Events computed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double align_start = realtime();
        align_db(core, db);
        double align_end = realtime();
        core->align_time += (align_end-align_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Banded alignment done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double est_scale_start = realtime();
        pthread_db(core,db,scaling_single);
        double est_scale_end = realtime();
        core->est_scale_time += (est_scale_end-est_scale_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Scaling calibration done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double meth_start = realtime();
        pthread_db(core, db, meth_single);
        double meth_end = realtime();
        core->meth_time += (meth_end-meth_start);

        fprintf(stderr, "[%s::%.3f*%.2f] HMM done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));


    }
    else{
        #ifdef FPGA_MOD

             double realtime0=core->realtime0;

             double event_start = realtime();
             pthread_db(core,db,event_single);
                double event_end = realtime();
                core->event_time += (event_end-event_start);

            fprintf(stderr, "[%s::%.3f*%.2f] Events computed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
            
            double align_start = realtime();
            if(db->batch_idx == 0)
                transfer_model_to_FPGA(core);
            align_db_FPGA_by_len(core,db);
            double align_end = realtime();
            core->align_FPGA_time += (align_end - align_start);
            fprintf(stderr, "[%s::%.3f*%.2f] Banded alignment on FPGA done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

            for(int i=0;i<db->n_bam_rec;i++)
            {
                int read_length = db->read_len[i];
                int event_length = db->et[i].n;
                if(read_length==33320&&event_length==66345)
                {
                    fprintf(stderr,"read_len %d event_len %d n_pairs %d \n", read_length, event_length, db->n_event_align_pairs[i]);
                    for(int j=0; j<db->n_event_align_pairs[i]; j++)
                        fprintf(stderr,"%d %d\n", db->event_align_pairs[i][j].ref_pos, db->event_align_pairs[i][j].read_pos);
                }
            }

            double est_scale_start = realtime();
            pthread_db(core,db,scaling_single);
            double est_scale_end = realtime();
            core->est_scale_time += (est_scale_end-est_scale_start);

            fprintf(stderr, "[%s::%.3f*%.2f] Scaling calibration done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

            double meth_start = realtime();
            pthread_db(core, db, meth_single);
            double meth_end = realtime();
            core->meth_time += (meth_end-meth_start);

            fprintf(stderr, "[%s::%.3f*%.2f] HMM done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
            

    	#else
        if (core->opt.num_thread == 1) {
            int32_t i=0;
            int sum_length=0;
            for (i = 0; i < db->n_bam_rec; i++) {
                process_single(core,db,i);
                sum_length+= db->et[i].n +1 + db->read_len[i] - core->kmer_size+ 1 + 1;
            }


        }
        else {
            pthread_db(core,db,process_single);
        }
       #endif 
    }

    double process_end= realtime();
    core->process_db_time += (process_end-process_start);

    return;
}

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db) {

    double output_start = realtime();

    if (core->opt.flag & F5C_PRINT_EVENTS) {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++) {
            printf(">%s\tLN:%d\tEVENTSTART:%d\tEVENTEND:%d\n",
                   bam_get_qname(db->bam_rec[i]), (int)db->et[i].n,
                   (int)db->et[i].start, (int)db->et[i].end);
            uint32_t j = 0;
            for (j = 0; j < db->et[i].n; j++) {
                printf("{%d,%f,%f,%f}\t", (int)db->et[i].event[j].start,
                       db->et[i].event[j].length, db->et[i].event[j].mean,
                       db->et[i].event[j].stdv);
            }
            printf("\n");
        }
    }
    if (core->opt.flag & F5C_PRINT_BANDED_ALN) {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++) {
            if((db->read_stat_flag[i]) & FAILED_ALIGNMENT){
                continue;
            }
            printf(">%s\tN_ALGN_PAIR:%d\t{ref_pos,read_pos}\n",
                   bam_get_qname(db->bam_rec[i]),
                   (int)db->n_event_align_pairs[i]);
            AlignedPair* event_align_pairs = db->event_align_pairs[i];
            int32_t j = 0;
            for (j = 0; j < db->n_event_align_pairs[i]; j++) {
                printf("{%d,%d}\t", event_align_pairs[j].ref_pos,
                       event_align_pairs[j].read_pos);
            }
            printf("\n");
        }
    }

    if (core->opt.flag & F5C_PRINT_SCALING) {
        int32_t i = 0;
        printf("read\tshift\tscale\tvar\n");

        for (i = 0; i < db->n_bam_rec; i++) {
            if((db->read_stat_flag[i])&(FAILED_ALIGNMENT|FAILED_CALIBRATION)){
                continue;
            }
            printf("%s\t%.2lf\t%.2lf\t%.2lf\n", bam_get_qname(db->bam_rec[i]),
                   db->scalings[i].shift, db->scalings[i].scale,
                   db->scalings[i].var);
        }
    }

    core->sum_bases += db->sum_bases;
    core->total_reads += db->total_reads;
    core->bad_fast5_file += db->bad_fast5_file;
    core->ultra_long_skipped += db->ultra_long_skipped;

    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; i++){
        if(!db->read_stat_flag[i]){
            char* qname = bam_get_qname(db->bam_rec[i]);
            char* contig = core->m_hdr->target_name[db->bam_rec[i]->core.tid];

            if(core->mode==0) {
                std::map<int, ScoredSite> *site_score_map = db->site_score_map[i];
                // write all sites for this read
                for(auto iter = site_score_map->begin(); iter != site_score_map->end(); ++iter) {

                    const ScoredSite& ss = iter->second;
                    double sum_ll_m = ss.ll_methylated[0]; //+ ss.ll_methylated[1];
                    double sum_ll_u = ss.ll_unmethylated[0]; //+ ss.ll_unmethylated[1];
                    double diff = sum_ll_m - sum_ll_u;

                    // fprintf(stderr, "%s\t%d\t%d\t", ss.chromosome.c_str(), ss.start_position, ss.end_position);
                    // fprintf(stderr, "%s\t%.2lf\t", qname, diff);
                    // fprintf(stderr, "%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
                    // fprintf(stderr, "%d\t%d\t%s\n", ss.strands_scored, ss.n_cpg, ss.sequence.c_str());

                    // output only if inside the window boundaries
                    if( !( (core->clip_start != -1 && ss.start_position < core->clip_start) ||
                        (core->clip_end != -1 && ss.end_position >= core->clip_end) ) ) {
                        if(core->opt.meth_out_version==1){
                            printf("%s\t%d\t%d\t", contig, ss.start_position, ss.end_position);
                        }
                        else if(core->opt.meth_out_version==2){
                            printf("%s\t%c\t%d\t%d\t", contig, bam_is_rev(db->bam_rec[i]) ? '-' : '+', ss.start_position, ss.end_position);
                        }
                        printf("%s\t%.2lf\t", qname, diff);
                        printf("%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
                        printf("%d\t%d\t%s\n", ss.strands_scored, ss.n_cpg, ss.sequence.c_str());
                    }

                }
            }

            else if(core->mode==1){
                FILE* summary_fp = core->event_summary_fp;
                EventalignSummary summary = db->eventalign_summary[i];
                scalings_t scalings = db->scalings[i];
                if(summary_fp != NULL && summary.num_events > 0) {
                    size_t strand_idx = 0;
                    std::string fast5_path_str = core->readbb->get_signal_path(qname);
                    fprintf(summary_fp, "%ld\t%s\t", (long)(db->read_idx[i]), qname);
                    fprintf(summary_fp, "%s\t%s\t%s\t",fast5_path_str.c_str(), (core->opt.flag & F5C_RNA) ? "rna" : "dna", strand_idx == 0 ? "template" : "complement" );
                    fprintf(summary_fp, "%d\t%d\t%d\t%d\t", summary.num_events, summary.num_steps, summary.num_skips, summary.num_stays);
                    fprintf(summary_fp, "%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", summary.sum_duration/(db->sig[i]->sample_rate), scalings.shift, scalings.scale, 0.0, scalings.var);
                }
                std::vector<event_alignment_t> *event_alignment_result = db->event_alignment_result[i];
                char *event_alignment_result_str = db->event_alignment_result_str[i];

                // int8_t print_read_names = (core->opt.flag & F5C_PRINT_RNAME) ? 1 : 0;
                // int8_t scale_events = (core->opt.flag & F5C_SCALE_EVENTS) ? 1 : 0;
                // int8_t write_samples = (core->opt.flag & F5C_PRINT_SAMPLES) ? 1 : 0;
                // int8_t write_signal_index = (core->opt.flag & F5C_PRINT_SIGNAL_INDEX) ? 1 : 0;
                int8_t sam_output = (core->opt.flag & F5C_SAM) ? 1 : 0;

                if(sam_output==0){
                    // emit_event_alignment_tsv(stdout,0,&(db->et[i]),core->model,db->scalings[i],*event_alignment_result, print_read_names, scale_events, write_samples, write_signal_index,
                    //           db->read_idx[i], qname, contig, db->sig[i]->sample_rate, db->sig[i]->rawptr);
                    fputs(event_alignment_result_str,stdout);
                }
                else{
                    emit_event_alignment_sam(core->sam_output , qname, core->m_hdr, db->bam_rec[i], *event_alignment_result);
                }
            }
        }
        else{
            if((db->read_stat_flag[i])&FAILED_CALIBRATION){
                core->failed_calibration_reads++;
            }
            else if ((db->read_stat_flag[i])&FAILED_ALIGNMENT){
                core->failed_alignment_reads++;
            }
            else if ((db->read_stat_flag[i])&FAILED_QUALITY_CHK){
                core->qc_fail_reads++;
            }
            else{
                assert(0);
            }
        }
    }
    //core->read_index = core->read_index + db->n_bam_rec;
    double output_end = realtime();
    core->output_time += (output_end-output_start);

}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; ++i) {
        bam_destroy1(db->bam_rec[i]);
        db->bam_rec[i] = bam_init1();
        free(db->fasta_cache[i]);
        free(db->read[i]);
        free(db->sig[i]->rawptr);
        free(db->sig[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);

        #ifndef FPGA_MOD
          free(db->base_to_event_map[i]);
        #endif

        delete db->site_score_map[i];
        db->site_score_map[i] = new std::map<int, ScoredSite>;

        if(db->event_alignment_result){ //eventalign related
            delete db->event_alignment_result[i];
            db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        }
        if(db->event_alignment_result_str){ //eventalign related
            free(db->event_alignment_result_str[i]);
            db->event_alignment_result_str[i]=NULL;
        }
    }
}

/* completely free a data batch */
void free_db(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        bam_destroy1(db->bam_rec[i]);
    }
    free(db->bam_rec);
    free(db->fasta_cache);
    free(db->read);
    free(db->read_len);
    free(db->read_idx);
    free(db->et);
    free(db->sig);
    free(db->scalings);
    free(db->event_align_pairs);
    free(db->n_event_align_pairs);
    free(db->event_alignment);
    free(db->n_event_alignment);
    free(db->events_per_base);
    free(db->base_to_event_map);
    free(db->read_stat_flag);
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        delete db->site_score_map[i];
    }
    free(db->site_score_map);
    //eventalign related
    if(db->eventalign_summary){
        free(db->eventalign_summary);
    }
    if(db->event_alignment_result){
        for (i = 0; i < db->capacity_bam_rec; ++i) {
            delete db->event_alignment_result[i];
        }
        free(db->event_alignment_result);
    }
    if(db->event_alignment_result_str){
        free(db->event_alignment_result_str);
    }

    free(db);
}

/* initialise user specified options */
void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->min_mapq = 20;
    opt->batch_size = 512;
    opt->batch_size_bases = 2*1000*1000;
    opt->num_thread = 8;
    opt->num_iop = 1;       //if changed, the SLOW5 mode must be handled by default in the arg parsing
    opt->region_str = NULL; //whole genome processing if null

    opt->min_num_events_to_rescale = 200;

#ifndef HAVE_CUDA
    opt->flag |= F5C_DISABLE_CUDA;
    opt->batch_size_bases = 5*1000*1000;
#endif

    opt->flag |= F5C_SKIP_UNREADABLE;
    opt->debug_break=-1;
    opt->ultra_thresh=100000;

    opt->meth_out_version=2;

    opt->cuda_block_size=64;
    opt->cuda_dev_id=0;
    opt->cuda_mem_frac=1.0f; //later set by cuda_init()

    //effective only if  CPU_GPU_PROC  is set
    opt->cuda_max_readlen=3.0f;
    opt->cuda_avg_events_per_kmer=2.0f; //only if CUDA_DYNAMIC_MALLOC is unset
    opt->cuda_max_avg_events_per_kmer=5.0f;
}

