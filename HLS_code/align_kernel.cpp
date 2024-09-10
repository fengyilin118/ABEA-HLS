#include "align_multi.h"
//keep band_score local, stream in events and model_levels

void band_aligner_stream_core_v2(int band_lower_left_event,int min_offset, int max_offset, int n_event,
event_info events[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], 
lp_info lps, FIXED_POINT band_score_pre2[BAND_WIDTH],FIXED_POINT band_score_pre[BAND_WIDTH],
FIXED_POINT band_score_current[BAND_WIDTH],
FIXED_TRACE trace_table[BAND_WIDTH],int trim_offset,FIXED_POINT trim_score
);
void band_aligner_stream_core_v3(int offset_left_init, int offset_diag_init,int band_lower_left_event,int min_offset, int max_offset, int n_event,
event_info events[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], 
lp_info lps, FIXED_POINT band_score_pre2[BAND_WIDTH],FIXED_POINT band_score_pre[BAND_WIDTH],
FIXED_POINT band_score_current[BAND_WIDTH],
FIXED_TRACE trace_table[BAND_WIDTH],int trim_offset
);
void align_kernel_stream_multi(
hls::stream<int>kmer_length_single[MULTI_READ],
hls::stream<int>event_length_single[MULTI_READ],
hls::stream<int>&kmer_length_post_single,
hls::stream<int>&event_length_post_single,
hls::stream<int>&band_length_single,
hls::stream<FIXED_POINT_MEAN> events_stream_single[MULTI_READ],
hls::stream<model_info>model_levels_stream[MULTI_READ],

hls::stream<FIXED_POINT>&max_cell_stream,
hls::stream<int>&max_cell_event_idx_stream,
hls::stream<TRACE_VEC>& trace_vec_kernel_stream,


hls::stream<int> &band_lower_left_event_stream_single,
lp_info lps_in[MULTI_READ],
hls::stream <int> &event_n_stream_max

)

{


    #pragma HLS pipeline off

int local_kmer_n[MULTI_READ];
int local_event_n[MULTI_READ];

int band_num_local[MULTI_READ];
model_info model_levels_cur[MULTI_READ][BAND_WIDTH];
event_info events[MULTI_READ][BAND_WIDTH+1];
int band_lower_left_event[MULTI_READ][2];
int band_lower_left_kmer[MULTI_READ][2];

FIXED_POINT band_score_work_pre[MULTI_READ][BAND_WIDTH];
FIXED_POINT band_score_work_pre2[MULTI_READ][BAND_WIDTH];
int band_lower_left_event_work_pre[MULTI_READ];
int band_lower_left_event_work_pre2[MULTI_READ];



#pragma HLS aggregate variable=model_levels_cur 


#pragma HLS array_partition variable=band_score_work_pre type=complete dim=2
#pragma HLS array_partition variable=band_score_work_pre2 type=complete dim=2

#pragma HLS array_partition variable=model_levels_cur type=complete dim=2
#pragma HLS array_partition variable=events type=complete dim=2



#pragma HLS array_partition variable=band_lower_left_event type=complete dim=2
#pragma HLS array_partition variable=band_lower_left_kmer type=complete dim=2


//#pragma HLS bind_storage variable=events type=RAM_T2P impl=auto
//#pragma HLS bind_storage variable=model_levels_cur type=RAM_T2P impl=auto


//#pragma HLS bind_storage variable=band_score_work_pre type=RAM_T2P impl=auto
//#pragma HLS bind_storage variable=band_score_work_pre2 type=RAM_T2P impl=auto

//send out two pre initialized band
//align_pre


lp_info lps[MULTI_READ];
int max_band=0;




FIXED_POINT lp_trim = -4.605170185988091;

int exceed_boundry_event=BAND_WIDTH / 2 - 1;//49
int exceed_boundry_kmer=-1 - BAND_WIDTH / 2;
int move_down_event=exceed_boundry_event+1;//50
int event_count[MULTI_READ];
int kmer_count[MULTI_READ];

int start_cell_offset = -1-exceed_boundry_kmer; 
int first_trim_sffset=move_down_event;//50




LOOP_INIT: 
for (int t = 0; t < 2; t++)
{
    for (int i = 0; i < MULTI_READ; i++)
    {
        #pragma HLS pipeline
        int temp_band_num;
        int temp_n_kmer;
        int temp_n_event;

        if(t==0)
        {


        kmer_length_single[i].read(temp_n_kmer);
        event_length_single[i].read(temp_n_event);

        kmer_length_post_single.write(temp_n_kmer);
        event_length_post_single.write(temp_n_event);

        local_kmer_n[i]=temp_n_kmer;
        local_event_n[i]=temp_n_event;

       
     
        temp_band_num=temp_n_kmer +temp_n_event+2;
        
        band_num_local[i]=temp_band_num;

        //printf("temp kmer %d \n",temp_n_kmer);
        //printf("temp n_event %d \n",temp_n_event);
        //printf("temp_band_num %d \n",temp_band_num);
        //printf("\n");

        if (temp_band_num>max_band)
        {
        max_band=temp_band_num;
        }
        
        //lps init
     
        lps[i]=lps_in[i];

        //band_num_out
        band_length_single.write(temp_band_num);
        event_n_stream_max.write(temp_n_event);

        }


        //left init
        band_lower_left_event[i][0]=exceed_boundry_event; // excced boundry, negative
        band_lower_left_kmer[i][0]=exceed_boundry_kmer;
        band_lower_left_event[i][1]=move_down_event;
        band_lower_left_kmer[i][1]=exceed_boundry_kmer;
 

        band_lower_left_event_stream_single.write(band_lower_left_event[i][t]);

        //counter init
        event_count[i]=0;
        kmer_count[i]=0;

    }
}

  
//printf("first trim %d \n",first_trim_sffset);
//printf("lp trim %f \n",lp_trim.to_double());

TRACE_VEC init_trace_vec_pre2;
TRACE_VEC init_trace_vec_pre;




LOOP_INIT_BAND_pre2:
for (int i = 0; i < MULTI_READ; i++)
{
    #pragma HLS pipeline
    for (int j = 0; j < BAND_WIDTH; j++)
    {
        FIXED_POINT pre2_score;
        if (j==start_cell_offset)
        {
            pre2_score=(FIXED_POINT)0;
        }
        else
        {         
            pre2_score=NEG_INF;
        }
         band_score_work_pre2[i][j]=pre2_score;        
    }   
}


LOOP_INIT_BAND_PRE:
for (int i = 0; i < MULTI_READ; i++)
{
    #pragma HLS pipeline
    for (int j = 0; j < BAND_WIDTH; j++)
    {
        FIXED_POINT pre_score;
        if (j==first_trim_sffset)
        {
            pre_score=lp_trim;
        }
        else
        {         
            pre_score=NEG_INF;
        }
         band_score_work_pre[i][j]=pre_score;     
    } 
}



LOOP_INIT_TRACE:
for (int j = 0; j < BAND_WIDTH; j++)
{
    #pragma HLS UNROLL
        FIXED_TRACE temp_trace;
        if (j==first_trim_sffset)
        {
            temp_trace=FROM_U;
        }
        else
        {         
            temp_trace=0;
        }
         init_trace_vec_pre2[j]=0;
         init_trace_vec_pre[j]=temp_trace;     
} 



for (int i = 0; i < MULTI_READ; i++)
{
    #pragma HLS pipeline
    trace_vec_kernel_stream.write(init_trace_vec_pre2);
    max_cell_stream.write(NEG_INF);
    max_cell_event_idx_stream.write(0);
}

for (int i = 0; i < MULTI_READ; i++)
{
        #pragma HLS pipeline

    trace_vec_kernel_stream.write(init_trace_vec_pre);
    max_cell_stream.write(NEG_INF);
    max_cell_event_idx_stream.write(0);
}



int SFR_PRE_MODEL_COUNT=BAND_WIDTH/2-1;
int SFR_PRE_EVENT_COUNT=SFR_PRE_MODEL_COUNT+2;

int band_idx_init=2-SFR_PRE_EVENT_COUNT;


LOOP_BAND_ITERATIVE: for (int band_idx = band_idx_init; band_idx<max_band; band_idx++)//question bound
{
    #pragma HLS loop_tripcount min=TEST_BAND_LEN max=TEST_BAND_LEN avg=TEST_BAND_LEN

 LOOP_MULTI_READ:   for (int i = 0; i < MULTI_READ; i++)
    {
	      #pragma HLS PIPELINE
#pragma HLS LOOP_FLATTEN
        #pragma HLS LATENCY MIN=10 MAX=10

        //model_info model_levels_local[BAND_WIDTH];
        //event_info events_local[BAND_WIDTH+1];

        //#pragma HLS aggregate variable=model_levels_local 

        
           // #pragma HLS array_partition variable=model_levels_local type=complete dim=0
          //  #pragma HLS array_partition variable=events_local type=complete dim=0
            
        bool right;
        right=false;
        
        int band_lower_left_kmer_current;
        int band_lower_left_event_current;

        int band_lower_left_kmer_pre;
        int band_lower_left_event_pre;


        int band_lower_left_kmer_pre2;
        int band_lower_left_event_pre2;

   
            FIXED_POINT ll = band_score_work_pre[i][0];
            FIXED_POINT ur = band_score_work_pre[i][BAND_WIDTH - 1];

            //FIXED_POINT ll = ll_align[i];
            //FIXED_POINT ur = ur_align[i];



    if (band_idx>1)
    {
            bool ll_ob = ll ==NEG_INF;
            bool ur_ob = ur ==NEG_INF;
           // right=false;
            



             band_lower_left_kmer_pre=band_lower_left_kmer[i][1];
             band_lower_left_event_pre=band_lower_left_event[i][1];


             band_lower_left_kmer_pre2=band_lower_left_kmer[i][0];
             band_lower_left_event_pre2=band_lower_left_event[i][0];

                
        if (ll_ob && ur_ob)
            {
                right = band_idx % 2 == 1;
            }
        else
            {
                right = ll < ur;
            }
    }  

    if ((right)||(band_idx<0))
        {

        model_info model_cur;
        if(kmer_count[i]<local_kmer_n[i])
        {

            model_levels_stream[i].read(model_cur);
    
        }


       // printf("kmer_count %d \n",  kmer_count[i]);

         LOOP_BAND_KMER_SHIFT:    for (int j = 0; j < BAND_WIDTH-1; j++)
        {
            #pragma HLS UNROLL
            model_levels_cur[i][j]= model_levels_cur[i][j+1];
        }
        model_levels_cur[i][BAND_WIDTH-1]=model_cur;   

        band_lower_left_kmer_current = band_lower_left_kmer_pre+1;
        band_lower_left_event_current = band_lower_left_event_pre;
        kmer_count[i]++;

    }



    if((!right)||(band_idx<2)) 
    {
       // printf("move down \n");

        FIXED_POINT_MEAN event_mean;
        event_info new_event;

        if (event_count[i]<(local_event_n[i]))
        {
            events_stream_single[i].read(event_mean);
            //printf("event_count %d %d \n",  event_count[i],local_event_n[i]);

        }
        new_event.mean=event_mean;

         LOOP_BAND_EVENT_SHIFT:     for (int j = 0; j < BAND_WIDTH; j++)
        {
            #pragma HLS UNROLL
            events[i][j]=events[i][j+1];            
        }
        
        events[i][BAND_WIDTH]=new_event;
        
        event_count[i]++;



        band_lower_left_event_current=band_lower_left_event_pre+1;
        band_lower_left_kmer_current=band_lower_left_kmer_pre;
    }

    
        FIXED_POINT band_score_out[BAND_WIDTH];
        FIXED_TRACE trace_table_out[BAND_WIDTH];


        #pragma HLS array_partition variable=band_score_out type=complete dim=0
        #pragma HLS array_partition variable=trace_table_out type=complete 

            
            if (band_idx>1&&band_idx<band_num_local[i])
            {

                int trim_offset=(-1)-band_lower_left_kmer_current;

                int event_idx_init=band_lower_left_event_current;
                int kmer_idx_init=band_lower_left_kmer_current;
                int offset_left_init = (kmer_idx_init-1)-band_lower_left_kmer_pre;
                int offset_diag_init = (kmer_idx_init-1)-band_lower_left_kmer_pre2;

                int kmer_min_offset=0-band_lower_left_kmer_current;
                int kmer_max_offset=local_kmer_n[i]-band_lower_left_kmer_current;
                int event_min_offset=band_lower_left_event_current-(local_event_n[i]-1);
                int event_max_offset=band_lower_left_event_current+1;

                 int min_offset;
                if (kmer_min_offset > event_min_offset)
                {
                    min_offset = kmer_min_offset;

                }
                else
                {
                    min_offset = event_min_offset;
                }

                if (min_offset > 0)
                {
                    min_offset = min_offset;//original no +1
                } 
                else
                {
                    min_offset = 0;
                }

                int max_offset;
                if (kmer_max_offset < event_max_offset)
                {
                    max_offset = kmer_max_offset;
                }
                    
                else
                {
                    max_offset = event_max_offset;
                }

                if (max_offset > BAND_WIDTH)
                {
                            max_offset = BAND_WIDTH;
                }
                else
                {
                            max_offset = max_offset;
                }

                FIXED_POINT band_aligner_pre[BAND_WIDTH];
                FIXED_POINT band_aligner_pre_up[BAND_WIDTH];
                FIXED_POINT band_aligner_pre2[BAND_WIDTH];

                #pragma HLS array_partition variable=band_aligner_pre type=complete
                #pragma HLS array_partition variable=band_aligner_pre_up type=complete

                #pragma HLS array_partition variable=band_aligner_pre2 type=complete

                 switch (offset_diag_init)
                {
                case 1:
                    for (int q = 1; q < BAND_WIDTH; q++)
                    {
                    band_aligner_pre2[q-1]=band_score_work_pre2[i][q];
                    }
                    band_aligner_pre2[BAND_WIDTH-1]=NEG_INF;

                    break;
                case -1:
                for (int q = 0; q < BAND_WIDTH-1; q++)
                {
                band_aligner_pre2[q+1]=band_score_work_pre2[i][q];
                }

                    band_aligner_pre2[0]=NEG_INF;
                    break;
                default:
                    for (int q = 0; q < BAND_WIDTH; q++)
                    {
                        band_aligner_pre2[q]=band_score_work_pre2[i][q];
                    }
                    
                    break;
                }


                switch (offset_left_init)
                {
                case 1:
                    for (int q = 1; q < BAND_WIDTH; q++)
                    {
                    band_aligner_pre[q-1]=band_score_work_pre[i][q];
                    
                    }
                    band_aligner_pre[BAND_WIDTH-1]=NEG_INF;


                    break;
                case -1:
                for (int q = 0; q < BAND_WIDTH-1; q++)
                {
                band_aligner_pre[q+1]=band_score_work_pre[i][q];
                }

                band_aligner_pre[0]=NEG_INF;

                    break;
                default:
                    for (int q = 0; q < BAND_WIDTH; q++)
                    {
                        band_aligner_pre[q]=band_score_work_pre[i][q];

                    }
                    
                    break;
                }
            FIXED_POINT trim_score=lp_trim_const*((FIXED_POINT)(band_lower_left_event_current-trim_offset+1));


            band_aligner_stream_core_v2( band_lower_left_event_current, min_offset,max_offset,local_event_n[i],events[i],model_levels_cur[i],lps[i],band_aligner_pre2,band_aligner_pre,band_score_out,trace_table_out,trim_offset,trim_score);



                int temp_lower_left_event_out;
                int temp_lower_left_kmer_out;
                
                band_lower_left_event[i][0]=band_lower_left_event[i][1];
                temp_lower_left_event_out=band_lower_left_event_current;
                band_lower_left_event_stream_single.write(temp_lower_left_event_out);
                band_lower_left_event[i][1]=temp_lower_left_event_out;



                band_lower_left_kmer[i][0]=band_lower_left_kmer[i][1];
                temp_lower_left_kmer_out=band_lower_left_kmer_current;
                band_lower_left_kmer[i][1]=temp_lower_left_kmer_out;


                TRACE_VEC trace_vec_stream_out;

                FIXED_POINT band_score_to_compare=NEG_INF;
                int temp_max_event_idx=0;                
          
                
            LOOP_RECEIVE:for (int j = 0; j < BAND_WIDTH; j++)
            {
                trace_vec_stream_out[j]=trace_table_out[j];                
                band_score_work_pre2[i][j]=band_score_work_pre[i][j];
                band_score_work_pre[i][j]=band_score_out[j];

            }

            int max_cell_offset=local_kmer_n[i]-1-band_lower_left_kmer_current;


            if (max_cell_offset>=0&&max_cell_offset<BAND_WIDTH)
            {
                band_score_to_compare=band_score_out[max_cell_offset];
                temp_max_event_idx=band_lower_left_event_current-max_cell_offset;

            }
            

            
            max_cell_stream.write(band_score_to_compare);
            max_cell_event_idx_stream.write(temp_max_event_idx);
            

            trace_vec_kernel_stream.write(trace_vec_stream_out);
                
            }
        
    }



}




for (int i = 0; i < MULTI_READ; i++)
{
    while (kmer_count[i]<local_kmer_n[i])
    {
            #pragma HLS pipeline

                        #pragma HLS loop_tripcount min=1 max=1 avg=1

        printf("left_over kernel kmer\n");
        model_levels_stream[i].read();
        kmer_count[i]++;
    }

    while(event_count[i]<local_event_n[i])
    {
        #pragma HLS pipeline
                        #pragma HLS loop_tripcount min=1 max=1 avg=1

                printf("left_over kernel event \n");

        events_stream_single[i].read();
        event_count[i]++;
    }
    
}


//printf("max band %d \n",max_band);
//printf("kmer_local %d \n",local_kmer_n[0]);
//printf("kmer_count %d \n", kmer_count[0] );
//printf("event_count %d \n",event_count[0]);



} 



void band_aligner_stream_core_v3(int offset_left_init, int offset_diag_init,int band_lower_left_event,int min_offset, int max_offset, int n_event,
event_info events[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], 
lp_info lps, FIXED_POINT band_score_pre2[BAND_WIDTH],FIXED_POINT band_score_pre[BAND_WIDTH],
FIXED_POINT band_score_current[BAND_WIDTH],
FIXED_TRACE trace_table[BAND_WIDTH],int trim_offset
)
{
#pragma HLS inline

    LOOP_BAND_CORE:for (int offset = 0; offset < BAND_WIDTH; offset++)
    {
        #pragma HLS UNROLL

            FIXED_POINT up, left, diag;
            up=NEG_INF;
            left=NEG_INF;
            diag=NEG_INF;
            #ifdef DEBUG_MOD

            //printf("band core \n");
       
            #endif
            
            if(offset >= min_offset && offset < max_offset)
            {

                #ifdef DEBUG_MOD
            printf("offset_left %d \n",offset);
            printf("offset_up %d \n",offset+1);

            printf("offset_diag %d \n",offset);
                #endif            


                    switch (offset_left_init)
                    {
                    case 1:
                        if(offset!=(BAND_WIDTH-1))
                        {
                            left=band_score_pre[offset+1];
                        }

                         if(offset+2<BAND_WIDTH)
                        {
                            up=band_score_pre[offset+2];
                        }
                 
                        break;

                    case -1:
                        if(offset!=0)
                        {
                            left=band_score_pre[offset-1];
                        }
                
                   
                        up=band_score_pre[offset];

                        break;


                    default:
                        left=band_score_pre[offset];
                        
                         if (offset!=(BAND_WIDTH-1))
                        {
                            up = band_score_pre[offset+1];// band_score_pre[offset+1].score;//replace offset_up with direct access
                        }
                                   
                        break;
                    }



                    switch (offset_diag_init)
                    {
                    case 1:
                        if (offset!=(BAND_WIDTH-1))
                        {
                           diag=band_score_pre2[offset+1];
                        }
                        break;

                    case -1:
                        if(offset!=0)
                        {
                            diag=band_score_pre2[offset-1];
                        }
                        break;
                    
                    default:
                        diag=band_score_pre2[offset];
                        break;
                    }
                                    
                                   
    #ifdef DEBUG_MOD
                
                    printf("left =%f \n",left.to_double());
                    printf("up =%f \n",up.to_double());
                    printf("offset %d diag= %f \n",offset,diag.to_double());
    #endif
                FIXED_POINT_MEAN unscaledLevel = events[BAND_WIDTH-offset].mean;//replace event_idx with direct access
                FIXED_POINT_MEAN scaledLevel = unscaledLevel;


                model_info curr_model=model_levels_cur[offset];
                #pragma HLS disaggregate variable=curr_model

                //FIXED_POINT gp_mean = scailings.scalings_scale * model_levels_cur[offset].mean + scailings.scalings_shift;
                FIXED_POINT gp_mean=curr_model.gp_mean;
                
                //

                //FIXED_POINT gp_mean=68.724602;
                
                //int kmer_idx=model_levels_cur[offset].id;
                
                //int read_idx=model_levels_cur[offset].read_id;
                //int read_idx_event=events[BAND_WIDTH-offset].read_id;

                FIXED_POINT_STDV gp_stdv = curr_model.stdv;
                FIXED_POINT_STDV gp_log_stdv = curr_model.log_stdv;


                FIXED_POINT_STDV log_inv_sqrt_2pi = -0.918938;
                FIXED_POINT a;


                if(gp_stdv==(FIXED_POINT_STDV)0)
                {
                    a=0;
                }
                else
                {
                    FIXED_POINT_STDV inverse_gp_stdv=curr_model.gp_stdv_inverse;
                    FIXED_POINT a_sub;
                    a_sub=(scaledLevel-gp_mean);
                    a=a_sub*inverse_gp_stdv;
                    //a = (scaledLevel - gp_mean) / gp_stdv;
                }
                //FIXED_POINT a;

                #ifdef DEBUG_MOD
                printf(" offset %d model_mean=%f \n", offset,model_levels_cur[offset].mean.to_double());
                printf(" offset%d scaliedLevel= %f \n",offset, scaledLevel.to_double());
                printf(" offset%d gp_mean= %f \n",offset,gp_mean.to_double());
                printf(" offset%d gp_stdv= %f \n",offset,gp_stdv.to_double());
                printf(" offset%d a= %f \n",offset,a.to_double());

                #endif

                



                ap_fixed<3,1> half = -0.5;
                FIXED_POINT_STDV lp_emission;
                FIXED_POINT_STDV lp_emission_sub;
                FIXED_POINT lp_emission_add;
                lp_emission =log_inv_sqrt_2pi;
                //printf("offset %d lp_emission= %f \n",lp_emission.to_double());
                lp_emission_sub =lp_emission-gp_log_stdv;
               // printf("offset %d lp_emission= %f \n",lp_emission.to_double());



                lp_emission_add =lp_emission_sub+half*a*a;


                //FIXED_POINT lp_emission = log_inv_sqrt_2pi - gp_log_stdv + (half * a * a);
                  //              printf("offset %d lp_emission= %f \n",lp_emission.to_double());


                FIXED_POINT score_d = NEG_INF;
                FIXED_POINT score_u = NEG_INF;
                FIXED_POINT score_l = NEG_INF;
              
                    score_d = diag + lps.lp_step + lp_emission_add;

             
                    score_u = up + lps.lp_stay + lp_emission_add;
                

            
                    score_l = left + lps.lp_skip;
                


    
            #ifdef DEBUG_MOD
                printf("offset %d lp_emission= %f \n",offset,lp_emission.to_double());
                printf("offset %d lp_step=%f \n",offset,lps.lp_step.to_double());
                printf("offset %d lp_skip=%f \n",offset,lps.lp_skip.to_double());

                printf("offset %d score d=%f score u=%f score l=%f \n",offset,score_d.to_double(),score_u.to_double(),score_l.to_double());
            #endif

                FIXED_POINT max_score = score_d;
                FIXED_TRACE from = FROM_D;



                align_kernel_score_calculate(score_d, score_u, score_l, &max_score, &from);
             
                
                band_score_current[offset] = max_score;

                trace_table[offset] = from;

            }
            else if(offset==trim_offset)
            {
                        int event_idx;
                        event_idx = band_lower_left_event-offset;

                        if (event_idx >= 0 && event_idx < n_event)
                        {


                            band_score_current[offset] = lp_trim_const * ((FIXED_POINT)event_idx + 1);
                            
                            trace_table[offset] = FROM_U;
                        }
                        else
                        {
                            band_score_current[offset] = NEG_INF;//
                            trace_table[offset]=0;
                        }
                }
             else
                {
                    band_score_current[offset] = NEG_INF;//
                    trace_table[offset]=0;
                }
       

            //printf("band score %f \n",max_score.to_double());    
    
    }

    



}




void band_aligner_stream_core_v2(int band_lower_left_event,int min_offset, int max_offset, int n_event,
event_info events[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], 
lp_info lps, FIXED_POINT band_score_pre2[BAND_WIDTH],FIXED_POINT band_score_pre[BAND_WIDTH],
FIXED_POINT band_score_current[BAND_WIDTH],
FIXED_TRACE trace_table[BAND_WIDTH],int trim_offset,FIXED_POINT trim_score
)
{
#pragma HLS inline

    LOOP_BAND_CORE:for (int offset = 0; offset < BAND_WIDTH; offset++)
    {
        #pragma HLS UNROLL

            FIXED_POINT up, left, diag;
            FIXED_POINT out_score;
            FIXED_TRACE out_trace;
            out_score=NEG_INF;
            out_trace=0;
            #ifdef DEBUG_MOD

            //printf("band core \n");
       
            #endif
            
            if(offset >= min_offset && offset < max_offset)
            {

                #ifdef DEBUG_MOD
            printf("offset_left %d \n",offset);
            printf("offset_up %d \n",offset+1);

            printf("offset_diag %d \n",offset);
                #endif            

                
                    left = band_score_pre[offset];//replace offset_levet with direct acess               
                    
                    
                    if (offset+1<BAND_WIDTH)
                    {
                        up = band_score_pre[offset+1];// band_score_pre[offset+1].score;//replace offset_up with direct access

                    }
                    else
                    {
                        up=NEG_INF;
                    }
                    

                    diag = band_score_pre2[offset];//replace offset_diag with direct access
                
    #ifdef DEBUG_MOD
                
                    printf("left =%f \n",left.to_double());
                    printf("up =%f \n",up.to_double());
                    printf("offset %d diag= %f \n",offset,diag.to_double());
    #endif
                FIXED_POINT_MEAN unscaledLevel = events[BAND_WIDTH-offset].mean;//replace event_idx with direct access
                FIXED_POINT_MEAN scaledLevel = unscaledLevel;


                model_info curr_model=model_levels_cur[offset];
                #pragma HLS disaggregate variable=curr_model

                //FIXED_POINT gp_mean = scailings.scalings_scale * model_levels_cur[offset].mean + scailings.scalings_shift;
                FIXED_POINT gp_mean=curr_model.gp_mean;
                
                //

                //FIXED_POINT gp_mean=68.724602;
                
                //int kmer_idx=model_levels_cur[offset].id;
                
                //int read_idx=model_levels_cur[offset].read_id;
                //int read_idx_event=events[BAND_WIDTH-offset].read_id;

                FIXED_POINT_STDV gp_stdv = curr_model.stdv;
                FIXED_POINT_STDV gp_log_stdv = curr_model.log_stdv;


                FIXED_POINT_STDV log_inv_sqrt_2pi = -0.918938;
                FIXED_POINT a;


                if(gp_stdv==(FIXED_POINT_STDV)0)
                {
                    a=0;
                }
                else
                {
                    FIXED_POINT_STDV inverse_gp_stdv=curr_model.gp_stdv_inverse;
                    FIXED_POINT a_sub;
                    a_sub=(scaledLevel-gp_mean);
                    a=a_sub*inverse_gp_stdv;
                    //a = (scaledLevel - gp_mean) / gp_stdv;
                }
                //FIXED_POINT a;

                #ifdef DEBUG_MOD
                printf(" offset %d model_mean=%f \n", offset,model_levels_cur[offset].mean.to_double());
                printf(" offset%d scaliedLevel= %f \n",offset, scaledLevel.to_double());
                printf(" offset%d gp_mean= %f \n",offset,gp_mean.to_double());
                printf(" offset%d gp_stdv= %f \n",offset,gp_stdv.to_double());
                printf(" offset%d a= %f \n",offset,a.to_double());

                #endif

                



                ap_fixed<3,1> half = -0.5;
                FIXED_POINT_STDV lp_emission;
                FIXED_POINT_STDV lp_emission_sub;
                FIXED_POINT lp_emission_add;
                lp_emission =log_inv_sqrt_2pi;
                //printf("offset %d lp_emission= %f \n",lp_emission.to_double());
                lp_emission_sub =lp_emission-gp_log_stdv;
               // printf("offset %d lp_emission= %f \n",lp_emission.to_double());



                lp_emission_add =lp_emission_sub+half*a*a;


                //FIXED_POINT lp_emission = log_inv_sqrt_2pi - gp_log_stdv + (half * a * a);
                  //              printf("offset %d lp_emission= %f \n",lp_emission.to_double());


                FIXED_POINT score_d = NEG_INF;
                FIXED_POINT score_u = NEG_INF;
                FIXED_POINT score_l = NEG_INF;
              
                    score_d = diag + lps.lp_step + lp_emission_add;

             
                    score_u = up + lps.lp_stay + lp_emission_add;
                

            
                    score_l = left + lps.lp_skip;
                


    
            #ifdef DEBUG_MOD
                printf("offset %d lp_emission= %f \n",offset,lp_emission.to_double());
                printf("offset %d lp_step=%f \n",offset,lps.lp_step.to_double());
                printf("offset %d lp_skip=%f \n",offset,lps.lp_skip.to_double());

                printf("offset %d score d=%f score u=%f score l=%f \n",offset,score_d.to_double(),score_u.to_double(),score_l.to_double());
            #endif

                FIXED_POINT max_score = score_d;
                FIXED_TRACE from = FROM_D;



                align_kernel_score_calculate(score_d, score_u, score_l, &max_score, &from);
             
                
                out_score = max_score;

                out_trace = from;

            }
            else if(offset==trim_offset)
            {
                    int event_idx;
                    event_idx = band_lower_left_event-offset;

                if (event_idx >= 0 && event_idx < n_event)
                {

                    out_score = trim_score;
                    out_trace = FROM_U;
                }
                 
            }
        
       
            band_score_current[offset]=out_score;
            trace_table[offset]=out_trace;

            //printf("band score %f \n",max_score.to_double());    
    
    }

    



}

