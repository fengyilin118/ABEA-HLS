#include "align_multi.h"




void align_dataflow_stream(int* sequence_len,int* event_t_length,
int event_offset, int sequence_offset,int reset,

EVENT_VEC *event_mean0,

READ_VEC*  read0,
scailing_info scailings[MULTI_TOTAL], 
lp_info lps[MULTI_TOTAL], 
FIXED_POINT_MEAN model_levels_mean[MULTI_TOTAL][MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_stdv[MULTI_TOTAL][MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_log_stdv[MULTI_TOTAL][MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_gp_stdv_inverse[MULTI_TOTAL][MAX_MODEL_NUM],

AP_TRACE_VEC* trace_out_0,
LOWER_LEFT_VEC* left_out_0,

POS_POST_VEC* aligned_ref_read_pos_0,

int* n_align
 )
{
#pragma HLS DATAFLOW 
#pragma HLS stable variable=event_mean0
#pragma HLS stable variable=read0
#pragma HLS stable variable=scailings
#pragma HLS stable variable=model_levels_mean
#pragma HLS stable variable=model_levels_stdv
#pragma HLS stable variable=model_levels_log_stdv
#pragma HLS stable variable=model_levels_gp_stdv_inverse

#pragma HLS stable variable=lps

#pragma HLS stable variable=sequence_len
#pragma HLS stable variable=event_t_length


         hls::stream <int> sequence_offset_stream_group;
         hls::stream <int> event_offset_stream_group;
        
         hls::stream <int> sequence_offset_stream_vec[LOAD_DATAFLOW];
         hls::stream <int> event_offset_stream_vec[LOAD_DATAFLOW];

         hls::stream <int> sequence_offset_stream_gen[MULTI_TOTAL];
         hls::stream <int> event_offset_stream_gen[MULTI_TOTAL];
      

         
        
         hls::stream <EVENT_VEC,8> event_vec_stream[LOAD_DATAFLOW];
         hls::stream <READ_VEC,8>read_vec_stream[LOAD_DATAFLOW];

         hls::stream <float,EVENT_FIFO_DEPTH>event_float_stream[MULTI_TOTAL];
         hls::stream <char,READ_FIFO_DEPTH>read_char_stream[MULTI_TOTAL];
        
        //#pragma HLS bind_storage type=FIFO impl=BRAM variable=read_char_stream
        #pragma HLS bind_storage type=FIFO impl=URAM variable=event_float_stream


         hls::stream <FIXED_POINT_MEAN,5> event_stream[MULTI_TOTAL];
         hls::stream <model_info,5>model_levels_stream[MULTI_TOTAL];

         hls::stream<int>event_t_length_single[MULTI_TOTAL];
         hls::stream<int>event_t_length_2kernel_single[MULTI_TOTAL];

         hls::stream<int>sequence_len_single[MULTI_TOTAL];
         hls::stream<int>n_kmers_2kernel_single[MULTI_TOTAL];



        hls::stream <FIXED_POINT,MULTI_READ> max_cell_stream;
        hls::stream <int,MULTI_READ>max_cell_event_idx_stream;
        hls::stream <TRACE_VEC,20> trace_vec_kernel_stream;
        hls::stream<int,20> band_lower_left_event_out_stream;


        hls::stream <int,MULTI_READ> band_length_single;
        hls::stream <int,MULTI_READ> event_n_stream_max;

        hls::stream <int,MULTI_TOTAL*2>band_length_max_out;


        hls::stream <int,MULTI_READ> max_event_idx_stream;

        hls::stream<bool>event_last[MULTI_TOTAL];
        hls::stream<bool>model_last[MULTI_TOTAL];


       //post_streams
        hls::stream<int,MULTI_READ>kmer_length_post_single;
        hls::stream<int,MULTI_READ>event_length_post_single;


        hls::stream <FIXED_POINT_MEAN,5> event_post_stream[MULTI_TOTAL];
        hls::stream <model_info,5>model_levels_post_stream[MULTI_TOTAL];

        hls::stream <int> max_band_post_out_stream;
        hls::stream <TRACE_VEC,5>trace_vec_post_stream[2];
        hls::stream <int,5>band_lower_left_event_post_stream[2];
        hls::stream <int,MULTI_TOTAL>max_event_idx_post_stream;

        hls::stream <int,3>ref_pos_post_stream[MULTI_TOTAL];
        hls::stream <int,3>read_pos_post_stream[MULTI_TOTAL];

        hls::stream<int,MULTI_TOTAL> n_pairs_post_stream;


        align_entry(reset,sequence_offset,event_offset,sequence_offset_stream_group,event_offset_stream_group);
        align_load_event_group(event_t_length,event_offset_stream_group,event_offset_stream_vec, event_t_length_single,event_mean0,event_vec_stream);
        align_load_read_group(sequence_len,sequence_offset_stream_group,sequence_offset_stream_vec,sequence_len_single,read0,read_vec_stream);



        align_stream_event_vec_load_v2<LOAD_DATAFLOW>(event_offset_stream_vec,event_offset_stream_gen,event_vec_stream,event_float_stream);
        align_stream_read_vec_load_v2<LOAD_DATAFLOW>(sequence_offset_stream_vec,sequence_offset_stream_gen,read_vec_stream,read_char_stream);
        
        
        align_stream_event_mean_gen<MULTI_TOTAL>(event_offset_stream_gen,event_t_length_single,event_t_length_2kernel_single,event_float_stream,event_stream,event_post_stream,event_last);
        align_stream_model_gen<MULTI_TOTAL>(sequence_offset_stream_gen,sequence_len_single,n_kmers_2kernel_single,read_char_stream,model_levels_mean,model_levels_stdv,model_levels_log_stdv,model_levels_gp_stdv_inverse,model_levels_stream,model_levels_post_stream,scailings,model_last);


        align_kernel_stream_multi(n_kmers_2kernel_single,event_t_length_2kernel_single,kmer_length_post_single,event_length_post_single,band_length_single,event_stream,model_levels_stream,max_cell_stream,max_cell_event_idx_stream,trace_vec_kernel_stream,band_lower_left_event_out_stream,lps,event_n_stream_max);

        align_max_v2(event_n_stream_max,band_length_single,band_length_max_out,max_cell_stream,max_cell_event_idx_stream,max_event_idx_stream);

        
        align_channel(band_length_max_out,trace_vec_kernel_stream,band_lower_left_event_out_stream,max_event_idx_stream,max_event_idx_post_stream,max_band_post_out_stream,trace_vec_post_stream,band_lower_left_event_post_stream,trace_out_0,left_out_0);

        align_post_multi(kmer_length_post_single,event_length_post_single,event_post_stream,model_levels_post_stream,trace_vec_post_stream,max_event_idx_post_stream,band_lower_left_event_post_stream,ref_pos_post_stream,read_pos_post_stream,n_pairs_post_stream);


        post_pos_out(event_last,model_last,max_band_post_out_stream,ref_pos_post_stream,read_pos_post_stream,n_pairs_post_stream,aligned_ref_read_pos_0, n_align);

  
}









