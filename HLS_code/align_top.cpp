#include "align_multi.h"


void align_top_stream (
int reset,

EVENT_VEC* event_mean0,

READ_VEC* read0,

int* n_align,
POS_POST_VEC* aligned_ref_read_pos_0,

float* model_data,
AP_TRACE_VEC* trace_out_0,
LOWER_LEFT_VEC* left_out_0,


int* length_int, float* scailings_double, float* lps_double,int max_event_t_lengh,int max_sequence_length,int batch_num

)
{


	#pragma HLS INTERFACE mode=m_axi port = model_data depth=MAX_MODEL_NUM*4 offset=slave bundle=gmem0
    
    #pragma HLS INTERFACE mode=m_axi port = length_int depth=MULTI_TOTAL*2*TEST_BATCH offset=slave bundle=gmem0
    #pragma HLS INTERFACE mode=m_axi port = scailings_double depth=MULTI_TOTAL*4*TEST_BATCH offset=slave bundle=gmem0
    #pragma HLS INTERFACE mode=m_axi port = lps_double depth=MULTI_TOTAL*3*TEST_BATCH offset=slave bundle=gmem0

    
    #pragma HLS INTERFACE mode=m_axi port = read0 depth=TEST_READ_OFFSET*MULTI_TOTAL/READ_VEC_WIDTH*TEST_BATCH offset=slave bundle=gmem0 max_widen_bitwidth=READ_VEC_WIDTH*8
    #pragma HLS INTERFACE mode=m_axi port = event_mean0 depth=TEST_EVENT_OFFSET*MULTI_TOTAL/EVENT_VEC_WIDTH*TEST_BATCH offset=slave bundle=gmem1 max_widen_bitwidth=EVENT_VEC_WIDTH*32
    #pragma HLS INTERFACE mode=m_axi port = trace_out_0 depth=TEST_BAND_LEN offset=slave bundle=gmem2 latency=64
    #pragma HLS INTERFACE mode=m_axi port = left_out_0 depth=TEST_BAND_LEN offset=slave bundle=gmem3 latency=64


    #pragma HLS INTERFACE mode=m_axi port = aligned_ref_read_pos_0 depth=TEST_BAND_LEN*TEST_BATCH offset=slave bundle=gmem4
    #pragma HLS INTERFACE mode=m_axi port = n_align depth=MULTI_TOTAL*TEST_BATCH offset=slave bundle=gmem4


	#pragma HLS INTERFACE s_axilite port=model_data
	#pragma HLS INTERFACE s_axilite port=length_int
	#pragma HLS INTERFACE s_axilite port=scailings_double
	#pragma HLS INTERFACE s_axilite port=lps_double


#pragma HLS INTERFACE s_axilite port=read0
#pragma HLS INTERFACE s_axilite port=event_mean0
#pragma HLS INTERFACE s_axilite port=trace_out_0
#pragma HLS INTERFACE s_axilite port=left_out_0


#pragma HLS INTERFACE s_axilite port=aligned_ref_read_pos_0


#pragma HLS INTERFACE s_axilite port=n_align

    #pragma HLS INTERFACE mode=s_axilite port=max_event_t_lengh
    #pragma HLS INTERFACE mode=s_axilite port=max_sequence_length

    #pragma HLS INTERFACE mode=s_axilite port=return
    #pragma HLS INTERFACE mode=s_axilite port=reset
   
    


    static FIXED_POINT_MEAN model_levels_mean[MULTI_TOTAL][MAX_MODEL_NUM];
    static FIXED_POINT_STDV model_levels_stdv[MULTI_TOTAL][MAX_MODEL_NUM];
    static FIXED_POINT_STDV model_levels_log_stdv[MULTI_TOTAL][MAX_MODEL_NUM];
    static FIXED_POINT_STDV model_levels_gp_stdv_inverse[MULTI_TOTAL][MAX_MODEL_NUM];


    #pragma HLS array_partition variable=model_levels_mean type=complete dim=1
    #pragma HLS array_partition variable=model_levels_stdv type=complete dim=1
    #pragma HLS array_partition variable=model_levels_log_stdv type=complete dim=1
    #pragma HLS array_partition variable=model_levels_gp_stdv_inverse type=complete dim=1
  




    FIXED_POINT scailings_local[512*MULTI_TOTAL*4];
    FIXED_POINT lps_local[512*3*MULTI_TOTAL];
    int length_local[512*2*MULTI_TOTAL];


    if(reset==1)

{
        float* model_data_offset=model_data;
        LOOP_LOAD_MODEL_MEAN: 
        
        for (int i = 0; i < MAX_MODEL_NUM; i++)
        {
            #pragma HLS pipeline
                FIXED_POINT_MEAN temp_model_data;
               temp_model_data=(FIXED_POINT_MEAN)model_data_offset[i];
                      for (int j = 0; j < MULTI_TOTAL; j++)
               {
                model_levels_mean[j][i]=temp_model_data;
               }

                                   
        }
        

        model_data_offset=model_data+MAX_MODEL_NUM;

        LOOP_LOAD_MODEL_STDV:for (int i = 0; i < MAX_MODEL_NUM; i++)
        {
            #pragma HLS pipeline
            FIXED_POINT_STDV temp_model_data;    
            temp_model_data=(FIXED_POINT_STDV)model_data_offset[i];
                   for (int j = 0; j < MULTI_TOTAL; j++)
               {
                model_levels_stdv[j][i]=temp_model_data;
               }

           
            
        }
        model_data_offset=model_data+MAX_MODEL_NUM*2;

        LOOP_LOAD_MODEL_LOG_STDV:
        for (int i = 0; i <  MAX_MODEL_NUM; i++)
        {
            #pragma HLS pipeline
          
            FIXED_POINT_STDV temp_model_data;    
            temp_model_data=(FIXED_POINT_STDV)model_data_offset[i];
                   for (int j = 0; j < MULTI_TOTAL; j++)
               {
                model_levels_log_stdv[j][i]=temp_model_data;
               }     

       
        }
        model_data_offset=model_data+MAX_MODEL_NUM*3;

         LOOP_LOAD_MODEL_GP_STDV_INVERSE:
        for (int i = 0; i <  MAX_MODEL_NUM; i++)
        {
            #pragma HLS pipeline
          
             FIXED_POINT_STDV temp_model_data;    
            temp_model_data=(FIXED_POINT_STDV)model_data_offset[i];
                for (int j = 0; j < MULTI_TOTAL; j++)
               {
                model_levels_gp_stdv_inverse[j][i]=temp_model_data;
               }
                 
        }
return;

}
else{
LOOP_LOAD_SCAILING:for (int i = 0; i < batch_num*4*MULTI_READ; i++)
    {
        #pragma HLS loop_tripcount min=4*MULTI_READ max=4*MULTI_READ avg=4*MULTI_READ
        #pragma HLS pipeline
        scailings_local[i]=(FIXED_POINT)scailings_double[i];
    }


LOOP_LOAD_LPS:for (int i = 0; i <batch_num*3*MULTI_READ; i++)
{
     #pragma HLS loop_tripcount min=3*MULTI_READ max=3*MULTI_READ avg=3*MULTI_READ
     #pragma HLS pipeline
     lps_local[i]=(FIXED_POINT)lps_double[i];

}

    
LOOP_LOAD_LENGTH:for (int i = 0; i < batch_num*2*MULTI_READ; i++)
{
    #pragma HLS loop_tripcount min=2*MULTI_READ max=2*MULTI_READ avg=2*MULTI_READ
     #pragma HLS pipeline
    length_local[i]=length_int[i];
}


    for (int batch = 0; batch < batch_num; batch++)
    {
        #pragma HLS loop_tripcount min=1 max=1 avg=1

            lp_info lps[MULTI_TOTAL];
    int sequence_len_dataflow[MULTI_TOTAL];
    int event_t_len_dataflow[MULTI_TOTAL];
        scailing_info scailings[MULTI_READ];
        #pragma HLS array_partition type=complete variable=scailings
    

        LOOP_LOAD_SCAILING_LOG_VAR:
        for (int j = 0; j < MULTI_READ; j++)
        {
            #pragma HLS pipeline
            scailings[j].scaling_log_var=scailings_local[batch*32+j];

        }

        LOOP_LOAD_SCAILING_SCALE:
            for (int j = 0; j < MULTI_READ; j++)
        {
            #pragma HLS pipeline
            scailings[j].scalings_scale=scailings_local[batch*32+j+8];
        }

        LOOP_LOAD_SCAILING_SHIFT:
            for (int j = 0; j < MULTI_READ; j++)
        {
            #pragma HLS pipeline
            scailings[j].scalings_shift=scailings_local[batch*32+j+16];
        }

            LOOP_LOAD_SCAILING_VAR:
            for (int j = 0; j < MULTI_READ; j++)
        {
            #pragma HLS pipeline
            scailings[j].scalings_var=scailings_local[batch*32+j+24];
        }
            
        LOOP_LOAD_LPS_SKIP:for (int i = 0; i < MULTI_TOTAL; i++)
        {        
                #pragma HLS pipeline

            lps[i].lp_skip = lps_local[batch*24+i];
        }

        LOOP_LOAD_LPS_STAY:for (int i = 0; i < MULTI_TOTAL; i++)
        {
                #pragma HLS pipeline

            lps[i].lp_stay = lps_local[batch*24+i+8];

        }
        LOOP_LOAD_LPS_STEP:for (int i = 0; i < MULTI_TOTAL; i++)
        {
            #pragma HLS pipeline
            lps[i].lp_step = lps_local[batch*24+i+16];
        }
        

        LOOP_LOAD_LENGTH_seq:
        for (int j = 0; j < MULTI_TOTAL; j++)
        {
            #pragma HLS pipeline
            sequence_len_dataflow[j]=length_local[batch*16+j];
                    
        }
        LOOP_LOAD_LENGTH_event:
        for (int j = 0; j < MULTI_TOTAL; j++)
        {
            #pragma HLS pipeline
            event_t_len_dataflow[j]=length_local[batch*16+j+8];
                    
        }
        


            int dataflow_max_event_offset;
            int dataflow_max_sequence_offset;

            int dataflow_band_offset;
        
            dataflow_max_event_offset=max_event_t_lengh;
            dataflow_max_sequence_offset=max_sequence_length;
                
        
            int dataflow_post_offset=dataflow_max_event_offset+dataflow_max_sequence_offset+2;


            int event_vec_offset=dataflow_max_event_offset/EVENT_VEC_WIDTH;
            int read_vec_offset=dataflow_max_sequence_offset/READ_VEC_WIDTH;

           align_dataflow_stream(sequence_len_dataflow,event_t_len_dataflow,dataflow_max_event_offset,dataflow_max_sequence_offset,reset,event_mean0+(event_vec_offset*batch*MULTI_TOTAL),read0+(read_vec_offset*batch*MULTI_TOTAL),scailings,lps,model_levels_mean,model_levels_stdv,model_levels_log_stdv,model_levels_gp_stdv_inverse,
           trace_out_0,
           left_out_0,
           aligned_ref_read_pos_0+(dataflow_post_offset*batch),
           n_align+(MULTI_TOTAL*batch));
    }
       
}

return;

    }
