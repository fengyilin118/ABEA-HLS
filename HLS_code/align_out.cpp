#include "align_multi.h"



void align_channel(
hls::stream <int>& band_length_single_out,

hls::stream <TRACE_VEC>&trace_vec_stream, 
hls::stream <int>&left_event_stream,

hls::stream<int>&max_event_idx_stream,

hls::stream <int> &max_event_idx_post_stream,
hls::stream <int> &max_band_post_out_stream,

hls::stream <TRACE_VEC>trace_vec_post_stream[2],
hls::stream <int>band_lower_left_event_post_stream[2],
AP_TRACE_VEC* trace_out_0,
LOWER_LEFT_VEC* left_out_0
)
{
 
int band_length[MULTI_TOTAL];

int max_length=0;


#pragma HLS array_partition variable=band_length type=complete


    for (int j = 0; j < MULTI_READ; j++)
    {

        
        #pragma HLS pipeline
            int temp_length;
            band_length_single_out.read(temp_length);
            band_length[j]=temp_length;
            if(temp_length>max_length)
            {
                max_length=temp_length;
            }

 

    }


for (int i = 0; i < max_length; i++)

    {
        #pragma HLS loop_tripcount min=TEST_BAND_LEN max=TEST_BAND_LEN avg=TEST_BAND_LEN
        #pragma HLS pipeline II=8

        TRACE_VEC arr_trace_vec[MULTI_READ];
        AP_TRACE_VEC ap_trace;
        LOWER_LEFT_VEC lower_left_out;

        for (int k = 0; k < MULTI_READ; k++)
        {               
                         
                    int temp_left_vec_event;
                    TRACE_VEC temp_trace_vec;

                    if (i<band_length[k])
                    {
                        
                    trace_vec_stream .read(temp_trace_vec);
                    left_event_stream.read(temp_left_vec_event);
                    
                    } 
                    arr_trace_vec[k]=temp_trace_vec;
                    lower_left_out[k]=temp_left_vec_event;

        
                
        }
    
            ap_trace_encode(arr_trace_vec,&ap_trace);


         trace_out_0[i]=ap_trace;
         left_out_0[i]=lower_left_out;

     
    }   

    int max_event_idx[MULTI_TOTAL];

 
        for (int k = 0; k < MULTI_READ; k++)
        {
       
            
                #pragma HLS pipeline
                int temp_max_event_idx;
                max_event_idx_stream.read(temp_max_event_idx);

                //printf("align_channel max _event %d \n",temp_max_event_idx);


                max_event_idx[k]=temp_max_event_idx;    
            


        }


   //post section

max_band_post_out_stream.write(max_length);

LOOP_META_POST:
    for (int i = 0; i < MULTI_TOTAL; i++)
    {
        #pragma HLS pipeline
        int temp_event_idx;
        temp_event_idx=max_event_idx[i];

        max_event_idx_post_stream.write(temp_event_idx);
    }


    
int max_offset=max_length/8+1;

LOOP_LOAD_TRACE_POST:
    for (int i = max_offset; i>=0; i--)
    {
            #pragma HLS loop_tripcount min=TEST_BAND_LEN max=TEST_BAND_LEN avg=TEST_BAND_LEN

        #pragma HLS pipeline II=32

        LOWER_LEFT_VEC lower_left_vec[8];
        AP_TRACE_VEC ap_trace[8];

        
        for (int j = 0; j < 8; j++)
        {
            ap_trace[j]=trace_out_0[i*8+j];
            lower_left_vec[j]=left_out_0[i*8+j];

        }
        
        for (int j = 7; j >=0; j--)
        {
            TRACE_VEC arr_trace_vec[MULTI_READ];
            ap_trace_decode(ap_trace[j],arr_trace_vec);

            if(i*8+j<max_length)
            {
                    for (int k = 0; k < MULTI_TOTAL; k+=2)
                {
                    for (int q = 0; q < 2; q++)
                    {
                    TRACE_VEC temp_vec;
                    int temp_lower_left_event=lower_left_vec[j][k+q];
                    temp_vec=arr_trace_vec[k+q];                      
                    trace_vec_post_stream[q].write(temp_vec);
                    band_lower_left_event_post_stream[q].write(temp_lower_left_event);             
                    }    
                }
                
            }

        }
        

       
                      
    }


    
}





void post_pos_out(
hls::stream<bool> event_last_stream[MULTI_TOTAL],
hls::stream<bool> model_last_stream[MULTI_TOTAL],
hls::stream<int> &max_band_num_stream,
hls::stream<int> ref_pos_post_stream[MULTI_TOTAL],
hls::stream<int> read_pos_post_stream[MULTI_TOTAL],
hls::stream<int> &n_pairs_post_stream,
POS_POST_VEC* aligned_ref_read_pos_0,

int* n_align
)
{


int max_band_num;
max_band_num_stream.read(max_band_num);
int max_band_offset=max_band_num/8+1;

for (int i=max_band_offset-1; i>=0; i--)
{
        #pragma HLS loop_tripcount min=TEST_BAND_LEN/8 max=TEST_BAND_LEN/8 avg=TEST_BAND_LEN/8
        POS_POST_VEC ref_read_pos_vec[8];
        #pragma HLS pipeline II=8


    for (int j = 7; j >=0; j--)
    {
        if(i*8+j<max_band_num)
        {
                for (int k = 0; k < MULTI_READ; k++)
            {
                int temp_ref_pos;
                int temp_read_pos;
    
                ref_pos_post_stream[k].read(temp_ref_pos);
                read_pos_post_stream[k].read(temp_read_pos);
        
                ref_read_pos_vec[j][k*2]=temp_ref_pos;
                ref_read_pos_vec[j][k*2+1]=temp_read_pos;
            }

        }
    }


    for (int j = 0; j < 8; j++)
    {
        aligned_ref_read_pos_0[i*8+j]=ref_read_pos_vec[j];
    }
    
    
    
  
}

for (int i = 0; i < MULTI_TOTAL; i++)
{

    #pragma HLS pipeline
    int temp_n_pairs;

    n_pairs_post_stream.read(temp_n_pairs);
    n_align[i]=temp_n_pairs;
}
  

for (int i = 0; i < MULTI_TOTAL; i++)
{
    event_last_stream[i].read();
    model_last_stream[i].read();
}


}



