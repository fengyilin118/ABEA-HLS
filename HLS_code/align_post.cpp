#include "align_multi.h"

void align_post_multi(

hls::stream<int> &kmer_n_post_stream,
hls::stream<int> &event_n_post_stream,
hls::stream<FIXED_POINT_MEAN>  events_post_stream_single[MULTI_TOTAL],

hls::stream<model_info>model_levels_post_stream[MULTI_TOTAL],
hls::stream<TRACE_VEC> trace_vec_post[2],

hls::stream<int> &max_event_idx_post_stream,
hls::stream<int> band_lower_left_event_post_stream[2],
hls::stream<int> ref_pos_post_stream[MULTI_TOTAL],
hls::stream<int> read_pos_post_stream[MULTI_TOTAL],
hls::stream<int> &n_pairs_post_stream
)


{
    //POST SECTION
#pragma HLS pipeline off
int max_band;
int local_kmer_n[MULTI_TOTAL];
int local_event_n[MULTI_TOTAL];
int band_num_local[MULTI_TOTAL];

int max_event_idx[MULTI_TOTAL];
int curr_event_idx[MULTI_TOTAL];



FIXED_POINT sum_emission[MULTI_TOTAL];
bool spanned[MULTI_TOTAL];
int last_ref_pos[MULTI_TOTAL];

bool forward_skip[MULTI_TOTAL];

int outIndex[MULTI_TOTAL];

int curr_kmer_idx[MULTI_TOTAL];
int curr_band_idx[MULTI_TOTAL];

int n_aligned_events[MULTI_TOTAL];
int max_gap[MULTI_TOTAL];

int curr_gap[MULTI_TOTAL];

bool finish_flag[MULTI_TOTAL];
max_band=0;

const int POST_UNROLL_FACTOR=MULTI_TOTAL;

#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=local_kmer_n
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=local_event_n
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=band_num_local
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=max_event_idx
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=curr_event_idx
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=sum_emission
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=spanned
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=last_ref_pos
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=forward_skip
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=outIndex
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=curr_kmer_idx
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=curr_band_idx
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=n_aligned_events
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=max_gap
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=curr_gap
#pragma HLS array_partition type=cyclic factor=POST_UNROLL_FACTOR variable=finish_flag



LOOP_POST_LENGTH_INIT: for (int i = 0; i < MULTI_READ; i++)
{
    for (int j = 0; j < MULTI_GROUP; j++)
    {
        
        
        #pragma HLS pipeline
        
        int temp_event_n;
        int temp_kmer_n;
        int temp_band_num;


        kmer_n_post_stream.read(temp_kmer_n);
        event_n_post_stream.read(temp_event_n);



        local_kmer_n[j*MULTI_READ+i]=temp_kmer_n;
        local_event_n[j*MULTI_READ+i]=temp_event_n;
        temp_band_num=temp_kmer_n+temp_event_n+2;
        band_num_local[j*MULTI_READ+i]=temp_band_num;


        if(max_band<temp_band_num)
        {
            max_band=temp_band_num;
        }


    }
    
    /* code */
}




LOOP_POST_INIT:for (int i = 0; i < MULTI_TOTAL; i++)
{
 #pragma HLS pipeline
        int temp_event_idx;

    
        max_event_idx_post_stream.read(temp_event_idx);
        if (temp_event_idx>local_event_n[i]||temp_event_idx<0)
        {
            temp_event_idx=local_event_n[i]-1;
        }
        
        max_event_idx[i]=temp_event_idx;



        curr_kmer_idx[i]=local_kmer_n[i];
        curr_event_idx[i]=local_event_n[i];
        curr_band_idx[i]=band_num_local[i];

        n_aligned_events[i]=0;

       
        finish_flag[i]=true;
        //printf("max_event_idx %d total band %d  \n",temp_event_idx,curr_band_idx[i]);

}

model_info curr_model[MULTI_TOTAL];
FIXED_POINT_MEAN curr_event_mean[MULTI_TOTAL];
#pragma HLS array_partition variable=curr_model type=complete
#pragma HLS array_partition variable=curr_event_mean type=complete

LOOP_POST: for (int i = max_band-1; i>=0; i--)
{
    #pragma HLS loop_tripcount min=TEST_BAND_LEN max=TEST_BAND_LEN avg=TEST_BAND_LEN
    #pragma HLS pipeline II=4

    TRACE_VEC temp_trace_vec[MULTI_TOTAL];
    int temp_band_lower_left_event_arr[MULTI_TOTAL];
    for (int j=0;j<MULTI_TOTAL;j+=MULTI_GROUP*2)
    {
        for (int k = 0; k < MULTI_GROUP*2; k++)
        {
            band_lower_left_event_post_stream[k].read(temp_band_lower_left_event_arr[j+k]);
            trace_vec_post[k].read(temp_trace_vec[j+k]);
      
        }
           
    }
    
    
    for (int j = 0; j < MULTI_TOTAL; j++)
    {
                        
                    FIXED_TRACE temp_trace_table[BAND_WIDTH];
                    #pragma HLS array_partition variable=temp_trace_table type=complete

                    int ref_pos=-999;
                    int read_pos=-999;

                

            
                    for (int q = 0; q < BAND_WIDTH; q++)
                        {
                            temp_trace_table[q]=temp_trace_vec[j][q];
                        }

                 int temp_band_lower_left_event=temp_band_lower_left_event_arr[j];
                if(finish_flag[j])
                {
                            
               
                    if(i<band_num_local[j]&&band_num_local[j]>2)
                    {
                       
                    
                        
                        if(curr_event_idx[j]>max_event_idx[j])
                        {
                        
                            events_post_stream_single[j].read(curr_event_mean[j]);
                            curr_event_idx[j]--;

                            if(curr_event_idx[j]==max_event_idx[j])
                            {
                                model_levels_post_stream[j].read(curr_model[j]);
                                curr_kmer_idx[j] --;
                                //post starts initalize variables
                                        outIndex[j]=0;
                                        sum_emission[j]=(FIXED_POINT)0;
                                        max_gap[j]=0;
                                        curr_gap[j]=0;
                                        spanned[j]=true;
                                        forward_skip[j]=false;
                            }
                        
                        }
                        else
                        {
                            if(!forward_skip[j])
                            {
                                    ref_pos=curr_kmer_idx[j];
                                    read_pos=curr_event_idx[j];

                                if(outIndex[j]==0&&ref_pos!=(local_kmer_n[j]-1))
                                {
                                    spanned[j]=false;
                                }

                                last_ref_pos[j] = curr_kmer_idx[j];


                                outIndex[j]++;

                                    FIXED_POINT_MEAN scaledLevel = curr_event_mean[j];

                                    FIXED_POINT gp_mean = curr_model[j].gp_mean;
                                    FIXED_POINT_STDV inverse_gp_stdv = curr_model[j].gp_stdv_inverse; //scaling.var = 1;
                                    FIXED_POINT_STDV gp_log_stdv = curr_model[j].log_stdv;


                                    FIXED_POINT_STDV log_inv_sqrt_2pi = (FIXED_POINT)-0.918938;
                                    FIXED_POINT a = (scaledLevel - gp_mean)*inverse_gp_stdv;
                                    FIXED_POINT tempLogProb  = log_inv_sqrt_2pi - gp_log_stdv + ((FIXED_POINT)-0.5 * a * a);

                                    sum_emission[j]+=tempLogProb;

                                    n_aligned_events[j]+=1;


                                    int offset = temp_band_lower_left_event - curr_event_idx[j];

                                    FIXED_TRACE from = temp_trace_table[offset];

                                

                                    if (from == FROM_D) {
                                    
                                        curr_gap[j] = 0;
                                        
                                        if(curr_kmer_idx[j]>0 &&curr_event_idx[j]>0)
                                        {
                                            events_post_stream_single[j].read(curr_event_mean[j]);
                                            model_levels_post_stream[j].read(curr_model[j]);
                                            curr_kmer_idx[j] --;
                                            curr_event_idx[j] --;

                                            forward_skip[j]=true;
                                        }
                                        else
                                        {
                                            finish_flag[j]=false;
                                        }                       

                                    } else if (from == FROM_U) {

                                            if(curr_event_idx[j]>0)

                                            {
                                            events_post_stream_single[j].read(curr_event_mean[j]);
                                            curr_event_idx[j] -= 1;   

                                        
                                            curr_gap[j] = 0;
                                            }
                                            else
                                            {
                                                finish_flag[j]=false;
                                            }
                                    
                                        } else {
                                            curr_gap[j] += 1;

                                            if(curr_gap[j] > max_gap[j])
                                                max_gap[j] = curr_gap[j];
                                
                                            if(curr_kmer_idx[j]>0)
                                            {
                                            model_levels_post_stream[j].read(curr_model[j]);
                                            curr_kmer_idx[j] --;
                                            }
                                            else
                                            {
                                                finish_flag[j]=false;
                                            }             
                                        
                                        }
                                }
                                else                    
                                {
                                    forward_skip[j]=false;
                            
                                }       
                            }
                        }  
                    }
                    ref_pos_post_stream[j].write(ref_pos);
                    read_pos_post_stream[j].write(read_pos);  
    
    
        }
    

}

 for (int i = 0; i < MULTI_TOTAL; i++)
            {

              while (curr_event_idx[i]>0)
              {
                #pragma HLS loop_tripcount min=1 max=1 avg=1

                #pragma HLS pipeline
                events_post_stream_single[i].read();
                curr_event_idx[i]--;
              }

              while (curr_kmer_idx[i]>0)
              {
                #pragma HLS loop_tripcount min=1 max=1 avg=1

                #pragma HLS pipeline
                model_levels_post_stream[i].read();
                curr_kmer_idx[i] --;
              }                   
                
            }


for (int i = 0; i < MULTI_TOTAL; i++)
{
    #pragma HLS pipeline II=4
    if(last_ref_pos[i]!=0)
    {
        spanned[i]=false;
    }


    if (n_aligned_events[i]==0)
    {
        outIndex[i]=0;
    }
    else
    {
            FIXED_POINT_STDV inverse_n_aligned;
        if(n_aligned_events[i]>32768)
        {
            inverse_n_aligned=0;
        }
        else
        {
            inverse_n_aligned=(ap_fixed<3,2>)1/(ap_fixed<18,17>)n_aligned_events[i];
        }
        FIXED_POINT avg_log_emission = sum_emission[i]*inverse_n_aligned;
        int  max_gap_threshold = 50;
        FIXED_POINT min_average_log_emission =  (FIXED_POINT)-5.0;
        if (avg_log_emission < min_average_log_emission || !spanned[i] ||max_gap[i] > max_gap_threshold)
        {
            outIndex[i] = 0;
        }	

    }

    
    int temp_outIndex=outIndex[i];
    n_pairs_post_stream.write(temp_outIndex);

    //printf("spanned %d last_ref_ post %d sum_emission %f n_aligned_events %d max_gap %d \n",spanned[i], last_ref_pos[i],sum_emission[i].to_float(),n_aligned_events[i],max_gap[i]);

}



}

