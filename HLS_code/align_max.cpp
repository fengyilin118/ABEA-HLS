#include "align_multi.h"



void align_max_v2(
hls::stream <int> &event_n_stream,
hls::stream <int> &band_length_single,
hls::stream <int> &band_length_single_out,
hls::stream <FIXED_POINT>&max_cell_stream,
hls::stream <int>& max_cell_event_idx_stream,
hls::stream <int> &max_event_idx_stream
)
{


    int band_length[MULTI_READ];
    int event_n[MULTI_READ];
    int max_length=0;
    int max_event_idx[MULTI_READ];
    FIXED_POINT max_cell[MULTI_READ];
 

    LOOP_LENGTH_INIT:for (int i = 0; i < MULTI_READ; i++)
    {

        #pragma HLS pipeline
        int temp_band_length;
        band_length_single.read(temp_band_length);
        band_length_single_out.write(temp_band_length);
        max_cell[i]=NEG_INF;
        max_event_idx[i]=0;

        band_length[i]=temp_band_length;
        if(temp_band_length>max_length)
        {
            max_length=temp_band_length;
        }
        
        int temp_event_n;

        event_n_stream.read(temp_event_n);
        event_n[i]=temp_event_n;

    }
       
    


LOOP_MAX:for (int i = 0; i <max_length; i++)
    {
        //printf("trace vec max %d \n",i);
    #pragma HLS loop_tripcount min=TEST_BAND_LEN max=TEST_BAND_LEN avg=TEST_BAND_LEN

    	LOOP_MAX_READ: for (int j = 0; j < MULTI_READ; j++)
        {
            #pragma HLS pipeline II=1
            #pragma HLS UNROLL factor=1

        


            //printf("align max \n");

            int current_max_event_idx;
            FIXED_POINT current_max_cell;

                     
            int temp_band_lower_left_event;
            int temp_band_lower_left_kmer;



            if (i<band_length[j])
           {     
                        
     
        



            FIXED_POINT temp_max_cell;
            int temp_max_event_idx;


            max_cell_stream.read(temp_max_cell);
            max_cell_event_idx_stream.read(temp_max_event_idx);




            compare_band_fix(event_n[j],max_cell[j],temp_max_cell,&current_max_cell,max_event_idx[j],temp_max_event_idx,&current_max_event_idx);
                
            max_cell[j]=current_max_cell;
            max_event_idx[j]=current_max_event_idx;

           }


        }   


    }


     LOOP_UPDATE_MAX: for (int i = 0; i < MULTI_READ; i++)
    {
        #pragma HLS pipeline
        int temp_event_idx;
        temp_event_idx=max_event_idx[i];
        

        //printf("align_max max _event %f %d \n",temp_max_score,temp_event_idx);
        max_event_idx_stream.write(temp_event_idx);
    }

}

