#include "align_multi.h"


void align_entry(int reset,int sequence_offset,int event_offset,
hls::stream<int>&sequence_offset_stream,
hls::stream<int>&event_offset_stream
)
{

    sequence_offset_stream.write(sequence_offset);

    event_offset_stream.write(event_offset);

    
}





void align_load_event_group(int event_t_length[MULTI_TOTAL],
hls::stream<int>&event_offset_stream_group,
hls::stream<int>event_offset_stream_vec[LOAD_DATAFLOW],
hls::stream<int>event_t_length_out[MULTI_TOTAL], EVENT_VEC* event_mean0,
hls::stream <EVENT_VEC>event_vec_stream[LOAD_DATAFLOW]
)

{

int temp_event_length[MULTI_TOTAL];
int event_offset;

event_offset_stream_group.read(event_offset);


  
    

    for (int i = 0; i < MULTI_TOTAL; i++)
    {     
        #pragma HLS pipeline    
        int temp_event_len_out;
        temp_event_len_out=event_t_length[i];
        temp_event_length[i]=temp_event_len_out;
        //printf("load_event %d length \n",temp_event_len_out);
        event_t_length_out[i].write(temp_event_len_out);  
    }


    for (int i = 0; i < LOAD_DATAFLOW; i++)
    {
        #pragma HLS pipeline
        event_offset_stream_vec[i].write(event_offset);        
    }
      



    int event_iteration=event_offset/EVENT_VEC_WIDTH;

    LOOP_EVENT_MAIN: for (int i = 0; i < event_iteration; i++)
    {
        #pragma HLS loop_tripcount min=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH max=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH avg=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH

         LOOP_LOAD_MEM:for (int j = 0; j<LOAD_DATAFLOW; j++)
            {
                for (int k = 0; k < HALF_READ; k++)
                {
                       #pragma HLS pipeline II=1
                        EVENT_VEC temp_vec;
                        temp_vec=event_mean0[i+event_iteration*(j*HALF_READ+k)]; 
                        event_vec_stream[j].write(temp_vec);    
                }                               
            }
    }     

    LOOP_POST_EVENT: for (int i = event_iteration-1; i >=0; i--)
    {     
                #pragma HLS loop_tripcount min=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH max=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH avg=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH

            LOOP_LOAD_POST_MEM:for (int j = 0; j<LOAD_DATAFLOW; j++)
            {
                for (int k = 0; k < HALF_READ; k++)
                {
                       #pragma HLS pipeline II=1
                        EVENT_VEC temp_vec;
                        temp_vec=event_mean0[i+event_iteration*(j*HALF_READ+k)]; 
                        event_vec_stream[j].write(temp_vec);    
                }                               
            }             
    }


 
}





void align_load_read_group( 
    
int sequence_length[MULTI_TOTAL],
hls::stream<int>&sequence_offset_stream_group,
hls::stream<int>sequence_offset_stream_vec[LOAD_DATAFLOW],

hls::stream<int>sequence_length_out[MULTI_TOTAL], 
READ_VEC* read0,
hls::stream <READ_VEC>read_vec_stream_single[LOAD_DATAFLOW]

)
{

    int temp_sequence_length[MULTI_TOTAL];

    int sequence_offset;
    sequence_offset_stream_group.read(sequence_offset);

        for (int i = 0; i < MULTI_TOTAL; i++)
        {
#pragma HLS pipeline II=1
         int temp_sequence_out;
         temp_sequence_out=sequence_length[i];
         temp_sequence_length[i]=temp_sequence_out;
         sequence_length_out[i].write(temp_sequence_out);
        }
        
   for (int i = 0; i < LOAD_DATAFLOW; i++)
   {
    #pragma HLS pipeline
             sequence_offset_stream_vec[i].write(sequence_offset);
   }
   

     
    int sequence_iteration;

    sequence_iteration=sequence_offset/READ_VEC_WIDTH;



    for (int i = 0; i < sequence_iteration; i++)
    {
            #pragma HLS loop_tripcount min=TEST_READ_OFFSET/READ_VEC_WIDTH max=TEST_READ_OFFSET/READ_VEC_WIDTH avg=TEST_READ_OFFSET/READ_VEC_WIDTH

           for (int j = 0; j<LOAD_DATAFLOW; j++)
            {
                for (int k = 0; k < HALF_READ; k++)
                {
                          
                    #pragma HLS pipeline II=1

                    READ_VEC temp_sequence_vec;
                 
                    temp_sequence_vec=read0[i+sequence_iteration*(j*HALF_READ+k)];
                    
                    read_vec_stream_single[j].write(temp_sequence_vec);
                }
            }

    }

    for (int i =sequence_iteration-1; i>=0; i--)    
    {
                    #pragma HLS loop_tripcount min=TEST_READ_OFFSET/READ_VEC_WIDTH max=TEST_READ_OFFSET/READ_VEC_WIDTH avg=TEST_READ_OFFSET/READ_VEC_WIDTH

           for (int j = 0; j<LOAD_DATAFLOW; j++)
            {
                for (int k = 0; k < HALF_READ; k++)
                {
                          
                    #pragma HLS pipeline II=1

                    READ_VEC temp_sequence_vec;
                 
                    temp_sequence_vec=read0[i+sequence_iteration*(j*HALF_READ+k)];
                    
                    read_vec_stream_single[j].write(temp_sequence_vec);
                }
            }
    }
    
}

