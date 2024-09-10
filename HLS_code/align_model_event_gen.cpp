#include "align_multi.h"


void align_event_vec_multi_v2(
hls::stream <int>&event_offset_stream_vec,
hls::stream <int>event_offset_stream_gen[HALF_READ],
hls::stream <EVENT_VEC>&event_vec_stream,
hls::stream <float>event_float_stream[HALF_READ])
{


    uint32_t event_offset;
    int event_offset_stream;
    event_offset_stream_vec.read(event_offset_stream);
    event_offset=event_offset_stream;
LOOP_EVENT_OFFSET_DIST:
for (int i = 0; i < HALF_READ; i++)
{
    #pragma HLS pipeline
    event_offset_stream_gen[i].write(event_offset_stream);
}

    uint32_t event_iteration=event_offset/EVENT_VEC_WIDTH;

LOOP_EVNT_VEC:
for (uint32_t i = 0; i < event_iteration*2; i++)
{
    #pragma HLS loop_tripcount min=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH*2 max=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH*2 avg=TEST_EVENT_OFFSET/EVENT_VEC_WIDTH*2

    for (uint32_t k = 0; k < HALF_READ; k++)
    {
        #pragma HLS pipeline II=EVENT_VEC_WIDTH
        EVENT_VEC temp_vec;
        event_vec_stream.read(temp_vec);
        for (uint32_t j = 0; j < EVENT_VEC_WIDTH; j++)
        {
            float event_to_stream;
            int event_idx;
            if(i<event_iteration)
            {
                event_to_stream=temp_vec[j];
            }
            else
            {
                event_to_stream=temp_vec[EVENT_VEC_WIDTH-j-1];
            }        
            event_float_stream[k].write(event_to_stream);
        
        }  

    }   
}

}

void align_read_vec_multi_v2(
hls::stream<int>&sequence_offset_stream_vec,
hls::stream<int>sequence_offset_stream_gen[HALF_READ],
hls::stream<READ_VEC>&read_vec_stream,
hls::stream<char>read_char_stream[HALF_READ])
{

uint32_t sequence_offset;
int sequence_offset_stream;
sequence_offset_stream_vec.read(sequence_offset_stream);
sequence_offset=sequence_offset_stream;
for (int i = 0; i < HALF_READ; i++)
{
    #pragma HLS pipeline
    sequence_offset_stream_gen[i].write(sequence_offset_stream);
}


    uint32_t kmer_iteration=sequence_offset/READ_VEC_WIDTH;

    LOOP_READ_VEC:
    for (uint32_t i = 0; i < kmer_iteration*2; i++)
    {
            #pragma HLS loop_tripcount min=TEST_READ_OFFSET/READ_VEC_WIDTH*2 max=TEST_READ_OFFSET/READ_VEC_WIDTH*2 avg=TEST_READ_OFFSET/READ_VEC_WIDTH*2
        for (uint32_t k = 0; k < HALF_READ; k++)
        {
         #pragma HLS pipeline II=READ_VEC_WIDTH
            READ_VEC temp_vec;
            read_vec_stream.read(temp_vec);


            for (uint32_t j = 0; j < READ_VEC_WIDTH; j++)
            {               
                char temp_read;

                if(i<kmer_iteration)
                {
                temp_read=temp_vec[j];

                }
                else
                {
                    temp_read=temp_vec[READ_VEC_WIDTH-j-1];
                }
                read_char_stream[k] .write(temp_read);

            }
        }


    }
    


}



void align_event_mean_gen(
hls::stream <int>&event_offset_stream_gen,
hls::stream<int>&event_t_length_in,
hls::stream<int>&event_t_length_out,
hls::stream <float>&event_float_stream,
hls::stream <FIXED_POINT_MEAN>&event_stream,
hls::stream <FIXED_POINT_MEAN>&event_stream_post,
hls::stream <bool> &event_last_stream)
{

uint32_t event_offset;
int event_offset_stream;
uint32_t temp_event_length;
int event_length_stream;
event_offset_stream_gen.read(event_offset_stream);
event_offset=event_offset_stream;

event_t_length_in.read(event_length_stream);
event_t_length_out.write(event_length_stream);
temp_event_length=event_length_stream;

LOOP_EVENT_MEAN_GEN:
for (uint32_t i = 0; i < event_offset*2; i++)
{
    #pragma HLS pipeline II=1
    float temp_event_float;
    FIXED_POINT_MEAN event_to_stream;
    uint32_t event_idx;
    event_float_stream.read(temp_event_float);  

    event_to_stream=(FIXED_POINT_MEAN)temp_event_float;

    #pragma HLS loop_tripcount min=TEST_EVENT_OFFSET*2 max=TEST_EVENT_OFFSET*2 avg=TEST_EVENT_OFFSET*2

    
        if(i<event_offset)
        {
            event_idx=i;

        }
        else
        {
            event_idx=event_offset*2-i-1;
        }

        if(event_idx<temp_event_length)
        {
            if(i<event_offset)
            {
                event_stream.write(event_to_stream);

            }
            else
            {
                event_stream_post.write(event_to_stream);
            }

        }

        if(i==event_offset*2-1)
        {
            event_last_stream.write(true);
        }     
      
}

}

void align_model_gen(
hls::stream<int>&sequence_offset_stream_gen,
hls::stream<int>&sequence_length_in,
hls::stream<int>&kmer_length_out,
hls::stream<char>&read_char_stream, 
FIXED_POINT_MEAN model_levels_mean[MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_stdv[MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_log_stdv[MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_gp_stdv_inverse[MAX_MODEL_NUM],

hls::stream<model_info>&model_levels_cur_stream,
hls::stream<model_info>&model_levels_cur_post_stream,
scailing_info scailings,
hls::stream<bool>&model_last_stream)
{
    scailing_info local_scailings;

    local_scailings=scailings;



        int temp_sequence_length;
        int temp_kmer_length;
        sequence_length_in.read(temp_sequence_length);

        if(temp_sequence_length==0)
        {
            temp_kmer_length=0;
        }
        else
        {
            temp_kmer_length=temp_sequence_length-KMER_SIZE+1;
        }
  
        kmer_length_out.write(temp_kmer_length); 



        int sequence_offset;
        sequence_offset_stream_gen.read(sequence_offset);

    char shift_register[KMER_SIZE];
    char shift_register_post[KMER_SIZE];

    
LOOP_MODEL_GEN:
for (int i = 0; i < sequence_offset*2; i++)
    {    
        #pragma HLS loop_tripcount min=TEST_READ_OFFSET max=TEST_READ_OFFSET avg=TEST_READ_OFFSET
        #pragma HLS pipeline II=1

        char temp_read;
        read_char_stream.read(temp_read);
  
        int kmer_idx;

                        
        uint32_t kmer_ranks;

            if(i<sequence_offset)
            {
               kmer_idx=i-KMER_SIZE+1;
               kmer_generator(temp_read,shift_register, &kmer_ranks);

            }
            else
            {
                kmer_idx=2*sequence_offset-i-1;
                kmer_generator_post(temp_read,shift_register_post, &kmer_ranks);
            }
                 model_info temp_model_levels_cur;
            if(kmer_ranks<4096)
            {
                
                FIXED_POINT_MEAN curr_mean=model_levels_mean[kmer_ranks];

                temp_model_levels_cur.stdv=model_levels_stdv[kmer_ranks];
                temp_model_levels_cur.log_stdv=model_levels_log_stdv[kmer_ranks];
                //temp_model_levels_cur.id=kmer_idx;

                temp_model_levels_cur.gp_mean=local_scailings.scalings_scale * curr_mean + local_scailings.scalings_shift;

       
                temp_model_levels_cur.gp_stdv_inverse=model_levels_gp_stdv_inverse[kmer_ranks];
            }
                                        
           


                   if(kmer_idx>=0&&kmer_idx<temp_kmer_length)
                {
                    
                    

                    if (i<sequence_offset)
                    {
                        model_levels_cur_stream.write(temp_model_levels_cur);
                    }
                    else
                    {
                        model_levels_cur_post_stream.write(temp_model_levels_cur);
                    }

                }

                if (i==(sequence_offset*2-1))
                {
                    model_last_stream.write(true);
                }
                

    }

    


}


