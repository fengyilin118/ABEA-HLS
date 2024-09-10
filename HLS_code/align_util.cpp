#include "align_multi.h"



void kmer_generator(char read, char shift_register[KMER_SIZE], uint32_t *kmer_ranks)
{
#pragma HLS PIPELINE

#pragma HLS array_partition variable=shift_register type=complete dim=0



SHIF_LOOP:
    for (int i = 0; i < KMER_SIZE - 1; i++)
    {
#pragma HLS unroll
        shift_register[i] = shift_register[i + 1];
    }
   
    shift_register[KMER_SIZE - 1] = read;
/*
printf("shifter \n");
    for (int i = 0; i < KMER_SIZE; i++)
    {
    printf(" %c  ",shift_register[i]); 
    }
    printf("shifter \n");
    */

    *kmer_ranks = get_kmer_rank(shift_register, KMER_SIZE);

}

void kmer_generator_post(char read, char shift_register[KMER_SIZE], uint32_t *kmer_ranks)
{
#pragma HLS PIPELINE

#pragma HLS array_partition variable=shift_register type=complete dim=0



SHIF_LOOP:
    for (int i = KMER_SIZE - 1; i>0; i--)
    {
#pragma HLS unroll
        shift_register[i] = shift_register[i- 1];
    }
   
    shift_register[0] = read;
/*
printf("shifter \n");
    for (int i = 0; i < KMER_SIZE; i++)
    {
    printf(" %c  ",shift_register[i]); 
    }
    printf("shifter \n");
    */

    *kmer_ranks = get_kmer_rank(shift_register, KMER_SIZE);

}


uint32_t get_kmer_rank(char *str, uint32_t k)
{
    uint32_t r = 0;
    for (uint32_t i = 0; i < k; i++)
    {
        char base = str[k - i - 1];
        uint32_t rank;
        
        
        switch (base)
        {
        case 'A':
            rank=0;
            break;
        case 'C':
            rank=1;
            break;

        case 'G':
            rank=2;
            break;

        case 'T':
            rank=3;

            break;
        default:
            rank=0;
            break;
        }
        r += rank << (i << 1);
    }
    return r;
}





void align_kernel_score_calculate(FIXED_POINT score_d, FIXED_POINT score_u, FIXED_POINT score_l, FIXED_POINT *max_score_out, FIXED_TRACE *from_score_out)
{

    
    FIXED_TRACE from;
    FIXED_POINT max_score;

    /*  float max_score = score_d;
        uint8_t from = FROM_D;

        max_score = score_u > max_score ? score_u : max_score;
        from = max_score == score_u ? FROM_U : from;
        max_score = score_l > max_score ? score_l : max_score;
        from = max_score == score_l ? FROM_L : from;
    */
    max_score = score_d;
    from = FROM_D;

    if (score_u > max_score)
        max_score = score_u;
    else
        max_score = max_score;

    if (max_score == score_u)
        from = FROM_U;
    else
        from = from;

    if (score_l > max_score)
        max_score = score_l;
    else
        max_score = max_score;

    if (max_score == score_l)
        from = FROM_L;
    else
        from = from;

    *max_score_out = max_score;
    *from_score_out = from;
}


void ap_trace_encode(TRACE_VEC input_vec[MULTI_READ],AP_TRACE_VEC* output_vec)
{
    #pragma HLS inline
    AP_TRACE_VEC temp_vec;

    for (int i = 0; i < MULTI_READ; i++)
    {            

        for (int j = 0; j < BAND_WIDTH; j++)
        {
            temp_vec[i](j*2+1,j*2)=input_vec[i][j];
        }
        
    }

    *output_vec=temp_vec;
}


void ap_trace_decode(AP_TRACE_VEC input_vec, TRACE_VEC output_vec[MULTI_READ])
{
        #pragma HLS inline

    for (int i = 0; i < MULTI_READ; i++)
    {
        for (int j = 0; j < BAND_WIDTH; j++)
        {

            output_vec[i][j]=input_vec[i](j*2+1,j*2);
        }
        
    }
    

}








void event_shift(event_info events[BAND_WIDTH+1],event_info events_local[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], model_info model_levels_local[BAND_WIDTH], event_info new_event)
{
#pragma HLS inline
     LOOP_BAND_EVENT_SHIFT:     for (int j = 0; j < BAND_WIDTH; j++)
        {
            #pragma HLS UNROLL
            events[j]=events[j+1];
            events_local[j]=events[j+1];


            model_levels_local[j]=model_levels_cur[j];
            
        }
        
        events[BAND_WIDTH]=new_event;
        events_local[BAND_WIDTH]=new_event;
}


void model_shift(event_info events[BAND_WIDTH+1],event_info events_local[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], model_info model_levels_local[BAND_WIDTH], model_info new_model)
{
#pragma HLS inline
     
    LOOP_BAND_KMER_SHIFT:    for (int j = 0; j < BAND_WIDTH-1; j++)
        {
            #pragma HLS UNROLL
            model_levels_cur[j]= model_levels_cur[j+1];
            model_levels_local[j]=model_levels_cur[j+1];

            events_local[j]=events[j];
        }
        model_levels_cur[BAND_WIDTH-1]=new_model;
        model_levels_local[BAND_WIDTH-1]=new_model;

        events_local[BAND_WIDTH-1]=events[BAND_WIDTH-1];
        events_local[BAND_WIDTH]=events[BAND_WIDTH];


    }

void compare_band(int n_event,band_score_info band_a,band_score_info band_b, band_score_info* band_result,int event_idx_a,int event_idx_b,int* event_idx_result)
{
    
    FIXED_POINT score_a,score_b;
    if (band_a.score==NEG_INF)
    {
        score_a=NEG_INF;
    }
    else
    {
             score_a=band_a.score+((FIXED_POINT)n_event-(FIXED_POINT)event_idx_a)*lp_trim_const;

    }
    
    if (band_b.score==NEG_INF)
    {
        score_b=NEG_INF;
    }
    else
    {
            score_b=band_b.score+((FIXED_POINT)n_event-(FIXED_POINT)event_idx_b)*lp_trim_const;

    }
    
    

#pragma HLS inline
    //printf("score a %d  %f, score b %d %f \n",band_a.event_idx, score_a.to_double(),band_b.event_idx,score_b.to_double());
    
    

    if(score_a>score_b)
    {
        *band_result=band_a;
        *event_idx_result=event_idx_a;
    }
    else
    {
        *band_result=band_b;
        *event_idx_result=event_idx_b;
    }
}


void compare_band_fix(int n_event,FIXED_POINT band_a,FIXED_POINT band_b, FIXED_POINT* band_result,int event_idx_a,int event_idx_b,int* event_idx_result)
{
    
    FIXED_POINT score_a,score_b;
    if (band_a==NEG_INF)
    {
        score_a=NEG_INF;
    }
    else
    {
             score_a=band_a+((FIXED_POINT)n_event-(FIXED_POINT)event_idx_a)*lp_trim_const;

    }
    
    if (band_b==NEG_INF)
    {
        score_b=NEG_INF;
    }
    else
    {
            score_b=band_b+((FIXED_POINT)n_event-(FIXED_POINT)event_idx_b)*lp_trim_const;

    }
    
    

#pragma HLS inline
    //printf("score a %d  %f, score b %d %f \n",band_a.event_idx, score_a.to_double(),band_b.event_idx,score_b.to_double());
    
    

    if(score_a>score_b)
    {
        *band_result=band_a;
        *event_idx_result=event_idx_a;
    }
    else
    {
        *band_result=band_b;
        *event_idx_result=event_idx_b;
    }
}

