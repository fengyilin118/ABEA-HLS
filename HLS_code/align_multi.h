#define AP_INT_MAX_W 2048
#include <ap_fixed.h>
#include <ap_int.h>
#include <math.h>

#include <hls_math.h>
#include <hls_stream.h>
#include <hls_vector.h>
#include "stdio.h"
#include "fstream"
//#include "iostream"

//#define DEBUG_MOD 1
const int TEST_BATCH=1;
const int TEST_READ_LEN=1414;
const int TEST_EVENT_LEN=2476;
const int TEST_READ_OFFSET=(TEST_READ_LEN/256+1)*256;
const int TEST_EVENT_OFFSET=(TEST_EVENT_LEN/256+1)*256;
const int TEST_BAND_LEN=TEST_READ_OFFSET+TEST_EVENT_OFFSET+2;




const int EVENT_FIFO_DEPTH=120000;// fit for largest read ever seen
const int READ_FIFO_DEPTH=60000;//
const int BAND_WIDTH=32;//user changeable, typically <100


const int READ_VEC_WIDTH=32;//< 64
const int EVENT_VEC_WIDTH=16;//< 64



#define EVENT_LENGTH 2476 //testing purposes
#define READ_LENGTH 1414 //testing purposes
// #define or function
#define MOVE_DOWN 0
#define MOVE_RIGHT 1
#define FROM_D 0
#define FROM_U 1
#define FROM_L 2
#define event_move_down(event_idx) (event_idx + 1)
#define kmer_move_right(kmer_id) (kmer_id + 1)


#define event_at_offset(bi, offset) band_lower_left_event[(bi)] - (offset)
#define kmer_at_offset(bi, offset) band_lower_left_kmer[(bi)] + (offset)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left_kmer[bi]
#define band_event_to_offset(bi, ei) band_lower_left_event[bi] - (ei)
#define is_offset_valid(offset) (offset) >= 0 && (offset) < BAND_WIDTH
#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)

#define event_at_offset_stream(offset) band_lower_left_event_current - (offset)
#define kmer_at_offset_stream(offset) band_lower_left_kmer_current + (offset)
#define event_at_offset_stream(offset) band_lower_left_event_current - (offset)
#define kmer_at_offset_stream(offset) band_lower_left_kmer_current + (offset)
#define band_kmer_to_offset_stream(ki) (ki) - band_lower_left_kmer_current
#define band_event_to_offset_stream(ei) band_lower_left_event_current - (ei)
#define is_offset_valid(offset) (offset) >= 0 && (offset) < BAND_WIDTH
#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)





#define TEST_BAND 3882
#define TEST_KMER 50
#define TEST_EVENT 50


typedef ap_fixed<34,24> FIXED_POINT;//balanced range and precision,
typedef ap_fixed<20,10> FIXED_POINT_MEAN;
typedef ap_fixed<20,10> FIXED_POINT_STDV;

typedef ap_uint <2> FIXED_TRACE;



typedef double EVENT_VEC_AP;
typedef char READ_VEC_AP; 


const FIXED_POINT lp_skip_const=-23.02585;
const FIXED_POINT lp_trim_const=-4.60517;


const int MULTI_READ=8;//number of pipeline reads
const int MULTI_GROUP=1;
const int MULTI_TOTAL=8;
const int HALF_READ=4;

const int LOAD_DATAFLOW=2;
const int KMER_SIZE=6;
const int MAX_MODEL_NUM=4096;//MODEL SIZE



const FIXED_POINT INF_SMALL=0.00003;
const FIXED_POINT INF_BIG= 4194304; //
const FIXED_POINT NEG_INF=-4194304;//




typedef struct band_score_info
{
   FIXED_POINT score;
   bool is_max_cell;
}
band_score_info;



typedef struct lp_info{
   FIXED_POINT lp_skip=lp_skip_const;
   FIXED_POINT lp_stay;
   FIXED_POINT lp_step;

}lp_info;



typedef struct model_info{
    //FIXED_POINT_MEAN mean;
    FIXED_POINT_STDV stdv;
    FIXED_POINT_STDV log_stdv;
    FIXED_POINT gp_mean;
    FIXED_POINT_STDV gp_stdv_inverse;
    //int read_id;
    //int id;
    //bool valid;
 } model_info;

 typedef struct event_info{
    FIXED_POINT_MEAN mean;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
 }event_info;

typedef struct scailing_info{
FIXED_POINT scalings_scale;
FIXED_POINT scalings_shift;
FIXED_POINT scalings_var;
FIXED_POINT scaling_log_var;
}scailing_info;

typedef struct band_aligner_info
{
int n_kmers;
int band_idx;
int n_events;
int band_lower_left_event_pre2;
int band_lower_left_kmer_pre2;
int band_lower_left_event_pre;
int band_lower_left_kmer_pre;
int band_lower_left_event_current;
int band_lower_left_kmer_current;
lp_info lps;
}band_aligner_info;


//const int MODEL_VEC_WIDTH=8;



typedef hls::vector<char,READ_VEC_WIDTH> READ_VEC;
typedef hls::vector<float,EVENT_VEC_WIDTH>EVENT_VEC;

typedef hls::vector<FIXED_TRACE,BAND_WIDTH> TRACE_VEC;

typedef ap_uint <256> ENCODE_AP_TRACE;
typedef hls::vector<ENCODE_AP_TRACE,MULTI_READ> AP_TRACE_VEC;

typedef hls::vector<int,MULTI_TOTAL*2> POS_POST_VEC;
typedef hls::vector<int,MULTI_TOTAL> LOWER_LEFT_VEC;



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

);


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
 );
//dataflow functions

void align_entry(int reset,int sequence_offset,int event_offset,
hls::stream<int>&sequence_offset_stream,
hls::stream<int>&event_offset_stream
);


void align_load_event_group(int event_t_length[MULTI_TOTAL],
hls::stream<int>&event_offset_stream_group,
hls::stream<int>event_offset_stream_vec[LOAD_DATAFLOW],
hls::stream<int>event_t_length_out[MULTI_TOTAL], EVENT_VEC* event_mean0,
hls::stream <EVENT_VEC>event_vec_stream[LOAD_DATAFLOW]
);



void align_load_read_group( 
    
int sequence_length[MULTI_TOTAL],
hls::stream<int>&sequence_offset_stream_group,
hls::stream<int>sequence_offset_stream_vec[LOAD_DATAFLOW],

hls::stream<int>sequence_length_out[MULTI_TOTAL], 
READ_VEC* read0,
hls::stream <READ_VEC>read_vec_stream_single[LOAD_DATAFLOW]

);


void align_event_vec_multi_v2(
hls::stream <int>&event_offset_stream_vec,
hls::stream <int>event_offset_stream_gen[HALF_READ],
hls::stream <EVENT_VEC>&event_vec_stream,
hls::stream <float>event_float_stream[HALF_READ]);


void align_read_vec_multi_v2(
hls::stream<int>&sequence_offset_stream_vec,
hls::stream<int>sequence_offset_stream_gen[HALF_READ],
hls::stream<READ_VEC>&read_vec_stream,
hls::stream<char>read_char_stream[HALF_READ]);

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
hls::stream<bool>&model_last_stream);



void align_event_mean_gen(
hls::stream <int>&event_offset_stream_gen,
hls::stream<int>&event_t_length_in,
hls::stream<int>&event_t_length_out,
hls::stream <float>&event_float_stream,
hls::stream <FIXED_POINT_MEAN>&event_stream,
hls::stream <FIXED_POINT_MEAN>&event_stream_post,
hls::stream <bool> &event_last_stream);



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

);




void align_max_v2(
hls::stream <int> &event_n_stream,
hls::stream <int> &band_length_single,
hls::stream <int> &band_length_single_out,
hls::stream <FIXED_POINT>&max_cell_stream,
hls::stream <int>& max_cell_event_idx_stream,
hls::stream <int> &max_event_idx_stream
);


void align_channel(
hls::stream <int> &band_length_single_out,

hls::stream <TRACE_VEC>&trace_vec_stream, 
hls::stream <int>&left_event_stream,

hls::stream<int> &max_event_idx_stream,

hls::stream <int> &max_event_idx_post_stream,
hls::stream <int> &max_band_post_out_stream,

hls::stream <TRACE_VEC>trace_vec_post_stream[2],
hls::stream <int>band_lower_left_event_post_stream[2],
AP_TRACE_VEC* trace_out_0,
LOWER_LEFT_VEC* left_out_0

);

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
);

void post_pos_out(
hls::stream<bool> event_last_stream[MULTI_TOTAL],
hls::stream<bool> model_last_stream[MULTI_TOTAL],
hls::stream<int> &max_band_num_stream,
hls::stream<int> ref_pos_post_stream[MULTI_TOTAL],
hls::stream<int> read_pos_post_stream[MULTI_TOTAL],
hls::stream<int> &n_pairs_post_stream,
POS_POST_VEC* aligned_ref_read_pos_0,

int* n_align
);
//align_utility functions

void kmer_generator(char read, char shift_register[KMER_SIZE], uint32_t *kmer_ranks);
void kmer_generator_post(char read, char shift_register[KMER_SIZE], uint32_t *kmer_ranks);
uint32_t get_kmer_rank(char *str, uint32_t k);
void ap_trace_encode(TRACE_VEC input_vec[MULTI_READ],AP_TRACE_VEC* output_vec);
void ap_trace_decode(AP_TRACE_VEC input_vec, TRACE_VEC output_vec[MULTI_READ]);

void align_kernel_score_calculate(FIXED_POINT score_d, FIXED_POINT score_u, FIXED_POINT score_l, FIXED_POINT *max_score_out, FIXED_TRACE *from_score_out);
void event_shift(event_info events[BAND_WIDTH+1],event_info events_local[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], model_info model_levels_local[BAND_WIDTH], event_info new_event);
void model_shift(event_info events[BAND_WIDTH+1],event_info events_local[BAND_WIDTH+1], model_info model_levels_cur[BAND_WIDTH], model_info model_levels_local[BAND_WIDTH], model_info new_model);
void compare_band(int n_event,band_score_info band_a,band_score_info band_b, band_score_info* band_result,int event_idx_a,int event_idx_b,int* event_idx_result);
void compare_band_fix(int n_event,FIXED_POINT band_a,FIXED_POINT band_b, FIXED_POINT* band_result,int event_idx_a,int event_idx_b,int* event_idx_result);








template <int MULTI>
void align_stream_event_mean_gen(hls::stream <int>event_offset_stream_gen[MULTI_TOTAL],hls::stream<int>event_t_length_in[MULTI_TOTAL],hls::stream<int>event_t_length_out[MULTI_TOTAL],hls::stream <float>event_float_stream[MULTI_TOTAL],hls::stream <FIXED_POINT_MEAN>event_stream[MULTI_TOTAL],hls::stream <FIXED_POINT_MEAN>event_stream_post[MULTI_TOTAL],hls::stream<bool>event_last_stream[MULTI_TOTAL])
{
    #pragma HLS INLINE

   for (int i = 0; i < MULTI; i++)
   {
      #pragma HLS UNROLL
      align_event_mean_gen(event_offset_stream_gen[i],event_t_length_in[i],event_t_length_out[i],event_float_stream[i],event_stream[i],event_stream_post[i],event_last_stream[i]);
   }

}


template <int MULTI>
void align_stream_model_gen(
hls::stream<int>sequence_offset_stream_gen[MULTI_TOTAL],
hls::stream<int>sequence_length_in[MULTI_TOTAL],
hls::stream<int>kmer_length_out[MULTI_TOTAL],
hls::stream<char>read_char_stream[MULTI_TOTAL], 
FIXED_POINT_MEAN model_levels_mean[MULTI_TOTAL][MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_stdv[MULTI_TOTAL][MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_log_stdv[MULTI_TOTAL][MAX_MODEL_NUM],
FIXED_POINT_STDV model_levels_gp_stdv_inverse[MULTI_TOTAL][MAX_MODEL_NUM],

hls::stream<model_info>model_levels_cur_stream[MULTI_TOTAL],
hls::stream<model_info>model_levels_cur_post_stream[MULTI_TOTAL],

scailing_info scailings[MULTI_TOTAL],
hls::stream<bool>model_last_stream[MULTI_TOTAL]
)
{

    #pragma HLS INLINE

   for (int i = 0; i < MULTI; i++)
   {
      #pragma HLS UNROLL
      align_model_gen(sequence_offset_stream_gen[i],sequence_length_in[i],kmer_length_out[i],read_char_stream[i],model_levels_mean[i],model_levels_stdv[i],model_levels_log_stdv[i],model_levels_gp_stdv_inverse[i],model_levels_cur_stream[i],model_levels_cur_post_stream[i],scailings[i],model_last_stream[i]);
   }
}



template <int MULTI>
void align_stream_event_vec_load_v2(
hls::stream <int>event_offset_stream_vec[LOAD_DATAFLOW],
hls::stream <int>event_offset_stream_gen[MULTI_TOTAL],
hls::stream <EVENT_VEC>event_vec_stream[LOAD_DATAFLOW],
hls::stream <float>event_float_stream[MULTI_TOTAL]
)
{
#pragma HLS INLINE
for (int i = 0; i < MULTI; i++)
{
   #pragma HLS UNROLL
   align_event_vec_multi_v2(event_offset_stream_vec[i],&event_offset_stream_gen[i*HALF_READ],event_vec_stream[i],&event_float_stream[i*HALF_READ]);
}

}



template <int MULTI>
void align_stream_read_vec_load_v2(
hls::stream<int>sequence_offset_stream_vec[LOAD_DATAFLOW],
hls::stream<int>sequence_offset_stream_gen[MULTI_TOTAL],
hls::stream<READ_VEC>read_vec_stream[LOAD_DATAFLOW],
hls::stream<char>read_char_stream[MULTI_TOTAL] 
)

{
   #pragma HLS INLINE
   for (int i = 0; i < MULTI; i++)
   {
      #pragma HLS UNROLL
      align_read_vec_multi_v2(sequence_offset_stream_vec[i],&sequence_offset_stream_gen[i*HALF_READ],read_vec_stream[i],&read_char_stream[i*HALF_READ]);

   }
   
}
