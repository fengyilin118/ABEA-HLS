# FPGA accelerator for ABEA using High-Level Synthesis

TThe Adaptive Banded Event Alignment (ABEA) is a crucial  component in sequence polishing and DNA methylation detection, utilizing a dynamic programming strategy to align raw Nanopore signals to a biological reference sequence to correct sequencing errors and detect the modified nucleotides. Originally introduced by [Nanopolish](https://github.com/jts/nanopolish), ABEA serves as the primary bottleneck in both the Nanopolish call-methylation and eventalign modules. The proposed accelerator, implemented using Xilinx Vivado High-Level Synthesis (HLS), targets the Xilinx VU9P FPGA board and is based on the GPU-accelerated [f5c](https://github.com/hasindu2008/f5c) framework. When tested on an Amazon AWS F1 instance, this FPGA-based accelerator achieves an average throughput speedup of 10.05X over the CPU-only implementation and an average 1.81X speedup over the state-of-the-art GPU acceleration on the NVIDIA V100, while consuming only 7.2% of the energy.

 ## Dataset

- 'Nanopore WGS Consortium' sequencing dataset was used to conduct our evaluation, which can downloaded from [Whole Human Genome Sequencing Project](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md). Specifically, we chose five individual flowcells FAB39088, FAF05869, FAF18554, FAB41174 and FAB42395.
 
- GRCh38/hg38 [hg38 FASTA Human Reference Genome](https://drive.google.com/file/d/1Ur3xybIzQGSxuqeByyp5OMrpaRJXCsMI/view?usp=sharing) was used as the reference genomes for alignment.
 
 ## Code and FPGA BIN file
 Our contribution mainly attributes to designing two type of Compute Unites on FPGA: Pipeline CU for regular reads and Ultra Long CU for ultra-long reads, as well as bucketing scheduling for load balance.
 - **The FPGA BIN file is in the [binary directory](https://github.com/fengyilin118/ABEA-HLS/tree/main/binary)**. The details of HLS code will be updated later.
 - Bucket scheduling is implemented as [*align_db_FPGA_by_len*](https://github.com/fengyilin118/ABEA-HLS/blob/main/src/f5c.c) function.

## Kernel

The host calls two kinds of kernels (Pipeline CU and Ultra Long CU). 

- `align_top_stream` for pipelined inter-read alignments; 2 instances are deployed named as _align_top_stream_1_ and _align_top_stream_2_.
- `align_top_stream_no_group` for ultra-long read alignments; 1 instance is deployed named as _align_top_stream_no_group_1_.

[Xilinx Runtime Library (XRT) API](https://docs.amd.com/r/en-US/ug1393-vitis-application-acceleration/Getting-Started-with-Vitis) is invoked to handle the data transfer between Host DRAM to FPGA DRAM and the launch of specific kernels.

The interface of the two kinds of kernel is as follows. Generally, one initialization launching is needed for all kernels; firstly. After that, regular reads are collected into buckets for being processed in pipelined fashion by launching   sucessive _align_top_stream_ kernels. Each ultra-long read is processed one by one by sucessive launching  _align_top_stream_no_group_ kernel.

- **reset** : flag to identify initialization or performing alignments; if launching kernel for initialzion, flag is set as 1. if launching kernel for alignments, flag is set as 0.
- **event_mean** : mean signal values of events
- **read** : bases of reads
- **n_align** : the output of the number of (events,k-mer) pairs
- **aligned_ref_read_pos** : the output of (events,k-mer) pairs
- **model_data** : pore-model table; model_table are transferred to FPGA DRAM in initialization launch.
- **trace_out** : allocated memory for trace table; Assigned values are not needed.
- **length_int** : number of events and size of read bases for each read; 0s are padded for unfulled bucket.
- **scailings** : scaling parameters for the signal
- **lps** : skip penalty
- **max_event_t_length** : maximum number of events in this bucket
- **max_sequence_length** : maximum size of read bases in this bucket
- **batch_num** : number of batch in this bucket; one batch of reads perform alignment in full pipeline and batch size is 8 reads in our implementaion. _batch_num_ = ceil(_bucket_size_/_batch_size_) for _align_top_stream_ kernel. _batch_num_ = 1 for  _align_top_stream_no_group_ kernel.

This is a simple example for launching kernel to perform inter-read alignments

```
# initialization
kernel_0 = xrt::kernel(core->device, uuid, "align_top_stream:{align_top_stream_1}");
auto run_0 = core->kernel_0(1 ,NULL,NULL,NULL,NULL, core->model_data_bo_0,NULL, NULL, NULL, NULL, NULL, NULL, NULL, 1);
run_0.wait();
# transfer input from FPGA DRAM to Host DRAM
run_0 = kernel_0(0, event_mean_bo_0,read_bo_0,n_align_bo_0,aligned_pairs_bo_0,NULL,trace_table_bo_0, left_out_bo_0, length_int_bo_0,scaling_info_bo_0,lp_info_bo_0,event_offset,reads_offset,batch_num_0);
run_0.wait();
# transfer output from FPGA DRAM to Host DRAM
```
