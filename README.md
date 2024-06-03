# FPGA accelerator for ABEA using HLS

Adaptive Banded Event Alignment (ABEA) stands as a crucial algorithmic component in sequence polishing and DNA methylation detection, employing a dynamic programming strategy to align raw Nanopore signal to reference reads. ABEA was first introduced by [Nanopolish](https://github.com/jts/nanopolish) and it is the bottleneck for both nanopolish call-methylation module and nanopolish eventalign module. The proposed accelerator with Xilinx Vivado High-Level Synthesis (HLS) targets Xilinx UV9P FPGA board and was implemented based on GPU accleration [f5c](https://github.com/hasindu2008/f5c). The proposed accelerator obtains an average throughput speedup of 10.9 $\times$ over the CPU-only implementation, and an average 1.82 $\times$ speedup over the-state-of-art GPU acceleration with 3.8\% of the energy.

 ## Dataset

 'Nanopore WGS Consortium' sequencing dataset was used to conduct our evaluation, which can downloaded from [Whole Human Genome Sequencing Project](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md). Specifically, we chose five individual flowcells FAB39088, FAB42395, FAB41174, FAF05869 and FAF18554.
 
 GRCh38/hg38 [hg38 FASTA Human Reference Genome](https://drive.google.com/file/d/1Ur3xybIzQGSxuqeByyp5OMrpaRJXCsMI/view?usp=sharing) was used as the reference genomes for alignment.
 
 ## Usage
  
In this work, we target Xilinx VU9P offered by Amazon AWS F1 instance. The FPGA configuration file is in the binary directory. The host calls two kinds of kernels implemented in the configuration file. 

- `align_top_stream` for pipelined inter-read alignments; 2 instances are deployed named as _align_top_stream_1_ and _align_top_stream_2_.
- `align_top_stream_no_group` for ultra-long read alignments; 1 instance is deployed named as _align_top_stream_no_group_1_.

[Xilinx Runtime Library (XRT) API](https://docs.amd.com/r/en-US/ug1393-vitis-application-acceleration/Getting-Started-with-Vitis) is invoked to handle the data transfer between Host DRAM to FPGA DRAM and the launch of specific kernels.

The interface of the two kinds of kernel is as follows. Generally, one initialization launching is needed for all kernels; firstly. After that, regular reads are collected into buckets for being processed in pipelined fashion by launching   sucessive _align_top_stream_ kernels. Each ultra-long read is processed one by one by sucessive launching  _align_top_stream_no_group_ kernel.

- **reset** : flag to identify initialization or performing alignments; if launching kernel for initialzion, flag is set as 1. if launching kernel for alignments, flag is set as 0.
- **event_mean** : mean signal values of events
- **read** : bases of reads
- **n_align** : the output of the number of (events,k-mer) pairs
- **aligned_ref_read_pos** : the output of (events,k-mer) pairs
- **model_data** : pore-model table
- **trace_out** : allocated memory for trace table; Assigned values are not needed.
- **length_int** : number of events and size of read bases for each read; 0s are padded for unfulled bucket.
- **scailings** : scaling parameters for the signal
- **lps** : skip penalty
- **max_event_t_length** : maximum number of events in this bucket
- **max_sequence_length** : maximum size of read bases in this bucket
- **batch_num** : number of batch in this bucket; one batch of reads perform alignment in full pipeline and batch size is 8 reads in our implementaion. _batch_num_ = ceil(_bucket_size_/_batch_size_) for _align_top_stream_ kernel. _batch_num_ = 1 for  _align_top_stream_no_group_ kernel.

This is a simple example for launching kernel to perform inter-read alignments

```
kernel_0 = xrt::kernel(core->device, uuid, "align_top_stream:{align_top_stream_1}");
auto run_0 = core->kernel_0(1 ,NULL,NULL,NULL,NULL, core->model_data_bo_0,NULL, NULL, NULL, NULL, NULL, NULL, NULL, 1);
run_0.wait();
# transfer input from FPGA DRAM to Host DRAM
run_0 = kernel_0(0, event_mean_bo_0,read_bo_0,n_align_bo_0,aligned_pairs_bo_0,NULL,trace_table_bo_0, left_out_bo_0, length_int_bo_0,scaling_info_bo_0,lp_info_bo_0,event_offset,reads_offset,batch_num_0);
run_0.wait();
# transfer output from FPGA DRAM to Host DRAM
```
