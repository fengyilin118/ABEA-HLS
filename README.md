# FPGA accelerator for ABEA using HLS

Adaptive Banded Event Alignment (ABEA) stands as a crucial algorithmic component in sequence polishing and DNA methylation detection, employing a dynamic programming strategy to align raw Nanopore signal to reference reads. ABEA was first introduced by [Nanopolish](https://github.com/jts/nanopolish) and it is the bottleneck for both nanopolish call-methylation module and nanopolish eventalign module. The proposed accelerator with Xilinx Vivado High-Level Synthesis (HLS) targets Xilinx UV9P FPGA board and was implemented based on GPU accleration [f5c](https://github.com/hasindu2008/f5c). The proposed accelerator obtains an average throughput speedup of 10.9 $\times$ over the CPU-only implementation, and an average 1.82 $\times$ speedup over the-state-of-art GPU acceleration with 3.8\% of the energy.

 ## Dataset

 `Nanopore WGS Consortium' sequencing dataset was used to conduct our evaluation, which can downloaded from [Whole Human Genome Sequencing Project](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md). Specifically, we chose five individual flowcells FAB39088, FAB42395, FAB41174, FAF05869 and FAF18554.
 
 GRCh38/hg38 [hg38 FASTA Human Reference Genome](https://drive.google.com/file/d/1Ur3xybIzQGSxuqeByyp5OMrpaRJXCsMI/view?usp=sharing) was used as the reference genomes for alignment.
 
 ## Usage
  
In this work, we target Xilinx VU9P offered by Amazon AWS F1 instance. The FPGA configuration file is in the binary directory. The host calls two kinds of kernels implemented in the configuration file. 

- `align_top_stream` for pipelined inter-read alignments; 2 instances are deployed named as _align_top_stream_1_ and _align_top_stream_2_.
- `align_top_stream_no_group` for ultra-long read alignments; 1 instance is deployed named as _align_top_stream_no_group_1_.

[Xilinx Runtime Library (XRT) API](https://docs.amd.com/r/en-US/ug1393-vitis-application-acceleration/Getting-Started-with-Vitis) is invoked to handle the data transfer between Host DRAM to FPGA DRAM and the launch of specific kernels.

The interface of the two kinds of kernel is as follows.
- **reset** :
- **event_mean** :
- **read** :
- **n_align** :
- **aligned_ref_read_pos** :
- **model_data** :
- **trace_out** :
- **length_int** :
- **scailings** :
- **lps** :
- **max_event_t_length** :
- **max_sequence_length** :
- **batch_num** :
