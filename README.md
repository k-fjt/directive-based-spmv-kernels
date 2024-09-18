# directive-based-spmv-kernels

This repository explains the sparse matrix-vector product kernels in the Artifact Description (AD)/ Artifact Evaluation (AE) appendix for paper "Heterogeneous computing in a strongly-connected CPU-GPU environment: fast multiple time-evolution equation-based modeling accelerated using data-driven approach" accepted for Eleventh Workshop on Accelerator Programming and Directives (WACCPD 2024).

## Description

In this paper, we proposed a CPU-GPU heterogeneous computing method for solving time-evolution partial differential equation problems many times with guaranteed accuracy, in short time-to-solution and low energy-tosolution. Here, a data-driven approach, which leverages large CPU memory capacity, and equation-based modeling, which leverages fast GPU computing, were used in a strongly connected CPU–GPU environment. The proposed method was realized using only directive-based parallel programming models, indicating that directives are highly effective in analyses in heterogeneous computing environments.

The findings of this paper were attained by a software implementation and evaluation on computer environments. The matrix-vector product kernels in this repository is used for the evaluation of Table II in the paper, which are the most computationally costly kernels in the baseline and proposed methods. Kernels not included in Table II are memory bandwidth bound and take less time than the matrix-vector product kernel in both the baseline and proposed methods; thus, the performance trends of the entire application programs can be estimated by the kernels in this repository.

Below we explain the details of the sparse matrix-vector product kernels.

## Requirement
* GH200 node with nvhpc complier

## Source codes
* `main.F90`: main code
* `kernelCRS.F`: CRS-based kernel code in OpenMP/OpenACC
* `kernelEBE.F`: EBE-based kernel code in OpenACC
* `kernelEBE_CUDA.cu`: EBE-based kernel code in CUDA
* `Makefile`: Makefile used for building program

## Compile
Code can be compiled by typing `make` in the source directory.

## Test data
Data is given in `./data`.
* `./data/setting.dat`: setting data of finite-element mesh (number of nodes and elements).
* `./data/conn.dat`: connectivity data of finite-element mesh.
* `./data/coor.dat`: coordinate data of finite-element mesh.

## Run
Kernels programs can be run as below.
```
export OMP_NUM_THREADS=72
./kernelCRS1_CPU_OPENMP
./kernelCRS1_GPU_OPENACC
./kernelEBE1_GPU_OPENACC
./kernelEBE4_GPU_OPENACC
./kernelEBE4_GPU_CUDA
```

Results measured on a Single-GH200 node with one 72-core ARMv9-a Grace CPU@3.1 GHz, 480 GB LPDDR5X memory, and one H100 96 GB GPU, with nvhpc/24.1, NVIDIA Driver Version 535.104.05, CUDA Version 12.2 on Ubuntu 22.04.4 LTS is shown below:
```
results...
```

## Publication

Tsuyoshi Ichimura, Kohei Fujita, Muneo Hori, Lalith Maddegedara, Jack Wells, Alan Gray, Ian Karlin, John Linford, "Heterogeneous computing in a strongly-connected CPU-GPU environment: fast multiple time-evolution equation-based modeling accelerated using data-driven approach", Eleventh Workshop on Accelerator Programming and Directives (WACCPD 2024) (Accepted).

## License

directive-based-spmv-kernels, version 1.0.0 (c) 2024 Tsuyoshi Ichimura et al. directive-based-spmv-kernels is freely distributable under the terms of an MIT-style license.
