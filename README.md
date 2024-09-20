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
* `set_CRS.F`: code used for generating CRS-based matrix
* `Makefile`: Makefile used for building program

## Compile
Code can be compiled by typing below in the source directory.
```
make clobber
make
```
Five executable files will be generated as below
```
./kernelCRS1_CPU_OPENMP.exe
./kernelCRS1_GPU_OPENACC.exe
./kernelEBE1_GPU_OPENACC.exe
./kernelEBE4_GPU_OPENACC.exe
./kernelEBE4_GPU_CUDA.exe
```

## Test data
Test data can be downloaded from https://www.eri.u-tokyo.ac.jp/cshpc/share/WACCPD2024/data.tar.gz (236 MB) and extracted as `tar xzf data.tar.gz`.
Contents of `./data` is as below.
* `./data/setting.dat`: setting data of finite-element model (number of nodes, number of elements, and number of materials).
* `./data/material.dat`: material properties of finite-element model.
* `./data/conn.dat`: connectivity data of finite-element model.
* `./data/coor.dat`: coordinate data of finite-element model.
* `./data/num.dat`: material number of each element of finite-element model.

## Run
Kernels programs can be run as below.
```
export OMP_NUM_THREADS=72
./kernelCRS1_CPU_OPENMP.exe
./kernelCRS1_GPU_OPENACC.exe
./kernelEBE1_GPU_OPENACC.exe
./kernelEBE4_GPU_OPENACC.exe
./kernelEBE4_GPU_CUDA.exe
```

Results measured on a Single-GH200 node with one 72-core ARMv9-a Grace CPU@3.1 GHz, 480 GB LPDDR5X memory, and one H100 96 GB GPU, with nvhpc/24.1, NVIDIA Driver Version 535.104.05, CUDA Version 12.2 on Ubuntu 22.04.4 LTS is shown below:

==== ./kernelCRS1_CPU_OPENMP.exe ====
```
 using CRS with NVEC            1
 not using OpenACC
 matvec took   0.1655118465423584
...
 matvec took   0.1645679473876953
 norm            0    6781.519443194073
 norm            1    314711506799008.2
```

==== ./kernelCRS1_GPU_OPENACC.exe ====
```
 using CRS with NVEC            1
 using OpenACC
 matvec took   1.6763925552368164E-002
...
 matvec took   1.6685962677001953E-002
 norm            0    6781.519443193030
 norm            1    314711506799008.4
```

==== ./kernelEBE1_GPU_OPENACC.exe ====
```
 using EBE with NVEC            1
 using OpenACC
 matvec took   4.4920444488525391E-003
 ...
 matvec took   4.4050216674804688E-003
 norm            0    6781.519443193030
 norm            1    314711506799008.4
```

==== ./kernelEBE4_GPU_OPENACC.exe ====
```
 using EBE with NVEC            4
 using OpenACC
 matvec took   9.6509456634521484E-003
...
 matvec took   9.5551013946533203E-003
 norm            0    6781.519443193030         6781.519443193030
    6781.519443193030         6781.519443193030
 norm            1    314711506799008.4         314711506799008.4
    314711506799008.4         314711506799008.4
```

==== ./kernelEBE4_GPU_CUDA.exe ====
```
 using EBE with NVEC            4
 using OpenACC
 matvec took   1.0213851928710938E-002
...
 matvec took   1.0118007659912109E-002
 norm            0    6781.519443193030         6781.519443193030
    6781.519443193030         6781.519443193030
 norm            1    314711506259667.2         314711506259667.2
    314711506259667.2         314711506259667.2
```

## Validation of results

The correctness of results can be checked by comparing the `norm` values.
In case of 4-case runs (`./kernelEBE4_GPU_OPENACC.exe` and `./kernelEBE4_GPU_CUDA.exe`) the four values of the norm should also be identical.

## Interpretation of results

Elapsed time for the matrix vector product is given by the number `matvec took` in seconds. 10 iterations are conducted to check the fluctuation in elapsed time.

In case of 4-case runs (`./kernelEBE4_GPU_OPENACC.exe` and `./kernelEBE4_GPU_CUDA.exe`) the elapsed time is total for the 4 cases.
Thus, the elapsed time in Table II of the paper (per case) is obtained by dividing elapsed time by four.

## Publication

Tsuyoshi Ichimura, Kohei Fujita, Muneo Hori, Lalith Maddegedara, Jack Wells, Alan Gray, Ian Karlin, John Linford, "Heterogeneous computing in a strongly-connected CPU-GPU environment: fast multiple time-evolution equation-based modeling accelerated using data-driven approach", Eleventh Workshop on Accelerator Programming and Directives (WACCPD 2024) (Accepted).

## License

directive-based-spmv-kernels, version 1.0.0 (c) 2024 Tsuyoshi Ichimura et al. directive-based-spmv-kernels is freely distributable under the terms of an MIT-style license.
