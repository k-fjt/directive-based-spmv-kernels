FC = nvfortran
FFLAGCPU = -O3 -mp=multicore -Minfo=all -llapack
FFLAGGPU = -O3 -acc=gpu -gpu=cc90,cuda12.3,fastmath,loadcache:L1,lineinfo,deepcopy -Minfo=all -mp=multicore

# options
FLAG1 = -DNVEC=1 -DCRS10
FLAG2 = -DNVEC=1 -DCRS10
FLAG3 = -DNVEC=1
FLAG4 = -DNVEC=4
FLAG5 = -DNVEC=4 -DUSECUDA

PROGRAM1 = kernelCRS1_CPU_OPENMP.exe
PROGRAM2 = kernelCRS1_GPU_OPENACC.exe
PROGRAM3 = kernelEBE1_GPU_OPENACC.exe
PROGRAM4 = kernelEBE4_GPU_OPENACC.exe
PROGRAM5 = kernelEBE4_GPU_CUDA.exe

SRCS = \
	main.F90 \
	kernelCRS.F \
	kernelEBE.F 

OBJS = $(SRCS:.F=.o)

.SUFFIXES: .o .F .F90

all: $(PROGRAM1) $(PROGRAM2) $(PROGRAM3) $(PROGRAM4) $(PROGRAM5)

$(PROGRAM1):
	make clean
	$(FC) $(FLAG1) $(FFLAGCPU) -c $(SRCS)
	$(FC) $(FLAG1) $(FFLAGCPU) $(OBJS) -o $@

$(PROGRAM2):
	make clean
	$(FC) $(FLAG2) $(FFLAGGPU) -c $(SRCS)
	$(FC) $(FLAG2) $(FFLAGGPU) $(OBJS) -o $@

$(PROGRAM3):
	make clean
	$(FC) $(FLAG3) $(FFLAGGPU) -c $(SRCS)
	$(FC) $(FLAG3) $(FFLAGGPU) $(OBJS) -o $@

$(PROGRAM4):
	make clean
	$(FC) $(FLAG4) $(FFLAGGPU) -c $(SRCS)
	$(FC) $(FLAG4) $(FFLAGGPU) $(OBJS) -o $@

$(PROGRAM5): 
	make clean
	$(FC) $(FLAG5) $(FFLAGGPU) -c $(SRCS)
	nvcc $(FLAG5) --resource-usage --generate-code arch=compute_90,code=sm_90 -O3 -c kernelEBE_CUDA.cu
	$(FC) $(FLAG5) $(FFLAGGPU) $(OBJS) kernelEBE_CUDA.o -cuda -o $@

###################################################################################################
clean:
	rm -f *.o *.out *.lst *.s
clobber:
	rm -f *.o *.out *.lst *.s *.exe

