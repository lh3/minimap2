GPU				?= 		AMD
GPUARCH 		?=		sm_86
CONFIG			+= $(if $(MAX_MICRO_BATCH),-DMICRO_BATCH=\($(MAX_MICRO_BATCH)\))

###################################################
############  	CPU Compile 	###################
###################################################
CU_SRC			= $(wildcard gpu/*.cu)
CU_OBJS			= $(CU_SRC:%.cu=%.o)
CU_PTX			= $(CU_SRC:%.cu=%.ptx)
C_SRC			= $(wildcard gpu/*.c)
OBJS			+= $(C_SRC:%.c=%.o)
INCLUDES		+= -I gpu

###################################################
############  	CUDA Compile 	###################
###################################################
COMPUTE_ARCH    = $(GPUARCH:sm_%=compute_%)
NVCC 			= nvcc
CUDAFLAGS		= -rdc=true -gencode arch=$(COMPUTE_ARCH),code=$(GPUARCH) -diag-suppress=177 -diag-suppress=1650 # supress unused variable / func warning
CUDANALYZEFLAG	= -Xptxas -v 
CUDATESTFLAG	= -G 

###################################################
############	HIP Compile		###################
###################################################
HIPCC			= hipcc
HIPFLAGS		= -DUSEHIP 
HIPANALYZEFLAG  = -Rpass-analysis=kernel-resource-usage
HIPTESTFLAGS	= -G -ggdb
HIPLIBS			= -L${ROCM_PATH}/lib -lroctx64 -lroctracer64

###################################################
############	DEBUG Options	###################
###################################################
ifeq ($(GPU), AMD)
	GPU_CC 		= $(HIPCC)
	GPU_FLAGS	= $(HIPFLAGS)
	GPU_TESTFL	= $(HIPTESTFLAGS)
	GPU_ANALYZE	= $(HIPANALYZEFLAG)
	LIBS		+= $(HIPLIBS)
else
	GPU_CC 		= $(NVCC)
	GPU_FLAGS	= $(CUDAFLAGS)
	GPU_ANALYZE = $(CUDANALYZEFLAG)
	GPU_TESTFL	= $(CUDATESTFLAG)
endif

ifeq ($(DEBUG),analyze)
	GPU_FLAGS	+= $(GPU_ANALYZE)
endif
ifeq ($(DEBUG),verbose)
	GPU_FLAGS	+= $(GPU_TESTFL)
endif


%.o: %.cu
	$(GPU_CC) -c $(GPU_FLAGS) $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(CONFIG) $< -o $@

%.ptx: %.cu
	$(GPU_CC) -ptx -src-in-ptx $(GPU_FLAGS) $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(CONFIG) $< -o $@

%.as: %.o
	cuobjdump -all $< > $@

cleangpu: 
	rm -f $(CU_OBJS) $(CU_PTX)

# profile:CFLAGS += -pg -g3
# profile:all
# 	perf record --call-graph=dwarf -e cycles:u time ./minimap2 -a test/MT-human.fa test/MT-orang.fa > test.sam

cudep: gpu/.depend

gpu/.depend: $(CU_SRC)
	rm -f gpu/.depend
	$(GPU_CC) -c $(GPU_FLAGS) $(CFLAGS)  $(CPPFLAGS) $(INCLUDES) -MM $^ > $@

include gpu/.depend