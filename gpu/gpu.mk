GPU				?= 		AMD
CONFIG 			= $(if $(GPU_CONFIG),-DGPU_CONFIG='"$(GPU_CONFIG)"')
CONFIG			+= $(if $(LONG_BLOCK_SIZE),-D__LONG_BLOCK_SIZE__=\($(LONG_BLOCK_SIZE)\))
CONFIG			+= $(if $(MID_BLOCK_SIZE),-D__MID_BLOCK_SIZE__=\($(MID_BLOCK_SIZE)\))
CONFIG			+= $(if $(SHORT_BLOCK_SIZE),-D__SHORT_BLOCK_SIZE__=\($(SHORT_BLOCK_SIZE)\))
CONFIG			+= $(if $(MID_CUT),-DMM_MID_SEG_CUTOFF=\($(MID_CUT)\))
CONFIG			+= $(if $(LONG_CUT),-DMM_LONG_SEG_CUTOFF=\($(LONG_CUT)\))
CONFIG			+= $(if $(MICRO_BATCH),-DMICRO_BATCH=\($(MICRO_BATCH)\))

###################################################
############  	CPU Compile 	###################
###################################################
CU_SRC			= $(wildcard gpu/*.cu)
CU_OBJS			= $(CU_SRC:%.cu=%.o)
C_SRC			= $(wildcard gpu/*.c)
OBJS			+= $(C_SRC:%.c=%.o)
INCLUDES		+= -I gpu

###################################################
############  	CUDA Compile 	###################
###################################################
NVCC 			= nvcc
CUDAFLAGS		= -rdc=true -lineinfo
CUDATESTFLAG	= -G

###################################################
############	HIP Compile		###################
###################################################
HIPCC			= hipcc
HIPFLAGS		= -DUSEHIP 
HIPTESTFLAGS	= -G -Rpass-analysis=kernel-resource-usage -ggdb
HIPLIBS			= -L${ROCM_PATH}/lib -lroctx64 -lroctracer64

###################################################
############	DEBUG Options	###################
###################################################
ifeq ($(GPU), AMD)
	GPU_CC 		= $(HIPCC)
	GPU_FLAGS	= $(HIPFLAGS)
	GPU_TESTFL	= $(HIPTESTFLAGS)
	LIBS		+= $(HIPLIBS)
else
	GPU_CC 		= $(NVCC)
	GPU_FLAGS	= $(CUDAFLAGS)
	GPU_TESTFL	= $(CUDATESTFLAG)
endif

ifeq ($(DEBUG),analyze)
	GPU_FLAGS	+= $(GPU_TESTFL)
endif

%.o: %.cu
	$(GPU_CC) -c $(GPU_FLAGS) $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(CONFIG) $< -o $@

cleangpu: 
	rm -f gpu/*.o

# profile:CFLAGS += -pg -g3
# profile:all
# 	perf record --call-graph=dwarf -e cycles:u time ./minimap2 -a test/MT-human.fa test/MT-orang.fa > test.sam

cudep: gpu/.depend

gpu/.depend: $(CU_SRC)
	rm -f gpu/.depend
	$(GPU_CC) -c $(GPU_FLAGS) $(CFLAGS)  $(CPPFLAGS) $(INCLUDES) -MM $^ > $@

include gpu/.depend