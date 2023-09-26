GPU				?= 		AMD
CONFIG 			= $(if $(GPU_CONFIG),-DGPU_CONFIG='"$(GPU_CONFIG)"')
CONFIG			+= $(if $(LONG_BLOCK_SIZE),-D__LONG_BLOCK_SIZE__=\($(LONG_BLOCK_SIZE)\))
CONFIG			+= $(if $(MID_BLOCK_SIZE),-D__MID_BLOCK_SIZE__=\($(MID_BLOCK_SIZE)\))
CONFIG			+= $(if $(SHORT_BLOCK_SIZE),-D__SHORT_BLOCK_SIZE__=\($(SHORT_BLOCK_SIZE)\))
CONFIG			+= $(if $(MID_CUT),-DMM_MID_SEG_CUTOFF=\($(MID_CUT)\))
CONFIG			+= $(if $(LONG_CUT),-DMM_LONG_SEG_CUTOFF=\($(LONG_CUT)\))

###################################################
############  	CPU Compile 	###################
###################################################
CU_SRC			= $(wildcard gpu/*.cu)
CU_OBJS			= $(CU_SRC:%.cu=%.o)
INCLUDES		+= -I gpu

###################################################
############  	CUDA Compile 	###################
###################################################
NVCC 			= nvcc
CUDAFLAGS		= -rdc=true -lineinfo
CUDATESTFLAG	= -G -DDEBUG_PRINT

###################################################
############	HIP Compile		###################
###################################################
HIPCC			= hipcc
HIPFLAGS		= -DUSEHIP -Rpass-analysis=kernel-resource-usage
HIPTESTFLAGS	= -g -DDEBUG_PRINT

###################################################
############	DEBUG Options	###################
###################################################
ifeq ($(GPU), AMD)
	GPU_CC 		= $(HIPCC)
	GPU_FLAGS	= $(HIPFLAGS)
	GPU_TESTFL	= $(HIPTESTFLAGS)
else
	GPU_CC 		= $(NVCC)
	GPU_FLAGS	= $(CUDAFLAGS)
	GPU_TESTFL	= $(CUDATESTFLAG)
endif

ifneq ($(DEBUG),)
	GPU_FLAGS	+= $(GPU_TESTFL)
endif

ifneq ($(DEBUG_ANALYSIS),)
	GPU_FLAGS	+= $(GPU_TESTFL)
	GPU_FLAGS	+= -DDEBUG_CHECK -DDEBUG_VERBOSE
endif 


%.o: %.cu
	$(GPU_CC) -c $(GPU_FLAGS) $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(CONFIG) $< -o $@

cleangpu: 
	rm -f $(CU_OBJS)

# profile:CFLAGS += -pg -g3
# profile:all
# 	perf record --call-graph=dwarf -e cycles:u time ./minimap2 -a test/MT-human.fa test/MT-orang.fa > test.sam

cudep: gpu/.depend

gpu/.depend: $(CU_SRC)
	rm -f gpu/.depend
	$(GPU_CC) -c $(GPU_FLAGS) $(CFLAGS)  $(CPPFLAGS) $(INCLUDES) -MM $^ > $@

include gpu/.depend