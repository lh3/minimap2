GPU				?= NVCC
CONFIG 			= $(if $(GPU_CONFIG),-DGPU_CONFIG='"$(GPU_CONFIG)"')

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
CUDAFLAGS		= -rdc=true -DNDEBUG -lineinfo ## turn off assert
CUDATESTFLAG	= -G

###################################################
############	HIP Compile		###################
###################################################
HIPCC			= hipcc
HIPFLAGS		= -DUSEHIP -DNDEBUG ## turn off assert
HIPTESTFLAGS	= -g

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

# check: CFLAGS += -DDEBUG_CHECK
# check: HIPFLAGS += -DDEBUG_CHECK
# check: CUDAFLAGS += -DDEBUG_CHECK
# check: HIPFLAGS += $(HIPTESTFLAGS)
# check: CUDAFLAGS += $(CUDATESTFLAG)

# verbose: CFLAGS += -DDEBUG_VERBOSE
# verbose: HIPFLAGS += -DDEBUG_VERBOSE
# verbose: CUDAFLAGS += -DDEBUG_VERBOSE



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
