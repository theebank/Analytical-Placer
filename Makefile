default: libs run

all: libs run hb

all32: libs run hb fortran

all64: libs run hb fortran64

include SuiteSparse-5.10.1/SuiteSparse_config/SuiteSparse_config.mk

#-------------------------------------------------------------------------------
# UMFPACK optionally uses the CHOLMOD Partition module
LIB_WITH_CHOLMOD =
ifeq (,$(findstring -DNCHOLMOD, $(UMFPACK_CONFIG)))
    LIB_WITH_CHOLMOD = $(LIB_WITH_PARTITION) $(CUBLAS_LIB) $(CUDART_LIB)
endif

#-------------------------------------------------------------------------------
CXXp = g++ -std=c++17
C = $(CC) $(CF) $(UMFPACK_CONFIG) $(CONFIG_PARTITION) \
    -I SuiteSparse-5.10.1/include
Cpp = $(CXXp) $(CF) $(UMFPACK_CONFIG) $(CONFIG_PARTITION) \
    -I SuiteSparse-5.10.1/include

# LINUX
# LIBS = $(LDLIBS) -L SuiteSparse-5.10.1/lib -lumfpack -lamd -lsuitesparseconfig \
# 	$(LIB_WITH_CHOLMOD) $(LAPACK)

# MAC
LIBS = $(LDLIBS) -L SuiteSparse-5.10.1/lib -lumfpack -lamd -lsuitesparseconfig \
    $(LIB_WITH_CHOLMOD) $(LAPACK) -Wl,-rpath,SuiteSparse-5.10.1/lib


libs: metis
	( cd SuiteSparse-5.10.1/SuiteSparse_config ; $(MAKE) )
	( cd SuiteSparse-5.10.1/AMD ; $(MAKE) library )
	( cd SuiteSparse-5.10.1/UMFPACK/Lib ; $(MAKE) )
	- ( cd SuiteSparse-5.10.1/CHOLMOD && $(MAKE) library )
	- ( cd SuiteSparse-5.10.1/COLAMD && $(MAKE) library )
	- ( cd SuiteSparse-5.10.1/CCOLAMD ; $(MAKE) library )
	- ( cd SuiteSparse-5.10.1/CAMD ; $(MAKE) library )

metis: SuiteSparse-5.10.1/include/metis.h

SuiteSparse-5.10.1/include/metis.h:
	- ( cd SuiteSparse-5.10.1 && $(MAKE) metis )


analyticalPlacer: analyticalPlacer.cpp 
	$(Cpp) -o analyticalPlacer analyticalPlacer.cpp $(LIBS)

umfpack_simple: umfpack_simple.c
	$(C) -o umfpack_simple umfpack_simple.c $(LIBS)
	
visualizer: visualizer.cpp
	$(ECHO) $(CXX) $(CXX_FLAGS) $(foreach D,$(INC_DIRS),-I$D) $(GTK_INCLUDE_DIRS) $(SRCS) $(GTK_LIBS) $(LIBS) -o $(TARGET_DIR)/$(TARGET)



run: umfpack_di_demo umfpack_zi_demo umfpack_dl_demo umfpack_zl_demo umfpack_simple
	./umfpack_simple



#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) umfpack_simple a.out

clean:
	- $(RM) -r $(CLEAN)
	
	
#-------------------------------------------------------------------------------
VERBOSE ?= 1
ifeq ($(VERBOSE),1)
	ECHO := 
else
	ECHO := @
endif

CONF ?= debug

TARGET_DIR = .
TARGET = visualizer

GTK_VERSION_NUM = 3.0
EZGL_DIR = ezgl

SRCS = $(wildcard ./*.cpp $(EZGL_DIR)/*.cpp)
HDRS = $(wildcard ./*.h $(EZGL_DIR)/*.hpp)

GTK_INCLUDE_DIRS := $(shell pkg-config --cflags gtk+-$(GTK_VERSION_NUM) x11)

GTK_LIBS := $(shell pkg-config --libs gtk+-$(GTK_VERSION_NUM) x11)

INC_DIRS = . $(EZGL_DIR) $(EZGL_DIR)/.. SuiteSparse-5.10.1/include

CXX_FLAGS = -g -Wall -std=c++17

ifeq (release, $(CONF))
	CXX_FLAGS += -O3
else ifeq (debug, $(CONF))
# Don't change anything
else
    $(error Invalid value for CONF: '$(CONF)', must be 'release' or 'debug'. Try 'make help' for usage)
endif

# create the exe
$(TARGET_DIR)/$(TARGET): Makefile $(SRCS)
	$(ECHO) $(CXX) $(CXX_FLAGS) $(foreach D,$(INC_DIRS),-I$D) $(GTK_INCLUDE_DIRS) $(SRCS) $(GTK_LIBS) $(LIBS) -o $(TARGET_DIR)/$(TARGET)


# clean the EXE 
clean:
	$(ECHO) rm -f $(TARGET_DIR)/$(TARGET)

help:
	@echo "Makefile for ezgl example program"
	@echo ""
	@echo "Usage: "
	@echo '    > make -j4'
	@echo "        Call the default make target (all)."
	@echo "        This builds the project executable: '$(TARGET)'."
	@echo "        Use -j4 option to do parallel builds."
	@echo "    > make clean"
	@echo "        Removes any generated files including exectuables "
	@echo "        and object files."
	@echo "    > make help"
	@echo "        Prints this help message."
	@echo ""
	@echo ""
	@echo "Configuration Variables: "
	@echo "    CONF={release | debug}"
	@echo "        Controls whether the build performs compiler optimizations"
	@echo "        to improve performance. Currently set to '$(CONF)'."
	@echo ""
	@echo "        With CONF=release compiler optimization is enabled."
	@echo ""
	@echo "        With CONF=debug compiler optimization is disabled to improve"
	@echo "        interactive debugging."


.PHONY: all clean help
	
