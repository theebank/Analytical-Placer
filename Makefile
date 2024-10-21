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
CXX = g++ -std=c++17
C = $(CC) $(CF) $(UMFPACK_CONFIG) $(CONFIG_PARTITION) \
    -I SuiteSparse-5.10.1/include
Cpp = $(CXX) $(CF) $(UMFPACK_CONFIG) $(CONFIG_PARTITION) \
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


run: umfpack_di_demo umfpack_zi_demo umfpack_dl_demo umfpack_zl_demo umfpack_simple
	./umfpack_simple



#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) umfpack_simple a.out

clean:
	- $(RM) -r $(CLEAN)

