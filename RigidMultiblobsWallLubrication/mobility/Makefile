# location of the Python header files
PYTHON_VERSION = 2.7
PYTHON_DIR = /usr
# PYTHON_DIR = /usr/local/pkg/python/$(PYTHON_VERSION)
PYTHON_INCLUDE = $(PYTHON_DIR)/include/python$(PYTHON_VERSION)
PYTHON_LIB = $(PYTHON_DIR)/lib/python$(PYTHON_VERSION)/config/

# location of the Boost Python include files and library
BOOST_INC = /usr/include/boost
BOOST_LIB = /usr/lib64

# Location of Fortran libraries,
# uncomment to use the Fast Multipole Metho library dstfmmlib3d-1.2
# FORTRAN_LIB = --link-fopenmp -L/path/to/stfmmlib3d-1.2/ -lstfmm3d 

# 
TARGET = mobility_ext
TARGET_FORTRAN = mobility_fmm

all: $(TARGET).so
ifneq ($(FORTRAN_LIB),)
all: $(TARGET_FORTRAN).so
endif

$(TARGET_FORTRAN).so: $(TARGET_FORTRAN).f90
	f2py $(FORTRAN_LIB) -c $(TARGET_FORTRAN).f90 -m $(TARGET_FORTRAN)

$(TARGET).so: $(TARGET).o
	g++ -O3 -shared  -dynamiclib $(TARGET).o -L$(BOOST_LIB)  -lboost_python -L$(PYTHON_LIB) -lpython$(PYTHON_VERSION) -o $(TARGET).so

$(TARGET).o: $(TARGET).cc
	g++ -O3 -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).cc

clean: 
	rm -f $(TARGET).o $(TARGET).so $(TARGET_FORTRAN).o $(TARGET_FORTRAN).so
