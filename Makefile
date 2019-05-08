#
# MNS_CD requires Sundials <= 2.7.0
#

SUNDIALS_DIR ?= /usr/local
SUNDIALS_INCLUDE ?= $(SUNDIALS_DIR)/include
SUNDIALS_LIB ?=$(SUNDIALS_DIR)/lib

CXX ?= clang++

all: MNS_CD

MNS_CD: MNS_CD.cpp Input.cpp Input.h Constants.h
	$(CXX) -O3 -o MNS_CD -I $(SUNDIALS_INCLUDE) MNS_CD.cpp Input.cpp -L $(SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial

clean:
	rm -f MNS_CD

cleanall: clean
	 rm -f Output Profile_[0123]
