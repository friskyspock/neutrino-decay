CC=g++
file_name := dtset.cpp
exec_name := ../bin/dtset
#wildcard, this finds all *.c file and make *.o files
objects = $(patsubst %.c,%.o,$(wildcard *.c))

#PYTHIA8=/home/akanksha/MG5_aMC_v2_5_5/HEPTools/pythia8
#fastjet=/usr/local
DELPHES_DIR := /home/blizzard/Research/tools/madGraph/Delphes/
FASTJET_DIR := /home/blizzard/Research/tools/fastjet-install/
#HDF5_DIR := /home/student03/tools/hdf5/

#DELPHES_INC := -I$(DELPHES_DIR) -I$(DELPHES_DIR)/classes -I$(DELPHES_DIR)/external
#DELPHES_LIB := -lDelphes -Wl -rpath $(DELPHES_DIR) -L$(DELPHES_DIR)

FJFLAGS := `$(FASTJET_DIR)/bin/fastjet-config --cxxflags`
FJLIBS := `$(FASTJET_DIR)/bin/fastjet-config --libs`

#ROOT_INC := `$(ROOT_DIR)/bin/root-config --incdir`
#ROOT_LIB := `$(ROOT_DIR)/bin/root-config --libs`

PLATFORM_TYPE := $(shell uname -s)
#FJFLAGS := $(shell fastjet-config --cxxflags)
#FJLIBS := $(shell fastjet-config --libs)
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTINC := $(shell root-config --incdir)

#include delphes exroot
INC := -I$(DELPHES_DIR) -I$(DELPHES_DIR)/classes -I$(DELPHES_DIR)/external -I$(HDF5_DIR)/include
LIB := -L$(DELPHES_DIR) -lDelphes -Wl,-rpath=$(DELPHES_DIR)

all: Test_project

Test_project: $(objects) $(file_name)
	$(CXX) -O3 $(objects) $(file_name) $(ROOTFLAGS) $(INC) $(ROOTLIBS) $(LIB) -lfastjet -lNsubjettiness -Wl,-rpath=$(DELPHES_DIR) -o $(exec_name) $(FJFLAGS) $(FJLIBS)  $(ROOTLIBS)

%.o: %.C %.h
	rm  lib/*; gcc -fPIC -O2 -c $(FJFLAGS) $< -o $@

clean:
	rm $(FILES) Test_project lib/*
