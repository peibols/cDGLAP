#!/bin/bash

g++ -g -std=c++14 -Wall \
	main.cc Kernels.cc no_OpenMP_MomentsMedEvol.cc cDGLAP.cc Util.cc Collimator.cc \
	MedMicrojet.cc VacMicrojet.cc MomentsVacEvol.cc \
	-o no_OpenMP_microjets \
	-lgsl -lgslcblas -lm
