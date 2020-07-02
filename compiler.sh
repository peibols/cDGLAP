#!/bin/bash

g++ -g -std=c++14 -fopenmp -Wall \
	main.cc Kernels.cc MomentsMedEvol.cc cDGLAP.cc Util.cc Collimator.cc \
	MedMicrojet.cc VacMicrojet.cc MomentsVacEvol.cc \
	-o microjets \
	-lgsl -lgslcblas -lm
