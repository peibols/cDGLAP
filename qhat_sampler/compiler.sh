#!/bin/bash

g++ -g -std=c++0x -mcmodel=medium \
	$1.cc Random.cc Glauber_$2.cc Hydro_$2.cc -o $1_$2 \
	-O2 -ansi -pedantic -W -Wall -Wshadow -Wcast-align -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winvalid-pch -Werror=missing-braces -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default -Wno-unused -Wvariadic-macros -Wwrite-strings -Werror=declaration-after-statement -Werror=implicit-function-declaration -Werror=nested-externs -Werror=old-style-definition -Werror=strict-prototypes -fPIC
