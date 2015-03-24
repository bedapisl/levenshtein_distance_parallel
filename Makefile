CPP=g++
#CFLAGS=-Wall -msse4.1 -O3 -std=c++11 -fopenmp			#release
CFLAGS=-Wall -msse4.1 -O3 -std=c++11 -fopenmp -ggdb		#debug
INCLUDE=./internal
LDFLAGS=
LIBS=
LIBDIRS=
SOURCE=levenshtein.cpp
HEADERS=$(shell find . -name '*.hpp')
EXECUTABLE=./levenshtein


.PHONY: all clear clean purge

all: $(EXECUTABLE)



# Building Targets

$(EXECUTABLE): $(SOURCE) $(HEADERS)
	@echo Compiling and linking executable "$@" ...
	@$(CPP) $(CFLAGS) $(addprefix -I,$(INCLUDE)) $(LDFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS)) $< -o $@



# Cleaning Stuff

clear:
	@echo Removing all generated files...
	-@rm -f $(EXECUTABLE)

clean: clear

purge: clear
