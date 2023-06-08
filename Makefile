# this is the makefile for the DSGF_C project

TARGET = DSGF
CC = icc
FLAGS = -Ofast -Wall -qopenmp -std=c99
INCLUDE_CUSTOM = include/
INCLUDE_LAPACK = ${MKLROOT}/include
SRC = src
LAPACK_DIR = ${MKLROOT}/lib/intel64
LAPACK_LIBRARIES = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
BUILD_DIR = build
OBJECT_FILES = $(wildcard $(BUILD_DIR)/*.o)

default: $(TARGET)

$(TARGET): $(OBJECT_FILES)
	$(CC) $(FLAGS) -I $(INCLUDE_CUSTOM) -I $(INCLUDE_LAPACK) -o $(TARGET) $(OBJECT_FILES) -L $(LAPACK_DIR) $(LAPACK_LIBRARIES)

DSGF_main.o: DSGF_main.c, 


$(BUILD_DIR):
	mkdir $(BUILD_DIR)

clean:
	rm -rf results
