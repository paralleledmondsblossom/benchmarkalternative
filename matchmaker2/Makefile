CC=g++
CXX=$(CC)
CFLAGS=-lm -O3
CXXFLAGS=$(CFLAGS)
LD=$(CXX)
NVCC=nvcc
LDFLAGS = $(CFLAGS)
NVCC=nvcc
CUDAFLAGS = -O3 -Xptxas -O3 -Xcompiler -O3 -w 


all: mmaker
mmaker: main.cu cudaBFS.h cudaBFS.cu matchgpu.cu MMArguments.cpp cheap.c matching.c MMArguments.h matchmaker.h
	$(NVCC) $(LDFLAGS) $(CUDAFLAGS)  -c -o main.o main.cu
	$(NVCC) $(LDFLAGS) $(CUDAFLAGS)  -c -o cudaBFS.o cudaBFS.cu
	$(NVCC) $(LDFLAGS) $(CUDAFLAGS)  -c -o matchgpu.o matchgpu.cu
	$(CC) $(LDFLAGS)   -c -o MMArguments.o MMArguments.cpp
	$(CC) $(LDFLAGS)   -c -o cheap.o cheap.c
	$(CC) $(LDFLAGS)   -c -o matching.o matching.c
	$(NVCC) *.o -o mmaker $(LDFLAGS)

mmaker_lib.a: main_lib.cu cudaBFS.h cudaBFS.cu matchgpu.cu MMArguments.cpp cheap.c matching.c MMArguments.h matchmaker.h
	$(NVCC) $(LDFLAGS) $(CUDAFLAGS)  -c -o main_lib.o main_lib.cu
	$(NVCC) $(LDFLAGS) $(CUDAFLAGS)  -c -o cudaBFS.o cudaBFS.cu
	$(NVCC) $(LDFLAGS) $(CUDAFLAGS)  -c -o matchgpu.o matchgpu.cu
	$(CC) $(LDFLAGS)   -c -o MMArguments.o MMArguments.cpp
	$(CC) $(LDFLAGS)   -c -o cheap.o cheap.c
	$(CC) $(LDFLAGS)   -c -o matching.o matching.c
	$(NVCC) -lib *.o -o mmaker_lib.a $(LDFLAGS)

clean: 
	-rm -rf *.o *.a mmaker
depend:
	makedepend -Y *.cpp *.c *.hpp
