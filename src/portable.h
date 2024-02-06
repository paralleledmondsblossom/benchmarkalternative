/*
 * portable.h -- Isolate system dependencies for code portability
 */

/*
 * Copyright 1989 by John Kececioglu
 */


#ifndef PortableInclude
#define PortableInclude


#define Ansi
   /*
    * Should not be defined if the compiler does not support the ANSI standard
    */


/*
 * Argument lists for function prototypes
 */
#ifdef Ansi
#define Proto(arguments) arguments
#else
#define Proto(arguments) ()
#endif


/*
 * Void data type
 */
#ifdef Ansi
typedef void Void;
#else
typedef int Void;
#endif


/*
 * Generic pointer data type
 */
#ifdef Ansi
typedef void *Pointer;
#else
typedef char *Pointer;
#endif


/*
 * Dynamic memory allocation
 */
#ifdef Ansi
#include <stdlib.h>
#define MallocArgumentType size_t
#else
#define MallocArgumentType unsigned int
#endif

#define Allocate(bytes) malloc((MallocArgumentType) (bytes))
#define Free(memory)    free((Pointer) (memory))


/*
 * Halting the program
 */
#define Halt() exit(1)


/*
 * Random number generation
 */
#define MaximumRandomInteger         ((long) 2147483647)
#define GenerateRandomInteger()      ((long) lrand48())
#define SeedRandomIntegerGenerator() srand48((long) 1)


#endif /* PortableInclude */
