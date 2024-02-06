/*
 * graph.h -- Directed graph definitions
 */

/*
 * Copyright 1989, 1992, 1996 by John Kececioglu
 */


#ifndef Edge_Vertex
#define Edge_Vertex


#include <stdio.h>
#include "portable.h"
#include "list.h"


/*
 * Graph attribute types
 */
typedef Pointer GraphData;
typedef Pointer VertexData;
typedef Pointer EdgeData;

/*
 * Graph vertex
 */
typedef struct {
   List *In;         /* In-edges */
   List *Out;        /* Out-edges */
   ListCell *Self;   /* Cell on graph vertex list */
   VertexData Label; /* Vertex attribute */
      /*
       * `Self' is reused for the pool of free vertices
       */
} Vertex;

/*
 * Directed edge
 */
typedef struct {
   Vertex *From;   /* Source vertex of edge */
   Vertex *To;     /* Target vertex of edge */
   ListCell *In;   /* Cell on in-edge list of target vertex */
   ListCell *Out;  /* Cell on out-edge list of source vertex */
   ListCell *Self; /* Cell on graph edge list */
   EdgeData Label; /* Edge attribute */
      /*
       * `Self' is reused for the pool of free edges
       */
} Edge;

                         
#endif /* Edge_Vertex */
