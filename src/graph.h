/*
 * graph.h -- Directed graph definitions
 */

/*
 * Copyright 1989, 1992, 1996 by John Kececioglu
 */


#ifndef GraphInclude
#define GraphInclude


#include <stdio.h>
#include "portable.h"
#include "list.h"
#include "edge_vertex.h"

/*
 * Graph attribute types
 */
typedef Pointer GraphData;
typedef Pointer VertexData;
typedef Pointer EdgeData;


/*
 * Edge list for CUDA
 */
typedef struct EdgeListSOA {
   int *Rows;
   int *Cols;
   int *Matching;
   int M;
   int N;
} EdgeList;

/*
 * Edge list for CUDA
 */
typedef struct Matchmaker2 {
        int nr, nc, nn;
        int *rxadj;
        int *radj;
        
        int *cxadj, *cadj;
        
        int *_rxadj;
        int *_radj;
        int *_cxadj;
        int *_cadj;
        int *_cmatch;
        int *_rmatch;
        int *_is_inserted;
        int *_bfs, *_preced;
        
        int *_non_matched_found;
        int *_is_inserted2;
        
        int *_root_array;
        
        int match_types[11];
} Matcher;


/*
 * Directed graph
 */
typedef struct GraphStruct {
   List *Vertices;
   List *Edges;
   GraphData Label; /* Graph attribute */
      /*
       * `Vertices' is reused for the pool of free graphs
       */
   
   //EdgeList EL;
   //Matcher mm;
   Vertex ** VertexArray;
   //HashTable * hash;
} Graph;


/*
 * Graphs
 */
extern Graph *CreateGraph  Proto(( GraphData D ));
extern Void   DestroyGraph Proto(( Graph *G ));
   
extern List     *GraphVertices Proto(( Graph *G ));
extern List     *GraphEdges    Proto(( Graph *G ));
extern GraphData GraphLabel    Proto(( Graph *G ));
extern GraphData GraphRelabel  Proto(( Graph *G, GraphData D ));


/*
 * Vertices
 */
extern Vertex *CreateVertex  Proto(( Graph *G, VertexData D ));
extern Void    DestroyVertex Proto(( Vertex *V));
   
extern List      *VertexIn      Proto(( Vertex *V ));
extern List      *VertexOut     Proto(( Vertex *V ));
extern VertexData VertexLabel   Proto(( Vertex *V ));
extern VertexData VertexRelabel Proto(( Vertex *V, VertexData D ));
   
extern Void VertexEject  Proto(( Vertex *V ));
extern Void VertexInject Proto(( Vertex *V ));


/*
 * Edges
 */
extern Edge *CreateEdge  Proto(( Graph *G, Vertex *V, Vertex *W, EdgeData D ));
extern Void  DestroyEdge Proto(( Edge *E ));
   
extern Vertex  *EdgeFrom    Proto(( Edge *E ));
extern Vertex  *EdgeTo      Proto(( Edge *E ));
extern EdgeData EdgeLabel   Proto(( Edge *E ));
extern EdgeData EdgeRelabel Proto(( Edge *E, EdgeData D ));
   
extern Void EdgeEject  Proto(( Edge *E ));
extern Void EdgeInject Proto(( Edge *E ));


/*
 * Reading and writing
 */
extern Void WriteGraph
   Proto(( Graph *G, FILE *stream ));

extern Graph *ReadGraph
   Proto(( FILE *stream , int * N, int * M));

extern Graph *CreateGraphFromCSC_MS_BFS_GRAFT
   Proto(( int *cxadj, int *cadj, int *matching, int *match_count_scalar, int nr_ptr, int nc_ptr, int nn_ptr, int just_read_file));

extern Graph *CreateGraphFromCSC
   Proto(( int *cxadj, int *cadj, int *matching, int *match_count_scalar, int nr_ptr, int nc_ptr, int nn_ptr, int just_read_file));

extern Void Initialize
   Proto(( Graph *, List * ));

extern Void Match
   Proto(( Edge * ));

extern Void WriteEdgeWeightedGraph
   Proto(( Graph *G, float (*Weight)(Edge *), FILE *stream ));

extern Graph *ReadEdgeWeightedGraph
   Proto(( FILE *stream, float (**Weight)(Edge *) ));


/*
 * Iteration
 */
#define ForAllVertices(V, L, P) \
        ForAllListElements(V, L, Vertex *, P)

#define ForAllGraphVertices(V, G, P) \
        ForAllVertices(V, GraphVertices(G), P)

#define ForAllInVertices(V, W, P) \
        for ((P) = ListHead(VertexIn(V)); \
             (W) = (ListItem(P) ? EdgeFrom((Edge *) ListNext(P)) : Nil); )

#define ForAllOutVertices(V, W, P) \
        for ((P) = ListHead(VertexOut(V)); \
             (W) = (ListItem(P) ? EdgeTo((Edge *) ListNext(P)) : Nil); )

#define ForAllEdges(E, L, P) \
        ForAllListElements(E, L, Edge *, P)

#define ForAllGraphEdges(E, G, P) \
        ForAllEdges(E, GraphEdges(G), P)

#define ForAllInEdges(V, E, P) \
        ForAllEdges(E, VertexIn(V), P)

#define ForAllOutEdges(V, E, P) \
        ForAllEdges(E, VertexOut(V), P)

#define ForAllIncidentEdges(V, E, P, Q) \
        for ((P) = ListHead(VertexIn(V)), \
             (Q) = ListHead(VertexOut(V)); \
             ((E) = (Edge *) ListItem(P)) || \
             ((E) = (Edge *) ListItem(Q)); \
             ListItem(P) ? ListNext(P) : ListNext(Q))

#define ForAllAdjacentVertices(V, W, P, Q) \
        for ((P) = ListHead(VertexIn(V)), \
             (Q) = ListHead(VertexOut(V)); \
             (ListItem(P) ? ((W) = EdgeFrom((Edge *) ListItem(P)), 1) : 0) || \
             (ListItem(Q) ? ((W) = EdgeTo(  (Edge *) ListItem(Q)), 1) : 0); \
             ListItem(P) ? ListNext(P) : ListNext(Q))

#define EdgeOther(E, V) \
        (EdgeFrom(E) != (V) ? EdgeFrom(E) : EdgeTo(E))

                         
#endif /* GraphInclude */
