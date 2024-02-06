/*
 * graph.c -- Directed graphs
 */

/*
 * Copyright 1989, 1992 by John Kececioglu
 */


/*
 * Synopsis
 *
 * This implementation of directed graphs uses the "forward star" and
 * "backward star" representation.  Undirected graphs may be encoded
 * by choosing an arbitrary direction for edges.
 *
 * Graphs, vertices, and edges may be attributed with an arbitrary label.
 * Subgraphs may be formed ejecting, and later injecting, vertices and edges.
 * Graphs may be written to a text file, and later read back in.
 *
 */

/*
 * Author
 *
 * John Kececioglu
 * kece@cs.uga.edu
 *
 * Department of Computer Science
 * The University of Georgia
 * Athens, GA 30602
 *
 */

/*
 * History
 *
 * 25 July 1993 JDK
 * Made the naming of functions and structures consistent with other libraries.
 * Added functions to read and write graphs.
 *
 */


#include <stdio.h>
#include "graph.h"
#include <assert.h>

#define VertexBlockSize 16 /* Number of vertices allocated per request */
#define EdgeBlockSize   32 /* Number of edges allocated per request */


typedef ListCell Cell;


static Graph  *GraphPool  = Nil; /* Pool of free graphs */
static Vertex *VertexPool = Nil; /* Pool of free vertices */
static Edge   *EdgePool   = Nil; /* Pool of free edges */


static Void Error Proto(( char *Message ));
static float EdgeWeight Proto(( Edge *E ));


#define FreeGraph(G)  (((G)->Vertices = (List *) GraphPool), GraphPool = (G))
#define FreeVertex(V) (((V)->Self = (Cell *) VertexPool), VertexPool = (V))
#define FreeEdge(E)   (((E)->Self = (Cell *) EdgePool), EdgePool = (E))

#define NewGraph(G)  (((G) = GraphPool), GraphPool = (Graph *) (G)->Vertices)
#define NewVertex(V) (((V) = VertexPool), VertexPool = (Vertex *) (V)->Self)
#define NewEdge(E)   (((E) = EdgePool), EdgePool = (Edge *) (E)->Self)
   /*
    * `NewGraph', `NewVertex', and `NewEdge' assume their argument is a
    *  variable
    */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Graphs
 *
 */


/*
 * CreateGraph -- Create an empty graph
 *
 */
Graph *CreateGraph

#ifdef Ansi
   (GraphData D)
#else
   (D) GraphData D;
#endif

{
   register Graph *G;
   
   if (GraphPool != Nil)
      NewGraph(G);
   else if ((G = (Graph *) Allocate(sizeof(Graph))) == NULL)
      Error("(CreateGraph) Memory allocation failed.");
   
   G->Vertices = CreateList();
   G->Edges = CreateList();
   G->Label = D;
   
   return G;
}


/*
 * DestroyGraph -- Destroy a graph
 *
 */
Void DestroyGraph

#ifdef Ansi
  (register Graph *G)
#else
   (G) register Graph *G;
#endif

{
   register Vertex *V;
   register Edge   *E;
   register Cell   *P;
   
   P = ListHead(G->Vertices);
   while ((V = (Vertex *) ListNext(P)))
   {
      DestroyList(V->In);
      DestroyList(V->Out);
      FreeVertex(V);
   }
   DestroyList(G->Vertices);

   P = ListHead(G->Edges);
   while ((E = (Edge *) ListNext(P)))
      FreeEdge(E);
   DestroyList(G->Edges);

   FreeGraph(G);
}


/*
 * GraphLabel -- Return the label of a graph
 *
 */
GraphData GraphLabel

#ifdef Ansi
   (Graph *G)
#else
   (G) Graph *G;
#endif

{
   return G->Label;
}


/*
 * GraphRelabel -- Relabel a graph
 *
 * The old label of the graph is returned.
 *
 */
GraphData GraphRelabel

#ifdef Ansi
   (register Graph *G, GraphData D)
#else
   (G, D) register Graph *G; GraphData D;
#endif

{
   register GraphData X;
   
   X = G->Label;
   G->Label = D;
   return X;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Vertices
 *
 */


/*
 * CreateVertex -- Create a graph vertex
 *
 */
Vertex *CreateVertex

#ifdef Ansi
   (Graph *G, VertexData D)
#else
   (G, D) Graph *G; VertexData D;
#endif

{
   register Vertex *V, *Block;

   if (VertexPool == Nil)
   {
      /*
       * Allocate a block of vertices
       */
      Block = (Vertex *) Allocate(VertexBlockSize * sizeof(Vertex));
      if (Block == NULL)
         Error("(CreateVertex) Memory allocation failed.");
      
      /*
       * Place the vertices in the block into the pool
       */
      for (V = Block; V - Block < VertexBlockSize; V++)
         FreeVertex(V);
   }
   
   NewVertex(V);
   V->In = CreateList();
   V->Out = CreateList();
   V->Self = ListPut((Pointer) V, G->Vertices);
   V->Label = D;
   
   return V;
}


/*
 * DestroyVertex -- Destroy a vertex of a graph
 *
 */
Void DestroyVertex

#ifdef Ansi
   (Vertex *V)
#else
   (V) Vertex *V;
#endif

{
   register Edge *E;
   register Cell *P;
   
   P = ListHead(V->In);
   while ((E = (Edge *) ListNext(P)))
   {
      ListDelete(E->Out);
      ListDelete(E->Self);
      FreeEdge(E);
   }
   DestroyList(V->In);
   
   P = ListHead(V->Out);
   while ((E = (Edge *) ListNext(P)))
   {
      ListDelete(E->In);
      ListDelete(E->Self);
      FreeEdge(E);
   }
   DestroyList(V->Out);
   
   ListDelete(V->Self);
   FreeVertex(V);
}


/*
 * GraphVertices -- Return the list of vertices of a graph
 *
 */
List *GraphVertices

#ifdef Ansi
   (Graph *G)
#else
   (G) Graph *G;
#endif

{
   return G->Vertices;
}


/*
 * VertexLabel -- Return the label of a vertex
 *
 */
VertexData VertexLabel

#ifdef Ansi
   (Vertex *V)
#else
   (V) Vertex *V;
#endif

{
   return V->Label;
}


/*
 * VertexRelabel -- Relabel a vertex
 *
 * The old label of the vertex is returned.
 *
 */
VertexData VertexRelabel

#ifdef Ansi
   (register Vertex *V, VertexData D)
#else
   (V, D) register Vertex *V; VertexData D;
#endif

{
   register VertexData X;

   X = V->Label;
   V->Label = D;
   return X;
}


/*
 * VertexIn -- Return the list of in-edges of a vertex
 *
 */
List *VertexIn

#ifdef Ansi
   (Vertex *V)
#else
   (V) Vertex *V;
#endif

{
   return V->In;
}


/*
 * VertexOut -- Return the list of out-edges of a vertex
 *
 */
List *VertexOut

#ifdef Ansi
   (Vertex *V)
#else
   (V) Vertex *V;
#endif

{
   return V->Out;
}


/*
 * VertexEject -- Eject a vertex from a graph
 *
 */
Void VertexEject

#ifdef Ansi
   (register Vertex *V)
#else
   (V) register Vertex *V;
#endif

{
   register Edge *E;
   register Cell *P;
   
   P = ListHead(V->In);
   while ((E = (Edge *) ListNext(P)))
   {
      ListEject(E->Out);
      ListEject(E->Self);
   }
   
   P = ListHead(V->Out);
   while ((E = (Edge *) ListNext(P)))
   {
      ListEject(E->In);
      ListEject(E->Self);
   }
   
   ListEject(V->Self);
}


/*
 * VertexInject -- Inject an ejected vertex back into a graph
 *
 * Vertices that are adjacent in the graph lists must be injected in reverse
 * order of ejection.
 *
 */
Void VertexInject

#ifdef Ansi
   (register Vertex *V)
#else
   (V) register Vertex *V;
#endif

{
   register Edge *E;
   register Cell *P;
   
   ListInject(V->Self);

   P = ListHead(V->Out);
   while ((E = (Edge *) ListPrev(P)))
   {
      ListInject(E->Self);
      ListInject(E->In);
   }

   P = ListHead(V->In);
   while ((E = (Edge *) ListPrev(P)))
   {
      ListInject(E->Self);
      ListInject(E->Out);
   }
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Edges
 *
 */


/*
 * CreateEdge -- Create a directed edge
 *
 */
Edge *CreateEdge

#ifdef Ansi
   (Graph *G, Vertex *V, Vertex *W, EdgeData D)
#else
   (G, V, W, D) Graph *G; Vertex *V, *W; EdgeData D;
#endif

{
   register Edge *E, *Block;

   if (EdgePool == Nil)
   {
      /*
       * Allocate a block of edges
       */
      Block = (Edge *) Allocate(EdgeBlockSize * sizeof(Edge));
      if (Block == NULL)
         Error("(CreateEdge) Memory allocation failed.");

      /*
       * Place the edges in the block into the pool
       */
      for (E = Block; E - Block < EdgeBlockSize; E++)
         FreeEdge(E);
   }
   
   NewEdge(E);
   E->From = V;
   E->To = W;
   E->In = ListPut((Pointer) E, W->In);
   E->Out = ListPut((Pointer) E, V->Out);
   E->Self = ListPut((Pointer) E, G->Edges);
   E->Label = D;
   
   return E;
}


/*
 * DestroyEdge -- Destroy an edge of a graph
 *
 */
Void DestroyEdge

#ifdef Ansi
   (register Edge *E)
#else
   (E) register Edge *E;
#endif

{
   ListDelete(E->In);
   ListDelete(E->Out);
   ListDelete(E->Self);
   FreeEdge(E);
}


/*
 * GraphEdges -- Return the list of edges of a graph
 *
 */
List *GraphEdges

#ifdef Ansi
   (Graph *G)
#else
   (G) Graph *G;
#endif

{
   return G->Edges;
}


/*
 * EdgeLabel -- Return the label of an edge
 *
 */
EdgeData EdgeLabel

#ifdef Ansi
   (Edge *E)
#else
   (E) Edge *E;
#endif

{
   return E->Label;
}


/*
 * EdgeRelabel -- Relabel an edge
 *
 * The old label of the edge is returned.
 *
 */
EdgeData EdgeRelabel

#ifdef Ansi
   (register Edge *E, EdgeData D)
#else
   (E, D) register Edge *E; EdgeData D;
#endif

{
   register EdgeData X;

   X = E->Label;
   E->Label = D;
   return X;
}


/*
 * EdgeFrom -- Return the source vertex of a directed edge
 *
 */
Vertex *EdgeFrom

#ifdef Ansi
   (Edge *E)
#else
   (E) Edge *E;
#endif

{
   return E->From;
}


/*
 * EdgeTo -- Return the target vertex of a directed edge
 *
 */
Vertex *EdgeTo

#ifdef Ansi
   (Edge *E)
#else
   (E) Edge *E;
#endif

{
   return E->To;
}


/*
 * EdgeEject -- Eject an edge from a graph
 *
 */
Void EdgeEject

#ifdef Ansi
   (register Edge *E)
#else
   (E) register Edge *E;
#endif

{
   ListEject(E->In);
   ListEject(E->Out);
   ListEject(E->Self);
}


/*
 * EdgeInject -- Inject an ejected edge back into a graph
 *
 * Edges that are adjacent in the graph lists must be injected in reverse order
 * of ejection.
 *
 */
Void EdgeInject

#ifdef Ansi
   (register Edge *E)
#else
   (E) register Edge *E;
#endif

{
   ListInject(E->Self);
   ListInject(E->Out);
   ListInject(E->In);
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Input and output
 *
 */


/*
 * WriteGraph -- Write out a graph to a text file
 *
 */
Void WriteGraph

#ifdef Ansi
   (Graph *G, FILE *stream)
#else
   (G, stream) Graph *G; FILE *stream;
#endif

{
   auto     int M, N;
   register int i;
   register Vertex *V;
   register Edge *E;
   register Cell *P;
   register VertexData *L;
   
   /*
    * Count the number of vertices and edges
    */
   N = ListSize(G->Vertices);
   M = ListSize(G->Edges);

   /*
    * Allocate an array to hold the original vertex labels
    */
   L = (VertexData *) Allocate(N * sizeof(VertexData));
   if (L == NULL)
      Error("(WriteGraph) Memory allocation failed.");
   
   /*
    * Relabel the vertices with pointers into the array
    */
   i = 0;
   P = G->Vertices;
   while ((V = (Vertex *) ListNext(P)))
   {
      L[i] = V->Label;
      V->Label = (VertexData) &L[i];
      i += 1;
   }
   
   /*
    * Write out the size of the graph
    */
   fprintf(stream, "vertices %d\n", N);
   fprintf(stream, "edges %d\n", M);
   fprintf(stream, "\n");
   
   /*
    * Write out the list of edges
    */
   P = ListHead(G->Edges);
   while ((E = (Edge *) ListNext(P)))
      fprintf(stream, "edge %ld %ld\n",
         ((VertexData *) E->From->Label) - L + 1,
         ((VertexData *) E->To->Label) - L + 1);
   
   /*
    * Restore the vertex labels
    */
   i = 0;
   P = ListHead(G->Vertices);
   while ((V = (Vertex *) ListNext(P)))
      V->Label = L[i++];
   
   /*
    * Free memory
    */
   Free(L);
}

/*
 * CreateGraphFromCSC -- Read in a graph from a CSC
 *
 * Assumes the graph is provided as a symmetric CSC.
 *
 */
Graph *CreateGraphFromCSC_MS_BFS_GRAFT

#ifdef Ansi
   (int *cxadj, int *cadj, int *matching, int *match_count_scalar, int nr, int nc, int nn, int just_read_file)
#else
   (stream) FILE *stream;
#endif

{
   //auto int M, N;
   register Graph *G;
   //register Vertex **V;
   auto int a, b;
   register int i;

   /*
    * Create an empty graph
    */
   G = CreateGraph(Nil);
   /*
   G->EL.Rows = (int *)malloc(2 * (*M) * sizeof(int));
   G->EL.Cols = (int *)malloc(2 * (*M) * sizeof(int));
   G->EL.Matching = (int *)calloc(*N, sizeof(int));
   G->EL.M = (*M);
   G->EL.N = (*N);
   G->hash = createHashTable(*M);
   //allocateGPUMatcher(G);
   if(G->mm._bfs==NULL){
      exit(0);
   }
   */
   /*
    * Allocate an array to hold onto vertices
    */
   G->VertexArray = (Vertex **) Allocate(nr * sizeof(Vertex *));
   if (G->VertexArray == NULL)
      Error("(ReadGraph) Memory allocation failed.");
   
   /*
    * Insert the vertices
    */
   for (i = 0; i < nr; i++)
      G->VertexArray[i] = CreateVertex(G, Nil);
   
   /*
    * Initialize graph
   */
   Initialize(G, CreateList());

   /*
    * Read the list of edges and insert them
    */
   int *match_count = (int *)calloc(nr, sizeof(int));
   register Edge *E;
   for (int r = 0; r < nr; ++r){
      int start = cxadj[r];
      int end = cxadj[r+1];
      //printf("col %d start %d end %d\n",r-nc,cadj[start],cadj[end]);
      for (;start<end;start++){
            //printf("u %d v %d\n",r-nc,cadj[start]);
            if (r<cadj[start]-nr){
               E = CreateEdge(G, G->VertexArray[r], G->VertexArray[cadj[start]-nr], Nil);
               //printf("u %d (%d) v %d\n",r,matching[r]-nr,cadj[start]-nr);
               if(!just_read_file && (r==matching[matching[r]]) && ((matching[r]-nr)==(cadj[start]-nr)) && match_count[r]<1 && match_count[(cadj[start]-nr)]<1){
                  //printf("MATCHED!!! u %d (%d) v %d\n",r,matching[r]-nr,cadj[start]-nr);
                  match_count[r]++;
                  match_count[cadj[start]-nr]++;
                  //(*match_count_scalar)+=1;
                  Match(E);
               }
            }
      }
   }
   /*
   for (int r = 0; r < nr; ++r){
      if (match_count[r]>1)
         printf("ERROR IN MATCHING %d match count %d\n",r,match_count[r]);
   }
   */
      //CreateEdge(G, V[a], V[b], Nil);
      /*
      G->EL.Rows[2*i] = a - 1;
      G->EL.Rows[2*i + 1] = b - 1;
      G->EL.Cols[2*i] = b - 1;
      G->EL.Cols[2*i + 1] = a - 1;
      if (a < b){
         OrderedPair key = {a - 1,b - 1};
         insert(G->hash, key, E);
         Edge *result1 = get(G->hash, key);
         assert(result1==E);
         Vertex * u = EdgeFrom(result1);
         Vertex * v = EdgeTo(result1);
         assert(u == G->VertexArray[a - 1]);
         assert(v == G->VertexArray[b - 1]);
      } else {
         OrderedPair key = {b - 1,a - 1};
         insert(G->hash, key, E);
         Edge *result1 = get(G->hash, key);
         assert(result1==E);
         Vertex * u = EdgeTo(result1);
         Vertex * v = EdgeFrom(result1);
         assert(u == G->VertexArray[b - 1]);
         assert(v == G->VertexArray[a - 1]); 
      }
      */

   /*
    * Free memory
    */
   //Free(V);
   
   /*
    * Return the graph
    */
   return G;
}


/*
 * CreateGraphFromCSC -- Read in a graph from a CSC - matchmaker2 and BFSHonestPath
 *
 * Assumes the graph is provided as a symmetric CSC.
 *
 */
Graph *CreateGraphFromCSC

#ifdef Ansi
   (int *cxadj, int *cadj, int *matching, int *match_count_scalar, int nr, int nc, int nn, int just_read_file)
#else
   (stream) FILE *stream;
#endif

{
   //auto int M, N;
   register Graph *G;
   //register Vertex **V;
   auto int a, b;
   register int i;

   /*
    * Create an empty graph
    */
   G = CreateGraph(Nil);
   /*
   G->EL.Rows = (int *)malloc(2 * (*M) * sizeof(int));
   G->EL.Cols = (int *)malloc(2 * (*M) * sizeof(int));
   G->EL.Matching = (int *)calloc(*N, sizeof(int));
   G->EL.M = (*M);
   G->EL.N = (*N);
   G->hash = createHashTable(*M);
   //allocateGPUMatcher(G);
   if(G->mm._bfs==NULL){
      exit(0);
   }
   */
   /*
    * Allocate an array to hold onto vertices
    */
   G->VertexArray = (Vertex **) Allocate(nr * sizeof(Vertex *));
   if (G->VertexArray == NULL)
      Error("(ReadGraph) Memory allocation failed.");
   
   /*
    * Insert the vertices
    */
   for (i = 0; i < nr; i++)
      G->VertexArray[i] = CreateVertex(G, Nil);
   
   /*
    * Initialize graph
   */
   Initialize(G, CreateList());

   /*
    * Read the list of edges and insert them
    */
   register Edge *E;
   for (int r = 0; r < nc; ++r){
      int start = cxadj[r];
      int end = cxadj[r+1];
      //printf("col %d start %d end %d\n",r,start,end);
      for (;start<end;start++)
         if (r<cadj[start]){
            E = CreateEdge(G, G->VertexArray[r], G->VertexArray[cadj[start]], Nil);
            if(!just_read_file && (r==matching[matching[r]]) && matching[r]==cadj[start]){
               //(*match_count_scalar)+=1;
               Match(E);
            }
         }
   }

      //CreateEdge(G, V[a], V[b], Nil);
      /*
      G->EL.Rows[2*i] = a - 1;
      G->EL.Rows[2*i + 1] = b - 1;
      G->EL.Cols[2*i] = b - 1;
      G->EL.Cols[2*i + 1] = a - 1;
      if (a < b){
         OrderedPair key = {a - 1,b - 1};
         insert(G->hash, key, E);
         Edge *result1 = get(G->hash, key);
         assert(result1==E);
         Vertex * u = EdgeFrom(result1);
         Vertex * v = EdgeTo(result1);
         assert(u == G->VertexArray[a - 1]);
         assert(v == G->VertexArray[b - 1]);
      } else {
         OrderedPair key = {b - 1,a - 1};
         insert(G->hash, key, E);
         Edge *result1 = get(G->hash, key);
         assert(result1==E);
         Vertex * u = EdgeTo(result1);
         Vertex * v = EdgeFrom(result1);
         assert(u == G->VertexArray[b - 1]);
         assert(v == G->VertexArray[a - 1]); 
      }
      */

   /*
    * Free memory
    */
   //Free(V);
   
   /*
    * Return the graph
    */
   return G;
}

/*
 * ReadGraph -- Read in a graph from a text file
 *
 * Assumes the graph was written using `WriteGraph'.
 *
 */
Graph *ReadGraph

#ifdef Ansi
   (FILE *stream, int * N, int * M)
#else
   (stream) FILE *stream;
#endif

{
   //auto int M, N;
   register Graph *G;
   //register Vertex **V;
   auto int a, b;
   register int i;
   
   /*
    * Read the size of the graph
    */
   if (fscanf(stream, " vertices %d", N) != 1)
      Error("(ReadGraph) Number of vertices not recognized.");
   if (fscanf(stream, " edges %d", M) != 1)
      Error("(ReadGraph) Number of edges not recognized.");
   

   /*
    * Create an empty graph
    */
   G = CreateGraph(Nil);
   /*
   G->EL.Rows = (int *)malloc(2 * (*M) * sizeof(int));
   G->EL.Cols = (int *)malloc(2 * (*M) * sizeof(int));
   G->EL.Matching = (int *)calloc(*N, sizeof(int));
   G->EL.M = (*M);
   G->EL.N = (*N);
   G->hash = createHashTable(*M);
   //allocateGPUMatcher(G);
   if(G->mm._bfs==NULL){
      exit(0);
   }
   */
   /*
    * Allocate an array to hold onto vertices
    */
   G->VertexArray = (Vertex **) Allocate(*N * sizeof(Vertex *));
   if (G->VertexArray == NULL)
      Error("(ReadGraph) Memory allocation failed.");
   
   /*
    * Insert the vertices
    */
   for (i = 0; i < *N; i++)
      G->VertexArray[i] = CreateVertex(G, Nil);
   
   /*
    * Read the list of edges and insert them
    */
   register Edge *E;
   for (i = 0; i < *M; i++)
   {
      if (fscanf(stream, " edge %d %d", &a, &b) != 2)
         Error("(ReadGraph) Edge not recognized.");
      // Instead of ignoring the return value, use it to make hash table.
      E = CreateEdge(G, G->VertexArray[a - 1], G->VertexArray[b - 1], Nil);
      //CreateEdge(G, V[a], V[b], Nil);
      /*
      G->EL.Rows[2*i] = a - 1;
      G->EL.Rows[2*i + 1] = b - 1;
      G->EL.Cols[2*i] = b - 1;
      G->EL.Cols[2*i + 1] = a - 1;
      if (a < b){
         OrderedPair key = {a - 1,b - 1};
         insert(G->hash, key, E);
         Edge *result1 = get(G->hash, key);
         assert(result1==E);
         Vertex * u = EdgeFrom(result1);
         Vertex * v = EdgeTo(result1);
         assert(u == G->VertexArray[a - 1]);
         assert(v == G->VertexArray[b - 1]);
      } else {
         OrderedPair key = {b - 1,a - 1};
         insert(G->hash, key, E);
         Edge *result1 = get(G->hash, key);
         assert(result1==E);
         Vertex * u = EdgeTo(result1);
         Vertex * v = EdgeFrom(result1);
         assert(u == G->VertexArray[b - 1]);
         assert(v == G->VertexArray[a - 1]); 
      }
      */
   }

   /*
    * Free memory
    */
   //Free(V);
   
   /*
    * Return the graph
    */
   return G;
}


/*
 * WriteEdgeWeightedGraph -- Write out an edge-weighted graph to a text file
 *
 */
Void WriteEdgeWeightedGraph

#ifdef Ansi
   (Graph *G, float (*Weight)(Edge *), FILE *stream)
#else
   (G, Weight, stream) Graph *G; float (*Weight)(); FILE *stream;
#endif

{
   auto     int M, N;
   register int i;
   register Vertex *V;
   register Edge *E;
   register Cell *P;
   register VertexData *L;
   
   /*
    * Count the number of vertices and edges
    */
   N = ListSize(G->Vertices);
   M = ListSize(G->Edges);

   /*
    * Allocate an array to hold the original vertex labels
    */
   L = (VertexData *) Allocate(N * sizeof(VertexData));
   if (L == NULL)
      Error("(WriteEdgeWeightedGraph) Memory allocation failed.");
   
   /*
    * Relabel the vertices with pointers into the array
    */
   i = 0;
   P = ListHead(G->Vertices);
   while ((V = (Vertex *) ListNext(P)))
   {
      L[i] = V->Label;
      V->Label = (VertexData) &L[i];
      i += 1;
   }
   
   /*
    * Write out the size of the graph
    */
   fprintf(stream, "vertices %d\n", N);
   fprintf(stream, "edges %d weighted\n", M);
   fprintf(stream, "\n");
   
   /*
    * Write out the list of edges
    */
   P = ListHead(GraphEdges(G));
   while ((E = (Edge *) ListNext(P)))
      fprintf(stream, "edge %ld %ld %g\n",
         ((VertexData *) E->From->Label) - L + 1,
         ((VertexData *) E->To->Label) - L + 1,
         Weight(E));
   
   /*
    * Restore the vertex labels
    */
   i = 0;
   P = ListHead(G->Vertices);
   while ((V = (Vertex *) ListNext(P)))
      V->Label = L[i++];
   
   /*
    * Free memory
    */
   Free(L);
}


/*
 * ReadEdgeWeightedGraph -- Read in an edge-weighted graph from a text file
 *
 * The graph is labelled with a pointer to an array of edge weights.
 * The edge weight function is returned through `Weight'.
 * This function assumes the graph was written using `WriteEdgeWeightedGraph'.
 *
 */
Graph *ReadEdgeWeightedGraph

#ifdef Ansi
   (FILE *stream, float (**Weight)(Edge *))
#else
   (stream, Weight) FILE *stream; float (**Weight)();
#endif

{
   auto int M, N;
   register Graph *G;
   register Vertex **V;
   register float *W;
   auto int a, b;
   register int i;
   
   /*
    * Read the size of the graph
    */
   if (fscanf(stream, " vertices %d", &N) != 1)
      Error("(ReadEdgeWeightedGraph) Number of vertices not recognized.");
   if (fscanf(stream, " edges %d weighted", &M) != 1)
      Error("(ReadEdgeWeightedGraph) Number of edges not recognized.");
   
   /*
    * Allocate arrays to hold onto vertices and edge weights
    */
   V = (Vertex **) Allocate(N * sizeof(Vertex *));
   W = (float *) Allocate(M * sizeof(float));
   if (V == NULL || W == NULL)
      Error("(ReadEdgeWeightedGraph) Memory allocation failed.");
   
   /*
    * Create an empty graph
    */
   G = CreateGraph((GraphData) W);
   
   /*
    * Insert the vertices
    */
   for (i = 0; i < N; i++)
      V[i] = CreateVertex(G, Nil);
   
   /*
    * Read the list of edges and insert them
    */
   for (i = 0; i < M; i++)
   {
      if (fscanf(stream, " edge %d %d %f", &a, &b, &W[i]) != 3)
         Error("(ReadEdgeWeightedGraph) Edge not recognized.");
      CreateEdge(G, V[a - 1], V[b - 1], (EdgeData) &W[i]);
   }

   /*
    * Free memory
    */
   Free(V);
   
   /*
    * Return the graph and its edge weight function
    */
   *Weight = EdgeWeight;
   return G;
}


/*
 * EdgeWeight -- Return the weight of an edge created by
 *               `ReadEdgeWeightedGraph'
 *
 */
static float EdgeWeight

#ifdef Ansi
   (Edge *E)
#else
   (E) Edge *E;
#endif

{
   return *((float *) EdgeLabel(E));
}


/*
 * Error -- Write an error message and halt
 *
 */
static Void Error

#ifdef Ansi
   (char *Message)
#else
   (Message) char *Message;
#endif

{
   fprintf(stderr, "%s\n", Message);
   Halt();
}
