/*
 * matching.c -- Maximum cardinality matching of a general graph
 */

/*
 * Copyright 1996 by John Kececioglu
 */


/*
 * Synopsis
 *
 * This is an implementation of Edmond's algorithm for computing a maximum
 * cardinality matching of a general graph.  For a graph of n vertices and
 * m edges, it runs in O(n m alpha(n,m)) time, where alpha(n,m) is an inverse
 * of Ackermann's function.
 *
 * To reduce the number of phases of Edmond's algorithm, the implementation
 * begins with a maximal matching obtained by a greedy heuristic in
 * O(m + n) time.  This initial matching always has at least half the number
 * of edges of a maximum matching.
 *
 * See:  Robert Endre Tarjan, Data Structures and Network Algorithms,
 * Society for Industrial and Applied Mathematics, Philadelphia, pp. 113-123,
 * 1983.
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
 * 6 November 1997 JDK
 * Implemented a more careful ordering of edges that are pushed on the
 * depth-first search stack.  Edges that form blossoms appear after edges
 * that do not.  Also, on encountering an edges that forms an augmenting
 * path, no more edges are pushed onto the search stack.
 *
 * 5 November 1997 JDK
 * Sped up the greedy maximal matching algorithm using a discrete bucketed heap
 * of vertices prioritized by degree, instead of a pairing heap as before.
 * This reduces the running time to O(n + m), down from O(n log n + m).
 *
 * Fixed a bug in Augment:  the procedure only expanded blossoms on the
 * augmenting path, which is incorrect.  The code now expands all blossoms
 * in the alternating tree that contains the augmenting path.
 *
 * 3 March 1997 JDK
 * Fixed a bug in MaximalMatching:  the reuse of the OriginalVertexLabel
 * macro from the maximum-cardinality matching code in the maximal matching
 * code, which cast the return value of OriginalLabel into a pointer to a
 * VertexAttribute, meant that the offset in the struct for the field
 * OriginalVertexLabelField was not correct.  A new OriginalMaximalVertexLabel
 * macro was added that correctly cast the return value of OriginalLabel into
 * a pointer to a MaximalVertexAttribute.
 *
 * 3 March 1997 JDK
 * Fixed a bug in Shrink:  in the loop that pushed onto the search stack the
 * edges incident to odd vertices on the odd-length cycle, the variable E,
 * which held the bridge edge along the cycle, was re-used as the loop
 * variable, destroying the original value of E.
 *
 * 28 February 1997 JDK
 * Fixed a bug in the implementation discovered by Justin Pecqueur.
 * The old implementation used a stack of vertices to grow the alternating
 * tree, and as a consequence was not performing a true depth-first search.
 * This caused the following invariant, that an edge between two even-labeled
 * vertices forms an odd-length alternating cycle, to not hold.  The new
 * implementation uses a stack of edges to grow the alternating tree, and
 * performs a true depth-first search.
 *
 * 17 December 1996 JDK
 * Completed the initial implementation with help from Justin Pecqueur.
 *
 */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Includes
 *
 */


#include <stdio.h>
#include "matching.h"
#include "set.h"
#include <time.h>
//#include "bipartite.h"
#include <assert.h>

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Types
 *
 */


typedef ListCell Cell;

typedef struct
{
   Vertex *BaseField;
}
SetAttribute;

typedef struct
{
   Element *BlossomField;
   Edge    *MatchField;
   Edge    *TreeField;
   Edge    *BridgeField;
   Vertex  *ShoreField;
   short    LabelField;
   int      AgeField;
   Cell    *SelfField;
   VertexData    OriginalVertexLabelField;
   SetAttribute *OriginalSetLabelField;

   
#ifdef Debug
   
   int   NameField;
   List *MembersField;
   List *ChildrenField;
   
#endif /* Debug */
   
   
}
VertexAttribute;


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Macros
 *
 */


#define False 0
#define True  1

#define EvenLabel      2
#define OddLabel       3
#define UnreachedLabel 4


#define Blossom(V) (((VertexAttribute *) VertexLabel(V))->BlossomField)
#define Match(V)   (((VertexAttribute *) VertexLabel(V))->MatchField)
#define Tree(V)    (((VertexAttribute *) VertexLabel(V))->TreeField)
#define Bridge(V)  (((VertexAttribute *) VertexLabel(V))->BridgeField)
#define Shore(V)   (((VertexAttribute *) VertexLabel(V))->ShoreField)
#define Label(V)   (((VertexAttribute *) VertexLabel(V))->LabelField)
#define Age(V)     (((VertexAttribute *) VertexLabel(V))->AgeField)
#define Self(V)    (((VertexAttribute *) VertexLabel(V))->SelfField)

#define OriginalVertexLabel(V) \
        (((VertexAttribute *) VertexLabel(V))->OriginalVertexLabelField)

#define OriginalSetLabel(V) \
        (((VertexAttribute *) VertexLabel(V))->OriginalSetLabelField)

#define Base(E) \
        (((SetAttribute *) SetLabel(E))->BaseField)

#define IsMatched(V) (Match(V) != Nil)
#define IsReached(V) (Label(V) != UnreachedLabel)
#define IsEven(V)    (Label(V) == EvenLabel)
#define IsOdd(V)     (Label(V) == OddLabel)

#define Other(E, V) (((E) == Nil) ? Nil : EdgeOther(E, V))


#ifdef Debug

#define Name(V) \
        ((V) != Nil ? ((VertexAttribute *) VertexLabel(V))->NameField : 0)

#define Members(V) \
        (((VertexAttribute *) VertexLabel(V))->MembersField)

#define Children(V) \
        (((VertexAttribute *) VertexLabel(V))->ChildrenField)

#endif /* Debug */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Function declarations
 *
 */


//static Void Initialize
//   Proto(( Graph *, List * ));

static Void Terminate
   Proto(( Graph * ));

static short Search
   Proto(( Vertex *, List **, List ** ));

static List *Recover
   Proto(( Vertex * ));

static Void Shrink
   Proto(( Edge *, List ** ));

static Void Augment
   Proto(( List *, List * ));

static Void Path
   Proto(( Vertex *, Vertex *, List * ));

static List *Matching
   Proto(( Graph * ));

static Void Error
   Proto(( char * ));


#ifdef Debug

static Void DumpAlternatingForest
   Proto(( Void ));

static Void Traverse
   Proto(( Vertex *, Edge *, int ));

#endif /* Debug */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Global variables
 *
 */


static VertexAttribute *VertexAttributes;

static SetAttribute *SetAttributes;

static int Time;


#ifdef Debug

static Graph *UnderlyingGraph;

#endif /* Debug */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Computing a maximum cardinality matching
 *
 */


/*
 * MaximumCardinalityMatching -- Compute a maximum cardinality matching of
 *                               a nonbipartite graph
 *
 * This implementation of Edmond's algorithm runs in O(n m alpha(n,m)) time
 * for a graph with n vertices and m edges, where alpha(n,m) is the amortized
 * time per operation for the disjoint-set data structure.
 *
 * A list of the edges in a maximum matching is returned.
 *
 */
List *MaximumCardinalityMatching

#ifdef Ansi
   (Graph *G)
#else
   (G) Graph *G;
#endif

{
   List   *Roots;
   Vertex *V;
   List   *M;
   List   *P, *T;
   Cell   *C;  
   
   //Initialize(G, CreateList());
   Roots = CreateList();
   ForAllGraphVertices(V, G, C)
      if (!IsMatched(V))
         Self(V) = ListPut(V, Roots);

   while ((V = ListGet(Roots)))
      if (Search(V, &P, &T))
         Augment(P, T);
   M = Matching(G);
      
   DestroyList(Roots);
   Terminate(G);
   
   return M;
}


List *MaximumCardinalityMatchingTrack

#ifdef Ansi
   (Graph *G, FILE * outputFileX,FILE * outputFileY,FILE * outputFileZ)
#else
   (G) Graph *G;
#endif

{
   List   *Roots;
   Vertex *V;
   List   *M;
   List   *P, *T;
   Cell   *C;  
   
   //Initialize(G, CreateList());
   Roots = CreateList();
   ForAllGraphVertices(V, G, C)
      if (!IsMatched(V))
         Self(V) = ListPut(V, Roots);
   #ifndef NDEBUG
   int numDead=0;
   #endif
   while ((V = ListGet(Roots))){
      if (Search(V, &P, &T)){
         
         #ifndef NDEBUG
         fprintf(outputFileX, "%d\n",2*ListSize(P)-1);
         fprintf(outputFileY, "%d\n",ListSize(T));
         #endif
         Augment(P, T);


      } 
      #ifndef NDEBUG
      else {
         numDead+=ListSize(T);
         DestroyList(T);
         P = Nil;
         T = Nil;
      }
      fprintf(outputFileZ, "%d\n",numDead);
      #endif
   }
   M = Matching(G);
      
   DestroyList(Roots);
   Terminate(G);
   
   return M;
}

/*
 * Initialize -- Given an approximate matching, initialize the vertex, edge,
 *               and disjoint set data structures for the maximum cardinality
 *               matching computation
 *
 */
//static Void Initialize
Void Initialize
   
#ifdef Ansi
   (Graph *G, List *M)
#else
   (G, M) Graph *G; List *M;
#endif

{
   Cell   *P;
   Vertex *V;
   Edge   *E;
   VertexAttribute *A;
   SetAttribute    *B;

   
#ifdef Debug
   
   int N;
   
   N = 1;
   UnderlyingGraph = G;
   
#endif /* Debug */

   int useInit = 0;
   //if (useInit)
   //   bipartite(G);
   
   Time = 1;
   VertexAttributes = (VertexAttribute *)
      Allocate(ListSize(GraphVertices(G)) * sizeof(VertexAttribute));
   SetAttributes = (SetAttribute *)
      Allocate(ListSize(GraphVertices(G)) * sizeof(SetAttribute));
   if (VertexAttributes == NULL || SetAttributes == NULL)
      Error("(MaximumCardinalityMatching) Memory allocation failed.");
  
   A = VertexAttributes;
   B = SetAttributes;
   ForAllGraphVertices(V, G, P)
   {
      Element   *X;
      VertexData D;

      D = VertexRelabel(V, A);
      X = CreateElement(B);
      Match(V) = Nil;
      Label(V) = UnreachedLabel;
      Age(V) = 0;
      Self(V) = Nil;
      Blossom(V) = X;
      Base(X) = V;
      OriginalVertexLabel(V) = D;
      OriginalSetLabel(V) = B;

      
#ifdef Debug
      
      A->NameField = N++;
      
#endif /* Debug */

      
      A++;
      B++;
   }
   /*
   if (useInit){
      clock_t start_time, end_time;
      double elapsed_time_ms;
      start_time = clock();
      int vertexID = 0;
      ForAllGraphVertices(V, G, P)
      {
         int matchPair = G->EL.Matching[vertexID];
         if (matchPair >= 0 && vertexID < matchPair){
            OrderedPair key = {vertexID,matchPair};
            Edge *result1 = get(G->hash, key);
            assert(result1 != NULL);
            Match(EdgeFrom(result1)) = result1;
            Match(EdgeTo(result1)) = result1;
         }
         vertexID++;
      }
      end_time = clock();
      // Calculate the elapsed time in milliseconds
      elapsed_time_ms = ((double)(end_time - start_time) / CLOCKS_PER_SEC) * 1000.0;
      // Print the elapsed time in milliseconds

      // Print the elapsed time in milliseconds
      printf("Load Matching from GPU : Elapsed Time: %.2f milliseconds\n", elapsed_time_ms);
      printf("Load Matching from GPU : Elapsed Time: %.2f seconds\n", elapsed_time_ms/1000.0);
   }
   */
   DestroyList(M);
}


/*
 * Match -- Given an approximate matching, initialize the vertex, edge,
 *               and disjoint set data structures for the maximum cardinality
 *               matching computation
 *
 */
Void Match
   
#ifdef Ansi
   (Edge *E)
#else
   (E) Edge *E;
#endif

{
   Match(EdgeFrom(E)) = E;
   Match(EdgeTo(E)) = E;
}

/*
 * Terminate -- Free the vertex, edge, and disjoint set data structures used
 *              by the matchings computation
 *
 */
static Void Terminate

#ifdef Ansi
   (Graph *G)
#else
   (G) Graph *G;
#endif

{
   Cell   *P;
   Vertex *V;
   
   ForAllGraphVertices(V, G, P)
   {
      DestroyElement(Blossom(V));
      VertexRelabel(V, OriginalVertexLabel(V));
   }
   
   Free(VertexAttributes);
   Free(SetAttributes);
}


/*
 * Search -- Explore an alternating tree rooted at V in depth-first order
 *
 * Returns true if an augmenting path starting from V exists.  If such a
 * path exists, the unmatched edges on the path are returned through P,
 * and all vertices in the alternating tree containing the path are returned
 * through Q.
 *
 */
static short Search

#ifdef Ansi
   (Vertex *V, List **P, List **Q)
#else
   (V, P, Q) Vertex *V; List **P, **Q;
#endif

{
   register short   Found;
   register Vertex *W;
   register Vertex *X, *Y, *Z;
   auto     Edge   *E, *F;
   auto     List   *S, *T, *U;
   register Cell   *C, *D;
  
   Label(V) = EvenLabel;
   Age(V) = Time++;
   Found = False;

   T = CreateList();
   ListPut(V, T);
   
   S = CreateList();
   ForAllIncidentEdges(V, E, C, D)
   {
      ListPush(E, S);

      W = Other(E, V);
      if (!IsReached(W) && !IsMatched(W))
         break;
   }
   
   while (!ListIsEmpty(S) && !Found)
   {
      E = ListPop(S);
      X = Base(Blossom(EdgeFrom(E)));
      Y = Base(Blossom(EdgeTo(E)));
      if (X == Y)
         continue;
      if (!IsEven(X))
      {
         Z = X;
         X = Y;
         Y = Z;
      }
         
      if (!IsReached(Y) && !IsMatched(Y))
      {
         Label(Y) = OddLabel;
         Tree(Y) = E;
         Age(Y) = Time++;
         ListPut(Y, T);
            
         U = Recover(Y);
         ListDelete(Self(Y));
               
         Found = True;
         break;
      }
         
      else if (!IsReached(Y) && IsMatched(Y))
      {
         Label(Y) = OddLabel;
         Tree(Y) = E;
         Age(Y) = Time++;
         ListPut(Y, T);
               
         F = Match(Y);
         Z = Other(F, Y);
         Label(Z) = EvenLabel;
         Age(Z) = Time++;
         ListPut(Z, T);

         ForAllIncidentEdges(Z, E, C, D)
            if (E != F)
            {
               ListPush(E, S);
               
               W = Other(E, Z);
               if (!IsReached(W) && !IsMatched(W))
                  break;
            }
      }
         
      else if (IsEven(Y))
         Shrink(E, &S);
   }
   DestroyList(S);

   if (!Found)
   {
      #ifndef NDEBUG
      *P = U;
      *Q = T;
      #else
      DestroyList(T);
      *P = Nil;
      *Q = Nil;
      #endif
   }
   else
   {
      *P = U;
      *Q = T;
   }
   
   return Found;
}


/*
 * Recover -- Recover an augmenting path ending at vertex V by walking
 *            up the tree back to the root.
 *
 * Returns a list of the unmatched edges on the path.
 *
 */
static List *Recover

#ifdef Ansi
   (register Vertex *V)
#else
   (V) register Vertex *V;
#endif

{
   register Vertex *W, *B;
   register List   *P;

   P = CreateList();
   
   do
   {
      ListPut(Tree(V), P);

      W = Other(Tree(V), V);
      B = Base(Blossom(W));
      Path(W, B, P);

      V = Other(Match(B), B);
   }
   while (V != Nil);

   return P;
}


/*
 * Path -- Recursively recover the even-length piece of an alternating path
 *         that begins at vertex V with a matched edge and ends at base B
 *         of its blossom
 *
 * The unmatched edges on the path are added to list P, and are in arbitrary
 * order.
 *
 */
static Void Path

#ifdef Ansi
   (register Vertex *V, Vertex *B, List *P)
#else
   (V, B, P) register Vertex *V, *B; List *P;
#endif

{
   register Vertex *W;
  
   if (V != B)
      if (IsOdd(V))
      {
         Path(Shore(V), Other(Match(V), V), P);
         ListPut(Bridge(V), P);
         Path(Other(Bridge(V), Shore(V)), B, P);
      }
      else if (IsEven(V))
      {
         W = Other(Match(V), V);
         ListPut(Tree(W), P);
         Path(Other(Tree(W), W), B, P);
      }
      else
         Error("(Path) Internal error.");
}


/*
 * Shrink -- Given an edge E between two even blossoms, shrink the implied
 *           cycle in the alternating tree into a superblossom
 *
 * Edges incident to odd vertices on the blossom are added to the stack S
 * of search edges.
 *
 */
static Void Shrink

#ifdef Ansi
   (Edge *E, List **S)
#else
   (E, S) Edge *E; List **S;
#endif

{
   auto     short    Found;
   register Vertex  *V, *W;
   register Vertex  *A, *B;
   register Element *X, *Y;
  
   V = EdgeFrom(E);
   W = EdgeTo(E);
   X = Blossom(V);
   Y = Blossom(W);
   B = Base(X);
   A = Base(Y);
   if (Age(A) > Age(B))
   {
      Vertex  *C;
      Element *Z;
      
      C = A;
      A = B;
      B = C;
      
      C = V;
      V = W;
      W = C;

      Z = X;
      X = Y;
      Y = Z;
   }
   
   /*
    * Walk up the alternating tree from vertex V to vertex A, shrinking
    * the blossoms into a superblossom.  Edges incident to the odd vertices
    * on the path from V to A are pushed onto stack S, to later search from.
    */
   Found = False;
   while (B != A)
   {
      Cell   *C, *D;
      Vertex *Z;
      Edge   *F, *M, *T;
  
      M = Match(B);
      W = Other(M, B);
      Bridge(W) = E;
      Shore(W) = V;
      
      T = Tree(W);
      if (!Found)
         ForAllIncidentEdges(W, F, C, D)
            if (F != M && F != T)
            {
               ListPush(F, *S);
               
               Z = Other(F, W);
               if (!IsReached(Z) && !IsMatched(Z))
               {
                  Found = True;
                  break;
               }
            }
      
      Y = Blossom(W);
      X = SetUnion(Y, X);
      E = T;
      V = Other(E, W);
      
      Y = Blossom(V);
      X = SetUnion(Y, X);
      B = Base(X);
   }
}


/*
 * Augment -- Augment the matching along augmenting path P, and expand
 *            into singleton sets all original vertices in T
 *
 * This assumes list P contains only the unmatched edges on the path,
 * and that list T contains all vertices in all blossoms in the alternating
 * tree containing the augmenting path.
 *
 */
static Void Augment

#ifdef Ansi
   (List *P, List *T)
#else
   (P, T) List *P, *T;
#endif

{
   register Vertex *V;
   register Edge   *E;
   register Cell   *C;
   
   ForAllEdges(E, P, C)
   {
      Match(EdgeFrom(E)) = E;
      Match(EdgeTo(E)) = E;
   }

   ForAllVertices(V, T, C)
   {
      Element      *B;
      SetAttribute *L;
      
      Label(V) = UnreachedLabel;
      
      L = OriginalSetLabel(V);
      DestroyElement(Blossom(V));
      B = CreateElement(L);
      Blossom(V) = B;
      Base(B) = V;
   }

   DestroyList(P);
   DestroyList(T);
}


/*
 * Matching -- Recover the final matching from the vertex match fields
 *
 */
static List *Matching

#ifdef Ansi
   (Graph *G)
#else
   (G) Graph *G;
#endif

{
   Vertex *V;
   Edge   *E;
   Cell   *P;
   List   *M;
  
   M = CreateList();
   
   ForAllGraphVertices(V, G, P)
   {
      E = Match(V);
      if (E != Nil && V == EdgeFrom(E))
         ListPut(E, M);
   }
   
   return M;
}


/*
 * Error -- Print an error message and halt
 *
 */
static Void Error

#ifdef Ansi
   (char *message)
#else
   (message) char *message;
#endif

{
   fprintf(stderr, "%s\n", message);
   Halt();
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Computing a maximal matching
 *
 */


typedef struct
{
   int    DegreeField;
   Cell  *EntryField;
   VertexData OriginalVertexLabelField;
}
MaximalVertexAttribute;
   

#define Degree(V) (((MaximalVertexAttribute *) VertexLabel(V))->DegreeField)
#define Entry(V)  (((MaximalVertexAttribute *) VertexLabel(V))->EntryField)

#define OriginalMaximalVertexLabel(V) \
        (((MaximalVertexAttribute *) VertexLabel(V))->OriginalVertexLabelField)


/*
 * MaximalMatching -- Compute a maximal matching of a nonbipartite graph
 *
 * A matching is maximal if no edge can be added to it.  This implementation
 * runs in O(m + n) time for a graph with m edges and n vertices.  A maximal
 * matching always has at least half the number of edges in a maximum
 * cardinality matching.
 *
 * The maximal matching that is returned is a greedy matching in the sense that
 * an edge is always added that is lexicographically minimum with respect to
 * the degrees of the two vertices touched by an edge, where degree is with
 * respect to the vertex subgraph induced by unmatched vertices.
 *
 */
List *MaximalMatching

#ifdef Ansi
   (Graph *G)
#else
   (G) Graph *G;
#endif

{
   MaximalVertexAttribute *VertexAttributes, *A;
   
   register int     I, D;
   auto     int     N;
   register Cell   *P, *Q;
   auto     List   *M;
   register List  **Heap;
   register Vertex *U, *V, *W;
   register Edge   *E, *F;

   
   N = ListSize(GraphVertices(G));
   VertexAttributes = (MaximalVertexAttribute *)
      Allocate(N * sizeof(MaximalVertexAttribute));
   Heap = (List **) Allocate((N - 1) * sizeof(List *));
   if (VertexAttributes == NULL || Heap == NULL)
      Error("(MaximalMatching) Memory allocation failed.");
   Heap -= 1;

   A = VertexAttributes;
   ForAllGraphVertices(V, G, P)
   {
      VertexData D;

      D = VertexRelabel(V, A);
      OriginalMaximalVertexLabel(V) = D;
      A++;
   }
   
   ForAllGraphVertices(V, G, P)
      Degree(V) = 0;
   ForAllGraphEdges(E, G, P)
   {
      Degree(EdgeFrom(E)) += 1;
      Degree(EdgeTo(E)) += 1;
   }
   
   for (I = 1; I < N; I++)
      Heap[I] = CreateList();
   D = N;
   ForAllGraphVertices(V, G, P)
   {
      I = Degree(V);
      if (I > 0)
      {
         Entry(V) = ListPut(V, Heap[I]);
         if (I < D)
            D = I;
      }
      else
         Entry(V) = Nil;
   }

   M = CreateList();
   for (;;)
   {
      V = Nil;
      for (; D < N; D++)
         if (!ListIsEmpty(Heap[D]))
         {
            V = (Vertex *) ListFront(Heap[D]);
            break;
         }
      if (!V)
         break;

      I = N;
      ForAllIncidentEdges(V, F, P, Q)
      {
         U = Other(F, V);
         if (Entry(U) && Degree(U) < I)
         {
            E = F;
            I = Degree(U);
         }
      }
      ListPut(E, M);

      W = Other(E, V);
      ListDelete(Entry(V));
      ListDelete(Entry(W));
      Entry(V) = Nil;
      Entry(W) = Nil;

      ForAllIncidentEdges(V, F, P, Q)
      {
         U = Other(F, V);
         if (Entry(U))
         {
            ListDelete(Entry(U));
            I = (Degree(U) -= 1);
            if (I > 0)
            {
               Entry(U) = ListPut(U, Heap[I]);
               if (I < D)
                  D = I;
            }
            else
               Entry(U) = Nil;
         }
      }
      ForAllIncidentEdges(W, F, P, Q)
      {
         U = Other(F, W);
         if (Entry(U))
         {
            ListDelete(Entry(U));
            I = (Degree(U) -= 1);
            if (I > 0)
            {
               Entry(U) = ListPut(U, Heap[I]);
               if (I < D)
                  D = I;
            }
            else
               Entry(U) = Nil;
         }
      }
   }

   for (I = 1; I < N; I++)
      DestroyList(Heap[I]);
   Heap += 1;
   Free(Heap);
   ForAllGraphVertices(V, G, P)
      VertexRelabel(V, OriginalMaximalVertexLabel(V));
   Free(VertexAttributes);
   
   return M;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Debugging scaffolding
 *
 */


#ifdef Debug


/*
 * DumpAlternatingForest -- Write to stderr a representation of the
 *                          alternating forest
 *
 */
static Void DumpAlternatingForest

#ifdef Ansi
   (Void)
#else
   ()
#endif

{
   Vertex *V;
   Cell   *P;
   
   ForAllGraphVertices(V, UnderlyingGraph, P)
   {
      if (V == Base(Blossom(V)))
         Members(V) = CreateList();
      Children(V) = CreateList();
   }
   
   ForAllGraphVertices(V, UnderlyingGraph, P)
   {
      ListPut(V, Members(Base(Blossom(V))));
      if (V == Base(Blossom(V)) && IsReached(V))
         if (IsEven(V))
         {
            if (IsMatched(V))
               ListPut(Match(V), Children(Other(Match(V), V)));
         }
         else /* IsOdd(V) */
         {
            ListPut(Tree(V),
                    Children(Base(Blossom(Other(Tree(V), V)))));
         }
   }       
   fprintf(stderr, "Alternating forest\n");
   ForAllGraphVertices(V, UnderlyingGraph, P)
      if (V == Base(Blossom(V)) && IsEven(V) && !IsMatched(V))
         Traverse(V, Nil, 0);
   fflush(stderr);

   ForAllGraphVertices(V, UnderlyingGraph, P)
   {
      if (V == Base(Blossom(V)))
         DestroyList(Members(V));
      DestroyList(Children(V));
   }
}

/*
 * Traverse -- Preorder traversal of a subtree of the alternating tree
 *
 */
static Void Traverse

#ifdef Ansi
   (Vertex *V, Edge *E, int D)
#else
   (V, E, D) Vertex *V; Edge *E; int D;
#endif

{
   Vertex *W, *A, *B;
   Cell   *P;
   int     i;
   
   for (i = 1; i <= D; i++)
      fprintf(stderr, "%s", " ");
   fprintf(stderr, "{");
   ForAllVertices(W, Members(Base(Blossom(V))), P)
      fprintf(stderr, "%s%d",
              W == ListFront(Members(Base(Blossom(V)))) ? "" : " ",
              Name(W));
   fprintf(stderr, "}");
   fprintf(stderr, " %s", IsEven(V) ? "even" : "odd");
   if (E)
      fprintf(stderr, " (%d, %d)", Name(EdgeFrom(E)), Name(EdgeTo(E)));
   fprintf(stderr, "\n");
   
   ForAllEdges(E, Children(V), P)
   {
      A = Base(Blossom(EdgeFrom(E)));
      B = Base(Blossom(EdgeTo(E)));
      if (A != V)
         Traverse(A, E, D + 1);
      else
         Traverse(B, E, D + 1);
   }
}


#endif /* Debug */
