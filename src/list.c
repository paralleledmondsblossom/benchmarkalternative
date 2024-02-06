/*
 * list.c -- Linear lists
 */

/*
 * Copyright 1989, 1992 by John Kececioglu
 */


/*
 * Synopsis
 *
 * The lists are circular, doubly-linked, and anchored.  This allows O(1) time
 * insertion, deletion, concatenation, and destruction.
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
 * 26 July 1991 JDK
 * Simplified the method of iterating over a list.  Replaced three iteration
 * functions with two macros.
 *
 * 22 June 1995 JDK
 * Fixed a bug in that `ListInsertAfter' did not return a pointer to the
 * inserted list cell.
 *
 */


#include <stdio.h>
#include "list.h"


typedef ListData Item;
typedef ListCell Cell;


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Creation and destruction
 *
 */


#define CellBlockSize 32 /* Number of cells allocated per memory request */


static Cell *CellPool = Nil; /* Pool of free cells */


/*
 * Cell pool maintenance
 */
#define FreeCell(C) (((C)->Next = CellPool), CellPool = (C))
#define NewCell(C)  (((C) = CellPool), CellPool = CellPool->Next)


/*
 * CreateList -- Create an empty list
 *
 */
List *CreateList
   
#ifdef Ansi
   (Void)
#else
   ()
#endif

{
   register List *L;
   register Cell *C, *Block;
   
   if (CellPool == Nil)
   {
      /*
       * Allocate a block of cells
       */
      Block = (Cell *) Allocate(sizeof(Cell) * CellBlockSize);
      if (Block == NULL)
      {
         fprintf(stderr, "(CreateList) Memory allocation failed.\n");
         Halt();
      }
      
      /*
       * Place the cells in the block into the pool
       */
      for (C = Block; C - Block < CellBlockSize; C++)
         FreeCell(C);
   }
   
   NewCell(L);
   L->Next = L;
   L->Prev = L;
   L->Item = Nil;
      /*
       * A list anchor has the `Nil' item so that `ListNext' and `ListPrev'
       * return `Nil' to signal completion of a list traversal.
       */
   
   return L;
}


/*
 * DestroyList -- Destroy a list
 *
 */
Void DestroyList
   
#ifdef Ansi
   (register List *L)
#else
   (L) register List *L;
#endif

{
   L->Prev->Next = CellPool;
   CellPool = L;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Arbitrary access
 *
 */


/*
 * ListInsertBefore -- Insert item into list before given cell
 *
 */
Cell *ListInsertBefore
   
#ifdef Ansi
   (Item I, Cell *C)
#else
   (I, C) Item I; Cell *C;
#endif

{
   register Cell *P, *Block;
   
   if (CellPool == Nil)
   {
      /*
       * Allocate a block of cells
       */
      Block = (Cell *) Allocate(sizeof(Cell) * CellBlockSize);
      if (Block == NULL)
      {
         fprintf(stderr, "(ListInsert) Memory allocation failed.\n");
         Halt();
      }
      
      /*
       * Place the cells in the block into the pool
       */
      for (P = Block; P - Block < CellBlockSize; P++)
         FreeCell(P);
   }

   NewCell(P);
   P->Item = I;
   P->Prev = C->Prev;
   P->Next = C;
   P->Prev->Next = P;
   P->Next->Prev = P;
   
   return P;
}


/*
 * ListInsertAfter -- Insert item into list after given cell
 *
 */
Cell *ListInsertAfter

#ifdef Ansi
   (Item I, Cell *C)
#else
   (I, C) Item I; Cell *C;
#endif

{
   return ListInsertBefore(I, C->Next);
}


/*
 * ListDelete -- Delete cell from list
 *
 * Returns the item of the deleted cell.
 *
 */
Item ListDelete
   
#ifdef Ansi
   (register Cell *C)
#else
   (C) register Cell *C;
#endif

{
   C->Prev->Next = C->Next;
   C->Next->Prev = C->Prev;
   FreeCell(C);
   
   return C->Item;
}


/*
 * ListEject -- Eject cell from list
 *
 */
Void ListEject
   
#ifdef Ansi
   (register Cell *C)
#else
   (C) register Cell *C;
#endif

{
   C->Prev->Next = C->Next;
   C->Next->Prev = C->Prev;
}


/*
 * ListInject -- Inject ejected cell back into list
 *
 * Adjacent ejected cells must be injected in reverse order.
 *
 */
Void ListInject
   
#ifdef Ansi
   (register Cell *C)
#else
   (C) register Cell *C;
#endif

{
   C->Prev->Next = C;
   C->Next->Prev = C;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Stacks, queues, and deques
 *
 */


/*
 * ListPush -- Push item onto front of list
 *
 */
Cell *ListPush
   
#ifdef Ansi
   (Item I, List *L)
#else
   (I, L) Item I; List *L;
#endif

{
   return ListInsertBefore(I, L->Next);
}


/*
 * ListPop -- Pop item off front of list
 *
 */
Item ListPop
   
#ifdef Ansi
   (List *L)
#else
   (L) List *L;
#endif

{
   return !ListIsEmpty(L) ? ListDelete(L->Next) : Nil;
}


/*
 * ListPut -- Put item onto rear of list
 *
 */
Cell *ListPut
   
#ifdef Ansi
   (Item I, List *L)
#else
   (I, L) Item I; List *L;
#endif

{
   return ListInsertAfter(I, L->Prev);
}


/*
 * ListPull -- Pull item off rear of list
 *
 */
Item ListPull
   
#ifdef Ansi
   (register List *L)
#else
   (L) register List *L;
#endif

{
   return !ListIsEmpty(L) ? ListDelete(L->Prev) : Nil;
}


/*
 * ListFront -- Item on front of list
 *
 */
Item ListFront
   
#ifdef Ansi
   (List *L)
#else
   (L) List *L;
#endif

{
   return !ListIsEmpty(L) ? L->Next->Item : Nil;
}


/*
 * ListRear -- Item at rear of list
 *
 */
Item ListRear
   
#ifdef Ansi
   (List *L)
#else
   (L) List *L;
#endif

{
   return !ListIsEmpty(L) ? L->Prev->Item : Nil;
}


/*
 * ListIsEmpty -- Is the list empty?
 *
 */
int ListIsEmpty
   
#ifdef Ansi
   (List *L)
#else
   (L) List *L;
#endif

{
   return L->Next == L;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Iteration
 *
 */


/*
 * The file `list.h' provides macros to iterate over lists.  A generic loop
 * using these macros looks like,
 *
 *    List *L; ListCell *P; ListData D;
 *
 *    P = ListHead(L);
 *    while (D = ListNext(P))
 *       ;
 *
 * `ListNext' returns the value `Nil' on termination.  Consequently, no list
 * can contain the value `Nil'.
 *
 */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Copying, concatenation, reversal, and size
 *
 */


/*
 * ListCopy -- Create a copy of a list
 *
 */
List *ListCopy
   
#ifdef Ansi
   (List *L)
#else
   (L) List *L;
#endif

{
   register List *C;
   register Cell *P;
   register Item X;

   C = CreateList();
   
   P = ListHead(L);
   while ((X = ListNext(P)))
      ListPut(X, C);
   
   return C;
}


/*
 * ListCat -- Concatenate the second list onto the first, destructively
 *
 */
List *ListCat
   
#ifdef Ansi
   (register List *A, register List *B)
#else
   (A, B) register List *A, *B;
#endif

{
   if (!ListIsEmpty(B))
   {
      A->Prev->Next = B->Next;
      B->Next->Prev = A->Prev;
      
      A->Prev = B->Prev;
      B->Prev->Next = A;
      
      B->Next = B->Prev = B;
   }
   
   return A;
}


/*
 * ListReverse -- Reverse a list destructively
 *
 * This is a linear time operation.
 *
 */
List *ListReverse
   
#ifdef Ansi
   (register List *A)
#else
   (A) register List *A;
#endif

{
   register List *B;
   register Item X;

   
   B = CreateList();
   while ((X = ListPop(A)))
      ListPush(X, B);
   
   ListCat(A, B);
   DestroyList(B);
   
   return A;
}


/*
 * ListSize -- Size of list
 *
 * This is a linear time operation.
 *
 */
int ListSize
   
#ifdef Ansi
   (register List *L)
#else
   (L) register List *L;
#endif

{
   register Cell *P;
   register int n;

   n = 0;
   for (P = ListHead(L); ListNext(P); )
      n += 1;
   return n;
}

