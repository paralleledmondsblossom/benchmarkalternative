/*
 * list.h -- Linear list definitions
 */

/*
 * Copyright 1989, 1992 by John Kececioglu
 */


#ifndef ListInclude
#define ListInclude


#include "portable.h"


/*
 * Nil pointer
 */
#ifndef Nil
#define Nil 0
#endif

/*
 * List item
 */
typedef Pointer ListData;

/*
 * List cell
 */
typedef struct ListStruct {
   ListData Item;
   struct ListStruct *Prev;
   struct ListStruct *Next;
      /*
       * `Next' is reused for the pool of free cells
       */
} List, ListCell;


/*
 * Creation and destruction
 */
extern List *CreateList Proto(( Void ));
extern Void DestroyList Proto(( List *L ));

   
/*
 * Arbitrary access
 */
extern ListCell *ListInsertBefore Proto(( ListData X, List *L ));
extern ListCell *ListInsertAfter  Proto(( ListData X, List *L ));
extern ListData ListDelete        Proto(( ListCell *P ));
extern Void     ListInject        Proto(( ListCell *P ));
extern Void     ListEject         Proto(( ListCell *P ));


/*
 * Stacks
 */
extern ListCell *ListPush Proto(( ListData X, List *L ));
extern ListData ListPop   Proto(( List *L ));
   
extern int ListIsEmpty Proto(( List *L ));

#define ListTop(L) ListFront(L)

   
/*
 * Queues and deques
 */
extern ListCell *ListPut  Proto(( ListData X, List *L ));
extern ListData ListPull  Proto(( List *L ));
extern ListData ListFront Proto(( List *L ));
extern ListData ListRear  Proto(( List *L ));

#define ListGet(L) ListPop(L)

   
/*
 * Copying, concatenation, reversal, and size
 */
extern List *ListCopy    Proto(( List *L ));
extern List *ListCat     Proto(( List *A, List *B ));
extern List *ListReverse Proto(( List *L ));
extern int  ListSize     Proto(( List *L ));

   
/*
 * Iteration
 */
#define ListHead(List) ((List)->Next)
#define ListTail(List) ((List)->Prev)
#define ListItem(Cell) ((Cell)->Item)
#define ListNext(Cell) (((Cell) = (Cell)->Next), (Cell)->Prev->Item)
#define ListPrev(Cell) (((Cell) = (Cell)->Prev), (Cell)->Next->Item)

#define ForAllListElements(Variable, List, Type, Cell) \
        for((Cell) = ListHead(List); \
            ((Variable) = (Type) ListItem(Cell)); \
            (Void) ListNext(Cell))


#endif /* ListInclude */
