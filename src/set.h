/*
 * disjointset.h -- Disjoint set definitions
 */

/*
 * Copyright 1989 by John Kececioglu
 */


#ifndef SetInclude
#define SetInclude


#include "portable.h"


/*
 * Nil pointer
 */
#ifndef Nil
#define Nil 0
#endif

/*
 * Data associated with a set
 */
typedef Pointer SetData;

/*
 * Element of a set
 */
typedef struct ElementStruct {
   SetData Label;
   int Rank;
   struct ElementStruct *Up;
      /*
       * `Up' is reused for the pool of free elements
       */
} Element;


extern Element *CreateElement Proto(( SetData D ));
extern Void    DestroyElement Proto(( Element *E ));

extern Element *SetFind  Proto(( Element *E ));
extern Element *SetUnion Proto(( Element *E, Element *F ));
   
extern SetData SetLabel Proto(( Element *E ));


#endif /* SetInclude */
