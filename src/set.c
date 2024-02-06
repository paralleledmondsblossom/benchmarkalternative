/*
 * set.c -- Disjoint sets
 */

/*
 * Copyright 1989 by John Kececioglu
 */


/*
 * Synopsis
 *
 * This is an implementation of Tarjan's disjoint set data structure,
 * with union by rank, and path halving.  A sequence of n union and
 * m find operations on a universe of n elements takes O(n + m A(m,n)) time,
 * where A(m,n) is an inverse of Ackermann's function.  For most conceivable
 * applications, A(m,n) is a constant.
 *
 * See the book, Robert Endre Tarjan, Data Structures and Network Algorithms,
 * Society for Industrial and Applied Mathematics, Philadelphia, 23-31, 1983.
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
 * 7 March 1997 JDK
 * Restructured the interface to use only an element type, and not a set type.
 * A call to Union now passes two arbitrary elements, and the sets containing
 * them are unioned if they are disjoint.
 *
 * 26 June 1994 JDK
 * Made the naming of functions and structures consistent with other libraries.
 *
 */


#include <stdio.h>
#include "set.h"


#define BlockSize 8 /* Number of elements allocated per memory request */


typedef SetData Data;


static Element *Pool = Nil; /* Pool of free elements */


/*
 * Element pool maintenance
 */
#define FreeElement(E) (((E)->Up = Pool), Pool = (E))
#define NewElement(E)  (((E) = Pool), Pool = Pool->Up)


/*
 * CreateElement -- Create a set element
 *
 * The element is the member of an implicit singleton set.  The data supplied
 * on creation labels this set.
 *
 */
Element *CreateElement
   
#ifdef Ansi
   (Data D)
#else
   (D) Data D;
#endif

{
   register Element *E;
	
   if (Pool == Nil)
   {
      register Element *Block;
      
      /*
       * Allocate a block of elements
       */
      Block = (Element *) Allocate(sizeof(Element) * BlockSize);
      if (Block == NULL)
      {
         fprintf(stderr, "(CreateElement) Memory allocation failed.\n");
         Halt();
      }

      /*
       * Place the elements in the block into the pool
       */
      for (E = Block; E - Block < BlockSize; E++)
         FreeElement(E);
   }
   
   NewElement(E);
   E->Up = E;
   E->Rank = 0;
   E->Label = D;
   
   return E;
}


/*
 * DestroyElement -- Destroy a set element
 *
 * Destruction of an element also destroys the implicit set containing it,
 * but not the other members of the set.
 *
 */
Void DestroyElement

#ifdef Ansi
   (register Element *E)
#else
   (E) register Element *E;
#endif

{
   FreeElement(E);
}


/*
 * SetUnion -- Destructive union of two disjoint sets
 *
 * The sets containing the two specified elements are unioned if they
 * are disjoint, and a representative element of the unioned set is returned.
 *
 * The data labeling the unioned set is the data labeling the set containing
 * the first argument (i.e. the second set is destructively merged into the
 * first).
 *
 */
Element *SetUnion

#ifdef Ansi
   (register Element *E, register Element *F)
#else
   (E, F) register Element *E, *F;
#endif

{
   register Data D;
   
   E = SetFind(E);
   F = SetFind(F);
   if (E == F)
      return E;
   
   D = E->Label;
   if (E->Rank < F->Rank)
   {
      E->Up = F;
      F->Label = D;
      return F;
   }
   else if (E->Rank > F->Rank)
   {
      F->Up = E;
      return E;
   }
   else
   {
      E->Rank += 1;
      F->Up = E;
      return E;
   }
}


/*
 * SetFind -- Find the set containing a given element
 *
 * A representative element of the containing set is returned.
 *
 */
Element *SetFind

#ifdef Ansi
   (register Element *E)
#else
   (E) register Element *E;
#endif

{
   while (E->Up->Up != E->Up)
      E = E->Up = E->Up->Up;
   return E->Up;
}


/*
 * SetLabel -- Return the data labeling a set
 *
 * The set is the one containing the specified element.
 *
 */
Data SetLabel

#ifdef Ansi
   (register Element *E)
#else
   (E) register Element *E;
#endif

{
   return SetFind(E)->Label;
}
