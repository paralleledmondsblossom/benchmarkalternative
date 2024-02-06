
/*
---------------------------------------------------------------------
 This file is part of matchmaker framework
 Copyright (c) 2012,
 By:    Mehmet Deveci,
        Kamer Kaya,
        Bora Ucar,
        Umit V. Catalyurek
---------------------------------------------------------------------
 For license info, please see the LICENSE.txt file in
 the main directory.
---------------------------------------------------------------------
*/



======
matchmaker
======

- To compile the code, make sure you have a version of gcc and nvcc
    $> make


- To run give the following arguments:
mmaker [argument_list]
Argument list
	IF: path to the input graph file (mandatory argument)
	IT: input file format. 0-chaco , 1-matrix market (default)
	MT: Maximum Matching Type
		0: Sequential DFS
		1: Sequential BFS
		2: Sequential PF
		3: Sequential PFP
		4: Sequential HK
		5: Sequential HK_DW
		6: Sequential ABMP
		7: Sequential ABMP_BFS
		8: GPU - APFB_GPUBFS
		9: GPU - APFB_GPUBFS_WR
		10: GPU - APsB_GPUBFS
		11: GPU - APsB_GPUBFS_WR
	IMT: Initial Matching Type
		0: 
		0: Cheap Matching
		1: SK
		2: SK_rand
		3: mind_cheap
		>=4: no initial matching
	R: Repetition number. Must be greater than 0.
	BD: Block Dimension (Grid Size). Must be greater than 0. 
	TD: Thread Dimension (Block Size). Must be greater than 0. 
	(if both BD and TD are not provided, then total thread number will be maximized) 

- For input graph, Chaco and Matrix Market formats are supported.

- Example run :
     $> ./mmaker IF=~/work/input/mm/Hamrle3.mtx R=1 MT=8 BD=256 TD=256

- For any question or problem, please contact:
    mdeveci@bmi.osu.edu
    kamer@bmi.osu.edu
    bora.ucar@ens-lyon.fr
    umit@bmi.osu.edu
