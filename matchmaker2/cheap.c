/*
 * File: cheap.c
 * Content: Contains jump-start heuristics
 * Author: Kamer Kaya, Johannes Langguth and Bora Ucar
 * Version: 0.3
 *
 * Please see the papers:
 *
 *  "Iain S. Duff, Kamer Kaya, and Bora Ucar, 'Design, analysis, and
    implementation of maximum transversal algorithms', ACM Transactions 
    on Mathematical Software, 38, 2, Article 13, (2011)."
 *
 * 	"Kamer Kaya, Johannes Langguth, Frederik Manne, and Bora Ucar, 
    'Push-relabel based algorithms for the maximum transversal', 
     Computers and Operations Research, 2013."
 *
 * for more details and cite them if you use the codes.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "matchmaker.h"

struct node {int id; int degree; struct node *next; struct node* prvs;};
typedef struct node Node;

void old_cheap(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m) {
	int ptr;
	int i = 0;
	for(; i < n; i++) {
		int s_ptr = col_ptrs[i];
		int e_ptr = col_ptrs[i + 1];
		for(ptr = s_ptr; ptr < e_ptr; ptr++) {
            int r_id = col_ids[ptr];
			if(row_match[r_id] == -1) {
				match[i] = r_id;
				row_match[r_id] = i;
				break;
			}
		}
	}
}

void sk_cheap(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m){
	int i;

	int* col_stack = (int*)malloc(n * sizeof(int));
	int* col_degrees = (int*)malloc(n * sizeof(int));
	memset(col_degrees, 0, n * sizeof(int));

	int no_of_d1_cols = 0;
	for(i = 0; i < n; i++) {
		col_degrees[i] = col_ptrs[i+1] - col_ptrs[i];
		if(col_degrees[i] == 1) {
			col_stack[no_of_d1_cols++] = i;
		}
	}

	int* row_stack = (int*)malloc(m * sizeof(int));
	int* row_degrees = (int*)malloc(m * sizeof(int));
	memset(row_degrees, 0, m * sizeof(int));

	int no_of_d1_rows = 0;
	for(i = 0; i < m; i++) {
		row_degrees[i] = row_ptrs[i+1] - row_ptrs[i];
		if(row_degrees[i] == 1) {
			row_stack[no_of_d1_rows++] = i;
		}
	}

	int stop = 0;
	int r_id = -1, c_id, r_id2, c_id2;
	int sptr, eptr, ptr;
	int sptr2, eptr2, ptr2;

	int remain = 0;
	int c_degree = 0;

	while(!stop) {
		while(no_of_d1_rows > 0 || no_of_d1_cols > 0) {
			if(no_of_d1_rows > 0) {
				r_id = row_stack[--no_of_d1_rows];
				if(row_degrees[r_id] == 1 && row_match[r_id] == -1) {
					sptr = row_ptrs[r_id];
					eptr = row_ptrs[r_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						c_id = row_ids[ptr];
						if(match[c_id] == -1) {
							match[c_id] = r_id;
							row_match[r_id] = c_id;

							sptr2 = col_ptrs[c_id];
							eptr2 = col_ptrs[c_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								r_id2 = col_ids[ptr2];
								if(row_match[r_id2] == -1) {
									if((--(row_degrees[r_id2])) == 1) {
										row_stack[no_of_d1_rows++] = r_id2;
									}
								}
							}
							break;
						}
					}
				}
			}

			if(no_of_d1_cols > 0) {
				c_id = col_stack[--no_of_d1_cols];
				if(col_degrees[c_id] == 1 && match[c_id] == -1) {
					sptr = col_ptrs[c_id];
					eptr = col_ptrs[c_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						r_id = col_ids[ptr];
						if(row_match[r_id] == -1) {
							row_match[r_id] = c_id;
							match[c_id] = r_id;

							sptr2 = row_ptrs[r_id];
							eptr2 = row_ptrs[r_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								c_id2 = row_ids[ptr2];
								if( match[c_id2] == -1) {
									if((--(col_degrees[c_id2])) == 1) {
										col_stack[no_of_d1_cols++] = c_id2;
									}
								}
							}
							break;
						}
					}
				}
			}
		}

		stop = 1;
		for(i = remain; i < n; i++) {
			c_id = i;
			c_degree = col_degrees[c_id];

			if(match[c_id] == -1 && c_degree != 0) {
				sptr = col_ptrs[c_id];
				eptr = col_ptrs[c_id + 1];

				for(ptr = sptr; ptr < eptr; ptr++) {
					r_id = col_ids[ptr];
					if(row_match[r_id] == -1) {
						match[c_id] = r_id;
						row_match[r_id] = c_id;
						stop = 0;
						break;
					}
				}
				ptr++;

				for(;ptr < eptr; ptr++) {
					r_id2 = col_ids[ptr];
					if(row_match[r_id2] == -1) {
						if((--(row_degrees[r_id2])) == 1) {
							row_stack[no_of_d1_rows++] = r_id2;
						}
					}
				}

				sptr = row_ptrs[r_id];
				eptr = row_ptrs[r_id + 1];
				int count = row_degrees[r_id];
				for(ptr = sptr;ptr < eptr && count > 0; ptr++) {
					c_id2 = row_ids[ptr];
					if(match[c_id2] == -1) {
						count--;
						if((--(col_degrees[c_id2])) == 1) {
							col_stack[no_of_d1_cols++] = c_id2;
						}
					}
				}
			}

			if(no_of_d1_cols + no_of_d1_rows > 0) {
				remain = i + 1;
				break;
			}

			if(i == n-1) {
				stop = 1;
			}
		}
	}

	free(row_degrees);
	free(row_stack);
	free(col_degrees);
	free(col_stack);
}

void sk_cheap_rand(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m) {
	int i;

	int* col_stack = (int*)malloc(n * sizeof(int));
	int* col_degrees = (int*)malloc(n * sizeof(int));
	memset(col_degrees, 0, n * sizeof(int));

	int no_of_d1_cols = 0;
	for(i = 0; i < n; i++) {
		col_degrees[i] = col_ptrs[i+1] - col_ptrs[i];
		if(col_degrees[i] == 1) {
			col_stack[no_of_d1_cols++] = i;
		}
	}

	int* row_stack = (int*)malloc(m * sizeof(int));
	int* row_degrees = (int*)malloc(m * sizeof(int));
	memset(row_degrees, 0, m * sizeof(int));

	int no_of_d1_rows = 0;
	for(i = 0; i < m; i++) {
		row_degrees[i] = row_ptrs[i+1] - row_ptrs[i];
		if(row_degrees[i] == 1) {
			row_stack[no_of_d1_rows++] = i;
		}
	}

	int* randarr = (int*)malloc(n * sizeof(int));
	for(i = 0; i < n; i++){randarr[i] = i;}
	int temp;
	for(i = n-1; i >= 0; i--) {
		int z = rand() % (i+1);
		temp = randarr[i]; randarr[i] = randarr[z]; randarr[z] = temp;
	}

	int stop = 0;
	int r_id = -1, c_id, r_id2, c_id2, e_id;
	int sptr, eptr, ptr;
	int sptr2, eptr2, ptr2;

	int remain = 0;
	int c_degree = 0;

	while(!stop) {
		while(no_of_d1_rows > 0 || no_of_d1_cols > 0) {
			if(no_of_d1_rows > 0) {
				r_id = row_stack[--no_of_d1_rows];
				if(row_degrees[r_id] == 1 && row_match[r_id] == -1) {
					sptr = row_ptrs[r_id];
					eptr = row_ptrs[r_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						c_id = row_ids[ptr];
						if(match[c_id] == -1) {
							match[c_id] = r_id;
							row_match[r_id] = c_id;

							sptr2 = col_ptrs[c_id];
							eptr2 = col_ptrs[c_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								r_id2 = col_ids[ptr2];
								if(row_match[r_id2] == -1) {
									if((--(row_degrees[r_id2])) == 1) {
										row_stack[no_of_d1_rows++] = r_id2;
									}
								}
							}
							break;
						}
					}
				}
			}

			if(no_of_d1_cols > 0) {
				c_id = col_stack[--no_of_d1_cols];
				if(col_degrees[c_id] == 1 && match[c_id] == -1) {
					sptr = col_ptrs[c_id];
					eptr = col_ptrs[c_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						r_id = col_ids[ptr];
						if(row_match[r_id] == -1) {
							row_match[r_id] = c_id;
							match[c_id] = r_id;

							sptr2 = row_ptrs[r_id];
							eptr2 = row_ptrs[r_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								c_id2 = row_ids[ptr2];
								if( match[c_id2] == -1) {
									if((--(col_degrees[c_id2])) == 1) {
										col_stack[no_of_d1_cols++] = c_id2;
									}
								}
							}
							break;
						}
					}
				}
			}
		}

		stop = 1;
		for(i = remain; i < n; i++) {
			c_id = randarr[i];
			c_degree = col_degrees[c_id];

			if(match[c_id] == -1 && c_degree != 0) {
				e_id = rand() % c_degree;

				sptr = col_ptrs[c_id];
				eptr = col_ptrs[c_id + 1];

				for(ptr = sptr; ptr < eptr; ptr++) {
					r_id = col_ids[ptr];
					if(row_match[r_id] == -1) {
						if(e_id == 0) {
							match[c_id] = r_id;
							row_match[r_id] = c_id;
							stop = 0;
							break;
						} else {
							if((--(row_degrees[r_id])) == 1) {
								row_stack[no_of_d1_rows++] = r_id;
							}
							e_id--;
						}
					}
				}
				ptr++;

				for(;ptr < eptr; ptr++) {
					r_id2 = col_ids[ptr];
					if(row_match[r_id2] == -1) {
						if((--(row_degrees[r_id2])) == 1) {
							row_stack[no_of_d1_rows++] = r_id2;
						}
					}
				}

				sptr = row_ptrs[r_id];
				eptr = row_ptrs[r_id + 1];
				int count = row_degrees[r_id];
				for(ptr = sptr;ptr < eptr && count > 0; ptr++) {
					c_id2 = row_ids[ptr];
					if(match[c_id2] == -1) {
						count--;
						if((--(col_degrees[c_id2])) == 1) {
							col_stack[no_of_d1_cols++] = c_id2;
						}
					}
				}
			}

			if(no_of_d1_cols + no_of_d1_rows > 0) {
				remain = i + 1;
				break;
			}

			if(i == n-1) {
				stop = 1;
			}
		}
	}

	free(randarr);
	free(row_degrees);
	free(row_stack);
	free(col_degrees);
	free(col_stack);
}


void mind_cheap(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m) {

	Node* rnodes = (Node*)malloc(sizeof(Node) * m);
	Node* cnodes = (Node*)malloc(sizeof(Node) * n);;
	Node* tptr;

	int i, deg, maxdeg = -1, cdeg, vtx, minnbr = -1, ptr, row ,col, temp;

	for(i = 0; i < n; i++) {
		deg = col_ptrs[i+1] - col_ptrs[i];
		cnodes[i].degree = deg;
		cnodes[i].id = i;
		if(deg > maxdeg) maxdeg = deg;
	}

	for(i = 0; i < m; i++) {
		deg = row_ptrs[i+1] - row_ptrs[i];
		rnodes[i].degree = deg;
		rnodes[i].id = i + n;
		if(deg > maxdeg) maxdeg = deg;
	}

	Node* lists = (Node*)malloc(sizeof(Node) * (maxdeg + 1));
	Node* listse = (Node*)malloc(sizeof(Node) * (maxdeg + 1));

	for(i = 0; i <= maxdeg; i++) {
		lists[i].next = &(listse[i]); lists[i].prvs = (Node*)0;
		listse[i].next = (Node*)0; listse[i].prvs = &(lists[i]);
		lists[i].id = -1; listse[i].id = -1;
		lists[i].degree = i; listse[i].degree = i;
	}

	for(i = 0; i < n; i++) {
		deg = cnodes[i].degree;
		if(deg > 0) {
			tptr = lists[deg].next;
			tptr->prvs = lists[deg].next = &(cnodes[i]);
			cnodes[i].next = tptr;
			cnodes[i].prvs = &(lists[deg]);
		}
	}
	for(i = 0; i < m; i++) {
		deg = rnodes[i].degree;
		if(deg > 0) {
			tptr = lists[deg].next;
			tptr->prvs = lists[deg].next = &(rnodes[i]);
			rnodes[i].next = tptr;
			rnodes[i].prvs = &(lists[deg]);
		}
	}

	cdeg = 1;
	while(cdeg <= maxdeg) {
		if(lists[cdeg].next == &(listse[cdeg])) {cdeg++; continue;}
		tptr = lists[cdeg].next;
		lists[cdeg].next = tptr->next;
		tptr->next->prvs = &(lists[cdeg]);
		vtx = tptr->id;

		if(vtx < n) {
			for(ptr = col_ptrs[vtx]; ptr < col_ptrs[vtx+1]; ptr++) {
				if(row_match[col_ids[ptr]] == -1) {
					minnbr = col_ids[ptr];
					break;
				}
			}

			for(ptr = ptr + 1; ptr < col_ptrs[vtx+1]; ptr++) {
				row = col_ids[ptr];
				if(row_match[row] == -1) {
					if(rnodes[row].degree < rnodes[minnbr].degree) {
						minnbr = col_ids[ptr];
					}
				}
			}

			match[vtx] = minnbr; row_match[minnbr] = vtx;
			rnodes[minnbr].next->prvs = rnodes[minnbr].prvs;
			rnodes[minnbr].prvs->next = rnodes[minnbr].next;
		} else {
			vtx = vtx - n;
			for(ptr = row_ptrs[vtx]; ptr < row_ptrs[vtx+1]; ptr++) {
				if(match[row_ids[ptr]] == -1) {
					minnbr = row_ids[ptr];
					break;
				}
			}

			for(ptr = ptr + 1; ptr < row_ptrs[vtx+1]; ptr++) {
				col = row_ids[ptr];
				if(match[col] == -1) {
					if(cnodes[col].degree < cnodes[minnbr].degree) {
						minnbr = row_ids[ptr];
					}
				}
			}

			row_match[vtx] = minnbr; match[minnbr] = vtx;
			cnodes[minnbr].next->prvs = cnodes[minnbr].prvs;
			cnodes[minnbr].prvs->next = cnodes[minnbr].next;
			temp = vtx; vtx = minnbr; minnbr = temp; /* swap */
		}

		for(ptr = col_ptrs[vtx]; ptr < col_ptrs[vtx+1]; ptr++) {
			row = col_ids[ptr];
			if(row_match[row] == -1) {
				deg = --(rnodes[row].degree);
				rnodes[row].next->prvs = rnodes[row].prvs;
				rnodes[row].prvs->next = rnodes[row].next;

				if(deg > 0) {
					tptr = lists[deg].next;
					tptr->prvs = lists[deg].next = &(rnodes[row]);
					rnodes[row].next = tptr;
					rnodes[row].prvs = &(lists[deg]);
				}
			}
		}

		for(ptr = row_ptrs[minnbr]; ptr < row_ptrs[minnbr+1]; ptr++) {
			col = row_ids[ptr];
			if(match[col] == -1) {
				deg = --(cnodes[col].degree);
				cnodes[col].next->prvs = cnodes[col].prvs;
				cnodes[col].prvs->next = cnodes[col].next;

				if(deg > 0) {
					tptr = lists[deg].next;
					tptr->prvs = lists[deg].next = &(cnodes[col]);
					cnodes[col].next = tptr;
					cnodes[col].prvs = &(lists[deg]);
				}
			}
		}
		cdeg--;
	}

	free(listse);
	free(lists);
	free(cnodes);
	free(rnodes);
}

void cheap_matching(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m, int cheap_id) {
	if(cheap_id == do_old_cheap) {
		old_cheap(col_ptrs, col_ids, match, row_match, n, m);
	} else if(cheap_id == do_sk_cheap) {
		sk_cheap(col_ptrs, col_ids, row_ptrs, row_ids, match, row_match, n, m);
	} else if(cheap_id == do_sk_cheap_rand) {
		sk_cheap_rand(col_ptrs, col_ids, row_ptrs, row_ids, match, row_match, n, m);
	} else if(cheap_id == do_mind_cheap) {
		mind_cheap(col_ptrs, col_ids, row_ptrs, row_ids, match, row_match, n, m);
	}
}
