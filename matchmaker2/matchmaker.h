/*
 * File: matchmaker.h
 * Content: Method descriptions and constants for matchmaker
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

 
#ifndef MATCHMAKER_H_
#define MATCHMAKER_H_

#define do_old_cheap 1
#define do_sk_cheap 2
#define do_sk_cheap_rand 3
#define do_mind_cheap 4

#define do_dfs 1
#define do_bfs 2
#define do_mc21 3
#define do_pf 4
#define do_pf_fair 5
#define do_hk 6
#define do_hk_dw 7
#define do_abmp 8
#define do_abmp_bfs 9
#define do_pr_fifo_fair 10

void old_cheap(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
void sk_cheap(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);
void sk_cheap_rand(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);
void mind_cheap(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);

void match_dfs(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
void match_bfs(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
void match_mc21(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
void match_pf(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
void match_pf_fair(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
void match_hk(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);
void match_hk_dw(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);
void match_abmp(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);
void match_abmp_bfs(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);
void match_pr_fifo_fair(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m, double relabel_period);

void pr_global_relabel(int* l_label, int* r_label, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m);

void cheap_matching(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, int* match, int* row_match, int n, int m, int cheap_id);
void matching(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m, int match_id, int cheap_id, double relabel_period);
#endif /* MATCHMAKER_H_ */
