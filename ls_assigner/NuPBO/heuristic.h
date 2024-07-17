#ifndef _HEURISTIC_H
#define _HEURISTIC_H
#include "basis_pms.h"
#include "my_class.h"
namespace lsa
{
void init_score_multi(externvariables &obj);

//extern void (*flip_ptr)(int flipvar);
void flip_with_neighbor(int flipvar,externvariables &obj);
void flip_no_neighbor(int flipvar,externvariables &obj);

void flip_update_score_multi(int flipvar,externvariables &obj);
void flip_update_score_no_neighbor_multi(int flipvar,externvariables &obj);

void update_weight_score_multi(int c,externvariables &obj);

//extern int (*select_var_after_update_weight_ptr)(externvariables &obj);
int select_var_after_update_weight_1(externvariables &obj);
int select_var_after_update_weight_2(externvariables &obj);

//extern double (*soft_var_greedy_ptr)(int v,externvariables &obj);
//extern double (*hard_var_greedy_ptr)(int v,externvariables &obj);
double var_greedy_hscore(int v,externvariables &obj);
double var_greedy_sscore(int v,externvariables &obj);
double var_greedy_score(int v,externvariables &obj);
} // namespace lsa
#endif