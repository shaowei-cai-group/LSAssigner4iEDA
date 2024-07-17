#include "heuristic.h"
#include "basis_pms.h"
#include "my_class.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace lsa
{
void init_score_multi(externvariables &obj)
{
    int sense, weight, v, c;
    for (v = 1; v <= obj.num_vars; v++)
        {
            obj.hscore[v] = 0;
            obj.sscore[v] = 0;
            for (int i = 0; i < obj.var_lit_count[v]; ++i)
                {
                    c = obj.var_lit[v][i].clause_num;
                    sense = obj.var_lit[v][i].sense;
                    weight = obj.var_lit[v][i].weight;

                    if (obj.org_clause_weight[c] == obj.top_clause_weight) // hard
                        {
                            if (obj.sat_count[c] < obj.clause_true_lit_thres[c]) // falsified
                                {
                                    if (sense != obj.cur_soln[v]) // flip better
                                        {
                                            obj.hscore[v] += double(obj.tuned_degree_unit_weight[c] *
                                                                    std::min(obj.clause_true_lit_thres[c] - obj.sat_count[c], weight));
                                        }
                                    else // flip worse
                                        {
                                            obj.hscore[v] -= double(obj.tuned_degree_unit_weight[c] * weight);
                                        }
                                }
                            else if (obj.sat_count[c] >= obj.clause_true_lit_thres[c]) // satisfied
                                {
                                    if (sense == obj.cur_soln[v]) // flip worse
                                        {
                                            obj.hscore[v] -= double(obj.tuned_degree_unit_weight[c] *
                                                                    std::max(0, obj.clause_true_lit_thres[c] - obj.sat_count[c] + weight));
                                        }
                                }
                        }
                    else // soft
                        {
                            if (obj.sat_count[c] < obj.clause_true_lit_thres[c]) // falsified
                                {
                                    if (sense != obj.cur_soln[v]) // flip better
                                        {
                                            obj.sscore[v] += obj.unit_weight[c] * obj.tune_soft_clause_weight[c];
                                        }
                                    else // flip worse
                                        {
                                            obj.sscore[v] -= obj.unit_weight[c] * obj.tune_soft_clause_weight[c];
                                        }
                                }
                            else if (obj.sat_count[c] >= obj.clause_true_lit_thres[c]) // satisfied
                                {
                                    if (sense == obj.cur_soln[v]) // flip worse
                                        {
                                            obj.sscore[v] -= obj.unit_weight[c] * obj.tune_soft_clause_weight[c];
                                        }
                                }
                        }
                }
        }
    return;
}

// void (*flip_ptr)(int flipvar);
void flip_update_score_multi(int flipvar, externvariables &obj)
{
    int i, v, c, j;
    int index;
    lit *clause_c;
    int weight;
    int gap = 0;
    double change = 0;
    for (i = 0; i < obj.var_lit_count[flipvar]; ++i)
        {
            c = obj.var_lit[flipvar][i].clause_num;
            clause_c = obj.clause_lit[c];
            weight = obj.var_lit[flipvar][i].weight;
            if (obj.org_clause_weight[c] == obj.top_clause_weight) // hard
                {
                    if (obj.cur_soln[flipvar] == obj.var_lit[flipvar][i].sense) // flip better
                        {
                            if (obj.sat_count[c] + weight < obj.clause_true_lit_thres[c]) // 1. falsified to
                                                                                          // falsified
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense != obj.cur_soln[v])
                                                {
                                                    change = double((obj.tuned_degree_unit_weight[c] *
                                                                     (std::min(gap, obj.clause_lit[c][j].weight) -
                                                                      std::min(gap - weight, obj.clause_lit[c][j].weight))));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] -= change;
                                                        }
                                                }
                                        }
                                }
                            else if (obj.sat_count[c] < obj.clause_true_lit_thres[c]) // 2. falsified to satisfied;
                                                                                      // //obj.sat_count[c]+weight >
                                                                                      // obj.clause_true_lit_thres[c]
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense != obj.cur_soln[v])
                                                {
                                                    obj.hscore[v] -=
                                                        double((obj.tuned_degree_unit_weight[c] * std::min(gap, obj.clause_lit[c][j].weight)));
                                                }
                                            else
                                                {
                                                    obj.hscore[v] += double(
                                                        obj.tuned_degree_unit_weight[c] *
                                                        (obj.clause_lit[c][j].weight - std::max(0, gap - weight + obj.clause_lit[c][j].weight)));
                                                }
                                        }
                                    sat(c, obj);
                                }
                            else // 3. satisfied to satisfied;
                                 // //obj.sat_count[c]+weight >
                                 // obj.clause_true_lit_thres[c],
                                 // obj.sat_count[c] >
                                 // obj.clause_true_lit_thres[c]
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense == obj.cur_soln[v])
                                                {
                                                    change = double(obj.tuned_degree_unit_weight[c] *
                                                                    (std::max(0, gap + obj.clause_lit[c][j].weight) -
                                                                     std::max(0, gap - weight + obj.clause_lit[c][j].weight)));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] += change;
                                                        }
                                                }
                                        }
                                }

                            obj.sat_count[c] += weight;
                        }
                    else // flip worse;// obj.cur_soln[flipvar] != cur_lit.sense
                        {
                            if (obj.sat_count[c] - weight >= obj.clause_true_lit_thres[c]) // 4. satisfied to
                                                                                           // satisfied
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense == obj.cur_soln[v])
                                                {
                                                    change = double(obj.tuned_degree_unit_weight[c] *
                                                                    (std::max(0, gap + weight + obj.clause_lit[c][j].weight) -
                                                                     std::max(0, gap + obj.clause_lit[c][j].weight)));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] -= change;
                                                        }
                                                }
                                        }
                                }
                            else if (obj.sat_count[c] >= obj.clause_true_lit_thres[c]) // 5. satisfied to
                                                                                       // falsified
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense == obj.cur_soln[v])
                                                {
                                                    obj.hscore[v] -=
                                                        double(obj.tuned_degree_unit_weight[c] *
                                                               (obj.clause_lit[c][j].weight - std::max(0, gap + obj.clause_lit[c][j].weight)));
                                                }
                                            else
                                                {
                                                    obj.hscore[v] +=
                                                        double(obj.tuned_degree_unit_weight[c] * std::min(obj.clause_lit[c][j].weight, gap + weight));
                                                }
                                        }
                                    unsat(c, obj);
                                }
                            else // 6.  falsified to falsified
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense != obj.cur_soln[v])
                                                {
                                                    change = double(obj.tuned_degree_unit_weight[c] *
                                                                    (std::min(obj.clause_lit[c][j].weight, gap + weight) -
                                                                     std::min(obj.clause_lit[c][j].weight, gap)));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] += change;
                                                        }
                                                }
                                        }
                                }

                            obj.sat_count[c] -= weight;

                        } // end else
                }
            else // soft
                {
                    if (obj.cur_soln[flipvar] == obj.var_lit[flipvar][i].sense) // flip better
                        {
                            obj.soft_unsat_weight -= obj.org_clause_weight[c];
                            sat(c, obj);
                            obj.sat_count[c] += weight;
                        }
                    else // flip worse
                        {
                            obj.soft_unsat_weight += obj.org_clause_weight[c];
                            unsat(c, obj);
                            obj.sat_count[c] -= weight;
                        } // end else
                }
        }
    return;
}
void flip_update_score_no_neighbor_multi(int flipvar, externvariables &obj)
{
    int i, v, c, j;
    int index;
    lit *clause_c;
    int weight;
    int gap = 0;
    double change = 0;
    for (i = 0; i < obj.var_lit_count[flipvar]; ++i)
        {
            c = obj.var_lit[flipvar][i].clause_num;
            clause_c = obj.clause_lit[c];
            weight = obj.var_lit[flipvar][i].weight;
            if (obj.org_clause_weight[c] == obj.top_clause_weight) // hard
                {
                    if (obj.cur_soln[flipvar] == obj.var_lit[flipvar][i].sense) // flip better
                        {
                            if (obj.sat_count[c] + weight < obj.clause_true_lit_thres[c]) // 1. falsified to
                                                                                          // falsified
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense != obj.cur_soln[v])
                                                {
                                                    change = double((obj.tuned_degree_unit_weight[c] *
                                                                     (std::min(gap, obj.clause_lit[c][j].weight) -
                                                                      std::min(gap - weight, obj.clause_lit[c][j].weight))));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] -= change;
                                                            if (obj.already_in_goodvar_stack[v] != -1 && obj.hscore[v] + obj.sscore[v] <= 0)
                                                                {
                                                                    int top_v = mypop(obj.goodvar_stack);
                                                                    obj.goodvar_stack[obj.already_in_goodvar_stack[v]] = top_v;
                                                                    obj.already_in_goodvar_stack[top_v] = obj.already_in_goodvar_stack[v];
                                                                    obj.already_in_goodvar_stack[v] = -1;
                                                                }
                                                        }
                                                }
                                        }
                                }
                            else if (obj.sat_count[c] < obj.clause_true_lit_thres[c]) // 2. falsified to satisfied;
                                                                                      // //obj.sat_count[c]+weight >
                                                                                      // obj.clause_true_lit_thres[c]
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense != obj.cur_soln[v])
                                                {
                                                    obj.hscore[v] -=
                                                        double((obj.tuned_degree_unit_weight[c] * std::min(gap, obj.clause_lit[c][j].weight)));
                                                    if (obj.already_in_goodvar_stack[v] != -1 && obj.hscore[v] + obj.sscore[v] <= 0)
                                                        {
                                                            int top_v = mypop(obj.goodvar_stack);
                                                            obj.goodvar_stack[obj.already_in_goodvar_stack[v]] = top_v;
                                                            obj.already_in_goodvar_stack[top_v] = obj.already_in_goodvar_stack[v];
                                                            obj.already_in_goodvar_stack[v] = -1;
                                                        }
                                                }
                                            else
                                                {
                                                    obj.hscore[v] += double(
                                                        obj.tuned_degree_unit_weight[c] *
                                                        (obj.clause_lit[c][j].weight - std::max(0, gap - weight + obj.clause_lit[c][j].weight)));
                                                    if (obj.hscore[v] + obj.sscore[v] > 0 && obj.already_in_goodvar_stack[v] == -1)
                                                        {
                                                            obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                                                            mypush(v, obj.goodvar_stack);
                                                        }
                                                }
                                        }
                                    sat(c, obj);
                                }
                            else // 3. satisfied to satisfied;
                                 // //obj.sat_count[c]+weight >
                                 // obj.clause_true_lit_thres[c],
                                 // obj.sat_count[c] >
                                 // obj.clause_true_lit_thres[c]
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense == obj.cur_soln[v])
                                                {
                                                    change = double(obj.tuned_degree_unit_weight[c] *
                                                                    (std::max(0, gap + obj.clause_lit[c][j].weight) -
                                                                     std::max(0, gap - weight + obj.clause_lit[c][j].weight)));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] += change;
                                                            if (obj.hscore[v] + obj.sscore[v] > 0 && obj.already_in_goodvar_stack[v] == -1)
                                                                {
                                                                    obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                                                                    mypush(v, obj.goodvar_stack);
                                                                }
                                                        }
                                                }
                                        }
                                }
                            obj.sat_count[c] += weight;
                        }
                    else // flip worse;// obj.cur_soln[flipvar] != cur_lit.sense
                        {
                            if (obj.sat_count[c] - weight >= obj.clause_true_lit_thres[c]) // 4. satisfied to
                                                                                           // satisfied
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense == obj.cur_soln[v])
                                                {
                                                    change = double(obj.tuned_degree_unit_weight[c] *
                                                                    (std::max(0, gap + weight + obj.clause_lit[c][j].weight) -
                                                                     std::max(0, gap + obj.clause_lit[c][j].weight)));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] -= change;
                                                            if (obj.already_in_goodvar_stack[v] != -1 && obj.hscore[v] + obj.sscore[v] <= 0)
                                                                {
                                                                    int top_v = mypop(obj.goodvar_stack);
                                                                    obj.goodvar_stack[obj.already_in_goodvar_stack[v]] = top_v;
                                                                    obj.already_in_goodvar_stack[top_v] = obj.already_in_goodvar_stack[v];
                                                                    obj.already_in_goodvar_stack[v] = -1;
                                                                }
                                                        }
                                                }
                                        }
                                }
                            else if (obj.sat_count[c] >= obj.clause_true_lit_thres[c]) // 5. satisfied to
                                                                                       // falsified
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense == obj.cur_soln[v])
                                                {
                                                    obj.hscore[v] -=
                                                        double(obj.tuned_degree_unit_weight[c] *
                                                               (obj.clause_lit[c][j].weight - std::max(0, gap + obj.clause_lit[c][j].weight)));
                                                    if (obj.already_in_goodvar_stack[v] != -1 && obj.hscore[v] + obj.sscore[v] <= 0)
                                                        {
                                                            int top_v = mypop(obj.goodvar_stack);
                                                            obj.goodvar_stack[obj.already_in_goodvar_stack[v]] = top_v;
                                                            obj.already_in_goodvar_stack[top_v] = obj.already_in_goodvar_stack[v];
                                                            obj.already_in_goodvar_stack[v] = -1;
                                                        }
                                                }
                                            else
                                                {
                                                    obj.hscore[v] +=
                                                        double(obj.tuned_degree_unit_weight[c] * std::min(obj.clause_lit[c][j].weight, gap + weight));
                                                    if (obj.hscore[v] + obj.sscore[v] > 0 && obj.already_in_goodvar_stack[v] == -1)
                                                        {
                                                            obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                                                            mypush(v, obj.goodvar_stack);
                                                        }
                                                }
                                        }
                                    unsat(c, obj);
                                }
                            else // 6.  falsified to falsified
                                {
                                    gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
                                    for (j = 0; j < obj.clause_lit_count[c]; j++)
                                        {
                                            v = obj.clause_lit[c][j].var_num;
                                            if (v == flipvar)
                                                {
                                                    continue;
                                                }
                                            if (obj.clause_lit[c][j].sense != obj.cur_soln[v])
                                                {
                                                    change = double(obj.tuned_degree_unit_weight[c] *
                                                                    (std::min(obj.clause_lit[c][j].weight, gap + weight) -
                                                                     std::min(obj.clause_lit[c][j].weight, gap)));
                                                    if (0 == change)
                                                        {
                                                            break;
                                                        }
                                                    else
                                                        {
                                                            obj.hscore[v] += change;
                                                            if (obj.hscore[v] + obj.sscore[v] > 0 && obj.already_in_goodvar_stack[v] == -1)
                                                                {
                                                                    obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                                                                    mypush(v, obj.goodvar_stack);
                                                                }
                                                        }
                                                }
                                        }
                                }

                            obj.sat_count[c] -= weight;

                        } // end else
                }
            else // soft
                {
                    if (obj.cur_soln[flipvar] == obj.var_lit[flipvar][i].sense) // flip better
                        {
                            obj.soft_unsat_weight -= obj.org_clause_weight[c];
                            sat(c, obj);
                            obj.sat_count[c] += weight;
                        }
                    else // flip worse
                        {
                            obj.soft_unsat_weight += obj.org_clause_weight[c];
                            unsat(c, obj);
                            obj.sat_count[c] -= weight;
                        } // end else
                }
        }
    return;
}

void update_weight_score_multi(int c, externvariables &obj)
{
    int i = 0, v = 0, weight;
    for (i = 0; i < obj.clause_lit_count[c]; i++)
        {
            weight = obj.clause_lit[c][i].weight;
            v = obj.clause_lit[c][i].var_num;
            if (obj.clause_lit[c][i].sense != obj.cur_soln[v])
                {
                    obj.hscore[v] += double(obj.h_inc * std::min(obj.clause_true_lit_thres[c] - obj.sat_count[c], weight)) / obj.avg_clause_coe[c];
                    if (obj.hscore[v] + obj.sscore[v] > 0 && obj.already_in_goodvar_stack[v] == -1)
                        {
                            obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                            mypush(v, obj.goodvar_stack);
                        }
                }
            else
                {
                    obj.hscore[v] -= double(obj.h_inc * weight) / obj.avg_clause_coe[c];
                    if (obj.already_in_goodvar_stack[v] != -1 && obj.hscore[v] + obj.sscore[v] <= 0)
                        {
                            int top_v = mypop(obj.goodvar_stack);
                            obj.goodvar_stack[obj.already_in_goodvar_stack[v]] = top_v;
                            obj.already_in_goodvar_stack[top_v] = obj.already_in_goodvar_stack[v];
                            obj.already_in_goodvar_stack[v] = -1;
                        }
                }
        }
    return;
}

// int (*select_var_after_update_weight_ptr)(externvariables &obj);

// double (*soft_var_greedy_ptr)(int v,externvariables &obj);
// double (*hard_var_greedy_ptr)(int v,externvariables &obj);
double var_greedy_hscore(int v, externvariables &obj)
{
    return obj.hscore[v];
}
double var_greedy_sscore(int v, externvariables &obj)
{
    return obj.sscore[v];
}
double var_greedy_score(int v, externvariables &obj)
{
    return obj.hscore[v] + obj.sscore[v];
}
int select_var_after_update_weight_1(externvariables &obj)
{
    int r, c, i, l, best_var;
    if (obj.hardunsat_stack_fill_pointer > 0)
        {
            int gap;
            c = obj.hardunsat_stack[rand() % obj.hardunsat_stack_fill_pointer];

            gap = obj.clause_true_lit_thres[c] - obj.sat_count[c];
            l = 0;
            for (i = 0; i < obj.clause_lit_count[c]; i++)
                {
                    if (obj.clause_lit[c][i].sense != obj.cur_soln[obj.clause_lit[c][i].var_num])
                        {
                            obj.temp_unsat[l].var_num = obj.clause_lit[c][i].var_num;
                            obj.temp_unsat[l].weight = obj.clause_lit[c][i].weight;
                            l++;
                        }
                }

            if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < obj.rwprob)
                return obj.temp_unsat[rand() % l].var_num;

            int best_weight;
            best_var = obj.temp_unsat[0].var_num;
            best_weight = obj.temp_unsat[0].weight;
            for (i = 1; i < l; i++)
                {
                    if (obj.hscore[obj.temp_unsat[i].var_num] + obj.sscore[obj.temp_unsat[i].var_num] > obj.hscore[best_var] + obj.sscore[best_var])
                        {
                            best_var = obj.temp_unsat[i].var_num;
                            best_weight = obj.temp_unsat[i].weight;
                        }
                    else if (obj.hscore[obj.temp_unsat[i].var_num] + obj.sscore[obj.temp_unsat[i].var_num] ==
                             obj.hscore[best_var] + obj.sscore[best_var])
                        {
                            if (obj.time_stamp[obj.temp_unsat[i].var_num] < obj.time_stamp[best_var])
                                {
                                    best_var = obj.temp_unsat[i].var_num;
                                    best_weight = obj.temp_unsat[i].weight;
                                }
                        }
                }
            return best_var;
        }
    else
        {
            if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < obj.rwprob)
                c = obj.softunsat_stack[rand() % obj.softunsat_stack_fill_pointer];
            else
                {
                    r = rand() % obj.softunsat_stack_fill_pointer;
                    c = obj.softunsat_stack[r];
                    for (i = 1; i < obj.hd_count_threshold; i++)
                        {
                            r = rand() % obj.softunsat_stack_fill_pointer;
                            if (obj.hscore[obj.clause_lit[obj.softunsat_stack[r]][0].var_num] +
                                    obj.sscore[obj.clause_lit[obj.softunsat_stack[r]][0].var_num] >
                                obj.hscore[obj.clause_lit[c][0].var_num] + obj.sscore[obj.clause_lit[c][0].var_num])
                                {
                                    c = obj.softunsat_stack[r];
                                }
                        }
                }
            // c = obj.softunsat_stack[rand() %
            // obj.softunsat_stack_fill_pointer];
            return obj.clause_lit[c][0].var_num;
        }
}
int select_var_after_update_weight_2(externvariables &obj)
{
    int r, c, i, l, best_var, best_w, temp_l, v, w;

    if (obj.hardunsat_stack_fill_pointer > 0)
        {
            c = obj.hardunsat_stack[rand() % obj.hardunsat_stack_fill_pointer];
            l = 0;
            for (i = 0; i < obj.clause_lit_count[c]; i++)
                {
                    if (obj.clause_lit[c][i].sense != obj.cur_soln[obj.clause_lit[c][i].var_num])
                        {
                            obj.temp_unsat[l].var_num = obj.clause_lit[c][i].var_num;
                            obj.temp_unsat[l].weight = obj.clause_lit[c][i].weight;
                            l++;
                        }
                }
            if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < obj.rwprob)
                {
                    return obj.temp_unsat[rand() % l].var_num;
                }
            else
                {
                    v = obj.temp_unsat[0].var_num;
                    best_w = obj.hscore[v] + obj.sscore[v];
                    temp_l = 0;
                    obj.temp_array[temp_l++] = v;
                    for (i = 1; i < l; i++)
                        {
                            v = obj.temp_unsat[i].var_num;
                            if (best_w < obj.hscore[v] + obj.sscore[v])
                                {
                                    temp_l = 0;
                                    obj.temp_array[temp_l++] = v;
                                    best_w = obj.hscore[v] + obj.sscore[v];
                                }
                            else if (best_w == obj.hscore[v] + obj.sscore[v] && obj.time_stamp[v] < obj.time_stamp[obj.temp_array[0]])
                                {
                                    temp_l = 0;
                                    obj.temp_array[temp_l++] = v;
                                }
                            else if (best_w == obj.hscore[v] + obj.sscore[v] && obj.time_stamp[v] == obj.time_stamp[obj.temp_array[0]])
                                {
                                    obj.temp_array[temp_l++] = v;
                                }
                        }
                    return obj.temp_array[rand() % temp_l];
                }
        }
    else
        {
            c = obj.softunsat_stack[rand() % obj.softunsat_stack_fill_pointer];
            return obj.clause_lit[c][0].var_num;
        }

    return rand() % obj.num_vars + 1;
}
} // namespace lsa