#include "basis_pms.h"
#include "MaxSATFormula.h"
#include "heuristic.h"
#include "init.h"
#include "my_class.h"
#include "struct_from_basis.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace lsa
{
long long total_step;

//控制local_search的执行时间
std::chrono::high_resolution_clock::time_point program_start_time;

bool compare(temp_var a, temp_var b)
{
    return a.weight > b.weight;
}

void settings(externvariables &obj)
{
    // steps
    total_step = 0;
    obj.max_tries = 100000000;
    obj.max_flips = 10000000;
    obj.max_non_improve_flip = 10000000;

    obj.rdprob = 0.01;
    obj.rwprob = 0.1;
    obj.hd_count_threshold = 15;
    obj.h_inc = 1;
    obj.s_inc = 1;

    obj.best_known = -2;
    if (obj.num_vars > 2000)
        {
            obj.rdprob = 0.01;
            obj.rwprob = 0.1;
            obj.hd_count_threshold = 50;
        }
    obj.hard_var_greedy_ptr = var_greedy_score;
    obj.soft_var_greedy_ptr = var_greedy_score;
}

void allocate_memory(externvariables &obj)
{
    int malloc_var_length = obj.num_vars + 10;
    int malloc_clause_length = obj.num_clauses + 10;

    obj.temp_unsat = new temp_var[malloc_var_length];
    obj.temp_array = new int[malloc_clause_length];

    obj.var_lit = new lit *[malloc_var_length];
    obj.var_lit_count = new int[malloc_var_length]();
    obj.clause_lit = new lit *[malloc_clause_length];
    obj.clause_lit_count = new int[malloc_clause_length]();
    obj.clause_true_lit_thres = new int[malloc_clause_length];
    obj.avg_clause_coe = new double[malloc_clause_length]();

    obj.hscore = new double[malloc_var_length];
    obj.sscore = new double[malloc_var_length];
    obj.var_neighbor = new int *[malloc_var_length];
    obj.var_neighbor_count = new int[malloc_var_length];
    obj.time_stamp = new long long[malloc_var_length];
    obj.neighbor_flag = new int[malloc_var_length];
    obj.temp_neighbor = new int[malloc_var_length];

    obj.soft_clause_num_index = new int[malloc_clause_length];
    obj.hard_clause_num_index = new int[malloc_clause_length];
    obj.org_clause_weight = new long long[malloc_clause_length];
    obj.tune_soft_clause_weight = new double[malloc_clause_length];
    obj.unit_weight = new double[malloc_clause_length];
    obj.tuned_degree_unit_weight = new double[malloc_clause_length];
    obj.sat_count = new int[malloc_clause_length];
    obj.sat_var = new int[malloc_clause_length];

    obj.hardunsat_stack = new int[malloc_clause_length];
    obj.index_in_hardunsat_stack = new int[malloc_clause_length];
    obj.softunsat_stack = new int[malloc_clause_length];
    obj.index_in_softunsat_stack = new int[malloc_clause_length];

    obj.goodvar_stack = new int[malloc_var_length];
    obj.already_in_goodvar_stack = new int[malloc_var_length];

    obj.cur_soln = new int[malloc_var_length];
    obj.best_soln = new int[malloc_var_length];
}

void free_memory(externvariables &obj)
{
    int i;
    for (i = 0; i < obj.num_clauses; i++)
        delete[] obj.clause_lit[i];

    for (i = 1; i <= obj.num_vars; ++i)
        {
            delete[] obj.var_lit[i];
            delete[] obj.var_neighbor[i];
        }
    delete[] obj.temp_array;
    delete[] obj.temp_unsat;
    delete[] obj.var_lit;
    delete[] obj.var_lit_count;
    delete[] obj.clause_lit;
    delete[] obj.clause_lit_count;
    delete[] obj.clause_true_lit_thres;
    delete[] obj.avg_clause_coe;

    delete[] obj.hscore;
    delete[] obj.sscore;
    delete[] obj.var_neighbor;
    delete[] obj.var_neighbor_count;
    delete[] obj.time_stamp;
    delete[] obj.neighbor_flag;
    delete[] obj.temp_neighbor;

    delete[] obj.soft_clause_num_index;
    delete[] obj.hard_clause_num_index;
    delete[] obj.org_clause_weight;
    delete[] obj.tune_soft_clause_weight;
    delete[] obj.unit_weight;
    delete[] obj.tuned_degree_unit_weight;
    delete[] obj.sat_count;
    delete[] obj.sat_var;

    delete[] obj.hardunsat_stack;
    delete[] obj.index_in_hardunsat_stack;
    delete[] obj.softunsat_stack;
    delete[] obj.index_in_softunsat_stack;

    delete[] obj.goodvar_stack;
    delete[] obj.already_in_goodvar_stack;

    delete[] obj.cur_soln;
    delete[] obj.best_soln;
}

void build_neighbor_relation(externvariables &obj)
{
    // obj.outputFile << "c start build neighbor" << std::endl;
    int i, j, count;
    int v, c, n;
    int temp_neighbor_count;

    for (v = 1; v <= obj.num_vars; ++v)
        {
            obj.neighbor_flag[v] = 1;
            temp_neighbor_count = 0;

            for (i = 0; i < obj.var_lit_count[v]; ++i)
                {
                    c = obj.var_lit[v][i].clause_num;
                    for (j = 0; j < obj.clause_lit_count[c]; ++j)
                        {
                            n = obj.clause_lit[c][j].var_num;
                            if (obj.neighbor_flag[n] != 1)
                                {
                                    obj.neighbor_flag[n] = 1;
                                    obj.temp_neighbor[temp_neighbor_count++] = n;
                                }
                        }
                }

            obj.neighbor_flag[v] = 0;
            obj.var_neighbor[v] = new int[temp_neighbor_count];
            obj.var_neighbor_count[v] = temp_neighbor_count;

            count = 0;
            for (i = 0; i < temp_neighbor_count; i++)
                {
                    obj.var_neighbor[v][count++] = obj.temp_neighbor[i];
                    obj.neighbor_flag[obj.temp_neighbor[i]] = 0;
                }
        }
    // obj.outputFile << "c end build neighbor" << std::endl;
}

void build_instance(const char *filename, externvariables &obj)
{
    int i = 0, v = 0, c = 0;
    long long j = 0, k = 0;
    MaxSATFormula maxsatformula;
    maxsatformula.build_ins(filename, maxsatformula);

    obj.num_vars = maxsatformula.num_vars;
    obj.num_clauses = maxsatformula.hclause.size() + maxsatformula.sclause.size();
    obj.top_clause_weight = maxsatformula.top_clause_weight + 1;

    allocate_memory(obj);

    obj.num_hclauses = obj.num_sclauses = 0;
    // Now, read clauses, one at a time.
    temp_var *temp_weight = new temp_var[obj.num_vars];
    int cur_weight;
    int cur_lit;

    c = 0;
    for (k = 0; k < maxsatformula.hclause.size(); ++k)
        {
            HardC &temc = maxsatformula.hclause[k];
            obj.clause_lit_count[c] = 0;
            obj.org_clause_weight[c] = obj.top_clause_weight;
            obj.clause_true_lit_thres[c] = temc.degree;
            if (obj.clause_true_lit_thres[c] <= 0)
                {
                    obj.num_clauses--;
                    continue;
                }
            obj.hard_clause_num_index[obj.num_hclauses++] = c;
            // obj.outputFile << c << " " << obj.top_clause_weight << " " <<
            // temc.degree << " ";
            for (j = 0; j < temc.clause.size(); ++j)
                {
                    cur_weight = temc.clause[j].weight;
                    cur_lit = temc.clause[j].var;
                    temp_weight[obj.clause_lit_count[c]].weight = cur_weight;
                    temp_weight[obj.clause_lit_count[c]].var_num = cur_lit;
                    obj.clause_lit_count[c]++;
                    // obj.outputFile  <<cur_weight << " " << cur_lit << " ";
                }
            // obj.outputFile << obj.clause_lit_count[c] <<  std::endl;
            std::sort(temp_weight, temp_weight + obj.clause_lit_count[c], compare);
            obj.clause_lit[c] = new lit[obj.clause_lit_count[c] + 1];

            for (i = 0; i < obj.clause_lit_count[c]; ++i)
                {
                    obj.clause_lit[c][i].clause_num = c;
                    obj.clause_lit[c][i].var_num = abs(temp_weight[i].var_num);
                    obj.clause_lit[c][i].weight = temp_weight[i].weight;
                    obj.avg_clause_coe[c] += double(obj.clause_lit[c][i].weight);

                    if (temp_weight[i].var_num > 0)
                        obj.clause_lit[c][i].sense = 1;
                    else
                        obj.clause_lit[c][i].sense = 0;

                    obj.var_lit_count[obj.clause_lit[c][i].var_num]++;
                }

            // round
            obj.avg_clause_coe[c] = round(double(obj.avg_clause_coe[c] / (double)obj.clause_lit_count[c]));
            if (obj.avg_clause_coe[c] < 1)
                obj.avg_clause_coe[c] = 1;

            obj.clause_lit[c][i].var_num = 0;
            obj.clause_lit[c][i].clause_num = -1;
            obj.clause_lit[c][i].weight = 0;

            c++;
        }

    for (k = 0; k < maxsatformula.sclause.size(); ++k)
        {
            SoftC &temc = maxsatformula.sclause[k];
            obj.clause_lit_count[c] = 0;
            obj.org_clause_weight[c] = temc.weight;
            obj.clause_true_lit_thres[c] = temc.degree;
            if (obj.clause_true_lit_thres[c] <= 0)
                {
                    obj.num_clauses--;
                    continue;
                }
            // obj.outputFile << c << " " << temc.weight << " " << temc.degree
            // << "
            // ";
            obj.soft_clause_num_index[obj.num_sclauses++] = c;

            cur_weight = temc.clause[0].weight;
            cur_lit = temc.clause[0].var;
            temp_weight[obj.clause_lit_count[c]].weight = cur_weight;
            temp_weight[obj.clause_lit_count[c]].var_num = cur_lit;
            obj.clause_lit_count[c] = 1;

            obj.clause_lit[c] = new lit[obj.clause_lit_count[c] + 1];

            obj.clause_lit[c][0].clause_num = c;
            obj.clause_lit[c][0].var_num = abs(temp_weight[0].var_num);
            obj.clause_lit[c][0].weight = temp_weight[0].weight;
            obj.avg_clause_coe[c] += double(obj.clause_lit[c][0].weight);

            // obj.outputFile  <<cur_weight << " " << cur_lit << " ";
            // obj.outputFile << std::endl;

            if (temp_weight[0].var_num > 0)
                obj.clause_lit[c][0].sense = 1;
            else
                obj.clause_lit[c][0].sense = 0;

            obj.var_lit_count[obj.clause_lit[c][0].var_num]++;

            // round
            obj.avg_clause_coe[c] = round(double(obj.avg_clause_coe[c] / (double)obj.clause_lit_count[c]));
            if (obj.avg_clause_coe[c] < 1)
                obj.avg_clause_coe[c] = 1;

            obj.clause_lit[c][1].var_num = 0;
            obj.clause_lit[c][1].clause_num = -1;
            obj.clause_lit[c][1].weight = 0;

            c++;
        }
    obj.outputFile << obj.num_clauses << " " << obj.num_hclauses << " " << obj.num_sclauses << std::endl;
    delete[] temp_weight;

    // creat var literal arrays
    long long tmp_lit_num = 0;
    double avg_neighbor_lit = 0;
    for (v = 1; v <= obj.num_vars; ++v)
        {
            obj.var_lit[v] = new lit[obj.var_lit_count[v] + 1];
            tmp_lit_num += obj.var_lit_count[v];
            obj.var_lit_count[v] = 0;
        }
    avg_neighbor_lit = double(tmp_lit_num - obj.num_sclauses) / (obj.num_vars - obj.num_sclauses + 1);
    // obj.outputFile << "c avg_neighbor_lit: " << avg_neighbor_lit<< std::endl;

    // scan all clauses to build up var literal arrays
    for (c = 0; c < obj.num_clauses; ++c)
        {
            for (i = 0; i < obj.clause_lit_count[c]; ++i)
                {
                    v = obj.clause_lit[c][i].var_num;
                    obj.var_lit[v][obj.var_lit_count[v]++] = obj.clause_lit[c][i];
                }
        }

    if (avg_neighbor_lit < 1e+7)
        {
            build_neighbor_relation(obj);
            obj.flip_ptr = flip_with_neighbor;
        }
    else
        {
            obj.flip_ptr = flip_no_neighbor;
        }
    obj.best_soln_feasible = 0;
    obj.opt_unsat_weight = obj.top_clause_weight;
}

void init_local_search(std::vector<int> &init_solution, externvariables &obj)
{
    int v, c;
    int i, j;

    obj.local_soln_feasible = 0;

    // Initialize clause information
    for (i = 0; i < obj.num_hclauses; i++)
        {
            c = obj.hard_clause_num_index[i];
            obj.unit_weight[c] = 1;
            obj.tuned_degree_unit_weight[c] = double(obj.unit_weight[c]) / obj.avg_clause_coe[c];
        }
    // round
    double tmp_avg_soft_clause_weight = 0.0;
    tmp_avg_soft_clause_weight = round(double(obj.top_clause_weight - 1) / obj.num_sclauses);
    if (tmp_avg_soft_clause_weight < 1)
        tmp_avg_soft_clause_weight = 1;
    for (i = 0; i < obj.num_sclauses; i++)
        {
            c = obj.soft_clause_num_index[i];
            obj.tune_soft_clause_weight[c] = double(obj.org_clause_weight[c] / tmp_avg_soft_clause_weight);
            obj.unit_weight[c] = obj.initsoftw;
        }

    // init solution
    obj.init_assignment_ptr(init_solution, obj);

    // init stacks
    obj.hard_unsat_nb = 0;
    obj.hardunsat_stack_fill_pointer = 0;
    obj.softunsat_stack_fill_pointer = 0;

    /* figure out sat_count, sat_var, soft_unsat_weight and init unsat_stack */
    obj.soft_unsat_weight = 0;

    for (i = 0; i < obj.num_hclauses; i++)
        {
            c = obj.hard_clause_num_index[i];
            obj.sat_count[c] = 0;
            for (j = 0; j < obj.clause_lit_count[c]; ++j)
                {
                    if (obj.cur_soln[obj.clause_lit[c][j].var_num] == obj.clause_lit[c][j].sense)
                        {
                            obj.sat_count[c] += obj.clause_lit[c][j].weight;
                            obj.sat_var[c] = obj.clause_lit[c][j].var_num;
                        }
                }
            if (obj.sat_count[c] < obj.clause_true_lit_thres[c])
                {
                    unsat(c, obj);
                }
        }
    for (i = 0; i < obj.num_sclauses; i++)
        {
            c = obj.soft_clause_num_index[i];
            obj.sat_count[c] = 0;

            if (obj.cur_soln[obj.clause_lit[c][0].var_num] == obj.clause_lit[c][0].sense)
                {
                    obj.sat_count[c] += obj.clause_lit[c][0].weight;
                    obj.sat_var[c] = obj.clause_lit[c][0].var_num;
                }
            else
                {
                    obj.soft_unsat_weight += (obj.clause_true_lit_thres[c] - obj.sat_count[c]) * obj.org_clause_weight[c];
                    unsat(c, obj);
                }
        }

    /*figure out score*/
    init_score_multi(obj);

    // init goodvars stack
    obj.goodvar_stack_fill_pointer = 0;
    for (v = 1; v <= obj.num_vars; v++)
        {
            if (obj.hscore[v] + obj.sscore[v] > 0)
                {
                    obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                    mypush(v, obj.goodvar_stack);
                }
            else
                obj.already_in_goodvar_stack[v] = -1;
        }
}

int pick_var(externvariables &obj)
{
    int i, v, r;
    int best_var;

    if (obj.goodvar_stack_fill_pointer > 0)
        {
            if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < obj.rdprob)
                return obj.goodvar_stack[rand() % obj.goodvar_stack_fill_pointer];

            if (obj.goodvar_stack_fill_pointer < obj.hd_count_threshold)
                {
                    best_var = obj.goodvar_stack[0];
                    for (i = 1; i < obj.goodvar_stack_fill_pointer; ++i)
                        {
                            v = obj.goodvar_stack[i];
                            if (obj.hscore[v] + obj.sscore[v] > obj.hscore[best_var] + obj.sscore[best_var])
                                best_var = v;
                            else if (obj.hscore[v] + obj.sscore[v] == obj.hscore[best_var] + obj.sscore[best_var])
                                {
                                    if (obj.time_stamp[v] < obj.time_stamp[best_var])
                                        best_var = v;
                                }
                        }
                    return best_var;
                }
            else
                {
                    r = rand() % obj.goodvar_stack_fill_pointer;
                    best_var = obj.goodvar_stack[r];

                    for (i = 1; i < obj.hd_count_threshold; ++i)
                        {
                            r = rand() % obj.goodvar_stack_fill_pointer;
                            v = obj.goodvar_stack[r];
                            if (obj.hscore[v] + obj.sscore[v] > obj.hscore[best_var] + obj.sscore[best_var])
                                best_var = v;
                            else if (obj.hscore[v] + obj.sscore[v] == obj.hscore[best_var] + obj.sscore[best_var])
                                {
                                    if (obj.time_stamp[v] < obj.time_stamp[best_var])
                                        best_var = v;
                                }
                        }
                    return best_var;
                }
        }
    update_clause_weights(obj);

    return obj.select_var_after_update_weight_ptr(obj);
}

void update_goodvarstack(int flipvar, externvariables &obj)
{
    int v;

    // remove the variables no longer goodvar in goodvar_stack
    for (int index = obj.goodvar_stack_fill_pointer - 1; index >= 0; index--)
        {
            v = obj.goodvar_stack[index];
            if (obj.hscore[v] + obj.sscore[v] <= 0)
                {
                    int top_v = mypop(obj.goodvar_stack);
                    obj.goodvar_stack[index] = top_v;
                    obj.already_in_goodvar_stack[top_v] = index;
                    obj.already_in_goodvar_stack[v] = -1;
                }
        }

    // add goodvar
    for (int i = 0; i < obj.var_neighbor_count[flipvar]; ++i)
        {
            v = obj.var_neighbor[flipvar][i];
            if (obj.hscore[v] + obj.sscore[v] > 0)
                {
                    if (obj.already_in_goodvar_stack[v] == -1)
                        {
                            obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                            mypush(v, obj.goodvar_stack);
                        }
                }
        }
}

void flip_with_neighbor(int flipvar, externvariables &obj)
{
    double org_flipvar_score = obj.hscore[flipvar];
    double org_sscore = obj.sscore[flipvar];
    obj.cur_soln[flipvar] = 1 - obj.cur_soln[flipvar];

    flip_update_score_multi(flipvar, obj);

    // update information of flipvar
    obj.hscore[flipvar] = -org_flipvar_score;
    obj.sscore[flipvar] = -org_sscore;
    update_goodvarstack(flipvar, obj);
}

void flip_no_neighbor(int flipvar, externvariables &obj)
{
    double org_flipvar_score = obj.hscore[flipvar];
    double org_sscore = obj.sscore[flipvar];
    obj.cur_soln[flipvar] = 1 - obj.cur_soln[flipvar];

    flip_update_score_no_neighbor_multi(flipvar, obj);

    // update information of flipvar
    obj.hscore[flipvar] = -org_flipvar_score;
    obj.sscore[flipvar] = -org_sscore;
    if (obj.already_in_goodvar_stack[flipvar] != -1)
        {
            int top_v = mypop(obj.goodvar_stack);
            obj.goodvar_stack[obj.already_in_goodvar_stack[flipvar]] = top_v;
            obj.already_in_goodvar_stack[top_v] = obj.already_in_goodvar_stack[flipvar];
            obj.already_in_goodvar_stack[flipvar] = -1;
        }
    return;
}

void print_best_solution(externvariables &obj)
{
    // obj.outputFile << "c total_step: " << total_step << " " << obj.tries <<
    // std::endl;
    if (0 == obj.is_print)
        {
            obj.is_print = 1;
            if (obj.best_soln_feasible == 1)
                {
                    if (verify_sol(obj))
                        {
                            obj.outputFile << "s ";
                            for (int i = 1; i <= obj.num_vars; i++)
                                {
                                    obj.outputFile << obj.best_soln[i] << " ";
                                }
                            obj.outputFile << std::endl;
                            // obj.outputFile << obj.opt_unsat_weight << '\t' <<
                            // obj.opt_time << '\t' << obj.tries << '\t' <<
                            // obj.step << std::endl;
                        }
                    else
                        obj.outputFile << "got an incorrect solution" << std::endl;
                }
            else
                obj.outputFile << "no feasible solution found" << std::endl;
        }
}

void local_search(std::vector<int> &init_solution, const char *inputfile, externvariables &obj)
{
    program_start_time = std::chrono::high_resolution_clock::now();
    for (obj.tries = 1; obj.tries < obj.max_tries; ++obj.tries)
        {
            init_local_search(init_solution, obj);
            obj.max_flips = 10000000;
            for (obj.step = 1; obj.step < obj.max_flips; ++obj.step)
                {
                    total_step++;
                    if (obj.hard_unsat_nb == 0)
                        {
                            obj.local_soln_feasible = 1;
                            obj.best_soln_feasible = 1;

                            if (obj.soft_unsat_weight < obj.opt_unsat_weight)
                                {
                                    obj.opt_unsat_weight = obj.soft_unsat_weight;
                                    obj.outputFile << "o " << obj.soft_unsat_weight << std::endl; 
                                    obj.opt_time = get_runtime(obj);
                                    for (int v = 1; v <= obj.num_vars; ++v)
                                        obj.best_soln[v] = obj.cur_soln[v];
                                    if (obj.opt_unsat_weight == 0 || obj.opt_unsat_weight <= obj.best_known)
                                        {
                                            return;
                                        }
                                    obj.max_flips = obj.step + 10000000;
                                }
                        }

                    int flipvar = pick_var(obj);

                    obj.flip_ptr(flipvar, obj);
                    // check_new_score();
                    obj.time_stamp[flipvar] = obj.step;
                    if (local_search_runtime() >= obj.time_limit)
                        {
                            // obj.outputFile << "Time limit reached. Exiting."
                            // << std::endl;
                            break;
                        }
                }
            if (local_search_runtime() >= obj.time_limit)
                {
                    // obj.outputFile << "Time limit reached. Exiting." <<
                    // std::endl;
                    return;
                }
            // obj.outputFile<<"外层for循环"<< fixed <<
            // setprecision(2)<<local_search_runtime()<<std::endl;
        }
}

void check_softunsat_weight(externvariables &obj)
{
    int c, j, flag;
    long long verify_unsat_weight = 0;

    for (c = 0; c < obj.num_clauses; ++c)
        {
            flag = 0;
            int tem_clause_true_lit_count = 0;
            for (j = 0; j < obj.clause_lit_count[c]; ++j)
                {
                    if (obj.cur_soln[obj.clause_lit[c][j].var_num] == obj.clause_lit[c][j].sense)
                        {
                            tem_clause_true_lit_count++;
                        }
                }
            if (tem_clause_true_lit_count >= obj.clause_true_lit_thres[c])
                flag = 1;

            if (flag == 0)
                {
                    if (obj.org_clause_weight[c] == obj.top_clause_weight) // verify hard clauses
                        {
                            continue;
                        }
                    else
                        {
                            verify_unsat_weight += obj.org_clause_weight[c] * (obj.clause_true_lit_thres[c] - tem_clause_true_lit_count);
                        }
                }
        }

    if (verify_unsat_weight != obj.soft_unsat_weight)
        {
            obj.outputFile << obj.step << std::endl;
            obj.outputFile << "verify unsat weight is" << verify_unsat_weight << " and soft unsat weight is " << obj.soft_unsat_weight << std::endl;
        }
    // return 0;
}

bool verify_sol(externvariables &obj)
{
    obj.outputFile << "c start verification" << std::endl;
    int c, j, flag;
    long long verify_unsat_weight = 0;
    long long real_min_weight = 0;

    for (c = 0; c < obj.num_clauses; ++c)
        {
            if (obj.org_clause_weight[c] != obj.top_clause_weight)
                {
                    if (obj.clause_lit[c][0].sense == 0)
                        {
                            real_min_weight += obj.org_clause_weight[c] * obj.best_soln[obj.clause_lit[c][0].var_num];
                        }
                    else
                        {
                            real_min_weight -= obj.org_clause_weight[c] * obj.best_soln[obj.clause_lit[c][0].var_num];
                        }
                }
            flag = 0;
            int tem_clause_true_lit_count = 0;
            for (j = 0; j < obj.clause_lit_count[c]; ++j)
                {
                    if (obj.best_soln[obj.clause_lit[c][j].var_num] == obj.clause_lit[c][j].sense)
                        {
                            tem_clause_true_lit_count += obj.clause_lit[c][j].weight;
                        }
                }
            if (tem_clause_true_lit_count >= obj.clause_true_lit_thres[c])
                flag = 1;

            if (flag == 0)
                {
                    if (obj.org_clause_weight[c] == obj.top_clause_weight) // verify hard clauses
                        {
                            // output the falsified clause under the assignment
                            obj.outputFile << "c Error: hard clause " << c << " is falsified" << std::endl;

                            obj.outputFile << "c ";
                            for (j = 0; j < obj.clause_lit_count[c]; ++j)
                                {
                                    if (obj.clause_lit[c][j].sense == 0)
                                        obj.outputFile << "-";
                                    obj.outputFile << obj.clause_lit[c][j].var_num << " ";
                                }
                            obj.outputFile << std::endl;
                            obj.outputFile << "c ";
                            for (j = 0; j < obj.clause_lit_count[c]; ++j)
                                obj.outputFile << obj.best_soln[obj.clause_lit[c][j].var_num] << " ";
                            obj.outputFile << std::endl;
                            return 0;
                        }
                    else
                        {
                            verify_unsat_weight += obj.org_clause_weight[c] * (obj.clause_true_lit_thres[c] - tem_clause_true_lit_count);
                        }
                }
        }

    if (verify_unsat_weight == obj.opt_unsat_weight)
        {
            obj.outputFile << "c realmin: " << real_min_weight << std::endl;
            return 1;
        }
    else
        {
            obj.outputFile << "c Error: find opt=" << obj.opt_unsat_weight << ", but verified opt=" << verify_unsat_weight << std::endl;
        }
    return 0;
}

void simple_print(externvariables &obj)
{
    if (obj.best_soln_feasible == 1)
        {
            if (verify_sol(obj) == 1)
                obj.outputFile << obj.opt_unsat_weight << '\t' << obj.opt_time << std::endl;
            else
                obj.outputFile << "solution is wrong " << std::endl;
        }
    else
        obj.outputFile << -1 << '\t' << -1 << std::endl;
}

void increase_weights(externvariables &obj)
{
    int i, c, v;
    int flag;

    for (i = 0; i < obj.hardunsat_stack_fill_pointer; ++i)
        {
            flag = 1;
            c = obj.hardunsat_stack[i];
            obj.unit_weight[c] += obj.h_inc;
            obj.tuned_degree_unit_weight[c] = double(obj.unit_weight[c]) / obj.avg_clause_coe[c];
            update_weight_score_multi(c, obj);
        }

    // increase all soft clause weights
    if (0 == obj.hard_unsat_nb)
        {
            for (i = 0; i < obj.num_sclauses; i++)
                {
                    c = obj.soft_clause_num_index[i];
                    obj.unit_weight[c] += obj.s_inc;
                    v = obj.clause_lit[c][0].var_num;

                    if (obj.clause_lit[c][0].sense != obj.cur_soln[v])
                        {
                            obj.sscore[v] += obj.s_inc * obj.tune_soft_clause_weight[c];
                            if (obj.hscore[v] + obj.sscore[v] > 0 && obj.already_in_goodvar_stack[v] == -1)
                                {
                                    obj.already_in_goodvar_stack[v] = obj.goodvar_stack_fill_pointer;
                                    mypush(v, obj.goodvar_stack);
                                }
                        }
                    else
                        {
                            obj.sscore[v] -= obj.s_inc * obj.tune_soft_clause_weight[c];
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

void update_clause_weights(externvariables &obj)
{
    increase_weights(obj);
}

void unsat(int clause, externvariables &obj)
{
    if (obj.org_clause_weight[clause] == obj.top_clause_weight) // hard
        {
            obj.index_in_hardunsat_stack[clause] = obj.hardunsat_stack_fill_pointer;
            mypush(clause, obj.hardunsat_stack);
            obj.hard_unsat_nb++;
        }
    else // soft
        {
            obj.index_in_softunsat_stack[clause] = obj.softunsat_stack_fill_pointer;
            mypush(clause, obj.softunsat_stack);
            // obj.soft_unsat_weight += obj.org_clause_weight[clause];
        }
}

void sat(int clause, externvariables &obj)
{
    int index, last_unsat_clause;

    if (obj.org_clause_weight[clause] == obj.top_clause_weight)
        {

            last_unsat_clause = mypop(obj.hardunsat_stack);
            index = obj.index_in_hardunsat_stack[clause];
            obj.hardunsat_stack[index] = last_unsat_clause;
            obj.index_in_hardunsat_stack[last_unsat_clause] = index;

            obj.hard_unsat_nb--;
        }
    else
        {
            last_unsat_clause = mypop(obj.softunsat_stack);
            index = obj.index_in_softunsat_stack[clause];
            obj.softunsat_stack[index] = last_unsat_clause;
            obj.index_in_softunsat_stack[last_unsat_clause] = index;

            // obj.soft_unsat_weight -= obj.org_clause_weight[clause];
        }
}

void check_new_score(externvariables &obj)
{
    long long tem_score = 0;
    long long tem_sscore = 0;
    int sense, c, v, i;
    int weight;
    for (v = 1; v <= obj.num_vars; v++)
        {
            tem_score = 0;
            tem_sscore = 0;
            for (i = 0; i < obj.var_lit_count[v]; ++i)
                {
                    c = obj.var_lit[v][i].clause_num;
                    sense = obj.var_lit[v][i].sense;
                    weight = obj.var_lit[v][i].weight;
                    if (obj.org_clause_weight[c] == obj.top_clause_weight)
                        {
                            if (obj.sat_count[c] < obj.clause_true_lit_thres[c])
                                {
                                    if (sense != obj.cur_soln[v])
                                        {
                                            tem_score += obj.unit_weight[c] * std::min(obj.clause_true_lit_thres[c] - obj.sat_count[c], weight);
                                        }
                                    else
                                        tem_score -= obj.unit_weight[c] * weight;
                                }
                            else if (obj.sat_count[c] >= obj.clause_true_lit_thres[c])
                                {
                                    if (sense == obj.cur_soln[v])
                                        {
                                            tem_score -= obj.unit_weight[c] * std::max(0, obj.clause_true_lit_thres[c] - obj.sat_count[c] + weight);
                                        }
                                }
                        }
                    else
                        {
                            if (obj.sat_count[c] < obj.clause_true_lit_thres[c])
                                {
                                    if (sense != obj.cur_soln[v])
                                        {
                                            tem_sscore += obj.unit_weight[c] * std::min(obj.clause_true_lit_thres[c] - obj.sat_count[c], weight);
                                        }
                                    else
                                        tem_sscore -= obj.unit_weight[c] * weight;
                                }
                            else if (obj.sat_count[c] >= obj.clause_true_lit_thres[c])
                                {
                                    if (sense == obj.cur_soln[v])
                                        {
                                            tem_sscore -= obj.unit_weight[c] * std::max(0, obj.clause_true_lit_thres[c] - obj.sat_count[c] + weight);
                                        }
                                }
                        }
                }
            if (tem_score != obj.hscore[v] || tem_sscore != obj.sscore[v])
                {

                    obj.outputFile << "score is worng in variable " << v << std::endl;
                    obj.outputFile << "tem_score is " << tem_score << std::endl;
                    obj.outputFile << "score function is " << obj.hscore[v] << std::endl;
                    obj.outputFile << "flip num is " << obj.step << std::endl;

                    for (i = 0; i < obj.var_lit_count[v]; ++i)
                        {
                            c = obj.var_lit[v][i].clause_num;
                            sense = obj.var_lit[v][i].sense;
                            weight = obj.var_lit[v][i].weight;
                            obj.outputFile << c << " ";
                        }
                    obj.outputFile << std::endl;
                    exit(0);
                    break;
                }
        }

    int tem_unsat_softweight = 0;
    for (int i = 0; i < obj.num_clauses; ++i)
        {
            if (obj.org_clause_weight[i] == obj.top_clause_weight)
                continue;
            if (obj.sat_count[i] < obj.clause_true_lit_thres[i])
                {
                    tem_unsat_softweight += (obj.clause_true_lit_thres[i] - obj.sat_count[i]);
                }
        }
    if (tem_unsat_softweight != obj.soft_unsat_weight)
        {
            obj.outputFile << "verify softunsat weight wrong " << std::endl;
            exit(0);
        }
}

void start_timing(externvariables &obj)
{
    times(&obj.start_time);
}

double get_runtime(externvariables &obj)
{
    struct tms stop;
    times(&stop);
    return (double)(stop.tms_utime - obj.start_time.tms_utime + stop.tms_stime - obj.start_time.tms_stime) / sysconf(_SC_CLK_TCK);
}
double local_search_runtime()
{
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(now - program_start_time);
    return duration.count();
}

} // namespace lsa