# pragma once
# include "struct_from_basis.h"
# include<vector>
#include <sys/times.h>
#include<fstream>

namespace lsa
{

class externvariables
{
	public:
	//basis_pms.h
	int* temp_array;
	temp_var* temp_unsat;
	struct tms start_time;
	
	int* soft_clause_num_index;
	int* hard_clause_num_index;

	int is_print;

	//控制local_search的执行时间
	int time_limit;

	// size of the instance
	int num_vars;	// var index from 1 to num_vars
	int num_clauses; // clause index from 0 to num_clauses-1
	int num_hclauses;
	int num_sclauses;

	// steps and time
	int tries;
	int max_tries;
	unsigned int max_flips;
	unsigned int max_non_improve_flip;
	unsigned int step;

	int cutoff_time;
	double opt_time;

	/* literal arrays */
	lit **var_lit;		  // var_lit[i][j] means the j'th literal of variable i.
	int *var_lit_count;	  // amount of literals of each variable
	lit **clause_lit;	  // clause_lit[i][j] means the j'th literal of clause i.
	int *clause_lit_count; // amount of literals in each clause
	int *clause_true_lit_thres;
	double *avg_clause_coe;

	/* Information about the variables. */
	double *hscore;
	double *sscore;
	long long *time_stamp;
	int **var_neighbor;
	int *var_neighbor_count;
	int *neighbor_flag;
	int *temp_neighbor;

	/* Information about the clauses */
	long long top_clause_weight;
	long long *org_clause_weight;
	double *tune_soft_clause_weight;
	double *unit_weight;
	double *tuned_degree_unit_weight;

	int *sat_count;
	int *sat_var;

	// unsat clauses stack
	int *hardunsat_stack;		  // store the falsified clause number
	int *index_in_hardunsat_stack; // which position is a clause in the unsat_stack
	int hardunsat_stack_fill_pointer;

	int *softunsat_stack;		  // store the falsified clause number
	int *index_in_softunsat_stack; // which position is a clause in the unsat_stack
	int softunsat_stack_fill_pointer;

	// good decreasing variables (dscore>0)
	int *goodvar_stack;
	int goodvar_stack_fill_pointer;
	int *already_in_goodvar_stack;

	/* Information about solution */
	int *cur_soln; // the current assignment, with 1's for True variables, and 0's for False variables
	int *best_soln;
	int best_soln_feasible; // when find a feasible solution, this is marked as 1.
	int local_soln_feasible;
	int hard_unsat_nb;
	long long soft_unsat_weight;
	long long opt_unsat_weight;
	long long best_known;

	// parameters used in algorithm
	float rwprob;
	float rdprob;
	int hd_count_threshold;
	int h_inc;
	int s_inc;
	float initsoftw;

	//heuristic.h
	void (*flip_ptr)(int flipvar,externvariables &obj);
	int (*select_var_after_update_weight_ptr)(externvariables &obj);
	double (*soft_var_greedy_ptr)(int v,externvariables &obj);
	double (*hard_var_greedy_ptr)(int v,externvariables &obj);

	//init.h
	void (*init_assignment_ptr)(std::vector<int> &init_solution,externvariables &obj);

	//控制输出到文件
	std::ofstream outputFile;

	externvariables();
    ~externvariables();
};
} // namespace lsa