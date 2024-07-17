#include "my_class.h"
#include "struct_from_basis.h"
namespace lsa
{
externvariables::externvariables()
{
    temp_array = nullptr;
    temp_unsat = nullptr;
    soft_clause_num_index = nullptr;
    hard_clause_num_index = nullptr;
    time_limit = 10;
    num_vars = 0;     // number of variables, var index from 1 to num_vars
    num_clauses = 0;  // number of clauses, clause index from 0 to num_clauses-1
    num_hclauses = 0; // number of hard clauses
    num_sclauses = 0; // number of soft clauses
    cutoff_time = 300;
    tuned_degree_unit_weight = nullptr;
    initsoftw = 0;
}

externvariables::~externvariables()
{
    // delete[] temp_array;
}
} // namespace lsa