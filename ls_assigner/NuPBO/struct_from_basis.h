# pragma once
namespace lsa
{
struct lit
{
	int clause_num; // clause num, begin with 0
	int var_num;	// variable num, begin with 1
	int weight;
	bool sense;		// 1 for true literals, 0 for false literals.
};
struct temp_var
{
	int var_num;
	int weight;
};
} // namespace lsa