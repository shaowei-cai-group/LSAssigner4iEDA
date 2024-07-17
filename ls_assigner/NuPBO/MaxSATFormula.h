#ifndef _MAXSATFORMULA_H_
#define _MAXSATFORMULA_H_

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

namespace lsa
{

class Lit
{
public:
  int weight;
  int var;
  // int sense;
};

class SoftC
{
public:
  std::vector<Lit> clause;
  long long weight;
  int degree;
};

class HardC
{
public:
  std::vector<Lit> clause;
  std::string str;
  long long weight;
  int degree;
};

class MaxSATFormula
{
public:
  int num_vars, num_clauses;
  long long top_clause_weight;
  std::vector<SoftC> sclause;
  std::vector<HardC> hclause;
  std::map<std::string, int> name2var;
  std::map<int, std::string> var2name;
  void build_ins(const char *ifilename, MaxSATFormula &maxsatformula);
};
} // namespace lsa
#endif