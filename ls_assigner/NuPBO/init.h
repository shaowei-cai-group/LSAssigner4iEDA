#ifndef _INIT_H
#define _INIT_H

#include "basis_pms.h"
namespace lsa
{
//extern void (*init_assignment_ptr)(vector<int> &init_solution);
void init_assignment_false(std::vector<int> &init_solution,externvariables &obj);
void init_assignment_true(std::vector<int> &init_solution,externvariables &obj);
void init_assignment_random(std::vector<int> &init_solution,externvariables &obj);
} // namespace lsa
#endif