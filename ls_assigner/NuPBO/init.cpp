#include "init.h"
#include "my_class.h"
// void (*init_assignment_ptr)(vector<int> &init_solution);
namespace lsa
{
void init_assignment_false(std::vector<int> &init_solution, externvariables &obj)
{
    int v = 0;
    for (v = 1; v <= obj.num_vars; v++)
        {
            obj.cur_soln[v] = 0;
            obj.time_stamp[v] = 0;
        }
    return;
}
void init_assignment_true(std::vector<int> &init_solution, externvariables &obj)
{
    int v = 0;
    if (init_solution.size() == 0)
        {
            for (v = 1; v <= obj.num_vars; v++)
                {
                    obj.cur_soln[v] = 1;
                    obj.time_stamp[v] = 0;
                }
        }
    else
        {
            for (v = 1; v <= obj.num_vars; v++)
                {
                    obj.cur_soln[v] = init_solution[v];
                    if (obj.cur_soln[v] == 2)
                        obj.cur_soln[v] = rand() % 2;
                    obj.time_stamp[v] = 0;
                }
        }
    return;
}
void init_assignment_random(std::vector<int> &init_solution, externvariables &obj)
{
    int v = 0;
    if (init_solution.size() == 0)
        {
            for (v = 1; v <= obj.num_vars; v++)
                {
                    obj.cur_soln[v] = rand() % 2;
                    obj.time_stamp[v] = 0;
                }
        }
    else
        {
            for (v = 1; v <= obj.num_vars; v++)
                {
                    obj.cur_soln[v] = init_solution[v];
                    if (obj.cur_soln[v] == 2)
                        obj.cur_soln[v] = rand() % 2;
                    obj.time_stamp[v] = 0;
                }
        }
    return;
}
} // namespace lsa