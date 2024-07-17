#include "settings.h"
#include "heuristic.h"
#include "init.h"
namespace lsa
{
void default_algorithm_settings(externvariables &obj)
{
    algorithm_settings_default(obj);
}
void algorithm_settings_default(externvariables &obj)
{
    obj.init_assignment_ptr = init_assignment_false;
    // init_assignment_ptr = init_assignment_true;
    // init_assignment_ptr = init_assignment_random;

    obj.flip_ptr = flip_with_neighbor;
    obj.select_var_after_update_weight_ptr = select_var_after_update_weight_2;
}
} // namespace lsa