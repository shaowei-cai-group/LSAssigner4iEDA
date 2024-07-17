#include "parse_arguments.h"
#include "basis_pms.h"
#include "cmdline.h"
#include "heuristic.h"
#include "init.h"
#include "my_class.h"
namespace lsa
{
bool parse_args(int argc, char *argv[], Config &cfg)
{
    bool ret = false;
    do
        {
            if (argc < 3 || NULL == argv)
                {
                    break;
                }

            cfg.input_file = argv[1];

            if (!StringUtil::String2Uint(argv[2], cfg.seed))
                {
                    break;
                }

            ret = true;
        }
    while (false);

    return ret;
}

bool parse_arguments(int argc, char **argv, externvariables &obj)
{
    /*cmdline:: parser a;
    a.add<int>("init_assignment", 'i', "init assignment", false, 1,
    cmdline::range(0,10));
    */
    if (sscanf(argv[3], "%d", &obj.time_limit) != 1)
        {
            std::cerr << "Invalid value for time_limit." << std::endl;
            return 1;
        }

    // 输出 time_limit 的值
    // std::obj.outputFile << "time_limit: " << obj.time_limit << std::endl;
    int i;
    for (i = 4; i < argc; i++)
        {
            if (0 == strcmp(argv[i], "-init_assign"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    if (0 == strcmp(argv[i], "1"))
                        {
                            obj.init_assignment_ptr = init_assignment_random;
                        }
                    else if (0 == strcmp(argv[i], "2"))
                        {
                            obj.init_assignment_ptr = init_assignment_true;
                        }
                    else if (0 == strcmp(argv[i], "3"))
                        {
                            obj.init_assignment_ptr = init_assignment_false;
                        }
                }
            else if (0 == strcmp(argv[i], "-best"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    sscanf(argv[i], "%lld", &obj.best_known);
                    obj.outputFile << "c best_known: " << obj.best_known << std::endl;
                    if (obj.best_known == -1)
                        {
                            obj.outputFile << "c no feasible solution" << std::endl;
                            exit(0);
                        }
                }
            else if (0 == strcmp(argv[i], "-rwprob"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    sscanf(argv[i], "%f", &obj.rwprob);
                    obj.outputFile << "c rwprob: " << obj.rwprob << std::endl;
                }
            else if (0 == strcmp(argv[i], "-initsoftw"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    sscanf(argv[i], "%f", &obj.initsoftw);
                    obj.outputFile << "c initsoftw: " << obj.initsoftw << std::endl;
                }
            else if (0 == strcmp(argv[i], "-hvargreedy"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    if (0 == strcmp(argv[i], "0"))
                        {
                            obj.hard_var_greedy_ptr = var_greedy_hscore;
                        }
                    else if (0 == strcmp(argv[i], "1"))
                        {
                            obj.hard_var_greedy_ptr = var_greedy_sscore;
                        }
                    else if (0 == strcmp(argv[i], "2"))
                        {
                            obj.hard_var_greedy_ptr = var_greedy_score;
                        }
                }
            else if (0 == strcmp(argv[i], "-svargreedy"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    if (0 == strcmp(argv[i], "0"))
                        {
                            obj.soft_var_greedy_ptr = var_greedy_hscore;
                        }
                    else if (0 == strcmp(argv[i], "1"))
                        {
                            obj.soft_var_greedy_ptr = var_greedy_sscore;
                        }
                    else if (0 == strcmp(argv[i], "2"))
                        {
                            obj.soft_var_greedy_ptr = var_greedy_score;
                        }
                }
            else if (0 == strcmp(argv[i], "-rdprob"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    sscanf(argv[i], "%f", &obj.rdprob);
                    obj.outputFile << "c rdprob: " << obj.rdprob << std::endl;
                }
            else if (0 == strcmp(argv[i], "-hinc"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    sscanf(argv[i], "%d", &obj.h_inc);
                    obj.outputFile << "c h_inc: " << obj.h_inc << std::endl;
                }
            else if (0 == strcmp(argv[i], "-sinc"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    sscanf(argv[i], "%d", &obj.s_inc);
                    obj.outputFile << "c s_inc: " << obj.s_inc << std::endl;
                }
            else if (0 == strcmp(argv[i], "-bms"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    sscanf(argv[i], "%d", &obj.hd_count_threshold);
                    obj.outputFile << "c hd_count_threshold: " << obj.hd_count_threshold << std::endl;
                }
            else if (0 == strcmp(argv[i], "-afterupdate"))
                {
                    i++;
                    if (i >= argc)
                        return false;
                    if (0 == strcmp(argv[i], "1"))
                        {
                            obj.select_var_after_update_weight_ptr = select_var_after_update_weight_1;
                        }
                    else if (0 == strcmp(argv[i], "2"))
                        {
                            obj.select_var_after_update_weight_ptr = select_var_after_update_weight_2;
                        }
                }
        }
    return true;
}
} // namespace lsa