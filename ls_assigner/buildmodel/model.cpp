#include "model.h"
#include "struct.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

//调试用
#include <iostream>
#include <ostream>

namespace lsa
{
LSVar LSModel::AddVar(const std::string var_name, const double var_obj)
{
    LSVar var;
    var.var_name = var_name;
    var.var_obj = var_obj;
    var.var_id = var_count;
    ++var_count;
    var_list.push_back(var);
    return var;
}

LSVar LSModel::AddVar(const std::string var_name)
{
    LSVar var;
    var.var_name = var_name;
    var.var_obj = 0;
    var.var_id = var_count;
    ++var_count;
    var_list.push_back(var);
    return var;
}

void LSModel::AddCon(const std::vector<std::pair<LSVar, double>> &con_var_list, const double rhs, const std::string relation)
{
    LSConstraint constraint;
    constraint.con_var_list = con_var_list;
    constraint.rhs = rhs;
    constraint.relation = relation;
    constraint_list.push_back(constraint);
}

void LSModel::SetSence(const std::string type)
{
    std::map<int, double> temp_map;
    for (int i = 0; i < var_list.size(); ++i)
        {
            temp_map[var_list[i].var_id] = var_list[i].var_obj;
        }
    objective.var_map = temp_map;
    objective.obj_type = type;
}

void LSModel::SetObj(const std::string type, const std::vector<std::pair<LSVar, double>> &var_list)
{
    std::map<int, double> temp_map;
    for (int i = 0; i < var_list.size(); ++i)
        {
            temp_map[var_list[i].first.var_id] = var_list[i].second;
        }
    objective.var_map = temp_map;
    objective.obj_type = type;
}

Solinfo LSModel::Solve(const std::string model_name, const int sol_time, const std::string out_name)
{
    Solinfo solinfo;
    std::string opb_name = model_name + ".opb";
    FILE *outfile = fopen(opb_name.c_str(), "w");

    if (outfile != NULL)
        {
            fprintf(outfile, "%s: ", objective.obj_type.c_str());
            for (auto v : objective.var_map)
                {
                    if (v.second >= 0)
                        {
                            fprintf(outfile, "+%.0f x%d ", v.second, v.first);
                        }
                    else
                        {
                            fprintf(outfile, "%.0f x%d ", v.second, v.first);
                        }
                }
            fprintf(outfile, ";\n");

            for (auto &con : constraint_list)
                {
                    for (auto &pair : con.con_var_list)
                        {
                            if (pair.second >= 0)
                                {
                                    fprintf(outfile, "+%.0f x%d ", pair.second, pair.first.var_id);
                                }
                            if (pair.second < 0)
                                {
                                    fprintf(outfile, "%.0f x%d ", pair.second, pair.first.var_id);
                                }
                        }

                    if (con.relation == "eq")
                        {
                            fprintf(outfile, "= %f;", con.rhs);
                        }
                    if (con.relation == "geq")
                        {
                            fprintf(outfile, ">= %f;", con.rhs);
                        }
                    if (con.relation == "leq")
                        {
                            fprintf(outfile, "<= %f;", con.rhs);
                        }

                    fprintf(outfile, ";\n");
                }

            fclose(outfile);

            //集成nupbo求解
            char **pboargv = new char *[5];
            const char *arg1 = "test";
            std::string cfg_input_file = opb_name;
            const char *arg2 = cfg_input_file.c_str();
            const char *arg3 = "2";
            const char *arg4 = std::to_string(sol_time).c_str();
            pboargv[0] = strdup(arg1);
            pboargv[1] = strdup(arg2);
            pboargv[2] = strdup(arg3);
            pboargv[3] = strdup(arg4);
            pboargv[4] = nullptr;
            nuPBOFunction(4, pboargv, out_name);
            for (int k = 0; pboargv[k] != nullptr; ++k)
                {
                    free(pboargv[k]);
                }
            delete[] pboargv;

            // 删除中间 opb文件
            const char *remove_opb = opb_name.c_str();
            remove(remove_opb);

            std::ifstream file(out_name);
            std::string last_line = "last_line";
            std::string second_last_line = "second_last_line";
            std::string line;
            while (std::getline(file, line))
                {
                    if (!line.empty())
                        {
                            second_last_line = last_line;
                            last_line = line;
                        }
                }
            file.close();

            //删除中间out文件
            const char *remove_out = out_name.c_str();
            remove(remove_out);

            // 如果求解器无解
            if (last_line.find("no") != std::string::npos)
                {
                    std::cout << model_name + "no solution" << std::endl;
                    solinfo.model_value = -1;
                }
            // 如果求解器成功求解
            else
                {
                    // 获取目标值
                    size_t value_found = second_last_line.find(":");
                    std::string value_str = second_last_line.substr(value_found + 1);
                    std::stringstream ss(value_str);
                    ss >> solinfo.model_value;
                    if (max_decimals != 0)
                        {
                            int divide = round(pow(10, max_decimals));
                            solinfo.model_value = solinfo.model_value / divide;
                        }

                    // 获取解方案
                    size_t sol_found = last_line.find("s");
                    std::string sol_str = last_line.substr(sol_found + 1);

                    std::istringstream iss(sol_str);
                    std::vector<int> sol_temp;
                    int temp;
                    while (iss >> temp)
                        {
                            sol_temp.push_back(temp);
                        }
                    for (std::vector<LSVar>::size_type i = 0; i < var_list.size(); ++i)
                        {
                            if (sol_temp[i] == 1)
                                {
                                    solinfo.vars_value.push_back(std::make_pair(var_list[i].var_name, sol_temp[i]));
                                }
                            else
                                {
                                    continue;
                                }
                        }
                }
        }
    else
        {
            std::cerr << "Unable to get opb_file";
        }
    return solinfo;
}

} // namespace lsa