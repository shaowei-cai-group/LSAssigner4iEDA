#ifndef MODEL_H
#define MODEL_H

#include "../NuPBO/nupbo.h"
#include <string>
#include "struct.h"
#include <vector>
# include<map>
namespace lsa
{
struct Trackinfo
{
    int track_begin;
    int track_length;
};

class LSVar{
private:
public:
    int var_id;
    std::string var_name = "";
    double var_obj = 0.0; // 在目标函数中的参数
};

class LSConstraint{
private:

public:
    std::string relation;
    double rhs;
    std::vector<std::pair<LSVar, double>> con_var_list;
};


class LSObjFunc{
private:
    
public:
    std::string obj_type;
    std::map<int, double> var_map;
};

class LSModel{
private:
    
public:
    std::vector<LSVar> var_list;
    std::vector<LSConstraint> constraint_list;
    LSObjFunc objective;
    int max_decimals = 0;
    int var_count = 0;
    
    static LSModel newModel()
    {
        LSModel model; 
        return model;
    }
    // 添加变量
    LSVar AddVar(const std::string var_name, const double var_obj);
    LSVar AddVar(const std::string var_name);
    // 添加约束
    void AddCon(const std::vector<std::pair<LSVar, double>> &con_var_list, const double rhs, const std::string relation);
    // 设置目标
    void SetSence(const std::string type);
    void SetObj(const std::string type, const std::vector<std::pair<LSVar, double>> &var_list);
    // 求解并显示结果
    Solinfo Solve(const std::string model_name, const int sol_time, const std::string out_name);
};
} // namespace lsa
#endif 
