#ifndef PROEDA_H
#define PROEDA_H

#include <string>
#include<vector>
#include<ctime>
#include"buildmodel/struct.h"
#include"buildmodel/model.h"
#include "input_struct.hpp"

namespace lsa
{
//建模和求解需要的参数
class TAParams
{
    public:
    // 需要的参数
    int sol_time = 10;
    double blockcoef = 500;
    double pincoef = 2;
    double routecoef = 1;
    double distcoef = 0.1;
    int solve_option = 1;
    int timeLimit = 60;
};

//单个panel的处理
class PanelInfo
{
    public:
    Panel input_panel;
    std::string unique_panel_name = "panel" + std::to_string(input_panel.ls_panel.layer_id) + "_" + std::to_string(input_panel.ls_panel.panel_id); 
    std::vector<int> assign_results;
    long long value_result;

    PanelInfo() {}
    PanelInfo(const LSPanel &ls_panel) 
    {
        Panel temp_panel(ls_panel);
        input_panel = temp_panel;
    }

    //建模相关函数
    double Max2(const double x, const double y);
    bool MayOverlap(const LSShape &rect1, const LSShape &rect2);
    LSShape WireToTrack(const LSShape &wire, const int track_id);
    double OverlapArea(const LSShape &rect1, const LSShape &rect2);
    Point MiddlePoint(const Point &a, const Point &b);
    void Shape2Rec(const LSShape &shape1,Rectangle &rec1);
    double Distance(const Point &p1, const Point &p2);
    double RPdistance(const Point &s1, const LSShape &s2);
    Cost ComputeCost();
    LSModel BuildModel(const TAParams &params);

    //轨道分配相关函数
    Trackinfo GetTrack();
    std::vector<int> GetAssignResults(const Solinfo &solinfo);

    // 更新panel的信息
    void UpdatePanel(const Trackinfo &trackinfo);
    // 对未解成功的panel随机分配
    long long RandomResult(const Trackinfo& trackinfo,const LSModel &model);

    int TA_panel(const TAParams &api_argu);

};
//panel数组的处理
class LSAssigner
{
    public:
    TAParams params;
    int nThreads =128;
    std::vector<PanelInfo> panel_infos;
    std::string output_file="output_file.txt";
    std::string value_file = "value_file.txt";
    std::vector<LSPanel> GetResult(const std::vector<LSPanel> &input,const std::string temp_directory_path);
};  

class Evaluator
{
    public:
    static bool ValidateSolution(const std::vector<LSPanel>& input);
    static void EvaluateSolution(const std::vector<LSPanel> &input);
    static void EvaluateOverallSolution(const std::vector<LSPanel> &input);
};
void ResultOutputFile(const std::string output_file, const std::vector<LSPanel> &output);
} // namespace lsa

#endif