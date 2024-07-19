#include "LSAssigner.h"
#include "NuPBO/nupbo.h"
#include "buildmodel/model.h"
#include "buildmodel/struct.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>

namespace lsa
{
// 预定义一个宏 用MAXDIST_i 访问 maxCDist[RRPairs[i].first][RRPairs[i].second]
#define MAXDIST_i maxCDist[RRPairs[i].first][RRPairs[i].second]

LSPanel LSAssigner::getResult(const LSPanel &input)
{
    PanelInfo panel_info(input);
    panel_info.unique_panel_name = "panel_" + std::to_string(panel_info.input_panel.ls_panel.layer_id) + "_" +
        std::to_string(panel_info.input_panel.ls_panel.panel_id);
    panel_info.TA_panel(params);
    LSPanel output = panel_info.input_panel.ls_panel;
    // Evaluator::EvaluateOverallSolution({output});
    return output;
}

std::vector<LSPanel> LSAssigner::GetResult(const std::vector<LSPanel> &input, const std::string temp_directory_path)
{
    std::vector<LSPanel> output(input.size());
    panel_infos.resize(input.size());
    omp_set_num_threads(nThreads);
#pragma omp parallel for
    for (int i = 0; i < input.size(); ++i)
        {
            PanelInfo temp_panel_info(input[i]);
            panel_infos[i] = temp_panel_info;
            panel_infos[i].unique_panel_name = "panel_" + std::to_string(panel_infos[i].input_panel.ls_panel.layer_id) + "_" +
                std::to_string(panel_infos[i].input_panel.ls_panel.panel_id);
            panel_infos[i].TA_panel(params);
            output[i] = panel_infos[i].input_panel.ls_panel;
        }
    std::string output_file_path = temp_directory_path + "output.txt";
    ResultOutputFile(output_file_path, output);
    Evaluator::EvaluateOverallSolution(output);
    std::cout << "Finished Track Assignment of LSAssigner" << std::endl;
    return output;
}

// 判断返回x和y中的较大值 condition ? value_if_true : value_if_false
double PanelInfo::Max2(const double x, const double y)
{
    return x > y ? x : y;
}

// 根据panel.prefer_direction判断特定方向上是否有可能重叠
bool PanelInfo::MayOverlap(const LSShape &rect1, const LSShape &rect2)
{
    double x_overlap = Max2(0, std::min(rect1.ur_x, rect2.ur_x) - std::max(rect1.ll_x, rect2.ll_x));
    double y_overlap = Max2(0, std::min(rect1.ur_y, rect2.ur_y) - std::max(rect1.ll_y, rect2.ll_y));
    // 如果panel是水平方向，判断x方向上是否有重叠
    if (input_panel.ls_panel.prefer_direction == "H")
        return x_overlap > 0;
    else
        return y_overlap > 0;
}

// 返回一个newRect，根据panel.prefer_direction修改wire的各属性值
LSShape PanelInfo::WireToTrack(const LSShape &wire, const int track_id)
{
    LSShape newRect = wire;
    int trackLine = 0;
    if (input_panel.ls_panel.prefer_direction == "H")
        {
            trackLine = input_panel.yTrack.start + track_id * input_panel.yTrack.step_length;
            int height = wire.ur_y - wire.ll_y;
            newRect.ll_y = trackLine - height / 2;
            newRect.ur_y = trackLine + height / 2;
        }
    else if (input_panel.ls_panel.prefer_direction == "V")
        {
            trackLine = input_panel.xTrack.start + track_id * input_panel.xTrack.step_length;
            int width = wire.ur_x - wire.ll_x;
            newRect.ll_x = trackLine - width / 2;
            newRect.ur_x = trackLine + width / 2;
        }
    return newRect;
}
// 计算重叠面积，分别计算x和y方向上的重叠长度再相乘
double PanelInfo::OverlapArea(const LSShape &rect1, const LSShape &rect2)
{
    double x_overlap = Max2(0, std::min(rect1.ur_x, rect2.ur_x) - std::max(rect1.ll_x, rect2.ll_x));
    double y_overlap = Max2(0, std::min(rect1.ur_y, rect2.ur_y) - std::max(rect1.ll_y, rect2.ll_y));
    return x_overlap * y_overlap;
}
// 返回两个点的中点
Point PanelInfo::MiddlePoint(const Point &a, const Point &b)
{
    return Point{(a.x + b.x) / 2.0, (a.y + b.y) / 2.0};
}
// 将Shape类型的对象转变为Rectangle类型的对象
void PanelInfo::Shape2Rec(const LSShape &shape1, Rectangle &rec1)
{
    rec1.lb.x = shape1.ll_x;
    rec1.lb.y = shape1.ll_y;
    rec1.rt.x = shape1.ur_x;
    rec1.rt.y = shape1.ur_y;
}

// 计算两点的曼哈顿距离
double PanelInfo::Distance(const Point &p1, const Point &p2)
{
    return abs(p1.x - p2.x) + abs(p1.y - p2.y);
}
// 返回s1和s2左下和右上中点之间的曼哈顿距离
double PanelInfo::RPdistance(const Point &s1, const LSShape &s2)
{
    Rectangle pinRect;
    Shape2Rec(s2, pinRect);

    Point pinCenter = MiddlePoint(pinRect.lb, pinRect.rt);

    return Distance(s1, pinCenter);
}
//建模准备函数
Cost PanelInfo::ComputeCost()
{
    int nTracks = input_panel.maxTrackId + 1;
    int nPins = input_panel.ls_panel.soft_shape_list.size();
    int nRoutes = input_panel.ls_panel.wire_list.size();
    int nBlocks = input_panel.ls_panel.hard_shape_list.size();
    Cost mycost(nRoutes, nTracks);

    // 大小为nRoutes ×
    // nRoutes的二维矩阵，元素初始化为fasle，存储每个路线对之间是否存在重叠
    std::vector<std::vector<bool>> mayOverlapRR(nRoutes, std::vector<bool>(nRoutes, false));
    // 大小为nRoutes ×
    // nRoutes的二维矩阵，元素初始化为0，记录每个路线对之间的最大距离
    std::vector<std::vector<int>> maxCDist(nRoutes, std::vector<int>(nRoutes, 0));

    for (int r1 = 0; r1 < nRoutes; r1++)
        {
            for (int r2 = 0; r2 < nRoutes; r2++)
                {
                    if (r1 != r2)
                        {
                            // 判断路径r1和r2是否有可能重叠
                            mayOverlapRR[r1][r2] = MayOverlap(input_panel.ls_panel.wire_list[r1], input_panel.ls_panel.wire_list[r2]);
                            // 如果没有就跳出比较下一个r2
                            if (!mayOverlapRR[r1][r2])
                                continue;

                            // 如果可能重叠修改r1对应的wire成为tmpR1，在随后的循环中，将r2对应的wire修改为tmpR2
                            LSShape tmpR1 = WireToTrack(input_panel.ls_panel.wire_list[r1], 0);
                            for (int d = 0; d < nTracks; ++d)
                                {
                                    LSShape tmpR2 = WireToTrack(input_panel.ls_panel.wire_list[r2], d);
                                    // 计算并确定tmpR1和tmpR2是否有重叠区域，并将结果保存到三维数组overlapRR中
                                    mycost.overlapRR[r1][r2][d] = OverlapArea(tmpR1, tmpR2);
                                    if (mycost.overlapRR[r1][r2][d] > 0.00001)
                                        {
                                            mycost.overlapRR[r1][r2][d] = 1;
                                        }
                                    // 如果两个形状没有重叠，则d-1保存为maxCDist[r1][r2]并跳出内层循环。
                                    if (mycost.overlapRR[r1][r2][d] == 0)
                                        {
                                            maxCDist[r1][r2] = d - 1;
                                            break;
                                        }
                                }
                        }
                }
        }

    // 大小为nRoutes × nTracks ×
    // nPins的三维矩阵，所有元素值都初始化为0，记录每个route和track之间每个pin重叠大小
    std::vector<std::vector<std::vector<int>>> overlapRP(nRoutes, std::vector<std::vector<int>>(nTracks, std::vector<int>(nPins, 0)));
    // 大小为nRoutes x
    // nPins的二维矩阵，元素初始化为false，判断route和pin之间有没重叠
    std::vector<std::vector<bool>> mayOverlapRP(nRoutes, std::vector<bool>(nPins, false));

    for (int r = 0; r < nRoutes; r++)
        {
            for (int p = 0; p < nPins; p++)
                {
                    // 判断路径r和p是否有可能重叠
                    mayOverlapRP[r][p] = MayOverlap(input_panel.ls_panel.wire_list[r], input_panel.ls_panel.soft_shape_list[p]);
                    if (!mayOverlapRP[r][p])
                        continue;
                    // 如果有可能重叠就计算每种route，track，pin对的重叠面积
                    for (int t = 0; t < nTracks; t++)
                        {
                            LSShape tmpRoute = WireToTrack(input_panel.ls_panel.wire_list[r], t);
                            overlapRP[r][t][p] = OverlapArea(tmpRoute, input_panel.ls_panel.soft_shape_list[p]);
                        }
                }
        }

    // 大小nRoutes x nTracks x
    // nBlock的三维数组，初始化为0，记录每个route，track，block之间重叠大小
    std::vector<std::vector<std::vector<int>>> overlapRB(nRoutes, std::vector<std::vector<int>>(nTracks, std::vector<int>(nBlocks, 0)));
    // 大小为nRoutes x
    // nBlocks的二维数组,初始化为false，判断每个route，block是否重叠
    std::vector<std::vector<bool>> mayOverlapRB(nRoutes, std::vector<bool>(nBlocks, false));

    for (int r = 0; r < nRoutes; r++)
        {
            for (int b = 0; b < nBlocks; b++)
                {
                    // 判断路径r和b是否有可能重叠
                    mayOverlapRB[r][b] = MayOverlap(input_panel.ls_panel.wire_list[r], input_panel.ls_panel.hard_shape_list[b]);
                    if (!mayOverlapRB[r][b])
                        continue;
                    // 如果有可能重叠就计算每种route，track，block对的重叠面积
                    for (int t = 0; t < nTracks; t++)
                        {
                            LSShape tmpRoute = WireToTrack(input_panel.ls_panel.wire_list[r], t);
                            overlapRB[r][t][b] = OverlapArea(tmpRoute, input_panel.ls_panel.hard_shape_list[b]);
                        }
                }
        }

    // 设置大小为nTracks x nRoutes 的二维数组存放assign的变量
    std::vector<std::pair<int, int>> assign;
    for (int t = 0; t < nTracks; t++)
        {
            for (int r = 0; r < nRoutes; r++)
                {
                    assign.push_back({t, r});
                }
        }

    // RRPairs 用于存储可能重叠并属于不同网络的路线对
    // homoRRPairs 用于存储可能重叠并且属于相同网络的路线对
    std::vector<std::pair<int, int>> homoRRPairs;
    for (int r1 = 0; r1 < nRoutes; r1++)
        {
            for (int r2 = r1 + 1; r2 < nRoutes; r2++)
                {
                    if (mayOverlapRR[r1][r2] && input_panel.ls_panel.wire_list[r1].net_id != input_panel.ls_panel.wire_list[r2].net_id)
                        {
                            mycost.RRPairs.push_back({r1, r2});
                        }
                    if (mayOverlapRR[r1][r2] && input_panel.ls_panel.wire_list[r1].net_id == input_panel.ls_panel.wire_list[r2].net_id)
                        {
                            homoRRPairs.push_back({r1, r2});
                        }
                }
        }

    // 将成本数据做成字典
    for (auto const &[t, r] : assign)
        {
            for (int b = 0; b < nBlocks; b++)
                {
                    if (overlapRB[r][t][b] > 0.00001)
                        {
                            mycost.blockCost[t][r] += 1;
                        }
                    // mycost.blockCost[t][r] += overlapRB[r][t][b];
                }
            for (int p = 0; p < nPins; p++)
                {
                    if (input_panel.ls_panel.soft_shape_list[p].net_id != input_panel.ls_panel.wire_list[r].net_id)
                        {
                            if (overlapRP[r][t][p] > 0.00001)
                                {
                                    mycost.pinCost[t][r] += 1;
                                }
                            // mycost.pinCost[t][r] += overlapRP[r][t][p];
                        }
                }
        }

    // 检查相对track和route的相对位置
    for (auto const &[t, r] : assign)
        {
            LSShape tmpRoute = WireToTrack(input_panel.ls_panel.wire_list[r], t);
            if (input_panel.ls_panel.soft_shape_list.size() == 0)
                {
                    mycost.distPR[t][r] = 0;
                    break;
                }

            // 根据预设的偏好获取左边和右边的中点坐标
            Point lend, rend;
            Point lt, lb, rt, rb;
            if (input_panel.ls_panel.prefer_direction == "H")
                {
                    lt = {tmpRoute.ll_x, tmpRoute.ur_y};
                    lb = {tmpRoute.ll_x, tmpRoute.ll_y};
                    lend = MiddlePoint(lt, lb);

                    rt = {tmpRoute.ur_x, tmpRoute.ur_y};
                    rb = {tmpRoute.ur_x, tmpRoute.ll_y};
                    rend = MiddlePoint(rt, rb);
                }
            else if (input_panel.ls_panel.prefer_direction == "V")
                {
                    lt = {tmpRoute.ur_x, tmpRoute.ll_y};
                    lb = {tmpRoute.ll_x, tmpRoute.ll_y};
                    lend = MiddlePoint(lt, lb);

                    rt = {tmpRoute.ur_x, tmpRoute.ur_y};
                    rb = {tmpRoute.ll_x, tmpRoute.ur_y};
                    rend = MiddlePoint(rt, rb);
                }

            // 计算panel坐标间曼哈顿距离
            double leftDist = abs(input_panel.ls_panel.ur_x - input_panel.ls_panel.ll_x) + abs(input_panel.ls_panel.ur_y - input_panel.ls_panel.ll_y);
            double rightDist =
                abs(input_panel.ls_panel.ur_x - input_panel.ls_panel.ll_x) + abs(input_panel.ls_panel.ur_y - input_panel.ls_panel.ll_y);
            bool existConn = false;
            // 循环检查 net.id是否一致
            for (int p = 0; p < nPins; p++)
                {
                    if (input_panel.ls_panel.soft_shape_list[p].net_id == input_panel.ls_panel.wire_list[r].net_id)
                        {
                            // 如果一致则设置existConn
                            // 为true，并检查是否更新left、rightDist
                            existConn = true;
                            if (leftDist > RPdistance(lend, input_panel.ls_panel.soft_shape_list[p]))
                                leftDist = RPdistance(lend, input_panel.ls_panel.soft_shape_list[p]);

                            if (rightDist > RPdistance(rend, input_panel.ls_panel.soft_shape_list[p]))
                                rightDist = RPdistance(rend, input_panel.ls_panel.soft_shape_list[p]);
                        }
                }
            // 如果existConn存在，就更新distPR为曼哈顿距离，否则为0
            if (existConn)
                mycost.distPR[t][r] = (leftDist + rightDist);
            else
                mycost.distPR[t][r] = 0;
        }
    return mycost;
}

//建模函数
LSModel PanelInfo::BuildModel(const TAParams &params)
{
    LSModel model = LSModel::newModel();

    int nTracks = input_panel.maxTrackId + 1;
    int nPins = input_panel.ls_panel.soft_shape_list.size();
    int nRoutes = input_panel.ls_panel.wire_list.size();
    int nBlocks = input_panel.ls_panel.hard_shape_list.size();
    int max_decimals = 0;

    Cost mycost = ComputeCost();

    std::vector<std::vector<int>> maxCDist(nRoutes, std::vector<int>(nRoutes, 0));

    std::vector<std::pair<LSVar, double>> obj_var_list;

    // 添加变量
    std::vector<std::vector<LSVar>> x(nTracks, std::vector<LSVar>(nRoutes));
    for (int t = 0; t < nTracks; ++t)
        {
            for (int r = 0; r < nRoutes; ++r)
                {
                    std::string var_name = "x " + std::to_string(t) + " " + std::to_string(r);
                    double obj_constant = params.blockcoef * mycost.blockCost[t][r] + params.pincoef * mycost.pinCost[t][r];
                    // double obj_constant = params.blockcoef *
                    // mycost.blockCost[t][r] + params.pincoef *
                    // mycost.pinCost[t][r] + params.distcoef *
                    // mycost.distPR[t][r]; std::string obj_str =
                    // std::to_string(obj_constant);
                    // obj_str.erase(obj_str.find_last_not_of('0') + 1,
                    // std::string::npos); size_t point = obj_str.find('.');
                    // if(point != std::string::npos)
                    // {
                    //     int decimals = obj_str.size() - point - 1;
                    //     max_decimals = std::max(max_decimals, decimals);
                    // }
                    x[t][r] = model.AddVar(var_name);
                    obj_var_list.push_back(std::make_pair(x[t][r], obj_constant));
                }
        }

    std::vector<std::vector<LSVar>> z;
    for (int i = 0; i < mycost.RRPairs.size(); ++i)
        {
            std::vector<LSVar> z_i;
            for (int d = 0; d <= maxCDist[mycost.RRPairs[i].first][mycost.RRPairs[i].second]; d++)
                {
                    std::string var_name = "z " + std::to_string(i) + " " + std::to_string(d);
                    double obj_constant = params.routecoef * mycost.overlapRR[mycost.RRPairs[i].first][mycost.RRPairs[i].second][d];
                    // std::string obj_str = std::to_string(obj_constant);
                    // obj_str.erase(obj_str.find_last_not_of('0') + 1,
                    // std::string::npos); size_t point = obj_str.find('.');
                    // if(point != std::string::npos)
                    // {
                    //     int decimals = obj_str.size() - point - 1;
                    //     max_decimals = std::max(max_decimals, decimals);
                    // }
                    LSVar z_i_d = model.AddVar(var_name, obj_constant);
                    z_i.push_back(z_i_d);
                }
            z.push_back(z_i);
        }

    // model.max_decimals = max_decimals;

    // // 目标函数设置
    // if (max_decimals != 0)
    // {
    //     for (int i = 0; i < model.var_list.size(); ++i)
    //     {
    //         model.var_list[i].var_obj *= pow(10,max_decimals);
    //     }
    // }
    std::string type = "min";
    model.SetObj(type, obj_var_list);

    // 设置约束
    for (int r = 0; r < nRoutes; r++)
        {
            std::vector<std::pair<LSVar, double>> var_list;
            std::string relation = "eq";
            double rhs = 1;
            for (int t = 0; t < nTracks; t++)
                {
                    var_list.push_back(std::make_pair(x[t][r], 1));
                }
            model.AddCon(var_list, rhs, relation);
        }

    for (int i = 0; i < mycost.RRPairs.size(); ++i)
        {
            for (int d = 0; d <= maxCDist[mycost.RRPairs[i].first][mycost.RRPairs[i].second]; d++)
                {
                    for (int t1 = 0; t1 < nTracks; t1++)
                        {
                            std::string relation = "leq";
                            double rhs = 1;
                            std::vector<std::pair<LSVar, double>> var_list;
                            if (t1 + d < nTracks)
                                {
                                    var_list.push_back(std::make_pair(z[i][d], -1));
                                    var_list.push_back(std::make_pair(x[t1][mycost.RRPairs[i].first], 1));
                                    var_list.push_back(std::make_pair(x[t1 + d][mycost.RRPairs[i].second], 1));
                                }
                            if (d != 0 && t1 - d >= 0)
                                {
                                    var_list.push_back(std::make_pair(z[i][d], -1));
                                    var_list.push_back(std::make_pair(x[t1][mycost.RRPairs[i].first], 1));
                                    var_list.push_back(std::make_pair(x[t1 - d][mycost.RRPairs[i].second], 1));
                                }
                            model.AddCon(var_list, rhs, relation);
                        }
                }
        }
    return model;
}

Trackinfo PanelInfo::GetTrack()
{
    Trackinfo trackinfo;
    for (int i = 0; i < input_panel.ls_panel.track_list.size(); ++i)
        {
            if (input_panel.ls_panel.track_list[i].axis == "X" && input_panel.ls_panel.prefer_direction == "V")
                {
                    trackinfo.track_begin = input_panel.ls_panel.track_list[i].start;
                    trackinfo.track_length = input_panel.ls_panel.track_list[i].step_length;
                }
            if (input_panel.ls_panel.track_list[i].axis == "Y" && input_panel.ls_panel.prefer_direction == "H")
                {
                    trackinfo.track_begin = input_panel.ls_panel.track_list[i].start;
                    trackinfo.track_length = input_panel.ls_panel.track_list[i].step_length;
                }
        }
    return trackinfo;
}

std::vector<int> PanelInfo::GetAssignResults(const Solinfo &solinfo)
{
    int wire_num = input_panel.ls_panel.wire_list.size();
    int track_num = input_panel.maxTrackId + 1;

    std::vector<int> assign(wire_num, 1);

    for (int i = 0; i < solinfo.vars_value.size(); ++i)
        {
            std::stringstream ss(solinfo.vars_value[i].first);
            std::string temp;
            std::vector<int> nums;
            while (ss >> temp)
                {
                    if (temp == "z")
                        {
                            break;
                        }
                    if (std::isdigit(temp[0]))
                        {
                            nums.push_back(std::stoi(temp));
                        }
                }
            if (temp == "z")
                {
                    break;
                }
            assign[nums[1]] = nums[0];
        }
    return assign;
}
// 更新panel的信息
void PanelInfo::UpdatePanel(const Trackinfo &trackinfo)
{
    int track_begin = trackinfo.track_begin;
    int track_length = trackinfo.track_length;

    int index = 0;

    if (input_panel.ls_panel.prefer_direction == "V")
        {
            for (int i = 0; i < input_panel.ls_panel.wire_list.size(); ++i)
                {
                    int d = assign_results[index];
                    int length = input_panel.ls_panel.wire_list[i].ll_x + input_panel.ls_panel.wire_list[i].ur_x;
                    input_panel.ls_panel.wire_list[i].ll_x = track_begin + track_length * d - length / 2;
                    input_panel.ls_panel.wire_list[i].ur_x = track_begin + track_length * d + length / 2;
                    index += 1;
                }
        }
    else if (input_panel.ls_panel.prefer_direction == "H")
        {
            for (int i = 0; i < input_panel.ls_panel.wire_list.size(); ++i)
                {
                    int d = assign_results[index];
                    int length = input_panel.ls_panel.wire_list[i].ll_y + input_panel.ls_panel.wire_list[i].ur_y;
                    input_panel.ls_panel.wire_list[i].ll_y = track_begin + track_length * d - length / 2;
                    input_panel.ls_panel.wire_list[i].ur_y = track_begin + track_length * d + length / 2;
                    index += 1;
                }
        }
}
long long PanelInfo::RandomResult(const Trackinfo &trackinfo, const LSModel &model)
{
    int track_begin = trackinfo.track_begin;
    int track_length = trackinfo.track_length;
    int max_d = input_panel.maxTrackId + 1;
    long long panel_value = 0;

    for (int i = 0; i < input_panel.ls_panel.wire_list.size(); ++i)
        {
            if (input_panel.ls_panel.prefer_direction == "V")
                {
                    int d = rand() % (max_d - 1);
                    int length = input_panel.ls_panel.wire_list[i].ll_x + input_panel.ls_panel.wire_list[i].ur_x;
                    input_panel.ls_panel.wire_list[i].ll_x = track_begin + track_length * d - length / 2;
                    input_panel.ls_panel.wire_list[i].ur_x = track_begin + track_length * d + length / 2;
                    panel_value += model.objective.var_map.at((d * i) + i) / pow(10, model.max_decimals);
                }
            if (input_panel.ls_panel.prefer_direction == "H")
                {
                    int d = rand() % (max_d - 1);
                    int length = input_panel.ls_panel.wire_list[i].ll_y + input_panel.ls_panel.wire_list[i].ur_y;
                    input_panel.ls_panel.wire_list[i].ll_y = track_begin + track_length * d - length / 2;
                    input_panel.ls_panel.wire_list[i].ur_y = track_begin + track_length * d + length / 2;
                    panel_value += model.objective.var_map.at((d * i) + i) / pow(10, model.max_decimals);
                }
        }
    return panel_value;
}

int PanelInfo::TA_panel(const TAParams &api_argu)
{
    std::string modelName = unique_panel_name;
    std::string outfile = unique_panel_name + ".out";
    LSModel model = BuildModel(api_argu);
    Trackinfo trackinfo = GetTrack();
    Solinfo solinfo = model.Solve(modelName, api_argu.sol_time, outfile);
    assign_results = GetAssignResults(solinfo);

    if (solinfo.model_value != -1)
        {
            value_result = solinfo.model_value;
        }
    else
        {
            // value_result = solinfo.model_value;
            value_result = RandomResult(trackinfo, model);
        }
    UpdatePanel(trackinfo);
    return 0;
}
void ResultOutputFile(const std::string output_file, const std::vector<LSPanel> &output)
{
    std::ofstream outFile(output_file);
    if (outFile.is_open())
        {
            for (int i = 0; i < output.size(); ++i)
                {
                    outFile << "panel"
                            << " " << std::fixed << std::setprecision(0) << output[i].layer_id << " " << output[i].panel_id << " " << output[i].ll_x
                            << " " << output[i].ll_y << " " << output[i].ur_x << " " << output[i].ur_y << " " << output[i].prefer_direction;
                    outFile << "\n{";
                    outFile << "\ntrack_list\n";
                    for (int j = 0; j < output[i].track_list.size(); ++j)
                        {
                            outFile << output[i].track_list[j].axis << " " << output[i].track_list[j].start << " "
                                    << output[i].track_list[j].step_length << " " << output[i].track_list[j].end << "\n";
                        }
                    outFile << "wire_list\n";
                    for (int j = 0; j < output[i].wire_list.size(); ++j)
                        {
                            outFile << output[i].wire_list[j].net_id << " " << output[i].wire_list[j].ll_x << " " << output[i].wire_list[j].ll_y
                                    << " " << output[i].wire_list[j].ur_x << " " << output[i].wire_list[j].ur_y << "\n";
                        }
                    outFile << "soft_shape_list\n";
                    for (int j = 0; j < output[i].soft_shape_list.size(); ++j)
                        {
                            outFile << output[i].soft_shape_list[j].net_id << " " << output[i].soft_shape_list[j].ll_x << " "
                                    << output[i].soft_shape_list[j].ll_y << " " << output[i].soft_shape_list[j].ur_x << " "
                                    << output[i].soft_shape_list[j].ur_y << "\n";
                        }
                    outFile << "hard_shape_list\n";
                    for (int j = 0; j < output[i].hard_shape_list.size(); ++j)
                        {
                            outFile << output[i].hard_shape_list[j].net_id << " " << output[i].hard_shape_list[j].ll_x << " "
                                    << output[i].hard_shape_list[j].ll_y << " " << output[i].hard_shape_list[j].ur_x << " "
                                    << output[i].hard_shape_list[j].ur_y << "\n";
                        }
                    outFile << "}\n";
                }
        }
    else
        {
            std::cerr << "Unable output";
        }
    outFile.close();
}
// check panel
bool Evaluator::ValidateSolution(const std::vector<LSPanel> &input)
{
    //目标：wire按照对应方向平移到了track上
    //平移：坐标的相对位置不变或者绝对位置的变化完全一样
    // 保证平移方向和对应轨道正确：h的话横坐标不能变，纵坐标的范围在ytrack范围之内
    // v的话纵坐标不能变，横坐标的范围在xtrack范围之内
    // TODO:现在只给出了得到的解的结果，没有给出初始的解，所以只能保证在track范围之内，不能保证平移方向和对应轨道正确
    int input_size = static_cast<int>(input.size());
    for (int i = 0; i < input_size; i++)
        {
            LSPanel temp_panel = input[i];
            if (input[i].prefer_direction == "H")
                {
                    for (auto const &track : temp_panel.track_list)
                        {
                            if (track.axis == "Y")
                                {
                                    for (auto const &wire : temp_panel.wire_list)
                                        {
                                            double wire_center = (wire.ll_y + wire.ur_y) / 2;
                                            if (wire_center < track.start || wire_center > track.end ||
                                                static_cast<int>(wire_center - track.start) % track.step_length != 0)
                                                {
                                                    std::cout << i << std::endl;
                                                    return false;
                                                }
                                        }
                                }
                        }
                }
            else
                {
                    for (auto const &track : temp_panel.track_list)
                        {
                            if (track.axis == "X")
                                {
                                    for (auto const &wire : temp_panel.wire_list)
                                        {
                                            double wire_center = (wire.ll_x + wire.ur_x) / 2;
                                            if (wire_center < track.start || wire_center > track.end ||
                                                static_cast<int>(wire_center - track.start) % track.step_length != 0)
                                                {
                                                    std::cout << i << std::endl;
                                                    return false;
                                                }
                                        }
                                }
                        }
                }
        }
    return true;
}
// eval panel
void Evaluator::EvaluateSolution(const std::vector<LSPanel> &out)
{
    std::vector<PanelInfo> panel_infos(out.size());
    for (int i = 0; i < out.size(); i++)
        {
            PanelInfo temp_panel_info(out[i]);
            panel_infos[i] = temp_panel_info;
            panel_infos[i].unique_panel_name = "panel_" + std::to_string(panel_infos[i].input_panel.ls_panel.layer_id) + "_" +
                std::to_string(panel_infos[i].input_panel.ls_panel.panel_id);
        }
    // std::vector<Cost> costs(panel_infos.size());
    double total_cost_rr = 0.0;
    double total_cost_rp = 0.0;
    double total_cost_rb = 0.0;
    for (int i = 0; i < panel_infos.size(); i++)
        {
            Cost cost = panel_infos[i].ComputeCost();
            int nTracks = panel_infos[i].input_panel.maxTrackId + 1;
            int nPins = panel_infos[i].input_panel.ls_panel.soft_shape_list.size();
            int nRoutes = panel_infos[i].input_panel.ls_panel.wire_list.size();
            int nBlocks = panel_infos[i].input_panel.ls_panel.hard_shape_list.size();

            std::string prefer_direction = panel_infos[i].input_panel.ls_panel.prefer_direction;
            std::vector<int> assign(nRoutes, 1);
            if (prefer_direction == "V")
                {
                    int track_start = panel_infos[i].input_panel.xTrack.start;
                    int track_length = panel_infos[i].input_panel.xTrack.step_length;
                    for (int r = 0; r < nRoutes; ++r)
                        {
                            int cur_lb = panel_infos[i].input_panel.ls_panel.wire_list[r].ll_x;
                            int cur_rt = panel_infos[i].input_panel.ls_panel.wire_list[r].ur_x;
                            int mid = (cur_lb + cur_rt) / 2;
                            int d = (mid - track_start) / track_length;
                            assign[r] = d;
                        }
                }
            else
                {
                    int track_start = panel_infos[i].input_panel.yTrack.start;
                    int track_length = panel_infos[i].input_panel.yTrack.step_length;
                    for (int r = 0; r < nRoutes; ++r)
                        {
                            int cur_lb = panel_infos[i].input_panel.ls_panel.wire_list[r].ll_y;
                            int cur_rt = panel_infos[i].input_panel.ls_panel.wire_list[r].ur_y;
                            int mid = (cur_lb + cur_rt) / 2;
                            int d = (mid - track_start) / track_length;
                            assign[r] = d;
                        }
                }

            for (int r = 0; r < nRoutes; r++)
                {
                    total_cost_rp += cost.pinCost[assign[r]][r];
                    total_cost_rb += cost.blockCost[assign[r]][r];
                }

            for (int i = 0; i < cost.RRPairs.size(); ++i)
                {
                    int route1 = cost.RRPairs[i].first;
                    int route2 = cost.RRPairs[i].second;
                    if (assign[route1] == assign[route2])
                        {
                            total_cost_rr += cost.overlapRR[route1][route2][assign[route1]];
                        }
                }
        }
    std::cout << "Total costs of wire overlaps: " << total_cost_rr << std::endl;
    std::cout << "Total costs of wire-pin overlaps: " << total_cost_rp << std::endl;
    std::cout << "Total costs of wire-block overlaps: " << total_cost_rb << std::endl;
}
void Evaluator::EvaluateOverallSolution(const std::vector<LSPanel> &input)
{
    if (ValidateSolution(input))
        {
            std::cout << "Valid solution" << std::endl;
            EvaluateSolution(input);
        }
    else
        {
            std::cout << "Invalid solution" << std::endl;
        }
}
} // namespace lsa