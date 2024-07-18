#include "../ls_assigner/LSAssigner.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <string.h>
#include <vector>
using namespace std;

bool parseInputFile(string input_file, std::vector<ids::LSPanel> &input)
{
    std::ifstream inputFile(input_file);
    // 文件不存在返回错误
    if (!inputFile)
        {
            std::cerr << "Failed to open input file: " << input_file << std::endl;
            return false;
        }

    std::string line;

    std::string current_begin;
    std::string currentState;

    // 逐行读取文件数据
    while (std::getline(inputFile, line))
        {
            if (line.empty() || line[0] == '#')
                {
                    continue;
                }
            std::istringstream iss1(line);
            // 获取该行第一个字符串内容
            iss1 >> current_begin;
            // 如果是panel，则记录以下信息
            if (current_begin == "panel")
                {
                    ids::LSPanel panel;
                    iss1 >> panel.layer_id >> panel.panel_id >> panel.lb_x >> panel.lb_y >> panel.rt_x >> panel.rt_y >> panel.prefer_direction;
                    input.push_back(panel);
                }
            // 如果是各种信息列表，则更新currentState
            else if (current_begin == "track_list" || current_begin == "wire_list" || current_begin == "soft_shape_list" ||
                     current_begin == "hard_shape_list")
                {
                    currentState = current_begin;
                }
            // 跳过{和}
            else if (current_begin == "{" || current_begin == "}")
                continue;
            // 如果panels非空且currentState为"track_list"
            else if (!input.empty() && currentState == "track_list")
                {
                    ids::LSTrack track;
                    std::istringstream iss(line);
                    // 更新track的信息，轴，开始点，每段长度和结束点
                    iss >> track.axis >> track.start >> track.step_length >> track.end;
                    {
                        input.back().track_list.push_back(track);
                    }
                }
            // 如果panels非空且currentState为"wire_list"
            else if (!input.empty() && currentState == "wire_list")
                {
                    ids::LSShape wire;
                    std::istringstream iss(line);
                    // 获取wire的信息，id，左下x，左下y，右上x，右上y
                    while (iss >> wire.net_id >> wire.lb_x >> wire.lb_y >> wire.rt_x >> wire.rt_y)
                        {
                            input.back().wire_list.push_back(wire);
                        }
                }
            // 如果panels非空且currentState为"soft_shape_list"
            else if (!input.empty() && currentState == "soft_shape_list")
                {
                    ids::LSShape shape;
                    std::istringstream iss(line);
                    // 获取shape的信息，net_id，左下x，左下y，右上x，右上y
                    while (iss >> shape.net_id >> shape.lb_x >> shape.lb_y >> shape.rt_x >> shape.rt_y)
                        {
                            input.back().soft_shape_list.push_back(shape);
                        }
                }
            // 如果panels非空且currentState为"hard_shape_list"
            else if (!input.empty() && currentState == "hard_shape_list")
                {
                    ids::LSShape shape;
                    std::istringstream iss(line);
                    // 获取shape的信息，net_id，左下x，左下y，右上x，右上y
                    while (iss >> shape.net_id >> shape.lb_x >> shape.lb_y >> shape.rt_x >> shape.rt_y)
                        {
                            input.back().hard_shape_list.push_back(shape);
                        }
                }
        }
    inputFile.close();
    return true;
}

int main(int argc, char *argv[])
{
    string input_file = argv[1];
    std::vector<ids::LSPanel> input;
    if (input_file != "")
        {
            if(parseInputFile(input_file, input))
            {
                lsa::LSAssigner ls_a;
                string output_file = "";
                std::vector<ids::LSPanel> output = ls_a.GetResult(input, output_file);
            }
        }
    else
        {
            cout << " Missing necessary parameter:input file" << endl;
        }
    return 0;
}