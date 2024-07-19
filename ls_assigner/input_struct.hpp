#pragma once

#include <map>
#include <string>
#include <vector>

namespace lsa {

struct LSTrack
{
    std::string axis;
    int start;
    int step_length;
    int end;
};

struct LSShape
{
    int net_id;
    int task_id;
    double ll_x; 
    double ll_y;
    double ur_x;
    double ur_y;
};

struct LSPanel{
    int layer_id;
    int panel_id;
    int ll_x;
    int ll_y;
    int ur_x;
    int ur_y;
    std::string prefer_direction;
    std::vector<LSTrack> track_list;
    std::vector<LSShape> wire_list;
    std::vector<LSShape> soft_shape_list;
    std::vector<LSShape> hard_shape_list;
};


}  // namespace lsa
