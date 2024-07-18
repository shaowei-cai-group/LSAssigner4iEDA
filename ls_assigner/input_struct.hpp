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
    double lb_x; 
    double lb_y;
    double rt_x;
    double rt_y;
};

struct LSPanel{
    int layer_id;
    int panel_id;
    int lb_x;
    int lb_y;
    int rt_x;
    int rt_y;
    std::string prefer_direction;
    std::vector<LSTrack> track_list;
    std::vector<LSShape> wire_list;
    std::vector<LSShape> soft_shape_list;
    std::vector<LSShape> hard_shape_list;
};


}  // namespace lsa
