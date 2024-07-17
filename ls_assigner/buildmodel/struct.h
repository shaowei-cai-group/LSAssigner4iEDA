#ifndef STRUCT_H
#define STRUCT_H
#include <vector>
#include <string>
#include <iostream>
#include "../ids.hpp"

namespace lsa
{
struct Panel
{
    int maxTrackId;//1
    ids::LSTrack xTrack;//2
    ids::LSTrack yTrack;//3
    ids::LSPanel ls_panel;

    Panel()
    {}

    Panel(const ids::LSPanel& input_ls_panel)
    {
        ls_panel = input_ls_panel;
        if (ls_panel.prefer_direction == "H")
        {
            for(int i = 0; i < static_cast<int>(ls_panel.track_list.size()); i++)
            {
                ids::LSTrack track = ls_panel.track_list[i];
                if (track.axis == "Y")
                {
                    maxTrackId = (track.end - track.start) / track.step_length;
                    yTrack = track;
                }
            }
        }
        else if (ls_panel.prefer_direction == "V")
        {
            for(int i = 0; i < ls_panel.track_list.size(); i++)
            {
                ids::LSTrack track = ls_panel.track_list[i];
                if (track.axis == "X")
                {
                    maxTrackId = (track.end - track.start) / track.step_length;
                    xTrack = track;
                }
            }
        }
    }

};

struct Point
{
    double x, y;
};

struct Rectangle
{
    Point lb; 
    Point rt; 
};

struct conflictPair
{
    int numPair;
    int numDist;
    double conflictCost;
    double pairCost;
};


struct Cost
{
    std::vector<std::vector<double>> blockCost;
    std::vector<std::vector<double>> pinCost;
    std::vector<std::vector<double>> distPR;
    std::vector<std::vector<std::vector<int>>> overlapRR;
    std::vector<std::pair<int, int>> RRPairs;
    Cost(const int nRoutes, const int nTracks) 
    {
        blockCost = std::vector<std::vector<double>>(nTracks, std::vector<double>(nRoutes, 0.0));
        pinCost = std::vector<std::vector<double>>(nTracks, std::vector<double>(nRoutes, 0.0));
        distPR = std::vector<std::vector<double>>(nTracks, std::vector<double>(nRoutes, 0.0));
        overlapRR = std::vector<std::vector<std::vector<int>>>(nRoutes, std::vector<std::vector<int>>(nRoutes, std::vector<int>(nTracks, 0)));
    }
};

struct Solinfo{
    double model_value;
    std::vector<std::pair<std::string, int>> vars_value;
};

} // namespace lsa
#endif