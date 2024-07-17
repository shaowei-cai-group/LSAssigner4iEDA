// ***************************************************************************************
// Copyright (c) 2023-2025 Peng Cheng Laboratory
// Copyright (c) 2023-2025 Institute of Computing Technology, Chinese Academy of Sciences
// Copyright (c) 2023-2025 Beijing Institute of Open Source Chip
//
// iEDA is licensed under Mulan PSL v2.
// You can use this software according to the terms and conditions of the Mulan PSL v2.
// You may obtain a copy of Mulan PSL v2 at:
// http://license.coscl.org.cn/MulanPSL2
//
// THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
// EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
// MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
//
// See the Mulan PSL v2 for more details.
// ***************************************************************************************
#pragma once

#include <map>
#include <string>
#include <vector>
namespace ieda {
class Log;
class Str;
class Time;
}  // namespace ieda
namespace icts {
struct PathInfo;
enum class LayerPattern;
class Node;
class Pin;
class Net;
class Inst;
enum class TopoType;
class CtsCellLib;
class CtsLibs;
class CtsReportTable;
class CtsLog;
class CtsInstance;
class CtsConfig;
class CtsDesign;
class CtsDBWrapper;
class CtsPin;
class CtsSignalWire;
class CtsNet;
class EvalNet;
class Evaluator;
class Endpoint;
enum class FitType;
class ModelBase;
class ModelFactory;
template <typename T>
class CtsPoint;
using Point = icts::CtsPoint<int>;
}  // namespace icts

namespace ista {
class TimingEngine;
class TimingIDBAdapter;
class DesignObject;
class RctNode;
class Net;
}  // namespace ista
namespace idb {
class IdbPin;
class IdbInstance;
class IdbNet;
}  // namespace idb

namespace ito {
class Tree;
}  // namespace ito
namespace idrc {

class RegionQuery;
class DrcRect;

}  // namespace idrc

namespace irt {

enum class Stage;
class PlanarRect;
class LayerRect;

template <typename T>
class Segment;

template <typename T>
class MTree;

class Net;
class AccessPoint;
class GRNode;

}  // namespace irt

namespace eval {

class TimingNet;
class TimingPin;
class TileGrid;

}  // namespace eval

namespace ids {

struct Segment
{
  int first_x;
  int first_y;
  std::string first_layer_name;
  int second_x;
  int second_y;
  std::string second_layer_name;
};

enum class PHYNodeType
{
  kNone = 0,
  kWire = 1,
  kVia = 2
};

struct Wire
{
  int first_x;
  int first_y;
  int second_x;
  int second_y;
  std::string layer_name;
};

struct Via
{
  std::string via_name;
  int x;
  int y;
};

struct PHYNode
{
  PHYNodeType type;
  Wire wire;
  Via via;
};

enum class AccessPointType
{
  kNone = 0,
  kPrefTrackGrid = 1,
  kCurrTrackGrid = 2,
  kCurrTrackCenter = 3,
  kCurrShapeCenter = 4
};

struct AccessPoint
{
  int x;
  int y;
  std::string layer_name;
  AccessPointType type;
  std::vector<std::string> via_name_list;
};

struct CellMaster
{
  std::string cell_master_name;
  std::map<std::string, std::vector<ids::AccessPoint>> pin_name_pa_list;
};

struct DRCTask
{
  idrc::RegionQuery* region_query;
  std::vector<idrc::DrcRect*> drc_rect_list;
};

struct DRCRect
{
  int32_t so_id = -1;  // Distinguish self and others

  int32_t lb_x = -1;
  int32_t lb_y = -1;

  int32_t rt_x = -1;
  int32_t rt_y = -1;

  std::string layer_name;

  bool is_artificial = false;
};


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


}  // namespace ids
