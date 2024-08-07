#include "RawClusterv1.h"

#include <phool/phool.h>

#include <cstdlib>
#include <limits>
#include <string>

void RawClusterv1::Reset()
{
  clusterid = 0;
  _z = std::numeric_limits<float>::signaling_NaN();
  _r = std::numeric_limits<float>::signaling_NaN();
  _phi = std::numeric_limits<float>::signaling_NaN();
  _energy = std::numeric_limits<float>::signaling_NaN();
  prop_map.clear();
  towermap.clear();
}

void RawClusterv1::addTower(const RawClusterDefs::keytype twrid, const float etower)
{
  if (towermap.find(twrid) != towermap.end())
  {
    std::cout << "tower 0x" << std::hex << twrid << ", dec: " << std::dec
              << twrid << " already exists, that is bad" << std::endl;
    exit(1);
  }
  towermap[twrid] = etower;
}

void RawClusterv1::identify(std::ostream& os) const
{
  os << "RawClusterv1 ID " << get_id() << " consist of " << getNTowers() << " towers with total energy of " << get_energy() << " GeV ";
  os << "@ (r,phi,z) = (" << get_r() << ", " << get_phi() << ", " << get_z() << "), (x,y,z) = (" << get_x() << ", " << get_y() << ", " << get_z() << ")";

  for (auto i : prop_map)
  {
    PROPERTY prop_id = static_cast<PROPERTY>(i.first);
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    os << "\t" << prop_id << ":\t" << property_info.first << " = \t";
    switch (property_info.second)
    {
    case type_int:
      os << get_property_int(prop_id);
      break;
    case type_uint:
      os << get_property_uint(prop_id);
      break;
    case type_float:
      os << get_property_float(prop_id);
      break;
    default:
      os << " unknown type ";
      break;
    }
    os << std::endl;
  }
}

////! convert cluster location to psuedo-rapidity given a user chosen z-location
// float RawClusterv1::get_eta(const float z) const
//{
//   if (get_r() <= 0) return numeric_limits<float>::signaling_NaN();
//   return asinh((get_z() - z) / get_r());
// }
//
////! convert cluster E_T given a user chosen z-location
// float RawClusterv1::get_et(const float z) const
//{
//   if (get_r() <= 0) return numeric_limits<float>::signaling_NaN();
//   return get_energy() * sin(atan2(get_r(), (get_z() - z)));
// }

bool RawClusterv1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i != prop_map.end();
}

float RawClusterv1::get_property_float(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_float))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_float) << std::endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end())
  {
    return u_property(i->second).fdata;
  }

  return std::numeric_limits<float>::signaling_NaN();
}

int RawClusterv1::get_property_int(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_int))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_int) << std::endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end())
  {
    return u_property(i->second).idata;
  }

  return std::numeric_limits<int>::min();
}

unsigned int
RawClusterv1::get_property_uint(const PROPERTY prop_id) const
{
  if (!check_property(prop_id, type_uint))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_uint) << std::endl;
    exit(1);
  }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i != prop_map.end())
  {
    return u_property(i->second).uidata;
  }

  return std::numeric_limits<unsigned int>::max();
}

void RawClusterv1::set_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id, type_float))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_float) << std::endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void RawClusterv1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id, type_int))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_int) << std::endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}

void RawClusterv1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id, type_uint))
  {
    std::pair<const std::string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    std::cout << PHWHERE << " Property " << property_info.first << " with id "
              << prop_id << " is of type " << get_property_type(property_info.second)
              << " not " << get_property_type(type_uint) << std::endl;
    exit(1);
  }
  prop_map[prop_id] = u_property(value).get_storage();
}
void RawClusterv1::set_et_iso(const float et_iso, const int radiusx10, bool subtracted, bool clusterTower = true)
{
  if (clusterTower)
  {
    if (subtracted)
    {
      switch (radiusx10)
      {
      case 1:
        set_property(prop_et_iso_calotower_sub_R01, et_iso);
        break;
      case 2:
        set_property(prop_et_iso_calotower_sub_R02, et_iso);
        break;
      case 3:
        set_property(prop_et_iso_calotower_sub_R03, et_iso);
        break;
      case 4:
        set_property(prop_et_iso_calotower_sub_R04, et_iso);
        break;
      default:
        std::string warning = "set_et_iso(const int radiusx10, bool subtracted, bool clusterTower) - radius:" + std::to_string(radiusx10) + " has not been defined";
        PHOOL_VIRTUAL_WARN(warning.c_str());
        break;
      }
    }
    else
    {
      switch (radiusx10)
      {
      case 1:
        set_property(prop_et_iso_calotower_R01, et_iso);
        break;
      case 2:
        set_property(prop_et_iso_calotower_R02, et_iso);
        break;
      case 3:
        set_property(prop_et_iso_calotower_R03, et_iso);
        break;
      case 4:
        set_property(prop_et_iso_calotower_R04, et_iso);
        break;
      default:
        std::string warning = "set_et_iso(const int radiusx10, bool subtracted, bool clusterTower) - radius:" + std::to_string(radiusx10) + " has not been defined";
        PHOOL_VIRTUAL_WARN(warning.c_str());
        break;
      }
    }
  }
  else
  {
    PHOOL_VIRTUAL_WARN("set_et_iso(const int radiusx10, bool subtracted, bool clusterTower) - nonclusterTower algorithms have not been defined");
  }
}

unsigned int
RawClusterv1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
  {
    return iter->second;
  }
  return std::numeric_limits<unsigned int>::max();
}

float RawClusterv1::get_et_iso(const int radiusx10 = 3, bool subtracted = false, bool clusterTower = true) const
{
  float r;
  if (clusterTower)
  {
    if (subtracted)
    {
      switch (radiusx10)
      {
      case 1:
        r = get_property_float(prop_et_iso_calotower_sub_R01);
        break;
      case 2:
        r = get_property_float(prop_et_iso_calotower_sub_R02);
        break;
      case 3:
        r = get_property_float(prop_et_iso_calotower_sub_R03);
        break;
      case 4:
        r = get_property_float(prop_et_iso_calotower_sub_R04);
        break;
      default:
        std::string warning = "get_et_iso(const int radiusx10, bool subtracted, bool clusterTower) - radius:" + std::to_string(radiusx10) + " has not been defined";
        PHOOL_VIRTUAL_WARN(warning.c_str());
        r = std::numeric_limits<float>::signaling_NaN();
        break;
      }
    }
    else
    {
      switch (radiusx10)
      {
      case 1:
        r = get_property_float(prop_et_iso_calotower_R01);
        break;
      case 2:
        r = get_property_float(prop_et_iso_calotower_R02);
        break;
      case 3:
        r = get_property_float(prop_et_iso_calotower_R03);
        break;
      case 4:
        r = get_property_float(prop_et_iso_calotower_R04);
        break;
      default:
        std::string warning = "get_et_iso(const int radiusx10, bool subtracted, bool clusterTower) - radius:" + std::to_string(radiusx10) + " has not been defined";
        PHOOL_VIRTUAL_WARN(warning.c_str());
        r = std::numeric_limits<float>::signaling_NaN();
        break;
      }
    }
  }
  else
  {
    PHOOL_VIRTUAL_WARN("get_et_iso(const int radiusx10, bool subtracted, bool clusterTower) - nonclusterTower algorithms have not been defined");
    r = std::numeric_limits<float>::signaling_NaN();
  }
  return r;
}
