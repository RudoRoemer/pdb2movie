#ifndef _UNIT_VALUE_H
#define _UNIT_VALUE_H

#include <string>

//
// Simple Unit + Value container
// 
////////////////////////////////////////////////////////////////////////////////
class UnitValue {
public: 
  std::string units; // for now, storing the units as a 
  float value;
  
  UnitValue & operator=(const UnitValue & unitValue);
  
/*  //
  // return true if this is comperable with another unitValue,
  //        false, otherwise
  // 
  //////////////////////////////////////////////////////////////////////////////
  bool isComprable(const UnitValue & unitValue);*/
};

bool operator<(const UnitValue & firstUnitValue,
               const UnitValue & secondUnitValue);

bool operator>(const UnitValue & firstUnitValue,
               const UnitValue & secondUnitValue);

#endif // _UNIT_VALUE_H
