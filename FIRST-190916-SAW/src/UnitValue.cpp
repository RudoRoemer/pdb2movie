#include "UnitValue.h"

UnitValue & UnitValue::operator=(const UnitValue & unitValue) {
  units = unitValue.units;
  value = unitValue.value;
  
  return *this;
}

bool operator<(const UnitValue & firstUnitValue,
               const UnitValue & secondUnitValue) {
  if (firstUnitValue.units == secondUnitValue.units) {
    // TODO - convert the units if they are compatible but aren't identical
    
    return (firstUnitValue.value < secondUnitValue.value);
  }
  
  // 
  
  // ideally, we would return non-comperable rather than false for
  // non-compatable units
  return false;
}


bool operator>(const UnitValue & firstUnitValue,
               const UnitValue & secondUnitValue) {
  if (firstUnitValue.units == secondUnitValue.units) {
    // TODO - convert the units if they are compatible but aren't identical
    
    return (firstUnitValue.value > secondUnitValue.value);
  }
  
  // 
  
  // ideally, we would return non-comperable rather than false for
  // non-compatable units
  return false;
}
