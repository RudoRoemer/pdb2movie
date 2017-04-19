#ifndef AMBERTRAJECTORY_H_
#define AMBERTRAJECTORY_H_

#include "Vec3.h"
#include <string>
#include <vector>
class AmberPrmtop;

class AmberTrajectory
{
public:
  AmberTrajectory();
  virtual ~AmberTrajectory();
  
  void initializeOutputFile( const std::string &filename_, const AmberPrmtop &prmtop );
  
  void append( const std::vector<Vec3> &coords );
private:
  std::string filename;
  bool readyToAppend;
  int ifbox;
};

#endif /*AMBERTRAJECTORY_H_*/
