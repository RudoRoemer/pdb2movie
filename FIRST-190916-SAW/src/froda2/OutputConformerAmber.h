#ifndef OUTPUTCONFORMERAMBER_H_
#define OUTPUTCONFORMERAMBER_H_

#include <string>
#include <vector>
#include "Vec3.h"
#include "AmberTrajectory.h"
class PerturbRelaxCycle;
class AmberPrmtop;

class OutputConformerAmber
{
public:
  OutputConformerAmber( const std::vector<Vec3> *meanPositions_,
      const PerturbRelaxCycle *cycle_,
      const AmberPrmtop *prmtop_ );
  virtual ~OutputConformerAmber();
  void setOutputPeriod( int outputPeriod_ ) { outputPeriod = outputPeriod_; }
  void write();

private:
  const std::vector<Vec3> *meanPositions;
  const PerturbRelaxCycle *cycle;
  const AmberPrmtop *prmtop;
  int outputPeriod;
  AmberTrajectory traj;
};

#endif /*OUTPUTCONFORMERAMBER_H_*/
