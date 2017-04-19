#include "OutputConformerAmber.h"
#include "PerturbRelaxCycle.h"
#include "AmberPrmtop.h"
#include <sstream>

using namespace std;

OutputConformerAmber::OutputConformerAmber(
    const std::vector<Vec3> *meanPositions_,
    const PerturbRelaxCycle *cycle_,
    const AmberPrmtop *prmtop_ ) : 
      meanPositions(meanPositions_),
      cycle(cycle_),
      prmtop(prmtop_),
      outputPeriod(0)
{
  string filename = "mdcrd";
  traj.initializeOutputFile( filename.c_str(), *prmtop );
}

OutputConformerAmber::~OutputConformerAmber()
{
}

void OutputConformerAmber::write() {
  if ( !outputPeriod || cycle->getCycleCount()%outputPeriod ) return;
  traj.append( *meanPositions );  
}
