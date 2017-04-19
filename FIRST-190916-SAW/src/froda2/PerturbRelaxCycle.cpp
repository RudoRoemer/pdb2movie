#include "PerturbRelaxCycle.h"
#include <algorithm>

using namespace std;

PerturbRelaxCycle::PerturbRelaxCycle(
    MinimizeSystem *minim_ ) :
  minim(minim_),
  enableObsMin( false ),
  count(0)
{
}

PerturbRelaxCycle::~PerturbRelaxCycle()
{
}

void PerturbRelaxCycle::receiveNotification( Observable *obs ) {
  mincycle = minim->getNumCompletedIterations();
  commands_minCycleReceived();
}

void PerturbRelaxCycle::doCycle()
{
  ++count;
  commands_cycleStart();
  commands_perturb();
  
  mincycle = 0;
  minim->minimize();
  mincycle = minim->getNumCompletedIterations();
  commands_cycleEnd();
}
