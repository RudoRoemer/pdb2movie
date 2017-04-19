#ifndef PERTURBRELAXCYCLE_H_
#define PERTURBRELAXCYCLE_H_

#include "Command.h"
#include "CommandList.h"
#include "CommandMemberFunction.h"
#include "Observable.h"
#include "MinimizeSystem.h"

class PerturbRelaxCycle : public Observer
{
public:
  PerturbRelaxCycle( MinimizeSystem *minim_ );
  virtual ~PerturbRelaxCycle();
  void doCycle();

  template <class T>
  void addCommand_perturb( T *t, void (T::*f)() ) {
    commands_perturb.addCommand( new CommandMemberFunction<T>( t, f ) );
  }
  
  template <class T>
  void addCommand_cycleStart( T *t, void (T::*f)() ) {
    commands_cycleStart.addCommand( new CommandMemberFunction<T>( t, f ) );
  }
  
  template <class T>
  void addCommand_cycleEnd( T *t, void (T::*f)() ) {
    commands_cycleEnd.addCommand( new CommandMemberFunction<T>( t, f ) );
  }
  
  template <class T>
  void addCommand_minCycleReceived( T *t, void (T::*f)() ) {
    if ( !enableObsMin ) {
      enableObsMin = true;
      minim->registerObserver( this );
    }
    commands_minCycleReceived.addCommand( new CommandMemberFunction<T>( t, f ) );
  }
  
  unsigned int getCycleCount() const { return count; }
  int getMinCycle() const { return mincycle; }
  void receiveNotification( Observable *obs );

private:
  CommandList commands_perturb;
  CommandList commands_cycleStart;
  CommandList commands_cycleEnd;
  CommandList commands_minCycleReceived;
  MinimizeSystem *minim;
  bool enableObsMin;
  unsigned int count;  
  int mincycle;
};

#endif /*PERTURBRELAXCYCLE_H_*/
