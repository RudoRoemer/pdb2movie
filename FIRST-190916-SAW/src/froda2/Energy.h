#ifndef ENERGY_H_
#define ENERGY_H_

#include <list>
#include "Observable.h"

class EnergyTerm
{
public:
  EnergyTerm() {}
  virtual ~EnergyTerm() {}  
  virtual double energy() = 0;
};

class Energy : public Observer {
public:
  Energy();
  virtual ~Energy();
  void addTerm( EnergyTerm* term ) { energyTerms.push_back(term); }
  double calc() {
    if ( !isUpdated ) update();
    return Etot;
  }
  double operator()() { return calc(); }
  
  void update();
  void receiveNotification( Observable *observable );

private:
  std::list<EnergyTerm*> energyTerms;
  bool isUpdated;
  double Etot;
};

#endif /*ENERGY_H_*/
