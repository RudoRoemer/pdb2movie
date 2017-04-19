#ifndef OUTPUTCONFORMERPDB_H_
#define OUTPUTCONFORMERPDB_H_

class PDB;
#include "PerturbRelaxCycle.h"
#include <vector>
#include "Vec3.h"

class OutputConformerPDB
{
public:
  OutputConformerPDB( 
      PDB *pdb_, 
      const std::vector<Vec3> *positions_, 
      const PerturbRelaxCycle *cycle_ );
  virtual ~OutputConformerPDB();
  void setPeriod( int outputPeriod_ ) { outputPeriod = outputPeriod_; }
  void write();
private:
  PDB *pdb;
  const std::vector<Vec3> *positions;
  const PerturbRelaxCycle *cycle;
  int outputPeriod;
};

#endif /*OUTPUTCONFORMERPDB_H_*/
