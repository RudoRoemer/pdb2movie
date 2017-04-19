#include "OutputConformerPDB.h"
#include "PDB.h"
#include "PerturbRelaxCycle.h"
#include <sstream>
#include <iostream>
#include <iomanip>

OutputConformerPDB::OutputConformerPDB( 
    PDB *pdb_,
    const std::vector<Vec3> *positions_,
    const PerturbRelaxCycle *cycle_ ) : 
      pdb(pdb_),
      positions(positions_),
      cycle(cycle_),
      outputPeriod(0)
{
}

OutputConformerPDB::~OutputConformerPDB()
{
}

void OutputConformerPDB::write()
{
  if ( !outputPeriod || cycle->getCycleCount()%outputPeriod  ) return;
  stringstream ss;
  string outfilename;
  
  pdb->setPositions( *positions );

  ss << "iter" << setfill('0') << setw(8) << cycle->getCycleCount() << ".pdb";
  outfilename = ss.str();

  pdb->write( outfilename );
}
