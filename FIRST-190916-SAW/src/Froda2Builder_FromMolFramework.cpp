#include "Froda2Builder_FromMolFramework.h"
#include "MolFramework.h"
#include "NeighborTable.h"
#include "SiteID.h"
#include <vector>
#include "Vec3.h"
#include "RigidUnitList.h"
#include <iostream>
#include "Parameters.h"
#include "RopesAssignment.h"
#include "FrodaRepulsion.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"
#include "ConstraintEnforcingPotentialBuilder.h"
#include "PDB.h"
#include "AmberPrmtop.h"
#include "AtomCategories.h"
#include "SymmetricCEPBuilder.h"

using namespace std;

Froda2Builder_FromMolFramework::
  Froda2Builder_FromMolFramework(
				 MolFramework &structure, Parameters &parameters, 
                                 SymmetryMatrices *symmetryMatrices )
{
  prmtop = NULL;
  
  //Neighbor Table
  NeighborTable *nt = new NeighborTable( structure.total_sites );
  for ( SiteID site1 = 1; site1 <= structure.total_sites; site1++ ) {
    SiteID site2;
    const vector<SiteID> *neighbors = &structure.site_info[site1].neighbor_list;
    for ( size_t i = 0; i < neighbors->size(); i++ ) {
      site2 = (*neighbors)[i];
      nt->insert( site1 - 1, site2 - 1 );
    }
  }
  nt->commit();
  
  //RigidUnitList
  vector<int> mapPtoRC( structure.total_sites );
  for ( SiteID site1 = 1; site1 <= structure.total_sites; site1++){
    mapPtoRC[site1 - 1] = structure.site_info[site1].rigid_label - 1;
  }
  RigidUnitList *rigidUnitList = 
    new RigidUnitList( mapPtoRC, nt );
  
  //Initial Points
  vector<Vec3> initialPoints( structure.total_sites );
  for ( SiteID site1 = 1; site1 <= structure.total_sites; site1++ ) {
    initialPoints[site1 - 1].x = structure.site_info[site1].coords[0];
    initialPoints[site1 - 1].y = structure.site_info[site1].coords[1];
    initialPoints[site1 - 1].z = structure.site_info[site1].coords[2];
  }
  
  rigidUnitSystem = new RigidUnitSystem( 
        rigidUnitList->getAssignment_RUtoPlist(),
        initialPoints );
  
  ConstraintEnforcingPotentialBuilder *cepBuilder;
  // if you want to add more terms to the constraint-enforcing potential
  // based on the information contained in parameters, this is the place
  // to do it
  if (parameters.useSymmetry) {
    cepBuilder = new SymmetricCEPBuilder(rigidUnitSystem,symmetryMatrices);
    // all setup for symmetry-related CEP term is done by constructor
  } else {
    cepBuilder = new ConstraintEnforcingPotentialBuilder(rigidUnitSystem);
    // setup for the default energy terms is below
  }
    
  cepBuilder->setSharedPointsEnergy();

  //Ropes
  RopesAssignment *ropesAssignment = new RopesAssignment;
  for ( size_t i = 0; i < structure.hydrophobic_tethers.size(); i++) {
    if ( structure.hydrophobic_tethers[i].energy >= 
            parameters.energy_cutoff ) continue;
    
    //it's a valid member
    SiteID site1 = structure.hydrophobic_tethers[i].site_1;
    SiteID site2 = structure.hydrophobic_tethers[i].site_2;      
    int p1 = site1 - 1;
    int p2 = site2 - 1;
    if ( rigidUnitSystem->doPointsBelongToSameRigidUnit( p1, p2 ) ) continue;

    double length = structure.site_info[site1].vdw_radius + 
                    structure.site_info[site2].vdw_radius +
                    parameters.ph_radius;
                    
    ropesAssignment->insert( p1, p2, length );
  }  
  cepBuilder->setRopesEnergy( ropesAssignment );
  
  string RCD_PDB_filename = structure.path_name + 
                            structure.base_name + 
                            "_RCD.pdb";
                            
  //Repulsion
  PDB *pdb = new PDB( RCD_PDB_filename );
  double offset;
  if ( parameters.froda2UseGhostTol ) {
    offset = parameters.ghostTol;
  }
  else {
    offset = 0.0;
  }
  
  //repulsion
  if ( parameters.froda2AmberPrmtopFile != "" ) {
    prmtop = new AmberPrmtop( parameters.froda2AmberPrmtopFile );
  } 

  Repulsion *repulsion = NULL;
  if ( parameters.froda2RepulsionType == "froda" ) {
    repulsion = new FrodaRepulsion( *pdb, *nt, offset );
  }
  else if ( parameters.froda2RepulsionType == "amber" && prmtop ) {
    repulsion = new AtomCategories( nt, prmtop );
  }
  else if ( parameters.froda2RepulsionType == "amber" && !prmtop ) {
    cout << "Prmtop not specified (required for Amber Repulsion settings)" << endl;
    exit(0);
  }
  else if ( parameters.froda2RepulsionType != "none" ) {
    cout << "\"" << parameters.froda2RepulsionType << "\" not a valid repulsion type for Froda2" << endl;
    exit(0);
  }
  
  if ( repulsion ) {
    cepBuilder->setOverlapEnergy( repulsion );
  }

  cep = cepBuilder->getConstraintEnforcingPotential();
  delete cepBuilder;
  
}

Froda2Builder_FromMolFramework::~Froda2Builder_FromMolFramework()
{
}

RigidUnitSystem *Froda2Builder_FromMolFramework::
    getRigidUnitSystem() {
  return rigidUnitSystem;
}
ConstraintEnforcingPotential *Froda2Builder_FromMolFramework::getConstraintEnforcingPotential() {
  return cep;
}

PDB *Froda2Builder_FromMolFramework::getPDB() {
  return pdb;
}

AmberPrmtop *Froda2Builder_FromMolFramework::getAmberPrmtop() {
  return prmtop;
}
