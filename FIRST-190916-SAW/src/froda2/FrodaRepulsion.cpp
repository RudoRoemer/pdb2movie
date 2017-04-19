#include "FrodaRepulsion.h"
#include "PDB.h"
#include "NeighborTable.h"
#include <cmath>

FrodaRepulsion::FrodaRepulsion( const PDB &pdb,
                                const NeighborTable &neighborTable_,
                                double mismatchTol_ ) :
  polarGeometry(true),
  polarHRadius(0.5),
  vdwTol(0.85),
  mismatchTol(mismatchTol_),
  neighborTable(neighborTable_)
{
  setupTypes();
  assignTypes( pdb );
  setupInteractionCutoffLookupTable();
}

FrodaRepulsion::~FrodaRepulsion()
{
}

void FrodaRepulsion::setupTypes() {
  nameOfType.push_back( "H" );
  radiusOfType.push_back( 1.00 ); 
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );  

  //The Hpolar type has a default radius of 1.00.
  //In certain situations it takes on the polarHRadius.
  nameOfType.push_back( "Hpolar" );
  radiusOfType.push_back( 1.00 );
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );  
  
  nameOfType.push_back( "C" );
  radiusOfType.push_back( 1.70 ); 
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );  
  
  nameOfType.push_back( "Cpolar" );
  radiusOfType.push_back( 1.70 ); 
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );  
  
  nameOfType.push_back( "N" );
  radiusOfType.push_back( 1.55 ); 
  chargeOfType.push_back( -1 );
  isPotentialDonorType.push_back( false );  
    
  nameOfType.push_back( "N_withAttachedH" );
  radiusOfType.push_back( 1.55 ); 
  chargeOfType.push_back( -1 );
  isPotentialDonorType.push_back( true );  
    
  nameOfType.push_back( "O" );
  radiusOfType.push_back( 1.40 ); 
  chargeOfType.push_back( -1 );
  
  nameOfType.push_back( "O_withAttachedH" );
  radiusOfType.push_back( 1.40 ); 
  chargeOfType.push_back( -1 );
  isPotentialDonorType.push_back( true );  
  
  nameOfType.push_back( "S" );
  radiusOfType.push_back( 1.80 ); 
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );  
  
  nameOfType.push_back( "P" );
  radiusOfType.push_back( 1.80 ); 
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );  

  nameOfType.push_back( "MG" );
  radiusOfType.push_back( 1.50 ); 
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );  

  nameOfType.push_back( "MN" );
  radiusOfType.push_back( 1.40 ); 
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );  

  nameOfType.push_back( "SI" );
  radiusOfType.push_back( 2.10 ); 
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );  

  nameOfType.push_back( "FE" );
  radiusOfType.push_back( 1.50 ); 
  chargeOfType.push_back( 1 );
  isPotentialDonorType.push_back( false );  

  nameOfType.push_back( "other" );
  radiusOfType.push_back( 1.50 ); 
  chargeOfType.push_back( 0 );
  isPotentialDonorType.push_back( false );  
  
  nTypes = nameOfType.size();
  for ( int i=0; i<nTypes; i++ ) {
    mapNameToType[nameOfType[i]] = i;
  }

}

void FrodaRepulsion::setupInteractionCutoffLookupTable() {
  //NOTE: rInner is a misnomer.
  //What is stored is the interaction cutoff distance for the 
  //repulsive potential energy
  rInnerLookup.resize(nTypes);
  rInnerLookupStorage.resize(nTypes*nTypes);
  for ( int i=0; i<nTypes; i++ ) {
    rInnerLookup[i] = &rInnerLookupStorage[i*nTypes];
  }

  rInner_ForThirdNeighbors_Lookup.resize(nTypes);
  rInner_ForThirdNeighbors_LookupStorage.resize(nTypes*nTypes);
  for ( int i=0; i<nTypes; i++ ) {
    rInner_ForThirdNeighbors_Lookup[i] = &rInner_ForThirdNeighbors_LookupStorage[i*nTypes];
  }

  for ( int t1 = 0; t1 < nTypes; t1++ ) {
    for ( int t2 = 0; t2 < nTypes; t2++ ) {
      double r;
      r = calcInteractionCutoffForTypePair( t1, t2 );
      
      rInnerLookup[t1][t2] = 
        rInnerLookup[t2][t1] = r + mismatchTol;
        
      rInner_ForThirdNeighbors_Lookup[t1][t2] = 
        rInner_ForThirdNeighbors_Lookup[t2][t1] = r * 0.8 + mismatchTol;
    }
  }
}

void FrodaRepulsion::assignTypes( const PDB &pdb ) {
  int nP = pdb.atomLines.size();
  atomTypeLookup.resize(nP);
  int atomType;
  isHydrogen.resize(nP, false);
  for ( int p=0; p<nP; p++) {
    if ( pdb.atomLines[p].element == "H" ) {
      isHydrogen[p] = true;
      // if any of p's neighbors is a N or O, atom is type Hpolar.
      // otherwise, it is type H
      atomType = mapNameToType["H"];
      vector<int>::const_iterator neighbor;
      for ( neighbor = neighborTable[p].begin(); 
            neighbor != neighborTable[p].end();
            neighbor++ ) {
        if ( pdb.atomLines[*neighbor].element == "N" ||
             pdb.atomLines[*neighbor].element == "O"){
          atomType = mapNameToType["Hpolar"];
          break;//break out of for loop.  Type has been set to Hpolar.
        }
      }
    }
    else if ( pdb.atomLines[p].element == "C" ) {
      if ( pdb.atomLines[p].name == " C  " ) atomType = mapNameToType["Cpolar"];
      else {
        int n_charged_neighbors = 0;
        // count the number of charged neighbors.
        vector<int>::const_iterator neighbor;
        for ( neighbor = neighborTable[p].begin(); 
              neighbor != neighborTable[p].end();
              neighbor++ ) {
          if ( pdb.atomLines[*neighbor].element == "N" ||
               pdb.atomLines[*neighbor].element == "O"){
              n_charged_neighbors++;
          }
        }
        if (n_charged_neighbors > 1){ // siteNumber carbon bonded to two Os or Ns gets siteNumber positive charge to balance them
          atomType = mapNameToType["Cpolar"];
        }
        else atomType = mapNameToType["C"];
      }
    }
    else if ( pdb.atomLines[p].element == "N" ) {
      atomType = mapNameToType["N"];
      vector<int>::const_iterator neighbor;
      for ( neighbor = neighborTable[p].begin(); 
            neighbor != neighborTable[p].end();
            neighbor++ ) {
        if ( pdb.atomLines[*neighbor].element == "H" ) {
          atomType = mapNameToType["N_withAttachedH"];
          break;
        }
      }
    }
    else if ( pdb.atomLines[p].element == "O" ) {
      atomType = mapNameToType["O"];
      vector<int>::const_iterator neighbor;
      for ( neighbor = neighborTable[p].begin(); 
            neighbor != neighborTable[p].end();
            neighbor++ ) {
        if ( pdb.atomLines[*neighbor].element == "H" ) {
          atomType = mapNameToType["O_withAttachedH"];
          break;
        }
      }
    }
    else if ( pdb.atomLines[p].element == "S" ) atomType = mapNameToType["S"];
    else if ( pdb.atomLines[p].element == "P" ) atomType = mapNameToType["P"];
    else if ( pdb.atomLines[p].element == "MG" ) atomType = mapNameToType["MG"];
    else if ( pdb.atomLines[p].element == "MN" ) atomType = mapNameToType["MN"];
    else if ( pdb.atomLines[p].element == "SI" ) atomType = mapNameToType["SI"];
    else if ( pdb.atomLines[p].element == "FE" ) atomType = mapNameToType["FE"];
    else {
      atomType = mapNameToType["other"];
      cout << "Sterics Warning: PDB Atom " << pdb.atomLines[p].serial << 
              " has unrecognized element \"" << pdb.atomLines[p].element << "\"" << endl;
    }
    
    atomTypeLookup[p] = atomType;
  }    
  
}

double FrodaRepulsion::calcInteractionCutoffForTypePair( int t1, int t2 ) const {
  double r1 = radiusOfType[t1];
  double r2 = radiusOfType[t2];
  double cutoff;
  if ( polarGeometry && (chargeOfType[t1] * chargeOfType[t2]) < 0  ) {
    if ( nameOfType[t1] == "Hpolar" ) cutoff = (polarHRadius + r2) * vdwTol;
    else if ( nameOfType[t2] == "Hpolar" ) cutoff = (r1 + polarHRadius) * vdwTol;
    else cutoff = ( r1 + r2 ) * vdwTol * 0.92; //tuned to backbone C O in alpha helix
  }
  else if ( polarGeometry &&
        (( isPotentialDonorType[t1] && chargeOfType[t2] == -1 ) ||
         ( isPotentialDonorType[t2] && chargeOfType[t1] == -1 )) ) {
    // we assume that the pair is a donor-acceptor pair
    // participating in a hydrogen bond.
    cutoff = ( r1 + r2 ) * vdwTol * 0.92;
  }
  else cutoff = ( r1 + r2 ) * vdwTol;
  
  return cutoff; 

}

void FrodaRepulsion::printinfo() const {
  int nP = atomTypeLookup.size();
  for ( int p=0; p<nP; p++ ) {
    int t = atomTypeLookup[p];
    cout << p << " " << nameOfType[t] << endl;
  }
}


double FrodaRepulsion::getInteractionCutoff( int p1, int p2 ) const {
  int t1 = atomTypeLookup[p1];
  int t2 = atomTypeLookup[p2];
  double r;
  if ( neighborTable.isThirdNeighbor(p1, p2) )
    r = rInner_ForThirdNeighbors_Lookup[t1][t2]; 
  else r = rInnerLookup[t1][t2];
  return r;  
}

double FrodaRepulsion::getMaxInteractionCutoff( int p ) const {
  int t1 = atomTypeLookup[p];
  double maxr = 0.0;
  for ( int t2 = 0; t2 < nTypes; t2++ ) { 
    double r = rInnerLookup[t1][t2];
    if ( r > maxr ) maxr = r;
  }
  return maxr;
}
