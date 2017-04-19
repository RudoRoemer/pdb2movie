#include "global_defs.h"
#include "interact.h"
#include "Froda.h"
#include "mt19937ar.h"
#include "flexweb.h"
#include "Spanner.h"
#include "Timer.h"
#include "Sterics.h"
#include "Grid.h"

#include <signal.h>

extern Parameters parameters;

bool Froda::allPairConstraintsSatisfied() {
  size_t nPairs = myPairs.size();
  int atom1;
  int atom2;
  double AB;
  double R;
  Vec3 diff;
  
  for ( size_t i = 0; i < nPairs; i++ ) {
    atom1 = myPairs[i].atomA;
    atom2 = myPairs[i].atomB;
    AB = myPairs[i].AB;
    diff = currentPos[atom2];
    diff -= currentPos[atom1];
    R = sqrt( diff.norm2() );
    if ( AB > 0 && R > AB ) return false;
    if ( AB < 0 && R < abs(AB) ) return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Description: returns true if two atoms are covalently bonded
//              by comparing their ghost membership
//              input is the lists of ghosts for each atom
////////////////////////////////////////////////////////////////////////////////
bool Froda::bonded( vector<unsigned int> &vec1, vector<unsigned int> &vec2 ){
  int size1 = vec1.size();
  int size2 = vec2.size();
  for ( int i = 0; i < size1; i++){
    for ( int j=0; j < size2; j++){
      if (vec1.at(i) == vec2.at(j)) return true;
    }
  }

  return false;
}

////////////////////////////////////////////////////////////////////////////////
// Description: says if two atoms might be donor/acceptor to an H bond
//              input is the two atoms and their neighbor lists
//              checks for intervening hydrogens
////////////////////////////////////////////////////////////////////////////////
bool Froda::hasHBondGeometry( unsigned int atom1, unsigned int atom2, vector< unsigned int > &nlist1, vector<unsigned int> &nlist2 ){


  Vector hcheck, abcheck;
  double hlength, ablength;
  double margin = 0.5; //margin in Angstroms
  int nlistsize;
  nlistsize = nlist1.size();
  double HAmax2 = 6.25; //2.5 squared

  abcheck.x = currentPos.at(atom1).x - currentPos.at(atom2).x;
  abcheck.y = currentPos.at(atom1).y - currentPos.at(atom2).y;
  abcheck.z = currentPos.at(atom1).z - currentPos.at(atom2).z;

  ablength = sqrt(dotProduct( abcheck, abcheck) );

  for ( int j = 0; j < nlistsize; j++){
    int possible = nlist1.at(j);
    if ( frodaAtom.at(possible).isHydrogen ){
      hcheck.x = currentPos.at(possible).x - currentPos.at(atom2).x;
      hcheck.y = currentPos.at(possible).y - currentPos.at(atom2).y;
      hcheck.z = currentPos.at(possible).z - currentPos.at(atom2).z;
      hlength = dotProduct(hcheck, hcheck);
      //if (hlength < HAmax2 ) return true;
      if (hlength < HAmax2 ){
        hlength = sqrt(hlength);
        if ( ablength - hlength > margin) return true;
      }
    }
  }
  nlistsize = nlist2.size();
  for ( int j = 0; j < nlistsize; j++){
    int possible = nlist2.at(j);
    if ( frodaAtom.at(possible).isHydrogen ){
      hcheck.x = currentPos.at(possible).x - currentPos.at(atom1).x;
      hcheck.y = currentPos.at(possible).y - currentPos.at(atom1).y;
      hcheck.z = currentPos.at(possible).z - currentPos.at(atom1).z;
      hlength = dotProduct(hcheck, hcheck);
      //if (hlength < HAmax2 ) return true;
      if (hlength < HAmax2 ){
        hlength = sqrt(hlength);
        if ( ablength - hlength > margin) return true;
      }
    }
  }
  return false;

}

////////////////////////////////////////////////////////////////////////////////
// Description: check for "blocking" of long-range interactions
//              atoms are 1,2 and n (neighbor of 1)
//              distances are 1-2 and n-2
//              if (1-2 minus n-2) > margin, then n blocks 1 and 2
////////////////////////////////////////////////////////////////////////////////
bool Froda::isBlocking( unsigned int atom1, unsigned int neighbor, unsigned int atom2, double margin ) {
  //checks if atom neighbor (bonded to atom1) blocks atom1 and atom2
  //used for polar/phobic routine; need LOTS of blocking

  Vector delta;
  double d12=0;
  double dn2=0;

  delta = currentPos.at(atom1);
  delta.x -= currentPos.at(atom2).x;
  delta.y -= currentPos.at(atom2).y;
  delta.z -= currentPos.at(atom2).z;

  d12 = sqrt(dotProduct( delta, delta ) );

  delta = currentPos.at(neighbor);
  delta.x -= currentPos.at(atom2).x;
  delta.y -= currentPos.at(atom2).y;
  delta.z -= currentPos.at(atom2).z;

  dn2 = sqrt(dotProduct( delta, delta ) );

  if ( (d12 - dn2) > margin ) {
    return true;
  }
  else return false;

}


////////////////////////////////////////////////////////////////////////////////
// Description: returns a unit vector normal to the plane of an sp2 center
//   if isFlat comes back false, not an sp2 center, vector is NULL_VEC
////////////////////////////////////////////////////////////////////////////////
Vector Froda::getNormalVector( unsigned int atom1, bool &isFlat, MolFramework &structure ) {
  int nCovN = structure.site_info[atom1].nCovalentNeighbors;
  if ( nCovN < 2 || nCovN > 3 ) {
    isFlat = false; //can't get plane
    return NULL_VEC;
  }
  isFlat = true;
  Vector pos1 = currentPos[atom1]; 
  Vector bondA, bondB, bondC;
  unsigned int atomA, atomB, atomC;
  atomA = structure.site_info[atom1].neighbor_list[0];
  bondA = currentPos[atomA];
  bondA.x -= pos1.x;
  bondA.y -= pos1.y;
  bondA.z -= pos1.z;
  atomB = structure.site_info[atom1].neighbor_list[1];
  bondB = currentPos[atomB];
  bondB.x -= pos1.x;
  bondB.y -= pos1.y;
  bondB.z -= pos1.z;
  if ( nCovN == 3 ) {
    atomC = structure.site_info[atom1].neighbor_list[2];
    bondC = currentPos[atomC];
    bondC.x -= pos1.x;
    bondC.y -= pos1.y;
    bondC.z -= pos1.z;
  }

  Vector myNormal = NULL_VEC;
  Vector temp;
  double dot;
  temp = crossProduct( bondA, bondB );
  dot = sqrt(dotProduct( temp, temp ) );
  if ( dot < epsilon ) {
    temp = NULL_VEC;
  }
  else {
    temp.x /= dot;
    temp.y /= dot;
    temp.z /= dot;
  }
  myNormal = temp;

  if ( nCovN == 3 ) {
    temp = crossProduct( bondB, bondC );
    dot = sqrt(dotProduct( temp, temp ) );
    if ( dot < epsilon ) {
      temp = NULL_VEC;
    }
    else {
      temp.x /= dot;
      temp.y /= dot;
      temp.z /= dot;
    }
    myNormal.x += temp.x;
    myNormal.y += temp.y;
    myNormal.z += temp.z;

    temp = crossProduct( bondC, bondA );
    dot = sqrt(dotProduct( temp, temp ) );
    if ( dot < epsilon ) {
      temp = NULL_VEC;
    }
    else {
      temp.x /= dot;
      temp.y /= dot;
      temp.z /= dot;
    }
    myNormal.x += temp.x;
    myNormal.y += temp.y;
    myNormal.z += temp.z;

    dot = sqrt(dotProduct( myNormal, myNormal ) );
    myNormal.x /= dot;
    myNormal.y /= dot;
    myNormal.z /= dot;
  }

  return myNormal;
}

double Froda::getDist2( SiteID atom1, SiteID atom2 ) {
  double dr;
  double value = 0;
  dr = currentPos[atom1].x - currentPos[atom2].x;
  value += dr * dr;
  dr = currentPos[atom1].y - currentPos[atom2].y;
  value += dr * dr;
  dr = currentPos[atom1].z - currentPos[atom2].z;
  value += dr * dr;
  return value;
}

double Froda::computeAngle( SiteID atom1, SiteID atom2, SiteID atom3 ) {
  double dist1_2 = getDist2( atom1, atom2 );
  double dist2_2 = getDist2( atom2, atom3 );
  double dist3_2 = getDist2( atom1, atom3 );

  double angle = ( - dist3_2 + dist1_2 + dist2_2 ) / ( 2 * sqrt( dist1_2 ) * sqrt( dist2_2 ) );
  return ( acos( angle ) );
}

bool Froda::isSaltBridge( SiteID hydrogen, SiteID acceptor, MolFramework &structure ) {
  SiteID donor = structure.site_info[hydrogen].neighbor_list[0];
  if ( structure.site_info[donor].hbond_status < 11 ) return false; //uncharged
  if ( structure.site_info[acceptor].hbond_status < 13 ) return false; //uncharged

  if ( getDist2( hydrogen, acceptor) > pow( parameters.cutoff_SB_hyd_accpt_dist, 2 ) ){
    return false; //too far away
  } 
  if ( getDist2( donor, acceptor) > pow( parameters.cutoff_SB_donor_accpt_dist, 2 ) ){
    return false; //too far away
  } 
  if (  RAD_TO_DEG * computeAngle( donor, hydrogen, acceptor )  <
        parameters.cutoff_SB_donor_hyd_accpt_angle ) {
    return false; //too tight angle
  }
 
  for ( int  b=0; b < structure.site_info[acceptor].nCovalentNeighbors; b++ ) {
    SiteID baseAtom = structure.site_info[acceptor].neighbor_list[b];
    if (  RAD_TO_DEG * computeAngle( hydrogen, acceptor, baseAtom ) <
          parameters.cutoff_SB_hyd_accpt_base_angle ) {
      return false; //too tight angle
    }
  } 

  return( true );
}

bool Froda::maybeHBond( SiteID hydrogen, SiteID acceptor, MolFramework &structure ) {
  SiteID donor = structure.site_info[hydrogen].neighbor_list[0];

  if ( getDist2( hydrogen, acceptor) > pow(parameters.cutoff_HB_hyd_accpt_dist,2) ) {
    return false; //too far
  }
  if ( getDist2( donor, acceptor) > pow(parameters.cutoff_HB_donor_accpt_dist,2) ) {
    return false; //too far
  }
  if (  RAD_TO_DEG * computeAngle( donor, hydrogen, acceptor )  <
        parameters.cutoff_HB_donor_hyd_accpt_angle ) {
    return false; //too tight angle
  }

  return( true );
}

int Froda::HBType( SiteID hydrogen, SiteID acceptor, MolFramework &structure ){
  //type 1 is sp2-sp2
  //type 2 is 3-2
  //type 3 is 2-3
  //type 4 is 3-3
  SiteID donor = structure.site_info[hydrogen].neighbor_list[0];
  int type = 0;

  if( structure.site_info[donor].hbond_status == DONOR_SP2 ||
      structure.site_info[donor].hbond_status == DONOR_AND_ACCEPTOR_SP2 ||
      structure.site_info[donor].hbond_status == DONOR_CHARGED_SP2 ){
    type = 1; //sp2 donor
  }
  else {
    type = 2; //sp3 donor
  }
  if( structure.site_info[acceptor].hbond_status == ACCEPTOR_SP2 ||
      structure.site_info[acceptor].hbond_status == DONOR_AND_ACCEPTOR_SP2 ||
      structure.site_info[acceptor].hbond_status == ACCEPTOR_CHARGED_SP2 ){
    type += 0; //sp2 acceptor
  }
  else {
    type += 2; //sp3 acceptor
  }

  return type;
}

void Froda::scanPhi( SiteID hydrogen, SiteID acceptor, MolFramework &structure, 
                     int hbType, double &phi, bool &isOK ) {
  isOK = true;
  //some types don't check phi?
  if ( hbType == 3 ) {
    phi = PI;
    isOK = true;
    return;
  }

  double myPhi;
  double bestPhi = 0;
  if ( hbType == 4 ) bestPhi = PI;  

  for ( int whichN = 0; whichN < structure.site_info[acceptor].nCovalentNeighbors;
        whichN++ ) {
    SiteID base = structure.site_info[acceptor].neighbor_list[whichN];
    myPhi = computeAngle( hydrogen, acceptor, base );

    if ( hbType < 3 ) {
      if ( RAD_TO_DEG * myPhi <= 90.0 ) {
        isOK = false;
        phi = myPhi;
        return;
      }  

      if ( myPhi > bestPhi ) {
        bestPhi = myPhi;
      }
    }
    //type 4 (sp3-sp3) is always acceptable but need phi for energy
    else if ( hbType == 4 ) { 
      if ( myPhi < bestPhi ) {
        bestPhi = myPhi;
      }
    }

  } 
  //if we got here, all phis are acceptable
  isOK = true;
  phi = bestPhi;
  return;
}

////////////////////////////////////////////////////////////////////////////////
bool Froda::hasMoved( unsigned int atom, MolFramework &structure ){
  if (  abs( currentPos[atom].x - structure.site_info[atom].coords[0]) > epsilon ) return true;
  if (  abs( currentPos[atom].y - structure.site_info[atom].coords[1]) > epsilon ) return true;
  if (  abs( currentPos[atom].z - structure.site_info[atom].coords[2]) > epsilon ) return true;

  return false; //currentPos and coords match each other
}

void Froda:: checkPairs(){
  for ( unsigned int whichPair = 0; whichPair < myPairs.size(); whichPair++ ) {
    myPairs[whichPair].tagGood = false;//indicates not successfully probed

    Vector AB = currentPos[ myPairs[whichPair].atomA ];
    AB.x -= currentPos[ myPairs[whichPair].atomB ].x;
    AB.y -= currentPos[ myPairs[whichPair].atomB ].y;
    AB.z -= currentPos[ myPairs[whichPair].atomB ].z;

    double dist = sqrt(dotProduct ( AB, AB ));
    myPairs[whichPair].ABcurrent = dist;

    if ( dist <  myPairs[whichPair].AB ) {
      myPairs[whichPair].tagGood = true;//indicates successfully probed
    }
  }

  return;
}




