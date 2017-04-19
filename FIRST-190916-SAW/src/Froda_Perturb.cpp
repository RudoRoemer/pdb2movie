#include "global_defs.h"
#include "interact.h"
#include "Froda.h"
#include "mt19937ar.h"
#include "flexweb.h"

#include <signal.h>

extern Parameters parameters;


////////////////////////////////////////////////////////////////////////////////
// Description: Randomly displace the atoms
////////////////////////////////////////////////////////////////////////////////
void Froda::newRandomMove(){ // moves all atoms by a random step
  int atoma; // atom number
  double stepsize;
  double myrad;
  double pmet; // temporary rescale
  Vector turnvec, radvec, movevec;

  if (persist || makeSingleMove) {
    for(  unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
      frodaAtom.at(siteNumber).randomMove = NULL_VEC;
    }
  }

  if (moveByCluster){

    stepsize = bounce;
    
    if ( makeSingleMove ) {
      //must pick one cluster number from the list
      //clusterNumber = 2 to myNClusters inclusive
      //genrand_real2 gives 0 to <1 range
      unsigned int clusterNumber = (unsigned int) ( genrand_real2() * (myNClusters -1 ) );
      //this gives a range of 0 to myNClusters -2
      clusterNumber += 2; //required range

      //cerr << "Throwing cluster " << clusterNumber << endl;

      tempvec = randomUnitVector();
      myrad = stepsize * genrand_real1(); // 0 to small_bounce

      tempvec.x *= myrad;
      tempvec.y *= myrad;
      tempvec.z *= myrad;

      //tempvec is the body displacement

      dummyvec = randomUnitVector();    

      //dummyvec is the unit axis

      if ( ghost[clusterNumber].bodysize == 1 ){ //traps single-atom bodies, was buggy before
        dummyvec = NULL_VEC;
      }
      else {
        pmet = 0.5* stepsize / ghost.at(clusterNumber).radRMS;
        dummyvec.x *= pmet;
        dummyvec.y *= pmet;
        dummyvec.z *= pmet;
      }

      //dummyvec is now a scaled axis for the cross product
      
      for(unsigned int member=0; member<ghost[clusterNumber].atoms.size();member++){

        atoma = ghost[clusterNumber].atoms.at(member);
        //cerr << "Throwing atom " << atoma << endl;
        if ( frodaAtom.at(atoma).isCore ) continue;

        radvec.x = currentPos.at(atoma).x - ghost.at(clusterNumber).posCentral.x;
        radvec.y = currentPos.at(atoma).y - ghost.at(clusterNumber).posCentral.y;
        radvec.z = currentPos.at(atoma).z - ghost.at(clusterNumber).posCentral.z;

        turnvec = crossProduct (radvec, dummyvec);

        movevec.x = tempvec.x + turnvec.x;
        movevec.y = tempvec.y + turnvec.y;
        movevec.z = tempvec.z + turnvec.z;

        currentPos.at(atoma).x += movevec.x;
        currentPos.at(atoma).y += movevec.y;
        currentPos.at(atoma).z += movevec.z;

        if (persist) {
          frodaAtom.at(atoma).randomMove.x += movevec.x;
          frodaAtom.at(atoma).randomMove.y += movevec.y;
          frodaAtom.at(atoma).randomMove.z += movevec.z;
        }

      }

      return; //done with single move
    }

    for( unsigned int clusterNumber = 2; clusterNumber <= myNClusters; clusterNumber++ ){

     //cerr << "Throwing cluster " << clusterNumber << endl;

      if (moveProb < 1.0){
        if ( genrand_real1() > moveProb) continue;
      }

      tempvec = randomUnitVector();
      myrad = stepsize * genrand_real1(); // 0 to small_bounce

      tempvec.x *= myrad;
      tempvec.y *= myrad;
      tempvec.z *= myrad;

      //tempvec is the body displacement

      dummyvec = randomUnitVector();    

      //dummyvec is the unit axis

      if ( ghost[clusterNumber].bodysize == 1 ){ //traps single-atom bodies, was buggy before
        dummyvec = NULL_VEC;
      }
      else {
        pmet = 0.5* stepsize / ghost.at(clusterNumber).radRMS;
        dummyvec.x *= pmet;
        dummyvec.y *= pmet;
        dummyvec.z *= pmet;
      }

      //dummyvec is now a scaled axis for the cross product
      
      for(unsigned int member=0; member<ghost[clusterNumber].atoms.size();member++){

        atoma = ghost[clusterNumber].atoms.at(member);
        //cerr << "Throwing atom " << atoma << endl;
        if ( frodaAtom.at(atoma).isCore ) continue;

        radvec.x = currentPos.at(atoma).x - ghost.at(clusterNumber).posCentral.x;
        radvec.y = currentPos.at(atoma).y - ghost.at(clusterNumber).posCentral.y;
        radvec.z = currentPos.at(atoma).z - ghost.at(clusterNumber).posCentral.z;

        turnvec = crossProduct (radvec, dummyvec);

        movevec.x = tempvec.x + turnvec.x;
        movevec.y = tempvec.y + turnvec.y;
        movevec.z = tempvec.z + turnvec.z;

        currentPos.at(atoma).x += movevec.x;
        currentPos.at(atoma).y += movevec.y;
        currentPos.at(atoma).z += movevec.z;

        if (persist) {
          frodaAtom.at(atoma).randomMove.x += movevec.x;
          frodaAtom.at(atoma).randomMove.y += movevec.y;
          frodaAtom.at(atoma).randomMove.z += movevec.z;
        }

      }
    }
  } else { // atom-move not body-move

    stepsize = bounce;
  
    if ( makeSingleMove ) {
      //must pick one atom number from the list
      //siteNumber = 1 to nTotalSites inclusive
      //genrand_real2 gives 0 to <1 range
      unsigned int siteNumber;

      bool notFound = true;

      while (notFound )  {
        siteNumber = (unsigned int) ( genrand_real2() * (nTotalSites) );
        siteNumber += 1; //required range
        if ( !frodaAtom.at(siteNumber).isCore ) notFound = false;
      } 

      //cerr << "Throwing atom " << siteNumber << endl;

      tempvec = randomUnitVector();
   
      myrad = stepsize*genrand_real1(); // 0 to stepsize
  
      tempvec.x *= myrad;
      tempvec.y *= myrad; 
      tempvec.z *= myrad;
      currentPos.at(siteNumber).x += tempvec.x;
      currentPos.at(siteNumber).y += tempvec.y;
      currentPos.at(siteNumber).z += tempvec.z;
      if (persist) {
        frodaAtom.at(siteNumber).randomMove.x += tempvec.x;
        frodaAtom.at(siteNumber).randomMove.y += tempvec.y;
        frodaAtom.at(siteNumber).randomMove.z += tempvec.z;
      }
      return; //done with single move
    }

    // Atoms in the core should not be moved
  
    for(  unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
      if ( frodaAtom.at(siteNumber).isCore ) continue;

        if (moveProb < 1.0){
          if ( genrand_real1() > moveProb) continue;
        }
       
        tempvec = randomUnitVector();
   
        myrad = stepsize*genrand_real1(); // 0 to stepsize
  
        tempvec.x *= myrad;
        tempvec.y *= myrad; 
        tempvec.z *= myrad;
        currentPos.at(siteNumber).x += tempvec.x;
        currentPos.at(siteNumber).y += tempvec.y;
        currentPos.at(siteNumber).z += tempvec.z;
        if (persist) {
          frodaAtom.at(siteNumber).randomMove.x += tempvec.x;
          frodaAtom.at(siteNumber).randomMove.y += tempvec.y;
          frodaAtom.at(siteNumber).randomMove.z += tempvec.z;
        }
      
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void Froda::oldRandomMove(){ // moves all atoms by a random step

  for(  unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
    if ( frodaAtom.at(siteNumber).isCore ) continue;
    currentPos.at(siteNumber).x += frodaAtom.at(siteNumber).randomMove.x;
    currentPos.at(siteNumber).y += frodaAtom.at(siteNumber).randomMove.y;
    currentPos.at(siteNumber).z += frodaAtom.at(siteNumber).randomMove.z;
  }
}


////////////////////////////////////////////////////////////////////////////////
// Description: Biases an atom to move away from the center of the rigid core
//              or the central point of the system if mobRC1 is true.
////////////////////////////////////////////////////////////////////////////////
void Froda::centrifuge(int an_atom){ 

  //cerr << "In centrifuge for atom " << an_atom << endl;

  Vector central_point;

  unsigned int whichLC = 0;

  if ( localCentrifuge ) {

    //cerr << "In local centrifuge." << endl;

    double dist2  = 1E10;

    for ( unsigned int LC = 0; LC < centerID.size(); LC++ ) {
      unsigned int thisLC = centerID[LC];
      dummyvec.x = currentPos.at(an_atom).x - currentPos.at(thisLC).x;
      dummyvec.y = currentPos.at(an_atom).y - currentPos.at(thisLC).y;
      dummyvec.z = currentPos.at(an_atom).z - currentPos.at(thisLC).z;

      dot = dotProduct( dummyvec, dummyvec );
      if ( dot < dist2 ) {
        whichLC = LC;
        dist2 = dot;
      }
    }
    //cerr << "Using local center " << whichLC << endl;

    central_point = currentPos.at( centerID[whichLC]);
  }
  else {

    //cerr << "Using global centrifuge." << endl;

    if (mobileRC1){ central_point = myCentralPoint;}
    else {central_point = ghost.at(1).posCentral;}
  }

  //cerr << "Central " << central_point.x << " ";
  //cerr << central_point.y << " ";
  //cerr << central_point.z << endl;

  //double small_bounce =bounce*bias_fraction;
  dummyvec.x = currentPos.at(an_atom).x - central_point.x;
  dummyvec.y = currentPos.at(an_atom).y - central_point.y;
  dummyvec.z = currentPos.at(an_atom).z - central_point.z;
  
  dot = dotProduct(dummyvec, dummyvec);
  veclength = sqrt(dot);
  
  if (veclength < epsilon){
    return; // bang on center, don't move it. Can this happen
  }
  else if ( veclength > centri_maxr ) {
    return; //too far, don't move it
  }

  dummyvec.x /= veclength;
  dummyvec.y /= veclength;
  dummyvec.z /= veclength;
  // dummyvec now a unit vector radial from the rigid core center
  // scale to bounce size
  
  dummyvec.x *= targetStep; // smaller then main bounce effect
  dummyvec.y *= targetStep; // smaller then main bounce effect
  dummyvec.z *= targetStep; // smaller then main bounce effect
 
  if ( localCentrifuge &&  ( centerDirection[whichLC] < 0 ) ) {
    dummyvec.x *= -1;
    dummyvec.y *= -1;
    dummyvec.z *= -1;
  }
 
  currentPos.at(an_atom).x += dummyvec.x;
  currentPos.at(an_atom).y += dummyvec.y;
  currentPos.at(an_atom).z += dummyvec.z;
}

////////////////////////////////////////////////////////////////////////////////
// Description: Biases an atom to move away from the z axis
//              Useful for membrane simulations
////////////////////////////////////////////////////////////////////////////////
void Froda::cylinderCentrifuge(int an_atom){ //moves things away from the z axis only

  Vector central_point = NULL_VEC;

  dummyvec.x = currentPos.at(an_atom).x - central_point.x;
  dummyvec.y = currentPos.at(an_atom).y - central_point.y;
  dummyvec.z = 0.0;
  
  dot = dotProduct(dummyvec, dummyvec);
  veclength = sqrt(dot);
  
  if (veclength < epsilon) return; // bang on center, don't move it. Can this happen? 
  dummyvec.x /= veclength;
  dummyvec.y /= veclength;
  dummyvec.z /= veclength;
  // dummyvec now a unit vector radial from the z axis
  // scale to bounce size
  
  dummyvec.x *= targetStep; // smaller then main bounce effect
  dummyvec.y *= targetStep; // smaller then main bounce effect
  
  currentPos.at(an_atom).x += dummyvec.x;
  currentPos.at(an_atom).y += dummyvec.y;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:  do a targeted move towards a particular structure (aim).
////////////////////////////////////////////////////////////////////////////////
void Froda::moveToTarget( MolFramework &aim){

 double tol, pmet;
 double largest_d = 0.0;
 double ild = 0.0; // equals inverse largest distance

 tol = tolTarget;

 // move all the atoms towards their targets 

 // if it's proportional-targetted, first find the maximum distance from the target
 //everyone else gets a smaller target step in proportion

 if (propTotarget){
   for (unsigned int siteNumber=1; siteNumber<=nTotalSites;siteNumber++){
     if (frodaAtom.at(siteNumber).isCore ) continue;
     // targeted- check myTarget
     unsigned int targatom = frodaAtom.at(siteNumber).myTarget;
     if (targatom == 0 ) continue;
  
     dummyvec.x = currentPos.at(siteNumber).x - aim.site_info[targatom].coords[0];
     dummyvec.y = currentPos.at(siteNumber).y - aim.site_info[targatom].coords[1];
     dummyvec.z = currentPos.at(siteNumber).z - aim.site_info[targatom].coords[2];
  
     dot = dotProduct(dummyvec, dummyvec);
     if ( dot > largest_d ) largest_d = dot; 
   }
 
   largest_d = sqrt(largest_d);
   ild = 1.0/largest_d;
 }
 

 for (unsigned int siteNumber=1; siteNumber<=nTotalSites;siteNumber++){

   if (frodaAtom.at(siteNumber).isCore ) continue;

   // targeted- check myTarget

   unsigned int targatom = frodaAtom.at(siteNumber).myTarget;
   if (targatom == 0 ) continue;

   dummyvec.x = currentPos.at(siteNumber).x - aim.site_info[targatom].coords[0];
   dummyvec.y = currentPos.at(siteNumber).y - aim.site_info[targatom].coords[1];
   dummyvec.z = currentPos.at(siteNumber).z - aim.site_info[targatom].coords[2];

   dot = sqrt(dotProduct(dummyvec, dummyvec));

   if (dot > abs(targetStep)){  // allow for negative steps!

     if (propTotarget){ pmet = targetStep * ild;} //step vector will have magnitude proportional to distance
     else{
       pmet = targetStep/dot; //step vector will have magnitude targetStep
     }
     tempvec.x = dummyvec.x * pmet;
     tempvec.y = dummyvec.y * pmet;
     tempvec.z = dummyvec.z * pmet;
     // tempvec is step vector from target to current
  
   } else{ //very close now, just apply remaining vector

     if (targetStep > 0.0){
       tempvec.x = dummyvec.x;
       tempvec.y = dummyvec.y;
       tempvec.z = dummyvec.z;
     } else{
       tempvec.x = -dummyvec.x;
       tempvec.y = -dummyvec.y;
       tempvec.z = -dummyvec.z;
       // negative step sizes give unfolding 
     }

   }

   currentPos.at(siteNumber).x -= tempvec.x;
   currentPos.at(siteNumber).y -= tempvec.y;
   currentPos.at(siteNumber).z -= tempvec.z;

 }

 return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: Go-type targeting towards target local contacts
////////////////////////////////////////////////////////////////////////////////
void Froda::moveToLocalTarget(){

  Vector dr;
  double dist, idist;

  for (unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++){
    if ( frodaAtom.at(siteNumber).myTarget == 0) continue;
    unsigned int nfriends = frodaAtom.at(siteNumber).friends.size();
    if (nfriends ==0) continue;

    double infriends = 1.0/nfriends;

    Vector centralpoint = NULL_VEC;

    for (unsigned int friendNumber = 0; friendNumber < nfriends; friendNumber++){
      int anatom = frodaAtom.at(siteNumber).friends.at(friendNumber);
      centralpoint.x += currentPos.at(anatom).x;
      centralpoint.y += currentPos.at(anatom).y;
      centralpoint.z += currentPos.at(anatom).z;
    }
    centralpoint.x *= infriends;
    centralpoint.y *= infriends;
    centralpoint.z *= infriends;

    dr.x = centralpoint.x - currentPos.at(siteNumber).x;
    dr.y = centralpoint.y - currentPos.at(siteNumber).y;
    dr.z = centralpoint.z - currentPos.at(siteNumber).z;

    dist = (dotProduct(dr,dr));

    if ( dist > frodaAtom.at(siteNumber).offset ){ //do some local targetting

      //cerr << "Atom " << siteNumber << " outside friend sphere, targetting..." << endl;

      idist = lstep/sqrt(dist);
      dr.x *= idist;
      dr.y *= idist;
      dr.z *= idist; //dr now siteNumber step vector from atom to centralpoint

      currentPos.at(siteNumber).x += dr.x;
      currentPos.at(siteNumber).y += dr.y;
      currentPos.at(siteNumber).z += dr.z;
    }

  }

  return;
}



////////////////////////////////////////////////////////////////////////////////
// Description: apply centrifuge (spherical or cylindrical) to the system
////////////////////////////////////////////////////////////////////////////////
void Froda::applyCentrifuge(){

  if (cylindrical){
    for (unsigned int siteNumber=1; siteNumber <= nTotalSites; siteNumber++){
      if( currentPos.at(siteNumber).z < lowZRange) continue;
      if( currentPos.at(siteNumber).z > highZRange) continue;
      cylinderCentrifuge(siteNumber);
    }
  } else{
    for (unsigned int siteNumber=1; siteNumber <= nTotalSites; siteNumber++){
      centrifuge(siteNumber);
    }
  }

  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: push a motion further (from Holger's conjmove idea)
//              by applying momentum to the ghost's translate and rotate
////////////////////////////////////////////////////////////////////////////////
void Froda::applyMomentum( int ghostid ){
  Vector dcpos;
  Vector drotor;

  Vector delta;

  dcpos.x = ghost.at(ghostid).posCentral.x - ghost.at(ghostid).olderPosCentral.x;
  dcpos.y = ghost.at(ghostid).posCentral.y - ghost.at(ghostid).olderPosCentral.y;
  dcpos.z = ghost.at(ghostid).posCentral.z - ghost.at(ghostid).olderPosCentral.z;

  drotor.x = ghost.at(ghostid).runningRotor.x + ghost.at(ghostid).baseRotor.x - ghost.at(ghostid).olderRotor.x;
  drotor.y = ghost.at(ghostid).runningRotor.y + ghost.at(ghostid).baseRotor.y - ghost.at(ghostid).olderRotor.y;
  drotor.z = ghost.at(ghostid).runningRotor.z + ghost.at(ghostid).baseRotor.z - ghost.at(ghostid).olderRotor.z;

  dcpos.x *= cmforwardfactor/(2*cacheEvery);
  dcpos.y *= cmforwardfactor/(2*cacheEvery);
  dcpos.z *= cmforwardfactor/(2*cacheEvery);
  drotor.x *= cmforwardfactor/(2*cacheEvery);
  drotor.y *= cmforwardfactor/(2*cacheEvery);
  drotor.z *= cmforwardfactor/(2*cacheEvery);

  //cerr << "Body " << ghostid << " has momentum: ";
  //cerr << dcpos.x << " " << dcpos.y << " " << dcpos.z;
  //cerr << " and rotation: ";
  //cerr << drotor.x << " " << drotor.y << " " << drotor.z << endl;

  for ( unsigned int ghostAtom = 0; ghostAtom < ghost.at(ghostid).bodysize; ghostAtom++){
    int atoma = ghost.at(ghostid).atoms.at(ghostAtom);
    Vector abond;
    abond.x = ghost.at(ghostid).pos.at(ghostAtom).x - ghost.at(ghostid).posCentral.x;
    abond.y = ghost.at(ghostid).pos.at(ghostAtom).y - ghost.at(ghostid).posCentral.y;
    abond.z = ghost.at(ghostid).pos.at(ghostAtom).z - ghost.at(ghostid).posCentral.z;
    
    delta = crossProduct ( drotor , abond );
    //cerr << "Atom " << atoma << " Turn " << delta.x << " " << delta.y << " " << delta.z << endl;
    delta.x += dcpos.x; 
    delta.y += dcpos.y; 
    delta.z += dcpos.z; 
    //c//err << "Atom " << atoma << " Move " << delta.x << " " << delta.y << " " << delta.z << endl;

    double dot = sqrt(dotProduct(delta, delta));

    if ( dot > cmmaxstep){
      //cerr << "Dot " << dot << " cmmaxstep " << cmmaxstep << " c/d " << cmmaxstep/dot << endl;
      delta.x *= cmmaxstep/dot;
      delta.y *= cmmaxstep/dot;
      delta.z *= cmmaxstep/dot;
    }
    //cerr << "Atom " << atoma << " Move " << delta.x << " " << delta.y << " " << delta.z << endl;

    //double momentumSize = sqrt(dotProduct(delta, delta));
    //cerr  << "Atom " << atoma << " Momentum step size: " << momentumSize << endl;

    currentPos.at(atoma).x += delta.x;
    currentPos.at(atoma).y += delta.y;
    currentPos.at(atoma).z += delta.z;
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: push a motion backwards (experimental)
////////////////////////////////////////////////////////////////////////////////
void Froda::reverseMomentum( int ghostid ){
  Vector dcpos;
  Vector drotor;

  Vector delta;

  dcpos.x = ghost.at(ghostid).posCentral.x - ghost.at(ghostid).olderPosCentral.x;
  dcpos.y = ghost.at(ghostid).posCentral.y - ghost.at(ghostid).olderPosCentral.y;
  dcpos.z = ghost.at(ghostid).posCentral.z - ghost.at(ghostid).olderPosCentral.z;

  drotor.x = ghost.at(ghostid).runningRotor.x + ghost.at(ghostid).baseRotor.x - ghost.at(ghostid).olderRotor.x;
  drotor.y = ghost.at(ghostid).runningRotor.y + ghost.at(ghostid).baseRotor.y - ghost.at(ghostid).olderRotor.y;
  drotor.z = ghost.at(ghostid).runningRotor.z + ghost.at(ghostid).baseRotor.z - ghost.at(ghostid).olderRotor.z;

  dcpos.x *= cmreversefactor/(2*cacheEvery);
  dcpos.y *= cmreversefactor/(2*cacheEvery);
  dcpos.z *= cmreversefactor/(2*cacheEvery);
  drotor.x *= cmreversefactor/(2*cacheEvery);
  drotor.y *= cmreversefactor/(2*cacheEvery);
  drotor.z *= cmreversefactor/(2*cacheEvery);

  for ( unsigned int ghostAtom = 0; ghostAtom < ghost.at(ghostid).bodysize; ghostAtom++){
    int atoma = ghost.at(ghostid).atoms.at(ghostAtom);
    
    delta = crossProduct ( drotor , ghost.at(ghostid).pos.at(ghostAtom) );
    delta.x += dcpos.x; 
    delta.y += dcpos.y; 
    delta.z += dcpos.z; 

    double dot = sqrt(dotProduct(delta, delta));

    if ( dot > cmmaxstep){
      delta.x *= cmmaxstep/dot;
      delta.y *= cmmaxstep/dot;
      delta.z *= cmmaxstep/dot;
    }

    currentPos.at(atoma).x -= delta.x;
    currentPos.at(atoma).y -= delta.y;
    currentPos.at(atoma).z -= delta.z;
  }
  return;

}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   apply the user-defined pair associations. Associated atoms move together. 
//   NEW August 21 2006 SAW: allowed for repulsive interactions
//   Defined by entry atom1 atom2 -R in pair.in file
//   atom1 and atom2 are biased to be at least R apart 
////////////////////////////////////////////////////////////////////////////////
void Froda::applyPairs() {
  int atom1;
  int atom2;
  double AB;
  Vector diff;
  double R;
  bool tagA;
  bool tagB;

  for ( int whichPair = 0; whichPair < nMyPairs; whichPair++ ) {
    //cerr << "Working on pair index " << whichPair << endl;
    atom1 = myPairs.at(whichPair).atomA;
    atom2 = myPairs.at(whichPair).atomB;
    AB = myPairs.at(whichPair).AB;
    tagA = myPairs.at(whichPair).tagA;
    tagB = myPairs.at(whichPair).tagB;

    diff = currentPos.at(atom2);
    diff.x -= currentPos.at(atom1).x;
    diff.y -= currentPos.at(atom1).y;
    diff.z -= currentPos.at(atom1).z;
    //diff is now a vector from A to B
    R = sqrt(dotProduct( diff, diff ) );

    if ( AB > 0 ) {
      if ( R > AB ) {
        //cerr << atom1 << " " << atom2 << " " << R << "/" << AB << endl;

        //apply motion to the atoms
        //make diff a unit vector
        diff.x /= R;
        diff.y /= R;
        diff.z /= R;
  
        //make diff a step vector
        diff.x *= targetStep;
        diff.y *= targetStep;
        diff.z *= targetStep;
  
        //move atom 1
        if ( !frodaAtom.at(atom1).isCore & tagA ) {
          currentPos.at(atom1).x += diff.x;
          currentPos.at(atom1).y += diff.y;
          currentPos.at(atom1).z += diff.z;
        }
  
        //move atom 2
        if ( !frodaAtom.at(atom2).isCore & tagB ) {
          currentPos.at(atom2).x -= diff.x;
          currentPos.at(atom2).y -= diff.y;
          currentPos.at(atom2).z -= diff.z;
        }

      }
    }
    else {
      AB *= -1.0;
      if ( R < AB ) {
        //cerr << atom1 << " " << atom2 << " " << R << "\\" << AB << endl;

        //apply motion to the atoms
        //make diff a unit vector
        diff.x /= R;
        diff.y /= R;
        diff.z /= R;
  
        //make diff a step vector
        diff.x *= targetStep;
        diff.y *= targetStep;
        diff.z *= targetStep;
  
        //move atom 1
        if ( !frodaAtom.at(atom1).isCore & tagA ) {
          currentPos.at(atom1).x -= diff.x;
          currentPos.at(atom1).y -= diff.y;
          currentPos.at(atom1).z -= diff.z;
        }
  
        //move atom 2
        if ( !frodaAtom.at(atom2).isCore & tagB ) {
          currentPos.at(atom2).x += diff.x;
          currentPos.at(atom2).y += diff.y;
          currentPos.at(atom2).z += diff.z;
        }

      }
    }

  }

  return;
}



////////////////////////////////////////////////////////////////////////////////
// Description:
//   builds an internal coordinate basis on a CA atom
//   used for adaptive eigenvetors.
////////////////////////////////////////////////////////////////////////////////
void Froda::makeResidueBasis() {
  Vector vec1, vec2, sumVector;
  Vector basisVector;
  double dot;

  for ( unsigned int carbon =0; carbon < nodeResidues.size(); carbon++ ) {
    unsigned int atomA, atomB, atomC;
    atomA = nodeResidues.at(carbon).cAlpha;
    atomB = nodeResidues.at(carbon).n1;
    atomC = nodeResidues.at(carbon).n2;

    vec1 = currentPos.at(atomB);
    vec1.x -= currentPos.at(atomA).x;
    vec1.y -= currentPos.at(atomA).y;
    vec1.z -= currentPos.at(atomA).z;

    vec2 = currentPos.at(atomC);
    vec2.x -= currentPos.at(atomA).x;
    vec2.y -= currentPos.at(atomA).y;
    vec2.z -= currentPos.at(atomA).z;

    sumVector.x = vec1.x - vec2.x;
    sumVector.y = vec1.y - vec2.y;
    sumVector.z = vec1.z - vec2.z;

    //make a basis out of the C-CA-N plane;
    basisVector = crossProduct( vec1, vec2 );
    dot = sqrt( dotProduct( basisVector, basisVector ) );
   
    basisVector.x /= dot;
    basisVector.y /= dot;
    basisVector.z /= dot;

    nodeResidues.at(carbon).e1 = basisVector;
    //first basis vector

    //now make the second unit basis
    basisVector = crossProduct( nodeResidues.at(carbon).e1, sumVector );
    dot = sqrt( dotProduct( basisVector, basisVector ) );
   
    basisVector.x /= dot;
    basisVector.y /= dot;
    basisVector.z /= dot;

    nodeResidues.at(carbon).e3 = basisVector;

    nodeResidues.at(carbon).e2 = crossProduct( nodeResidues.at(carbon).e3, nodeResidues.at(carbon).e1 );
    //completes formation of unit basis e123

  }
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   applies an elastic network model eigenvector
//   the biasing is applied per-residue 
////////////////////////////////////////////////////////////////////////////////
void Froda::addModeBias() {
   for ( unsigned int node =0; node < nodeResidues.size(); node++ ) {
     if ( !nodeResidues.at(node).isBiased ) continue;
     //what's the mode bias?
     Vector bias;
     if ( useInternalBias ) {
       //you must have run makeResidueBasis for this to work!
       nodeResidues.at(node).makeExternalBias();
       bias = nodeResidues.at( node ).externalBias;
     }
     else {
       bias = nodeResidues.at( node ).originalBias;
     }

     //shift the alpha carbon and all the others too
     unsigned int catom = nodeResidues.at(node).cAlpha;
     currentPos.at(catom).x += bias.x;
     currentPos.at(catom).y += bias.y;
     currentPos.at(catom).z += bias.z;

     //cerr << "Biasing alpha-carbon " << catom << " by vector " << bias.x << " " << bias.y << " " << bias.z << endl;

     //the other atoms belonging to a residue are
     //sameResidue.at( nodeResidue.at(carbon).cAlpha )

     for ( unsigned int k = 0; k < nodeResidues.at(node).sameResidue.size(); k++ ) {
       unsigned int theAtom = nodeResidues.at(node).sameResidue.at(k);
       currentPos.at(theAtom).x += bias.x;
       currentPos.at(theAtom).y += bias.y;
       currentPos.at(theAtom).z += bias.z;
     }
   }
}


void Froda::pickPeptideDihedral( SiteID &CA, SiteID &base, MolFramework &structure ) {

  bool stillLooking = true;
  SiteID tryCA = 0;
  SiteID tryBase = 0;

  while( stillLooking ) {

    unsigned int whichTry = (unsigned int ) ( genrand_real2() * alphaCarbon.size() ) ;
    tryCA = alphaCarbon[whichTry];

    //cerr << "Trying CA atom " << structure.site_info[tryCA].orig_atom_number << endl;

    bool seekingBase = true;
    while( seekingBase ) {
      unsigned int whichOther = (unsigned int) (genrand_real2() * frodaAtom.at(tryCA).neighbors.size() );
      
      if ( frodaAtom.at( frodaAtom.at(tryCA).neighbors[whichOther] ).isMain ) {
        seekingBase = false;
        tryBase = frodaAtom.at(tryCA).neighbors[whichOther];
      }
    }

    if ( structure.site_info[tryCA].rigid_label != structure.site_info[tryBase].rigid_label ) {
      stillLooking = false;
    } 

  }

  CA = tryCA;
  base = tryBase;

  return;
}

void Froda::perturbBackboneDihedral( SiteID CA, SiteID base, MolFramework &structure ){
  double stepsize;
  double largestCross;
  double scaler;

  Vector axis;
  Vector arm;
  Vector crossvec;
  Vector movevec;

  for(  unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
    frodaAtom.at(siteNumber).randomMove = NULL_VEC;
    frodaAtom.at(siteNumber).moved = false;
  }

  vector< SiteID > listToMove;
  vector< SiteID > nextListToMove;

  //i/ no step can exceed stepsize
  //ii/ atoms near the axis get smaller moves
  //the magnitude of arm cross axis measures the radius of the rotation

  //first, get the axis for the rotation
  axis = currentPos[base];
  axis.x -= currentPos[CA].x;
  axis.y -= currentPos[CA].y;
  axis.z -= currentPos[CA].z;
  double mag = sqrt( dotProduct( axis, axis ) );
  axis.x /= mag;
  axis.y /= mag;
  axis.z /= mag;

  //second, get the step size
  stepsize = bounce;
  stepsize *= ( 2* genrand_real1() -1 );
  //now positive or negative
  //interpretation; no atom moves further than stepsize

  //third, set up an initial set of atoms to move
  listToMove.clear();

  //cerr << "Working from CA atom " << structure.site_info[CA].orig_atom_number;
  //cerr << " and base atom " << structure.site_info[base].orig_atom_number << endl;

  //check to avoid perturbing a core
  if ( structure.site_info[base].rigid_label == 1 &&
       frodaAtom.at(base).isCore ) {
    SiteID temp = base; base = CA; CA= temp;
  }

  for ( int whichN =0; whichN < structure.site_info[base].nCovalentNeighbors;
        whichN++ ) {
    SiteID thisN = structure.site_info[base].neighbor_list[whichN];
    if (thisN == CA) continue;
    if ( frodaAtom.at(thisN).isCore ) continue; 
    listToMove.push_back(thisN);
  } 

  //enter loop of perturbing atoms

  largestCross = 0;

  bool stillWorking = true;
  while(stillWorking){
    //perturb the atoms on the list
    for ( unsigned int whichA = 0; whichA < listToMove.size(); whichA++ ) {
      SiteID thisAtom = listToMove[whichA];
      if ( frodaAtom.at(thisAtom).isCore ) continue;    
   
      arm = currentPos[thisAtom];
      arm.x -= currentPos[base].x;
      arm.y -= currentPos[base].y;
      arm.z -= currentPos[base].z;

      crossvec = crossProduct(arm, axis);
      double mag2 = dotProduct( crossvec, crossvec ) ;

      if ( mag2 > largestCross  ){
        largestCross = mag2 ;
      } 

      frodaAtom.at(thisAtom).randomMove = crossvec;
      frodaAtom.at(thisAtom).moved = true;
    }

    //fill up a new list
    nextListToMove.clear();
    for ( unsigned int whichA = 0; whichA < listToMove.size(); whichA++ ) {
      SiteID thisAtom = listToMove[whichA];
      for ( int whichN = 0; whichN < structure.site_info[thisAtom].nCovalentNeighbors;
            whichN++ ) {
        SiteID thisN = structure.site_info[thisAtom].neighbor_list.at(whichN);
        if ( frodaAtom.at(thisN).moved ) continue;
        if ( frodaAtom.at(thisN).isCore ) continue;    
        if ( thisN == base ) continue;
        if ( thisN == CA ) continue;
        nextListToMove.push_back(thisN);
      }

    }
    if (nextListToMove.size() == 0 ) {
     stillWorking = false; //end loop
    }
    else{
      listToMove.clear();
      listToMove = nextListToMove;
    }

  }

  //finally, normalise and apply the rotation

  largestCross = sqrt(largestCross);
  scaler = stepsize / largestCross;

  for ( SiteID whichA = 1; whichA < frodaAtom.size(); whichA++ ) {
    if ( frodaAtom[whichA].moved ) {
      frodaAtom[whichA].randomMove.x *= scaler;
      frodaAtom[whichA].randomMove.y *= scaler;
      frodaAtom[whichA].randomMove.z *= scaler;

      currentPos[whichA].x += frodaAtom[whichA].randomMove.x;
      currentPos[whichA].y += frodaAtom[whichA].randomMove.y;
      currentPos[whichA].z += frodaAtom[whichA].randomMove.z;
    }

  }

  return;
}


void Froda::pickPeptidePlane( SiteID &CA1, SiteID &b1, SiteID &b2, SiteID &CA2, MolFramework &structure ) {
  //first, use the pickDihedral function to get a suitable CA1 and b1;
  SiteID tryCA1, tryb1, tryb2, tryCA2;
  bool looking = true;

  while (looking) {
    tryCA1 = tryb1 = tryb2 = tryCA2 = 0;
    pickPeptideDihedral( tryCA1, tryb1, structure);

    if ( frodaAtom.at(tryb1).isCore && structure.site_info[tryb1].rigid_label == 1 ) continue;
    bool gotb2 = false;
    for ( int whichN = 0; whichN < structure.site_info[tryb1].nCovalentNeighbors;
          whichN++ ) {
      SiteID thisN = structure.site_info[tryb1].neighbor_list[whichN];
      if ( thisN == tryCA1 ) continue;
      if ( structure.site_info[thisN].nCovalentNeighbors > 1 ) {
        tryb2 = thisN;
        gotb2 = true;
      }
    }
    if ( !gotb2 ) continue; //try again

    bool gotCA2 = false;
    for ( int whichN = 0; whichN < structure.site_info[tryb2].nCovalentNeighbors;
          whichN++ ) {
      SiteID thisN = structure.site_info[tryb2].neighbor_list[whichN];
      if ( thisN == tryb1 ) continue;
      if ( structure.site_info[thisN].nCovalentNeighbors > 1 ) {
        if ( structure.isMainchain( thisN ) == 2 ) {
          tryCA2 = thisN;
          gotCA2 = true;
        } 
      }

    }

    if (!gotCA2 ) continue;

    //if ( structure.site_info[tryCA2].rigid_label != 
    //     structure.site_info[tryb2].rigid_label ) {
    looking = false;
    //}
           
  }

  CA1 = tryCA1;
  b1 = tryb1;
  b2 = tryb2;
  CA2 = tryCA2;

  return;
}


void Froda::perturbPeptidePlane( SiteID CA1, SiteID b1, SiteID b2, SiteID CA2, MolFramework &structure ) {

  double stepsize;
  Vector axis;
  Vector arm;
  Vector crossvec;
  Vector movevec;

  for(  unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
    frodaAtom.at(siteNumber).randomMove = NULL_VEC;
    frodaAtom.at(siteNumber).moved = false;
  }

  vector< SiteID > listToMove;

  //first, get the axis for the rotation
  axis = currentPos[CA1];
  axis.x -= currentPos[CA2].x;
  axis.y -= currentPos[CA2].y;
  axis.z -= currentPos[CA2].z;

  double mag = sqrt( dotProduct( axis, axis ) );
  axis.x /= mag;
  axis.y /= mag;
  axis.z /= mag;

  //second, get the step size
  stepsize = bounce;
  stepsize *= ( 2* genrand_real1() -1 );
  //now positive or negative
  //interpretation; no atom moves further than stepsize

  //third, set up an initial set of atoms to move
  listToMove.clear();

  listToMove.push_back( b1 );
  listToMove.push_back( b2 );

  for ( int whichN =0; whichN < structure.site_info[b1].nCovalentNeighbors;
        whichN++ ) {
    SiteID thisN = structure.site_info[b1].neighbor_list[whichN];
    if ( structure.site_info[thisN].nCovalentNeighbors > 1 ) continue;
    listToMove.push_back(thisN);
  } 
  for ( int whichN =0; whichN < structure.site_info[b2].nCovalentNeighbors;
        whichN++ ) {
    SiteID thisN = structure.site_info[b2].neighbor_list[whichN];
    if ( structure.site_info[thisN].nCovalentNeighbors > 1 ) continue;
    listToMove.push_back(thisN);
  } 

  //finally, perturb these atoms

  for ( unsigned int whichA = 0; whichA < listToMove.size(); whichA++ ) {
    SiteID thisAtom = listToMove[whichA];
    if ( frodaAtom.at(thisAtom).isCore ) continue;
       
    arm = currentPos[thisAtom];
    arm.x -= currentPos[CA1].x;
    arm.y -= currentPos[CA1].y;
    arm.z -= currentPos[CA1].z;

    crossvec = crossProduct(arm, axis);

    double mag2 = dotProduct( crossvec, crossvec );
    if ( mag2 > 1 ) {
      mag2 = sqrt(mag2);
      crossvec.x /= mag2;
      crossvec.y /= mag2;
      crossvec.z /= mag2;
    }
    movevec.x = crossvec.x * stepsize;
    movevec.y = crossvec.y * stepsize;
    movevec.z = crossvec.z * stepsize;

    currentPos.at(thisAtom).x += movevec.x;
    currentPos.at(thisAtom).y += movevec.y;
    currentPos.at(thisAtom).z += movevec.z;

    frodaAtom.at(thisAtom).randomMove = movevec;
    frodaAtom.at(thisAtom).moved = true;
  }

  return;
}

void Froda::twitchSidechain( MolFramework &structure ) {
  SiteID pick1, pick2;

  double largestCross;
  double scaler;

  SiteID base, other;
  Vector axis;
  Vector arm;
  Vector crossvec;
  Vector movevec;
  double stepsize;

  vector< SiteID > listToMove, nextListToMove;

  for(  unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
    frodaAtom.at(siteNumber).randomMove = NULL_VEC;
    frodaAtom.at(siteNumber).moved = false;
  }
  
  bool looking = true;
  pick2 = 0;

  while ( looking ){
    pick1 = (unsigned int ) ( genrand_real2() * frodaAtom.size() ) ;
    if ( pick1 == 0 ) continue;   

    unsigned int whichN = (unsigned int) ( genrand_real2() *
     structure.site_info[pick1].nCovalentNeighbors );
    pick2 = structure.site_info[pick1].neighbor_list[whichN]; 

    if ( structure.site_info[pick1].rigid_label ==
         structure.site_info[pick2].rigid_label ) continue;
 
    if ( frodaAtom[pick1].isMain && frodaAtom[pick2].isMain ) continue;
    if ( frodaAtom[pick1].isMain && frodaAtom[pick2].isCore ) continue;
    if ( frodaAtom[pick1].isCore && frodaAtom[pick2].isMain ) continue;

    //if we got here, there's at least one non-main atom, so we're good
    looking = false; 
  }

  //now go by cases to build the set of atoms to move

  if ( frodaAtom[pick1].isMain ) {
    other = pick1;
    base = pick2;
  }
  else {
    other = pick2;
    base = pick1;
  }

  //now twitch everything aroung the axis
  //remember to:
  //set the moved flag of moved atoms and skip them
  //stop at the Core
  //stop at Main
  //stop at base and other

  //first, get the axis for the rotation
  axis = currentPos[base];
  axis.x -= currentPos[other].x;
  axis.y -= currentPos[other].y;
  axis.z -= currentPos[other].z;

  double mag = sqrt( dotProduct( axis, axis ) );
  axis.x /= mag;
  axis.y /= mag;
  axis.z /= mag;

  //second, get the step size
  stepsize = bounce;
  stepsize *= ( 2* genrand_real1() -1 );
  //now positive or negative
  //interpretation; no atom moves further than stepsize

  //third, set up an initial set of atoms to move
  listToMove.clear();

  for ( int whichN =0; whichN < structure.site_info[base].nCovalentNeighbors;
        whichN++ ) {
    SiteID thisN = structure.site_info[base].neighbor_list[whichN];
    if (thisN == other) continue;
    if ( frodaAtom.at(thisN).isCore ) continue; 
    listToMove.push_back(thisN);
  } 

  largestCross = 0;

  //enter loop of perturbing atoms
  bool stillWorking = true;
  while(stillWorking){
    //perturb the atoms on the list
    for ( unsigned int whichA = 0; whichA < listToMove.size(); whichA++ ) {
      SiteID thisAtom = listToMove[whichA];
      if ( frodaAtom.at(thisAtom).isCore ) continue;    
      if ( frodaAtom.at(thisAtom).isMain ) continue;    
      if ( frodaAtom.at(thisAtom).moved ) continue;    
   
      arm = currentPos[thisAtom];
      arm.x -= currentPos[base].x;
      arm.y -= currentPos[base].y;
      arm.z -= currentPos[base].z;

      crossvec = crossProduct(arm, axis);

      double mag2 = dotProduct( crossvec, crossvec );

      if ( mag2 > largestCross ) {
        largestCross = mag2;
      }

      frodaAtom.at(thisAtom).randomMove = crossvec;
      frodaAtom.at(thisAtom).moved = true;
    }

    //fill up a new list
    nextListToMove.clear();
    for ( unsigned int whichA = 0; whichA < listToMove.size(); whichA++ ) {
      SiteID thisAtom = listToMove[whichA];
      for ( int whichN = 0; whichN < structure.site_info[thisAtom].nCovalentNeighbors;
            whichN++ ) {
        SiteID thisN = structure.site_info[thisAtom].neighbor_list.at(whichN);
        if ( frodaAtom.at(thisN).moved ) continue;
        if ( frodaAtom.at(thisN).isCore ) continue;    
        if ( frodaAtom.at(thisN).isMain ) continue;    
        if ( thisN == base ) continue;
        if ( thisN == other ) continue;
        nextListToMove.push_back(thisN);
      }

    }
    if (nextListToMove.size() == 0 ) {
     stillWorking = false; //end loop
    }
    else{
      listToMove.clear();
      listToMove = nextListToMove;
    }

  }

  //finally, normalise and apply the rotation

  largestCross = sqrt(largestCross);
  scaler = stepsize / largestCross;

  for ( SiteID whichA = 1; whichA < frodaAtom.size(); whichA++ ) {
    if ( frodaAtom[whichA].moved ) {
      frodaAtom[whichA].randomMove.x *= scaler;
      frodaAtom[whichA].randomMove.y *= scaler;
      frodaAtom[whichA].randomMove.z *= scaler;

      currentPos[whichA].x += frodaAtom[whichA].randomMove.x;
      currentPos[whichA].y += frodaAtom[whichA].randomMove.y;
      currentPos[whichA].z += frodaAtom[whichA].randomMove.z;
    }

  }

  return;
}


void Froda::newPerturbation( MolFramework &structure){
  if (!moveDihedral ) {
    newRandomMove();
  }
  else {
    //pick between angle perturbation, peptide perturbation, and sidechain perturbation
    //for now, only the first two are running, test on alanine
    double random = genrand_real1();
    if ( random < 0.333) {
      SiteID CA, base;
      pickPeptideDihedral( CA, base, structure );
      perturbBackboneDihedral( CA, base, structure );
    }
    else if (random < 0.666){
      SiteID CA1, b1, b2, CA2;
      pickPeptidePlane( CA1, b1, b2, CA2, structure );
      perturbPeptidePlane( CA1, b1, b2, CA2, structure ); 
    }
    else {
      twitchSidechain( structure );
    }
  }


  return;
}


void Froda::applyCMbias() {
      CMWorstFitCycles = cyclesToFit;
      for ( unsigned int clusterNumber = 2; clusterNumber <= myNClusters; clusterNumber++){
        //update the old and older Pos at this point
        if ( i%cacheEvery == 0 ) {
          ghost.at(clusterNumber).olderPosCentral = ghost.at(clusterNumber).oldPosCentral;
          ghost.at(clusterNumber).oldPosCentral = ghost.at(clusterNumber).posCentral;
          ghost.at(clusterNumber).olderRotor = ghost.at(clusterNumber).oldRotor;
          ghost.at(clusterNumber).oldRotor = ghost.at(clusterNumber).runningRotor;
          ghost.at(clusterNumber).oldRotor.x += ghost.at(clusterNumber).baseRotor.x;
          ghost.at(clusterNumber).oldRotor.y += ghost.at(clusterNumber).baseRotor.y;
          ghost.at(clusterNumber).oldRotor.z += ghost.at(clusterNumber).baseRotor.z;
        }
      }
      if (useCMforward && CMWorstFitCycles < cmlowerFitcycles ) {
        for ( unsigned int clusterNumber = 2; clusterNumber <= myNClusters; clusterNumber++){
          //cerr << "Applying momentum to body " << clusterNumber << endl;
          applyMomentum( clusterNumber);
        }
      }
      else if (useCMreverse && CMWorstFitCycles > cmupperFitcycles  ){
        for ( unsigned int clusterNumber = 2; clusterNumber <= myNClusters; clusterNumber++){
          reverseMomentum( clusterNumber);
        }
      }  


  return;
}


