#include "global_defs.h"
#include "generalUtils.h"
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

////////////////////////////////////////////////////////////////////////////////
// Description: The geometric fitting of ghosts to atoms. Fits ghost bodies over
//   real atoms using rotor equations
////////////////////////////////////////////////////////////////////////////////
void Froda::fitGhosts(){ 

  //LABELFIT
  double step;
  int atoma;
  double small_grad_square =pow(smallgrad,2);
  bool update_base=false;
  register double pmet; // one over temp

  Vector grad0, grad1; //gradients here and at  x+sigma d
  double g02, g12; // moduli squared
  double alpha; // step length from line minim
  Vector trial; // trial rotor value

  // Turn atoms positions into a new posCentralit and target bonds 
  // rotor-fit the ghost bonds to the target bonds 
  // apply the rotor so that current ghost is in the right place and orientation
  
  // the new posCentral comes from the currentPos of the atoms
  // the ideal bonds are in ghost_bondVector
  // the real bonds are made from currentPos
  // the rotated ideal bonds go into ghost_pos
  
  unsigned int ghost_size;

  double aweight = 0.0;
  
  for( unsigned int clusterNumber = 2; clusterNumber <= myNClusters; clusterNumber++ ){

    //cerr << "Cluster number: " << clusterNumber << endl;

    ghost.at(clusterNumber).posCentral.x=0;
    ghost.at(clusterNumber).posCentral.y=0;
    ghost.at(clusterNumber).posCentral.z=0;
 
    ghost_size = ghost[clusterNumber].bodysize;

    //redefine the bond vectors if the weights have changed
    //the change of central pos is sum ( new weight * bondvector ) / sum ( new weight )
    //////////////////////////////////////////////////////////////////////
    if ( ghost[clusterNumber].changedWeight ) {
      if ( ghost[clusterNumber].bodysize == 1 ) 
	continue;
      Vector newC = NULL_VEC;
      ghost.at(clusterNumber).bodyweight = 0;
      ghost.at(clusterNumber).radMS = 0;

      for( unsigned int ghostVertex = 0; ghostVertex < ghost_size; ghostVertex++ ){
        atoma = ghost[clusterNumber].atoms.at(ghostVertex);
        //double aweight = frodaAtom.at(atoma).weight;
	aweight = frodaAtom.at(atoma).weight;
        ghost.at(clusterNumber).bodyweight += aweight;

        newC.x += aweight * ghost[clusterNumber].bondVector[ghostVertex].x;
        newC.y += aweight * ghost[clusterNumber].bondVector[ghostVertex].y;
        newC.z += aweight * ghost[clusterNumber].bondVector[ghostVertex].z;
      }
      newC.x /= ghost.at(clusterNumber).bodyweight;
      newC.y /= ghost.at(clusterNumber).bodyweight;
      newC.z /= ghost.at(clusterNumber).bodyweight;

      //this newC change is due to reweighting, not motion
      //so pass into CM function
      ghost.at(clusterNumber).oldPosCentral.x += newC.x;
      ghost.at(clusterNumber).oldPosCentral.y += newC.y;
      ghost.at(clusterNumber).oldPosCentral.z += newC.z;
      ghost.at(clusterNumber).olderPosCentral.x += newC.x;
      ghost.at(clusterNumber).olderPosCentral.y += newC.y;
      ghost.at(clusterNumber).olderPosCentral.z += newC.z;

      //cerr << "Center of body " << clusterNumber << " reweighted by ";
      //cerr << newC.x << " " << newC.y << " " << newC.z << endl;

      for( unsigned int ghostVertex = 0; ghostVertex < ghost_size; ghostVertex++ ){
        ghost[clusterNumber].bondVector[ghostVertex].x -= newC.x;
        ghost[clusterNumber].bondVector[ghostVertex].y -= newC.y;
        ghost[clusterNumber].bondVector[ghostVertex].z -= newC.z;

        double rad2 = dotProduct( ghost[clusterNumber].bondVector.at(ghostVertex),
                                  ghost[clusterNumber].bondVector.at(ghostVertex) );

        atoma = ghost[clusterNumber].atoms.at(ghostVertex);
	    aweight = frodaAtom.at(atoma).weight;
        ghost.at(clusterNumber).radMS += aweight * rad2;
      }
      ghost.at(clusterNumber).radMS /= ghost.at(clusterNumber).bodyweight;
      ghost.at(clusterNumber).radRMS = sqrt(ghost.at(clusterNumber).radMS);
    }

    for (unsigned int ghostAtomNumber=0; ghostAtomNumber < ghost_size; ghostAtomNumber++){
      atoma = ghost[clusterNumber].atoms.at(ghostAtomNumber);
      aweight = frodaAtom.at(atoma).weight;
      ghost.at(clusterNumber).posCentral.x += aweight * currentPos.at(atoma).x;
      ghost.at(clusterNumber).posCentral.y += aweight * currentPos.at(atoma).y;
      ghost.at(clusterNumber).posCentral.z += aweight * currentPos.at(atoma).z;
    }
 
    pmet = 1.0/ghost.at(clusterNumber).bodyweight;
    ghost.at(clusterNumber).posCentral.x *= pmet;
    ghost.at(clusterNumber).posCentral.y *= pmet;
    ghost.at(clusterNumber).posCentral.z *= pmet;

    for( unsigned int ghostAtomNumber = 0; ghostAtomNumber < ghost_size; ghostAtomNumber++ ){
      atoma = ghost[clusterNumber].atoms.at(ghostAtomNumber);
      ghost[clusterNumber].posAtom.at(ghostAtomNumber).x = currentPos.at(atoma).x - ghost.at(clusterNumber).posCentral.x;
      ghost[clusterNumber].posAtom.at(ghostAtomNumber).y = currentPos.at(atoma).y - ghost.at(clusterNumber).posCentral.y;
      ghost[clusterNumber].posAtom.at(ghostAtomNumber).z = currentPos.at(atoma).z - ghost.at(clusterNumber).posCentral.z;
    }


  }

  // Now find the rotor that matches the ghost_bondVector bonds to posAtom
  // use a line minimiser (secant method)
  
  for( unsigned int clusterNumber = 2;clusterNumber <= myNClusters ; clusterNumber++ ){
    ghost_size = ghost[clusterNumber].bodysize;

    if ( ghost_size == 1){ // single-site ghost
      ghost.at(clusterNumber).pos.at(0).x = currentPos.at( ghost[clusterNumber].atoms.at(0) ).x;
      ghost.at(clusterNumber).pos.at(0).y = currentPos.at( ghost[clusterNumber].atoms.at(0) ).y;
      ghost.at(clusterNumber).pos.at(0).z = currentPos.at( ghost[clusterNumber].atoms.at(0) ).z;

      continue;
    }
    rotor = ghost.at(clusterNumber).runningRotor;
    update_base = false;

    step = 1/ (2.0*ghost.at(clusterNumber).bodyweight * ghost.at(clusterNumber).radMS);

    for(unsigned int rotorCycle=0; rotorCycle < maxRotorCycles ; rotorCycle++){

      bool leave = false;

      grad0 = NULL_VEC;
      grad1 = NULL_VEC;

      for(unsigned int ghostAtomNumber = 0; ghostAtomNumber < ghost_size; ghostAtomNumber++ ){

        atoma = ghost[clusterNumber].atoms.at(ghostAtomNumber);

	// BMH 8/30/06 - replaced next 7 lines with getGradient function
	//////////////////////////////////////////////////////////////////////
        //aweight = frodaAtom.at(atoma).weight;
        //dot = deps2dbx(ghost[clusterNumber].bondVector.at(ghostAtomNumber), ghost[clusterNumber].posAtom.at(ghostAtomNumber), rotor);
        //grad0.x += aweight *dot;
        //dot = deps2dby(ghost[clusterNumber].bondVector.at(ghostAtomNumber), ghost[clusterNumber].posAtom.at(ghostAtomNumber), rotor);
        //grad0.y += aweight *dot;
        //dot = deps2dbz(ghost[clusterNumber].bondVector.at(ghostAtomNumber), ghost[clusterNumber].posAtom.at(ghostAtomNumber), rotor);
        //grad0.z += aweight *dot;

	grad0 += getGradient( frodaAtom.at(atoma).weight,
			      ghost[clusterNumber].bondVector.at(ghostAtomNumber), 
			      ghost[clusterNumber].posAtom.at(ghostAtomNumber), 
			      rotor );
      } // end loop over atoms

      //grad0 is the gradient at this point

      g02 = dotProduct(grad0, grad0);
      if(g02 < small_grad_square){
        leave = true;
      }

      if (leave) break; //already good enough

      //form an experimental step and get the new gradient;
      tempvec.x = step*grad0.x;
      tempvec.y = step*grad0.y;
      tempvec.z = step*grad0.z;

      dot = ( dotProduct(tempvec,tempvec));
      if (dot > large_rotor_2){ // too big clusterNumber step
        tempvec.x *= large_rotor_2/dot;
        tempvec.y *= large_rotor_2/dot;
        tempvec.z *= large_rotor_2/dot;
      }
      trial.x = rotor.x - tempvec.x;
      trial.y = rotor.y - tempvec.y;
      trial.z = rotor.z - tempvec.z;

      //calculate the gradient at the trial position 

      //for( unsigned int ghostAtomNumber = 0; ghostAtomNumber < ghost_size; ghostAtomNumber++ ){
      for( unsigned int ghostAtomNumber = 0; ghostAtomNumber < ghost[clusterNumber].atoms.size(); ghostAtomNumber++ ){
        atoma = ghost[clusterNumber].atoms.at(ghostAtomNumber);

	// BMH 8/30/06 - replaced next 7 lines with getGradient function
	//////////////////////////////////////////////////////////////////////
	//aweight = frodaAtom.at(atoma).weight;
	//dot = deps2dbx(ghost[clusterNumber].bondVector.at(ghostAtomNumber), ghost[clusterNumber].posAtom.at(ghostAtomNumber), trial);
        //grad1.x += aweight *dot;
        //dot = deps2dby(ghost[clusterNumber].bondVector.at(ghostAtomNumber), ghost[clusterNumber].posAtom.at(ghostAtomNumber), trial);
        //grad1.y += aweight *dot;
        //dot = deps2dbz(ghost[clusterNumber].bondVector.at(ghostAtomNumber), ghost[clusterNumber].posAtom.at(ghostAtomNumber), trial);
        //grad1.z += aweight *dot;

	grad1 += getGradient( frodaAtom.at(atoma).weight, 
			      ghost[clusterNumber].bondVector.at(ghostAtomNumber), 
			      ghost[clusterNumber].posAtom.at(ghostAtomNumber), 
			      trial );

      } // end loop over atoms

      //grad1 is the gradient at the trial position;
      g12 = dotProduct(grad1, grad1);
      if(g12 < small_grad_square){
        leave = true;
        rotor = trial;
      }

      if (leave) break; //the trial position is good- carry on

      //if we get this far, then we need to calculate clusterNumber step and iterate again
      //we are either in the quadratic regime or not
      //so either steepest-descent or secant
      //if using SD, we already have our tempvec step
      //if using secant, recalculate

      if ( (rotorCycle > 2) && (g12 < g02)) { //use secant
        alpha = 1.0 / (1.0 - (dotProduct(grad0, grad1)/(g02)));

        tempvec.x *= alpha;
        tempvec.y *= alpha;
        tempvec.z *= alpha;
      }

      dot = ( dotProduct(tempvec,tempvec));
      if (dot > large_rotor_2){ // too big clusterNumber step
        //cout << "Body " << clusterNumber << " has rotor step " << dot <<  "; Rotor " << rotor.x << ", " << rotor.y << ", " << rotor.z << endl;
        tempvec.x *= large_rotor_2/dot;
        tempvec.y *= large_rotor_2/dot;
        tempvec.z *= large_rotor_2/dot;
      }

      //Finally - rotor update

      rotor.x -= tempvec.x;
      rotor.y -= tempvec.y;
      rotor.z -= tempvec.z;
     
      dot = dotProduct (rotor, rotor);
      if (dot > large_rotor_2) update_base=true;
      

      // if the rotor is large, redefine the ghost_pos_base to the new orientation

      if (update_base){
        //cerr << "Updating ghost for body " << clusterNumber <<  " at " << rotorCycle+1 << " " << fitCycles << endl;
        myX = Xval(rotor);
        for( unsigned int ghostSiteNumber = 0; ghostSiteNumber < ghost_size; ghostSiteNumber++ ){
          //tempvec.x = epsx(ghost[clusterNumber].bondVector.at(ghostSiteNumber), 0.0, rotor, myX);
          //tempvec.y = epsy(ghost[clusterNumber].bondVector.at(ghostSiteNumber), 0.0, rotor, myX);
          //tempvec.z = epsz(ghost[clusterNumber].bondVector.at(ghostSiteNumber), 0.0, rotor, myX);
          tempvec.x = epsx(ghost[clusterNumber].bondVector[ghostSiteNumber], 0.0, rotor, myX);
          tempvec.y = epsy(ghost[clusterNumber].bondVector[ghostSiteNumber], 0.0, rotor, myX);
          tempvec.z = epsz(ghost[clusterNumber].bondVector[ghostSiteNumber], 0.0, rotor, myX);
        
          //ghost[clusterNumber].bondVector.at(ghostSiteNumber).x = tempvec.x;
          //ghost[clusterNumber].bondVector.at(ghostSiteNumber).y = tempvec.y;
          //ghost[clusterNumber].bondVector.at(ghostSiteNumber).z = tempvec.z;
          ghost[clusterNumber].bondVector[ghostSiteNumber].x = tempvec.x;
          ghost[clusterNumber].bondVector[ghostSiteNumber].y = tempvec.y;
          ghost[clusterNumber].bondVector[ghostSiteNumber].z = tempvec.z;
        }
        ghost[clusterNumber].baseRotor.x += rotor.x;
        ghost[clusterNumber].baseRotor.y += rotor.y;
        ghost[clusterNumber].baseRotor.z += rotor.z;
        rotor.x = rotor.y = rotor.z =0.0;
        update_base =false;
      }
 
    } // end loop over iterations

    // create ghost_pos with the new target positions for the atoms
    
    myX = Xval(rotor);
    ghost.at(clusterNumber).runningRotor = rotor;
    for( unsigned int ghostNumber = 0; ghostNumber < ghost_size; ghostNumber++ ){
      //ghost[clusterNumber].pos.at(ghostNumber).x = epsx(ghost[clusterNumber].bondVector.at(ghostNumber), 0.0, rotor, myX) + ghost.at(clusterNumber).posCentral.x;
      //ghost[clusterNumber].pos.at(ghostNumber).y = epsy(ghost[clusterNumber].bondVector.at(ghostNumber), 0.0, rotor, myX) + ghost.at(clusterNumber).posCentral.y;
      //ghost[clusterNumber].pos.at(ghostNumber).z = epsz(ghost[clusterNumber].bondVector.at(ghostNumber), 0.0, rotor, myX) + ghost.at(clusterNumber).posCentral.z;  
      ghost[clusterNumber].pos.at(ghostNumber).x = epsx(ghost[clusterNumber].bondVector[ghostNumber], 0.0, rotor, myX) + ghost.at(clusterNumber).posCentral.x;
      ghost[clusterNumber].pos.at(ghostNumber).y = epsy(ghost[clusterNumber].bondVector[ghostNumber], 0.0, rotor, myX) + ghost.at(clusterNumber).posCentral.y;
      ghost[clusterNumber].pos.at(ghostNumber).z = epsz(ghost[clusterNumber].bondVector[ghostNumber], 0.0, rotor, myX) + ghost.at(clusterNumber).posCentral.z;  
    }
    
  } // end loop over clusters  
  
  return;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description: The matching of atoms to ghosts. Calculates the mismatch relative
//   to ghosts.
////////////////////////////////////////////////////////////////////////////////
int Froda::findGhostMismatch(){ 
  //LABELGHOSTMISMATCH

  int return_val = 0;

  worstGhostMismatch = 0.0;
  worstGhostAtom = 0;
  worstGhost = 0;

  unsigned int myatom;
  double ghostTol_square=pow(ghostTol,2);
  // for each ghost vertex, give its atom a mismatch vector
  // If all mismatches are less than tolerance, return 0 (all is good)
  // otherwise, return n mismatches - still working

  for( unsigned int clusterNumber = 2;clusterNumber <= myNClusters ; clusterNumber++ ){
    ghost[clusterNumber].changedWeight = false;
  }
  
  for( unsigned int clusterNumber = 2;clusterNumber <= myNClusters ; clusterNumber++ ){

    for( unsigned int ghostSiteNumber = 0; ghostSiteNumber < ghost[clusterNumber].bodysize; ghostSiteNumber++ ){
      unsigned int atoma = ghost[clusterNumber].atoms.at(ghostSiteNumber);
      
      frodaAtom.at(atoma).nMainMismatches++;
      
      dummyvec.x = ghost[clusterNumber].pos.at(ghostSiteNumber).x - currentPos.at(atoma).x;
      dummyvec.y = ghost[clusterNumber].pos.at(ghostSiteNumber).y - currentPos.at(atoma).y;
      dummyvec.z = ghost[clusterNumber].pos.at(ghostSiteNumber).z - currentPos.at(atoma).z;
      
      dot = dotProduct(dummyvec,dummyvec);

      if ( dot > worstGhostMismatch ) {
        worstGhostMismatch = dot;
        worstGhostAtom = atoma;
        worstGhost = clusterNumber;
      }

      if (checkingWeights) {
        int factor = (int) (2*dot/ghostTol_square);
        if ( factor > 0 ) {
          frodaAtom.at(atoma).hasWeight = true;
          frodaAtom.at(atoma).weightIncreasing = true;
          frodaAtom.at(atoma).weightDecreasing = false;
          //if ( frodaAtom.at(atoma).weight > 200 ) continue;
          frodaAtom.at(atoma).weight += factor;
          //cerr << "Atom " << atoma << " gained weight (" << frodaAtom.at(atoma).weight << ") due to mismatch in body " << clusterNumber << endl;

          unsigned int gSize = frodaAtom.at(atoma).ghostlist.size(); 
          for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
            ghost[ frodaAtom.at(atoma).ghostlist.at(whichGhost) ].changedWeight = true;
            //cerr << "Ghost " << frodaAtom.at(atoma).ghostlist.at(whichGhost) << " gained weight from atom " << atoma << endl;
          }
        }
      }

      if (dot > ghostTol_square){
        return_val++; 

        if (fitCycles > cyclesToFit) cyclesToFit = fitCycles;

        //pass information to body-response
        if (do_bodyResponse){
          //cluster is clusterNumber
          //atom is atoma
          //mismatch is dummyvec

          //put a mismatch on the rest of cluster
          //carrying them along with us
          //apply to atoms in the same half of the ghost
          //mismatch is -dummyvec

            Vector passvec;
            if ( frodaAtom.at(atoma).isCore ) continue;
            passvec.x = - dummyvec.x;
            passvec.y = - dummyvec.y;
            passvec.z = - dummyvec.z;
            

              Vector abond;
              abond.x = currentPos.at(atoma).x - ghost[clusterNumber].posCentral.x;
              abond.y = currentPos.at(atoma).y - ghost[clusterNumber].posCentral.y;
              abond.z = currentPos.at(atoma).z - ghost[clusterNumber].posCentral.z;
  
              unsigned int csize = ghost[clusterNumber].bodysize;
              for (unsigned int secondGhostSiteNumber=0; secondGhostSiteNumber < csize; secondGhostSiteNumber++){
                if (ghost[clusterNumber].atoms.at(secondGhostSiteNumber)==atoma) continue;
                myatom = ghost[clusterNumber].atoms.at(secondGhostSiteNumber);
                Vector bbond;
                bbond.x = ghost[clusterNumber].pos.at(secondGhostSiteNumber).x - ghost[clusterNumber].posCentral.x;
                bbond.y = ghost[clusterNumber].pos.at(secondGhostSiteNumber).y - ghost[clusterNumber].posCentral.y;
                bbond.z = ghost[clusterNumber].pos.at(secondGhostSiteNumber).z - ghost[clusterNumber].posCentral.z;
                if ( dotProduct( abond,bbond) < 0.0 ) continue; //wrong half of the body 
                frodaAtom.at(myatom).nMainMismatches ++;
                frodaAtom.at(myatom).mainMismatch.x += passvec.x;
                frodaAtom.at(myatom).mainMismatch.y += passvec.y;
                frodaAtom.at(myatom).mainMismatch.z += passvec.z;
                
              }
        } //end of body response clause

        stringstream probLine;
        probLine << "G" << clusterNumber << "A" << atoma << " Bsize" << ghost[clusterNumber].bodysize;
        probLine << " Aweight" << frodaAtom.at(atoma).weight;
        probLine << " Alist ";
        for (unsigned int ghostlistNumber =0; ghostlistNumber< frodaAtom.at(atoma).ghostlist.size(); ghostlistNumber++){
          probLine << " " << frodaAtom.at(atoma).ghostlist.at(ghostlistNumber);
          if ( frodaAtom.at(atoma).ghostlist.at(ghostlistNumber) > nTotalClusters) probLine << "(H)";
        }
        problemReport.push_back(probLine.str() );
        probLine.clear();

        if(debug_ghost && (fitCycles > 1)){
          cerr << "Body " << clusterNumber << ", Atom " << atoma << " off ghostpos by: " << sqrt(dot) << endl;
          cerr << "(Body has " << ghost[clusterNumber].bodysize << " members)";
          cerr << "(Atom is in " <<  frodaAtom.at(atoma).ghostlist.size() << " ghosts:" ;
          for (unsigned int ghostlistNumber =0; ghostlistNumber< frodaAtom.at(atoma).ghostlist.size(); ghostlistNumber++){
            cerr << " " << frodaAtom.at(atoma).ghostlist.at(ghostlistNumber);
            if ( frodaAtom.at(atoma).ghostlist.at(ghostlistNumber) > nTotalClusters) cerr << "(H)";
          }
          cerr << endl;
        }
      }
       if(dot != dot){ // NaN!
        cerr << "FATAL NaN ERROR" << endl;
        cerr << "Body " << clusterNumber << ", Atom " << atoma << " off ghostpos by NaN!" << endl;
        cerr << "(Body has " << ghost[clusterNumber].bodysize << " members)";
        cerr << "(Atom is in " << frodaAtom.at(atoma).ghostlist.size() << " ghosts)" << endl;
        isBad=true;
        dummyvec=NULL_VEC;
      } 
      
      frodaAtom.at(atoma).mainMismatch.x += dummyvec.x;
      frodaAtom.at(atoma).mainMismatch.y += dummyvec.y;
      frodaAtom.at(atoma).mainMismatch.z += dummyvec.z;
             
    }
    
  }
 
  worstGhostMismatch = sqrt( worstGhostMismatch ); 
 
  if (verbose && chatty_fitting) cerr << "GHOSTS" << setw(5) << return_val << "#";
  return(return_val);
}
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////
// Description: Prevent stretching of hydrophobic tethers
////////////////////////////////////////////////////////////////////////////////
int Froda::phobicMismatch(){

  double ghostTol_square=pow(ghostTol,2); 

  worstTetherStretch = 0.0;
  worstStretchAtom1 = 0;
  worstStretchAtom2 = 0;
  
  unsigned int atom1,atom2, myatom;
  double rad;
  register double pmet;
  int return_val = 0;
  // during initialise, we got a list of active ph constraints
  // The two atoms are not mutually rigid but have a distance tether
  // For each such pair, check distance and apply correction if need by 
  // If all distances are tolerable, return 0
  
  // There are n_myPhobes myPhobe objects in the vector Phobes 
  // each one has two atoms A,B and a distance AB 
  // find the distance AB', compare to AB, apply correction to A and B
  // if AB' exceeds a threshold, return 1
  
  for (int p=0; p<n_myPhobes;p++){
    atom1 = Phobes.at(p).atomA;
    atom2 = Phobes.at(p).atomB; 
    rad = Phobes.at(p).AB;

    double rad1 = rad + phobeTol ;
    double rad2 = rad + phobeTol + ghostTol;
    
    dummyvec.x = currentPos.at(atom1).x - currentPos.at(atom2).x;
    dummyvec.y = currentPos.at(atom1).y - currentPos.at(atom2).y;
    dummyvec.z = currentPos.at(atom1).z - currentPos.at(atom2).z;
    
    dot = dotProduct(dummyvec,dummyvec);
    veclength = sqrt(dot);
   
    
 
    if (veclength > (rad2) ){
      return_val++;

      //warn the ghosts;
      cyclesToFit = fitCycles;

      stringstream probLine;
      probLine << "R(" << atom1 << "--" << atom2 << ") = " << veclength;
      problemReport.push_back(probLine.str() );
      probLine.clear();

      if(debug_phobe){
        cerr << "R(" << atom1 << "--" << atom2 << ") = " << veclength << endl;
      }
    }
    
    if (veclength > rad1){ 
      pmet = 1.0/veclength;
      tempvec.x = dummyvec.x*pmet;
      tempvec.y = dummyvec.y*pmet;
      tempvec.z = dummyvec.z*pmet;
      //tempvec is now a unit vector from atom2 to atom1
      pmet = rad1-veclength; // negative as vec > rad; 

      if ( pmet < worstTetherStretch ) {
        worstTetherStretch = pmet;
        worstStretchAtom1 = atom1;
        worstStretchAtom2 = atom2;
      }

      if (checkingWeights) {
        int factor = (int) (2*pmet*pmet/ghostTol_square);
        if ( factor > 0 ) {
          unsigned int atoma = atom1;
          frodaAtom.at(atoma).hasWeight = true;
          frodaAtom.at(atoma).weightIncreasing = true;
          frodaAtom.at(atoma).weightDecreasing = false;
          frodaAtom.at(atoma).weight += factor;
          //cerr << "Atom " << atoma << " gained weight (" << frodaAtom.at(atoma).weight << ") due to mismatch in body " << clusterNumber << endl;
          unsigned int gSize = frodaAtom.at(atoma).ghostlist.size(); 
          for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
            ghost[ frodaAtom.at(atoma).ghostlist.at(whichGhost) ].changedWeight = true;
            //cerr << "Ghost " << frodaAtom.at(atoma).ghostlist.at(whichGhost) << " gained weight from atom " << atoma << endl;
          }
          atoma = atom2;
          frodaAtom.at(atoma).hasWeight = true;
          frodaAtom.at(atoma).weightIncreasing = true;
          frodaAtom.at(atoma).weightDecreasing = false;
          frodaAtom.at(atoma).weight += factor;
          //cerr << "Atom " << atoma << " gained weight (" << frodaAtom.at(atoma).weight << ") due to mismatch in body " << clusterNumber << endl;
          gSize = frodaAtom.at(atoma).ghostlist.size(); 
          for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
            ghost[ frodaAtom.at(atoma).ghostlist.at(whichGhost) ].changedWeight = true;
            //cerr << "Ghost " << frodaAtom.at(atoma).ghostlist.at(whichGhost) << " gained weight from atom " << atoma << endl;
          }
        }
      }

      //if one atom is the core, the other must take responsibility for the correction
      //otherwise they get half a mismatch vector each
      if ( frodaAtom.at(atom1).isCore || frodaAtom.at(atom2).isCore ) {
      } else {
        pmet *= 0.5;
      }

      tempvec.x *= pmet;
      tempvec.y *= pmet;
      tempvec.z *= pmet; //tempvec now mismatch vector
      
      frodaAtom.at(atom1).nMainMismatches++;
      frodaAtom.at(atom2).nMainMismatches++;
            
      frodaAtom.at(atom1).mainMismatch.x += tempvec.x;
      frodaAtom.at(atom1).mainMismatch.y += tempvec.y;
      frodaAtom.at(atom1).mainMismatch.z += tempvec.z;
      frodaAtom.at(atom2).mainMismatch.x -= tempvec.x;      
      frodaAtom.at(atom2).mainMismatch.y -= tempvec.y;
      frodaAtom.at(atom2).mainMismatch.z -= tempvec.z;                        

      //if it's body-response time, drag the clusters along, too.
      //do each atom in turn; do each cluster for each atom

      if ( do_bodyResponse && (veclength > rad2)){
          for( int which = 0; which < 2; which++){
            unsigned int atoma;
            Vector passvec;
            if ( which ==0){
              atoma = atom1; 
              if ( frodaAtom.at(atoma).isCore ) continue;
              passvec = tempvec;
            }
            else{
              atoma = atom2;
              if ( frodaAtom.at(atoma).isCore ) continue;
              passvec.x = - tempvec.x;
              passvec.y = - tempvec.y;
              passvec.z = - tempvec.z;
            }
  
            int gsize = frodaAtom.at(atoma).ghostlist.size();
            for ( int whichc = 0; whichc < gsize; whichc++){
              int clusterNumber = frodaAtom.at(atoma).ghostlist.at(whichc);
              if ( clusterNumber == 1) continue; //don't push the rigid core
  
              Vector abond;
              abond.x = currentPos.at(atoma).x - ghost[clusterNumber].posCentral.x;
              abond.y = currentPos.at(atoma).y - ghost[clusterNumber].posCentral.y;
              abond.z = currentPos.at(atoma).z - ghost[clusterNumber].posCentral.z;
  
              unsigned int csize = ghost[clusterNumber].bodysize;
              for (unsigned int secondGhostSiteNumber=0; secondGhostSiteNumber < csize; secondGhostSiteNumber++){
                if (ghost[clusterNumber].atoms.at(secondGhostSiteNumber)==atoma) continue;
                myatom = ghost[clusterNumber].atoms.at(secondGhostSiteNumber);
                Vector bbond;
                bbond.x = ghost[clusterNumber].pos.at(secondGhostSiteNumber).x - ghost[clusterNumber].posCentral.x;
                bbond.y = ghost[clusterNumber].pos.at(secondGhostSiteNumber).y - ghost[clusterNumber].posCentral.y;
                bbond.z = ghost[clusterNumber].pos.at(secondGhostSiteNumber).z - ghost[clusterNumber].posCentral.z;
                if ( dotProduct( abond,bbond) < 0.0 ) continue; //wrong half of the body 
                frodaAtom.at(myatom).nMainMismatches ++;
                frodaAtom.at(myatom).mainMismatch.x += passvec.x;
                frodaAtom.at(myatom).mainMismatch.y += passvec.y;
                frodaAtom.at(myatom).mainMismatch.z += passvec.z;
                
              }
            }
          }
      } //end of body response clause

    }
    
  }

  worstTetherStretch *= -1;
  if (verbose && chatty_fitting) cerr << "HYDROPHOBIC" << setw(3) << return_val << "#";
  return(return_val);
}
////////////////////////////////////////////////////////////////////////////////

void Froda::getStericParms( const MolFramework& structure, SiteID current_atom, SiteID other_atom, double& rInner, double& r0ptimum, double& rOuter ) {
  double R0 = 0.0; //"optimum" distance
  double Rout = 0.0;//begin corrections at this distance
  double Rin = 0.0;//reject conformer at this distance
  int myclass = 0; //classifies the "type" of interaction
  // myclass 1 covers interactions with a defined R0
  //myclass 2 covers interactions that use the "notice" parameter
  
  bool attractive = false;

  //here R0 and myclass are set
  if (polarGeometry ){ //might be in polar attractive regime- check for hbonds
    attractive = false;
    if (frodaAtom.at(current_atom).charge * frodaAtom.at(other_atom).charge < 0 ) attractive = true;
  }

  if (attractive){
    if ( frodaAtom.at(current_atom).isHydrogen ){
      R0 = polarHRadius + structure.site_info[other_atom].vdw_radius;
      myclass = 1;
    }
    else if (frodaAtom.at(other_atom).isHydrogen ){
      R0 = polarHRadius + structure.site_info[current_atom].vdw_radius;
      myclass = 1;
    } else{
      R0 = structure.site_info[current_atom].vdw_radius + structure.site_info[other_atom].vdw_radius;
      R0 *= 0.92; //tuned to backbone C O in alpha helix
      myclass = 1;
    }
  }
  else{
    R0 = structure.site_info[current_atom].vdw_radius + structure.site_info[other_atom].vdw_radius;
    myclass = 2;

    //check for unscreened repulsions
    if (polarGeometry && frodaAtom.at(current_atom).charge < 0 && frodaAtom.at(other_atom).charge < 0 ){
      if ( !hasHBondGeometry(current_atom, other_atom, structure.site_info[current_atom].neighbor_list,
         structure.site_info[other_atom].neighbor_list) ){
         myclass =1; //repulsive interaction
      }
      else{ //hbond relation
        R0 *= 0.92; // C 1.7 + O 1.4 = 3.1; * 0.92 ~ 2.85 Angstroms for hbond
        myclass = 1;
      }
    }

    //don't let hydrogens get too close together
    if ( frodaAtom.at(current_atom).isHydrogen && frodaAtom.at(other_atom).isHydrogen ) {
      myclass = 1;
    }

  }

  if ( isThirdNeighbor(current_atom, other_atom) ){
      R0 *= 0.8;
      myclass = 1;
  }

  //here Rout and Rin are set, using myclass and R0
  if ( myclass == 1){ //fixed outer radius- esp. hbond geometry
    Rout = R0;
    Rin = vdwTol * R0;
  }
  else if ( myclass == 2){ //general case, use tolerance vdwNotice
    Rout = vdwNotice * R0;
    Rin = vdwTol * R0;
  }
  
  //assign output values
  rInner = Rin;
  r0ptimum = R0;
  rOuter = Rout;
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: Detection and correction for collisions
////////////////////////////////////////////////////////////////////////////////
int Froda::stericMismatch( MolFramework& structure ){
  //Timer timer;
  //timer.reset();
  //timer.start();
  int return_val = 0;
  int myatom;
  //was previously register double pmet
  double pmet;
  bool myflag;
  double overlap;

  worstClash =0.0;
  worstClashAtom1 = 0;
  worstClashAtom2 = 0;

  for (int bin=0; bin < 7; bin++){
    overlapDegrees[bin] = 0;
  }  
  nOverlaps = 0;  

  SiteID current_atom, other_atom;
  Vector delta;
  vector<SiteID> possibleContacts;
  
  //int totPairs=0;
  //int totFilteredPairs=0;
  //int totFiltered2=0;
  //int totClashingPairs=0;

  for ( current_atom = 1; current_atom <= nTotalSites; current_atom++ ) {

    double rInner;
    double rOptimum;
    double rOuter;
    double rad2;
    double distToCheck = structure.site_info[current_atom].vdw_radius + sterics->getMaxVdwRadius();
    possibleContacts = sterics->getProximityList( current_atom );

    for ( size_t i = 0; i < possibleContacts.size(); i++ ) {

      other_atom = possibleContacts[i];

      if ( other_atom <= current_atom ) continue; //sterics doesn't filter this

      delta = currentPos[current_atom] - currentPos[other_atom];
      if ( abs(delta.x) > distToCheck || 
	   abs(delta.y) > distToCheck || 
	   abs(delta.z) > distToCheck ) 
	continue;
     
      getStericParms( structure, current_atom, other_atom, rInner, rOptimum, rOuter );
      if ( abs(delta.x) > rOuter || 
	   abs(delta.y) > rOuter || 
	   abs(delta.z) > rOuter ) 
	continue;

      rad2 = rOuter * rOuter;
      double dot = delta.norm2();

      //here we determine whether this pair is clashing at all
      //by comparing the square distance between the atoms to rad2
      myflag =false;
      if(dot < (rad2)){      
 

        //here the code determines if the pair exceeds tolerance, and cumulates the sum of such pairs
        //if the pair exceeds tolerance it also sets some field in the ghosts of these atoms
        //myflag is set to true if the pair exceeds tolerance
        veclength = sqrt(dot);

          //if ( checkingWeights ) {
          if ( false ) { //no weight on clashes?
            unsigned int atoma = current_atom;
            frodaAtom.at(atoma).hasWeight = true;
            frodaAtom.at(atoma).weightIncreasing = true;
            frodaAtom.at(atoma).weightDecreasing = false;
            //if ( frodaAtom.at(atoma).weight > 200 ) continue;
            frodaAtom.at(atoma).weight += 1;
            unsigned int gSize = frodaAtom.at(atoma).ghostlist.size(); 
            for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
              ghost[ frodaAtom.at(atoma).ghostlist.at(whichGhost) ].changedWeight = true;
            }
            atoma = other_atom;
            frodaAtom.at(atoma).hasWeight = true;
            frodaAtom.at(atoma).weightIncreasing = true;
            frodaAtom.at(atoma).weightDecreasing = false;
            //if ( frodaAtom.at(atoma).weight > 200 ) continue;
            frodaAtom.at(atoma).weight += 1;
            gSize = frodaAtom.at(atoma).ghostlist.size(); 
            for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
              ghost[ frodaAtom.at(atoma).ghostlist.at(whichGhost) ].changedWeight = true;
            }
          }

        if (veclength < (rInner) ){
          return_val++;
          myflag =true;
          stringstream probLine;
          probLine << "Contact exceeds critical limit " << current_atom << " " << other_atom;
          problemReport.push_back(probLine.str() );
          probLine.clear();
          if ( debug_vdw ) cerr << "Contact exceeds critical limit, " << current_atom << ", " << other_atom << endl;
          cyclesToFit = fitCycles;
        }
      
 
        //calculates the vector response to apply to the atoms
        overlap = rOuter - veclength;
        pmet = overlap/veclength;

        if ( overlap > worstClash ) {
          worstClash = overlap;
          worstClashAtom1 = current_atom;
          worstClashAtom2 = other_atom;
        } 
  
        //statistics
        int whichbin = (int) ((rOptimum - veclength)/overlapStep);
        if (whichbin > 6) whichbin = 6;
        overlapDegrees.at(whichbin)++;
        nOverlaps++;

        if ( checkingWeights ) {
        //if ( false ) { //no weight on clashes?
          unsigned int atoma = current_atom;
          frodaAtom.at(atoma).hasWeight = true;
          frodaAtom.at(atoma).weightIncreasing = true;
          frodaAtom.at(atoma).weightDecreasing = false;
          frodaAtom.at(atoma).weight += whichbin;
          unsigned int gSize = frodaAtom.at(atoma).ghostlist.size(); 
          for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
            ghost[ frodaAtom.at(atoma).ghostlist.at(whichGhost) ].changedWeight = true;
          }
          atoma = other_atom;
          frodaAtom.at(atoma).hasWeight = true;
          frodaAtom.at(atoma).weightIncreasing = true;
          frodaAtom.at(atoma).weightDecreasing = false;
          frodaAtom.at(atoma).weight += whichbin;
          gSize = frodaAtom.at(atoma).ghostlist.size(); 
          for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
            ghost[ frodaAtom.at(atoma).ghostlist.at(whichGhost) ].changedWeight = true;
          }
        }
  
        dummyvec.x = delta.x*pmet;
        dummyvec.y = delta.y*pmet;
        dummyvec.z = delta.z*pmet;
        
        // mismatch vector for current; negative for other
  
        //if one of the atoms is a core atom, the other must take responsibility for the motion
        //otherwise they get half a mismatch vector each
        if ( frodaAtom.at(current_atom).isCore || frodaAtom.at(other_atom).isCore ){
          //dummyvec.x *= 2.0;
          //dummyvec.y *= 2.0;
          //dummyvec.z *= 2.0;
        } else {
          dummyvec.x *= 0.5;
          dummyvec.y *= 0.5;
          dummyvec.z *= 0.5;
        }
        
        
        //store the response  
        frodaAtom.at(current_atom).nMainMismatches++;
        frodaAtom.at(other_atom).nMainMismatches++;
        frodaAtom.at(current_atom).mainMismatch.x += dummyvec.x;
        frodaAtom.at(current_atom).mainMismatch.y += dummyvec.y;
        frodaAtom.at(current_atom).mainMismatch.z += dummyvec.z;
        frodaAtom.at(other_atom).mainMismatch.x -= dummyvec.x;        
        frodaAtom.at(other_atom).mainMismatch.y -= dummyvec.y;        
        frodaAtom.at(other_atom).mainMismatch.z -= dummyvec.z;                        
  
        //also possibly calculate a bodyResponse, and store it
        if (myflag && do_bodyResponse){
        //if ( false ){
          for( int which = 0; which < 2; which++){
            unsigned int atoma;
            Vector passvec;
            if ( which ==0){
              atoma = current_atom; 
              if ( frodaAtom.at(atoma).isCore ) continue;
              passvec = dummyvec;
            }
            else{
              atoma = other_atom;
              if ( frodaAtom.at(atoma).isCore ) continue;
              passvec.x = - dummyvec.x;
              passvec.y = - dummyvec.y;
              passvec.z = - dummyvec.z;
            }
  
            int gsize = frodaAtom.at(atoma).ghostlist.size();
            for ( int whichc = 0; whichc < gsize; whichc++){
              int clusterNumber = frodaAtom.at(atoma).ghostlist.at(whichc);
              if ( clusterNumber == 1) continue; //don't push the rigid core
  
              Vector abond;
              abond.x = currentPos.at(atoma).x - ghost[clusterNumber].posCentral.x;
              abond.y = currentPos.at(atoma).y - ghost[clusterNumber].posCentral.y;
              abond.z = currentPos.at(atoma).z - ghost[clusterNumber].posCentral.z;
  
              unsigned int csize = ghost[clusterNumber].bodysize;
              for (unsigned int secondGhostSiteNumber=0; secondGhostSiteNumber < csize; secondGhostSiteNumber++){
                if (ghost[clusterNumber].atoms.at(secondGhostSiteNumber)==atoma) continue;
                myatom = ghost[clusterNumber].atoms.at(secondGhostSiteNumber);
                Vector bbond;
                bbond.x = ghost[clusterNumber].pos.at(secondGhostSiteNumber).x - ghost[clusterNumber].posCentral.x;
                bbond.y = ghost[clusterNumber].pos.at(secondGhostSiteNumber).y - ghost[clusterNumber].posCentral.y;
                bbond.z = ghost[clusterNumber].pos.at(secondGhostSiteNumber).z - ghost[clusterNumber].posCentral.z;
                if ( dotProduct( abond,bbond) < 0.0 ) continue; //wrong half of the body 
                frodaAtom.at(myatom).nMainMismatches ++;
                frodaAtom.at(myatom).mainMismatch.x += passvec.x;
                frodaAtom.at(myatom).mainMismatch.y += passvec.y;
                frodaAtom.at(myatom).mainMismatch.z += passvec.z;
                
              }
            }
          }
        }        
      }
    }
  }
  //timer.stop();
  //cout << "  stericMismatch time " << timer.gettime() << " ms. " << endl; 
 
  if (verbose && chatty_fitting) cerr << "CONTACT" << setw(3) << return_val << "#";
  return(return_val);
}



////////////////////////////////////////////////////////////////////////////////
// Description: reset the mismatches to zero
////////////////////////////////////////////////////////////////////////////////
void Froda::clearMismatchArrays(){

      for( unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
        frodaAtom.at(siteNumber).mainMismatch = NULL_VEC;
        frodaAtom.at(siteNumber).nMainMismatches = 0;
      }
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: keep the system inside a defined cylindrical block of space
////////////////////////////////////////////////////////////////////////////////
void Froda::applyHatbox() {

  //cerr << "Running Hatbox constraint." << endl;
  Vector radial;
  Vector correction;
  double radius;
  double over;

  double za = parameters.hatboxZa;
  double zb = parameters.hatboxZb;
  double xy = parameters.hatboxXY;
  
  for ( unsigned int atom = 1; atom <= nTotalSites; atom++ ) {
    if ( frodaAtom.at(atom).isCore ) continue;

    radial = currentPos.at(atom);
    radial.z = 0;
    radius = sqrt( dotProduct ( radial, radial ) );

    correction = NULL_VEC;

    //first check the Z factor

    if ( currentPos.at(atom).z < za ) {
      over = za - currentPos.at(atom).z;
      correction.z += over;
    }
    else if ( currentPos.at(atom).z > zb ) {
      over = zb - currentPos.at(atom).z;
      correction.z += over;
    }

    //now check the xy factor

    if ( radius > xy ) {
      radial.x /= radius; 
      radial.y /= radius; 
      over = xy - radius;

      correction.x += ( over * radial.x );
      correction.y += ( over * radial.y );
    }

    //cerr << "Z value " << currentPos.at(atom).z << " ";
    //cerr << "Radial value " << radius << " ";
    //cerr << "Correction vector " << correction.x << " " << correction.y << " " << correction.z << endl;

    frodaAtom.at(atom).nMainMismatches++;
    frodaAtom.at(atom).mainMismatch.x += correction.x;
    frodaAtom.at(atom).mainMismatch.y += correction.y;
    frodaAtom.at(atom).mainMismatch.z += correction.z;
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: Updating all the atom positions. Moves all the atoms to reduce
//   their mismatches.
////////////////////////////////////////////////////////////////////////////////
void Froda::updateAtoms( MolFramework &structure ){
  
  register double pmet;
  Vector misvec;

  // Every atom has a total-mismatch vector from n interactions
  // hence a resultant average mismatch vector
  // Every atom gets its coordinates changed

    for( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ){
     
      if( !frodaAtom.at(siteNumber).isCore ) {

        if ( frodaAtom.at(siteNumber).nMainMismatches ==0 ) {
          //cerr << "NO MISMATCHES for atom " << siteNumber << endl;
          continue;
        }

        misvec = frodaAtom.at(siteNumber).mainMismatch;


        pmet = 1.0/frodaAtom.at(siteNumber).nMainMismatches;
        misvec.x *= pmet;
        misvec.y *= pmet;
        misvec.z *= pmet;
        
        //cerr << "Atom " << siteNumber << " has " << frodaAtom.at(siteNumber).nMainMismatches;
        //cerr << " mismatches and a final vector of ";
        //cerr << misvec.x << " ";
        //cerr << misvec.y << " ";
        //cerr << misvec.z << endl;        
                
        currentPos.at(siteNumber).x += misvec.x;
        currentPos.at(siteNumber).y += misvec.y;
        currentPos.at(siteNumber).z += misvec.z;         
      }
      else if( frodaAtom.at(siteNumber).isCore ) {
        currentPos.at(siteNumber) = initialPos.at(siteNumber);
      }
    }

  // pass over into structure information
  for( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ){
    structure.site_info[siteNumber].coords[0] = currentPos.at(siteNumber).x;
    structure.site_info[siteNumber].coords[1] = currentPos.at(siteNumber).y;
    structure.site_info[siteNumber].coords[2] = currentPos.at(siteNumber).z;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
Vector Froda::getGradient( double &weight, Vector &pq, Vector &pqdash, Vector &b ){
  
  Vector gradient;
  double myX = Xval(b);
  double myEPSX = epsx( pq, pqdash.x, b, myX);
  double myEPSY = epsy( pq, pqdash.y, b, myX);
  double myEPSZ = epsz( pq, pqdash.z, b, myX);

  double BXBZ = b.x*b.z;

  double myDEPSXDBX = pq.y*( (0.25*BXBZ/myX) + 0.5*b.y) + pq.z*( (-0.25*b.x*b.y/myX) + 0.5*b.z);
  double myDEPSYDBX = pq.x*( (-0.25*BXBZ/myX) + 0.5*b.y) + pq.y *(-b.x) + pq.z*(-myX + (0.25*b.x*b.x/myX)); 
  double myDEPSZDBX = pq.x*( (0.25*b.x*b.y/myX) + 0.5*b.z)  + pq.y*(myX - (0.25*b.x*b.x/myX)) + pq.z*(-b.x);

  double myDEPSXDBY = pq.x*(-b.y) + pq.y*( (0.25*b.y*b.z/myX) + 0.5*b.x) + pq.z*(myX - (0.25*b.y*b.y/myX));
  double myDEPSYDBY = pq.x*((-0.25*b.y*b.z/myX) + 0.5*b.x) + pq.z*((0.25*b.x*b.y/myX) + 0.5*b.z);
  double myDEPSZDBY = pq.x*(-myX + (0.25*b.y*b.y/myX)) + pq.y*((-0.25*b.x*b.y/myX) + 0.5*b.z) + pq.z*(-b.y);

  double myDEPSXDBZ = pq.x*(-b.z) + pq.z*( (-0.25*b.y*b.z/myX) + 0.5*b.x) + pq.y*(-myX + (0.25*b.z*b.z/myX));
  double myDEPSYDBZ = pq.x*(myX - (0.25*b.z*b.z/myX)) + pq.y*(-b.z) + pq.z*((0.25*BXBZ/myX) + 0.5*b.y);
  double myDEPSZDBZ = pq.x*((0.25*b.y*b.z/myX) + 0.5*b.x) + pq.y*((-0.25*BXBZ/myX) + 0.5*b.y); 

  gradient.x = weight * 2 * ( (myEPSX * myDEPSXDBX) + (myEPSY * myDEPSYDBX) + (myEPSZ * myDEPSZDBX) );
  gradient.y = weight * 2 * ( (myEPSX * myDEPSXDBY) + (myEPSY * myDEPSYDBY) + (myEPSZ * myDEPSZDBY) );
  gradient.z = weight * 2 * ( (myEPSX * myDEPSXDBZ) + (myEPSY * myDEPSYDBZ) + (myEPSZ * myDEPSZDBZ) );

  return gradient;
}


////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void Froda::resetWeights( ){
      for ( unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ) {
        if ( frodaAtom.at(siteNumber).hasWeight ) {
              unsigned int gSize = frodaAtom.at(siteNumber).ghostlist.size();
              for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
                ghost.at( frodaAtom.at(siteNumber).ghostlist.at(whichGhost) ).changedWeight = true;
              }
        }
        frodaAtom.at(siteNumber).weightIncreasing = false;
        frodaAtom.at(siteNumber).weightDecreasing = false;
        frodaAtom.at(siteNumber).hasWeight = false;
        frodaAtom.at(siteNumber).weight = 1.0;
      }

  return;
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void Froda::checkWeights( ){
        for ( unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ) {
          if ( frodaAtom.at(siteNumber).weightDecreasing ) {
            frodaAtom.at(siteNumber).weight -= 1.0;
            //warn the ghosts
            unsigned int gSize = frodaAtom.at(siteNumber).ghostlist.size();
            for ( unsigned int whichGhost = 0; whichGhost < gSize; whichGhost++ ) {
              ghost.at( frodaAtom.at(siteNumber).ghostlist.at(whichGhost) ).changedWeight = true;
            }

            if ( frodaAtom.at(siteNumber).weight < 1.01 ) {
              frodaAtom.at(siteNumber).weight = 1.0;
              frodaAtom.at(siteNumber).weightDecreasing = false;
              frodaAtom.at(siteNumber).hasWeight = false;
            }
          }
          else if ( frodaAtom.at(siteNumber).hasWeight and !frodaAtom.at(siteNumber).weightIncreasing ) {
            frodaAtom.at(siteNumber).weightDecreasing = true;
          }
        }

  return;
}


