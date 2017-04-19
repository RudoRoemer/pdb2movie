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
#include "hybrid_36_c.h"

#include <signal.h>
#include <cstring>

extern Parameters parameters;


////////////////////////////////////////////////////////////////////////////////
// Description: Initial setup of ghost templates
////////////////////////////////////////////////////////////////////////////////
//void Froda::initialiseGhosts( MolFramework &structure ){
void Froda::stripNonCovalentNeighbors( MolFramework &structure ){
  structure.stripNonCovalentNeighbors();
}


void Froda::getGhostMembershipFromFile( MolFramework &structure ){

  ifstream gfile;
  string filename = structure.path_name + "g.in";
  gfile.open( filename.c_str() );

  string linein;
  vector< string > splitLine;  

  unsigned int maxGhostID = 0;

  //first pass; get highest ghost id

  while ( !gfile.eof() ) {
    getline ( gfile, linein);
    splitLine = tokenize_string( linein );
    if (splitLine.size() < 2 ) continue;
    if ( splitLine[0] != "GHOST" ) continue; //ignore nonghost lines
   
    unsigned int thisGhostID = atoi( splitLine[1].c_str() ); 
    //cerr << "Found GID " << thisGhostID << endl;
    if ( thisGhostID > maxGhostID ) {
      maxGhostID = thisGhostID;
    } 
  }

  myNClusters = maxGhostID;
  structure.total_clusters = myNClusters;
  cerr << "myNClusters = " << myNClusters << endl; 

  ghost.resize(myNClusters + 1 );

  //rewind the file:
  //gfile.close();
  //gfile.open( filename.c_str() );
  gfile.clear();
  gfile.seekg( 0, ios::beg );

  //now actually get the data
  while ( !gfile.eof() ) {
    getline ( gfile, linein);
    splitLine = tokenize_string( linein );
    if (splitLine.size() < 4 ) continue; //ignore blank entries
    if ( splitLine[0] != "GHOST" ) continue; //ignore nonghost lines
   
    unsigned int thisGhostID = atoi( splitLine[1].c_str() ); 
    //cerr << "Reading Ghost ID " << thisGhostID << endl;
    unsigned int thisAtomID; 

    for ( unsigned int whichEntry = 3; whichEntry < splitLine.size(); whichEntry++ ) {
      thisAtomID = atoi( splitLine[whichEntry].c_str() );
      thisAtomID = structure.orig_2_FIRST_atom_number[thisAtomID];

      ghost[thisGhostID].atoms.push_back( thisAtomID );
      frodaAtom.at(thisAtomID).ghostlist.push_back( thisGhostID );

    }
  }

  gfile.close();

  for ( unsigned int whichGhost = 1; whichGhost <= myNClusters; whichGhost++ ) {
    if ( ghost[whichGhost].atoms.size() == 0 ) {
      cerr << "Warning: ghost " << whichGhost << " has no members!" << endl;
    }
  }
  bool badAtom = false;

  for ( unsigned int whichAtom = 1; whichAtom < frodaAtom.size(); whichAtom++ ) {
    if ( frodaAtom[whichAtom].ghostlist.size() == 0 ) {
      cerr << "ERROR: atom " << whichAtom << " has no ghosts!" << endl;
      badAtom = true;
    }
  }


  if ( badAtom ) exit(1);
  return;
}

void Froda:: getRCDFromGhosts( MolFramework &structure ) {
  //we've read the ghost files
  //now back-build the RCD
  //this should run after we've stripped the non-covalent interactions
  //so loop over neighbor_lists

  for ( unsigned int whichAtom = 1; whichAtom < frodaAtom.size(); whichAtom++ ) {
    unsigned int thisGhost;
    unsigned int otherGhost;
    unsigned int otherAtom;
    unsigned int thisRC = 0;

    for ( unsigned int whichG = 0; whichG < frodaAtom[whichAtom].ghostlist.size(); whichG++ ) {
      thisGhost = frodaAtom[whichAtom].ghostlist[whichG];
      bool matchesAll = true; //check if this ghost is shared with all neighbors

      for ( unsigned int whichN = 0; whichN < structure.site_info[whichAtom].neighbor_list.size(); whichN++ ) {
        otherAtom = structure.site_info[whichAtom].neighbor_list[whichN];
        bool matchesOne = false; //true if my neighbor shares thisGhost

        for ( unsigned int whichOG = 0; whichOG < frodaAtom[otherAtom].ghostlist.size(); whichOG++ ) {
          otherGhost = frodaAtom[otherAtom].ghostlist[whichOG];
          if ( otherGhost == thisGhost ) { matchesOne = true; break; }
        }

        if ( !matchesOne ) {
          //thisGhost is not shared with all neighbors; so not RC
          matchesAll = false;
          break;
        }

      }   

      if ( matchesAll ) {
        thisRC = thisGhost; //found a RC label
        break;
      } 
    }
    if ( thisRC == 0 ) {
      cerr << "WARNING: could not reconstruct RC label for atom " << whichAtom << endl;
    }
    else {
      //cerr << "Atom " << whichAtom << " has RC label " << thisRC << endl;
    }

    structure.site_info[whichAtom].rigid_label = thisRC;
    structure.site_info[whichAtom].coll_mode_label = 0;
  }

  structure.total_clusters_larger_than_one = 0;
  structure.RC_atom_list.resize( structure.total_clusters + 1);
  for ( unsigned int whichAtom = 1; whichAtom < frodaAtom.size(); whichAtom++ ) {
    structure.RC_atom_list.at( structure.site_info[whichAtom].rigid_label ).push_back( whichAtom );
  }
  for ( unsigned int whichRC = 1; whichRC < structure.RC_atom_list.size(); whichRC ++ ) {
    if ( structure.RC_atom_list.at(whichRC).size() < 2 ) continue;
    structure.total_clusters_larger_than_one++;
  }

  return;
}

//void Froda::initialiseGhosts( MolFramework &structure ){
void Froda::getGhostMembershipFromRCD( MolFramework &structure ){
 
  /* Ghost body n includes rigid cluster n */
  /* and the neighbours bonded to it */
  /* so ghosts overlap, giving angular constraints */
  /* ghost_atoms recored ths atoms that are in the nth ghost body */
  /* ghostlist records the ghost bodies that the nth atom is in */

  //IMPORTANT NOTE
  //when mobileRC1
  //cluster 1 is empty
  //and ghost body n is rigid cluster (n-1)
  //beware!
  if( parameters.verbose >= 2 ){
    cout << "  Setting up the ghost bodies..." << endl;
  }

  //now renumber the clusters
  //and quote the real number of clusters
  myNClusters = 0;
  if (mobileRC1) myNClusters = 1;

  if( parameters.verbose >= 2 ){
    cout << "  Eliminating rigid clusters of size 1..." << endl;
  }

  for ( unsigned int rcAtomNumber =1; rcAtomNumber < structure.RC_atom_list.size(); rcAtomNumber++){
    if ( structure.RC_atom_list[rcAtomNumber].size() ==1){
      int whichatom = structure.RC_atom_list[rcAtomNumber].at(0);
      if (structure.site_info[whichatom].neighbor_list.size() ==1){
        int otheratom = structure.site_info[whichatom].neighbor_list.at(0);
        int othercluster = structure.site_info[otheratom].rigid_label;
        structure.RC_atom_list[rcAtomNumber].clear();
        structure.RC_atom_list[othercluster].push_back(whichatom);
        structure.site_info[whichatom].rigid_label = othercluster;
      }
      else{
        myNClusters +=1;
      }
    } else{
      myNClusters +=1;
    }
  } 

  // now we've eliminated fictitious single clusters
  if( parameters.verbose >= 2 ){
    cout << "  New total number of rigid clusters: " << myNClusters << endl;
  }

  ghost.resize(myNClusters + 1 );

  //build the real ghosts now; 
  int whichghost = 0;
  if (mobileRC1) whichghost = 1;

  for ( unsigned int rcAtomNumber =1; rcAtomNumber < structure.RC_atom_list.size(); rcAtomNumber++){
    if ( structure.RC_atom_list.at(rcAtomNumber).size() == 0){
      continue;
    } else{
      whichghost +=1;
      //cerr << "Now doing ghost: " << whichghost << endl;
      for (unsigned int secondRcAtomNumber=0; secondRcAtomNumber < structure.RC_atom_list.at(rcAtomNumber).size(); secondRcAtomNumber++){
        unsigned int whichatom = structure.RC_atom_list.at(rcAtomNumber).at(secondRcAtomNumber);
        ghost.at(whichghost).atoms.push_back(whichatom);
        frodaAtom.at(whichatom).ghostlist.push_back(whichghost);
        for (unsigned int neighborNumber=0; neighborNumber < structure.site_info[whichatom].neighbor_list.size(); neighborNumber++){
          int otheratom = structure.site_info[whichatom].neighbor_list.at(neighborNumber);
          if (structure.site_info[otheratom].rigid_label != rcAtomNumber){ // FIXME - have to change rigid_label to unsigned int (difficult since it may be assigned -1 for unknown; check carefuly that it isn't)
            ghost.at(whichghost).atoms.push_back(otheratom);
            frodaAtom.at(otheratom).ghostlist.push_back(whichghost);
          }
        }   
      } 
    }
  }

  if (debug_rigid){
    for( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ){
      cout << "ATOM " << siteNumber << " IS IN: ";
      for (unsigned int ghostlistNumber=0; ghostlistNumber < frodaAtom.at(siteNumber).ghostlist.size();ghostlistNumber++){
        cout << frodaAtom.at(siteNumber).ghostlist.at(ghostlistNumber) << " ";
      }
      cout << endl;
    }
  }

  if ( parameters.outputGhosts ) {
    outputGhosts( structure );
  }
 
  return;
}

//void Froda::initialiseGhosts( MolFramework &structure ){
void Froda::buildGhostGeometry( MolFramework &structure ){
 
  // put the ghost positions into the initial ghost body
  if( parameters.verbose >= 2 ){
    cout << "  Putting ghost positions into bodies..." << endl;
  }

  for (unsigned int clusterNumber=1; clusterNumber<=myNClusters;clusterNumber++) {
    ghost.at(clusterNumber).runningRotor = NULL_VEC;
      ghost[clusterNumber].bondVector.clear();
      ghost[clusterNumber].pos.clear();
      ghost[clusterNumber].posAtom.clear();
      ghost[clusterNumber].posCentral = NULL_VEC;
      ghost[clusterNumber].radMS = 0;
      ghost[clusterNumber].radRMS = 0;
      ghost[clusterNumber].posCentral = NULL_VEC;
      ghost[clusterNumber].oldPosCentral = NULL_VEC;
      ghost[clusterNumber].olderPosCentral = NULL_VEC;
      ghost[clusterNumber].oldRotor = NULL_VEC;
      ghost[clusterNumber].olderRotor = NULL_VEC;
  }


  for( unsigned int clusterNumber = 1; clusterNumber <= myNClusters; clusterNumber++ ){
    ghost.at(clusterNumber).bodysize = ghost.at(clusterNumber).atoms.size();
    
    for( unsigned int ghostAtomNumber = 0; ghostAtomNumber < ghost[clusterNumber].bodysize; ghostAtomNumber++ ){ // this line picks up if entry 1 has no size...
      if (clusterNumber==1) frodaAtom.at( ghost[clusterNumber].atoms.at(ghostAtomNumber)).isCore =true; // tag core atoms
      tempvec = initialPos.at(ghost[clusterNumber].atoms.at(ghostAtomNumber));
      ghost[clusterNumber].bondVector.push_back(tempvec);
      ghost[clusterNumber].pos.push_back(tempvec);
      ghost[clusterNumber].posAtom.push_back(tempvec);
    }
    if (ghost[clusterNumber].bodysize == 1) cerr << "WARNING: cluster " << clusterNumber << " is a single-site ghost!" << endl;
  }


  //setting all atom weights to 1
  for (unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++){
    frodaAtom.at(siteNumber).weight = 1;
  }

  // for each ghost body, give it a central position, and convert its positions into bonds (relative vectors)
  // posCentral is the average of the positions of the members 
    
  for (unsigned int clusterNumber=1; clusterNumber<= myNClusters; clusterNumber++) {
    ghost[clusterNumber].posCentral = NULL_VEC;
    ghost[clusterNumber].runningRotor = NULL_VEC;
    ghost[clusterNumber].posCentral = NULL_VEC;
    ghost[clusterNumber].radMS = 0;
    ghost[clusterNumber].radRMS = 0;
  }

  for (unsigned int clusterNumber=1; clusterNumber<= myNClusters; clusterNumber++) {
    ghost.at( clusterNumber ).bodyweight = 0;
    for( unsigned int ghostAtomNumber = 0; ghostAtomNumber < ghost[clusterNumber].bodysize; ghostAtomNumber++ ){
      double aweight = frodaAtom.at( ghost[clusterNumber].atoms.at(ghostAtomNumber) ).weight;
      ghost.at(clusterNumber).posCentral.x += aweight * ghost[clusterNumber].bondVector.at(ghostAtomNumber).x;
      ghost.at(clusterNumber).posCentral.y += aweight * ghost[clusterNumber].bondVector.at(ghostAtomNumber).y;
      ghost.at(clusterNumber).posCentral.z += aweight * ghost[clusterNumber].bondVector.at(ghostAtomNumber).z;
      ghost.at( clusterNumber ).bodyweight += aweight;
    }
    ghost.at(clusterNumber).posCentral.x /= ghost.at( clusterNumber ).bodyweight;
    ghost.at(clusterNumber).posCentral.y /= ghost.at( clusterNumber ).bodyweight;
    ghost.at(clusterNumber).posCentral.z /= ghost.at( clusterNumber ).bodyweight;

    ghost.at(clusterNumber).oldPosCentral.x = ghost.at(clusterNumber).posCentral.x;
    ghost.at(clusterNumber).oldPosCentral.y = ghost.at(clusterNumber).posCentral.y;
    ghost.at(clusterNumber).oldPosCentral.z = ghost.at(clusterNumber).posCentral.z;
    ghost.at(clusterNumber).olderPosCentral.x = ghost.at(clusterNumber).posCentral.x;
    ghost.at(clusterNumber).olderPosCentral.y = ghost.at(clusterNumber).posCentral.y;
    ghost.at(clusterNumber).olderPosCentral.z = ghost.at(clusterNumber).posCentral.z;
    
    // get bonds and calculate radii
    
    for( unsigned int ghostAtomNumber = 0; ghostAtomNumber < ghost[clusterNumber].bodysize; ghostAtomNumber++ ){
      ghost[clusterNumber].bondVector.at(ghostAtomNumber).x -= ghost.at(clusterNumber).posCentral.x;
      ghost[clusterNumber].bondVector.at(ghostAtomNumber).y -= ghost.at(clusterNumber).posCentral.y;
      ghost[clusterNumber].bondVector.at(ghostAtomNumber).z -= ghost.at(clusterNumber).posCentral.z;
      if (ghost[clusterNumber].bodysize == 1){ ghost.at(clusterNumber).radMS = -1000.0; } else{
        ghost.at(clusterNumber).radMS += frodaAtom.at( ghost[clusterNumber].atoms.at(ghostAtomNumber) ).weight * dotProduct(ghost[clusterNumber].bondVector.at(ghostAtomNumber), ghost[clusterNumber].bondVector.at(ghostAtomNumber));
      }
    }
    ghost.at(clusterNumber).radMS /= ghost.at( clusterNumber ).bodyweight;
    if ( ghost.at(clusterNumber).radMS < 0.0 ){ radius = 0.0; ghost.at(clusterNumber).radMS = 0.0; } else{
      temp = ghost.at(clusterNumber).radMS;
      radius = sqrt(temp);
    }
    ghost.at(clusterNumber).radRMS = radius;
    
  }        
}


////////////////////////////////////////////////////////////////////////////////
// Description: Make a list of "active" hydrophobic tethers
// These are tethers between different rigid clusters
////////////////////////////////////////////////////////////////////////////////
void Froda::initialisePhobic( MolFramework &structure ){ 
  
  if( parameters.verbose >= 2 ){
    cout << "  Setting up the hydrophobic interactions..." << endl;
  }
  
  double l1;
  
  for (unsigned int hydrophobicTetherNumber =0; hydrophobicTetherNumber < structure.hydrophobic_tethers.size();hydrophobicTetherNumber++){
    if (structure.hydrophobic_tethers.at(hydrophobicTetherNumber).energy < energy_cutoff){ //it's a valid member
      
      dummyPhobe.atomA = structure.hydrophobic_tethers.at(hydrophobicTetherNumber).site_1;
      dummyPhobe.atomB = structure.hydrophobic_tethers.at(hydrophobicTetherNumber).site_2;
      
      // don't count mutually rigid atoms
      if( !parameters.inputGhosts ) {
        if(structure.site_info[dummyPhobe.atomA].rigid_label == structure.site_info[dummyPhobe.atomB].rigid_label) continue;
      }
      else { 
        if ( bonded( frodaAtom.at(dummyPhobe.atomA).ghostlist,
                     frodaAtom.at(dummyPhobe.atomB).ghostlist )
                         ) continue;
      }
      
      n_myPhobes++;
      l1= structure.site_info[dummyPhobe.atomA].vdw_radius + structure.site_info[dummyPhobe.atomB].vdw_radius;

      dummyPhobe.AB = l1;
 
      Phobes.push_back(dummyPhobe);

    }
  }

  //cout << " Found " << Phobes.size() << " active hydrophobic tethers." << endl;
}
////////////////////////////////////////////////////////////////////////////////
// Make Backbone list for per-residue RMSD calculations
void Froda::makeBackboneList( MolFramework &structure ){

  for (unsigned int siteNumber =1; siteNumber<= structure.total_sites; siteNumber++){
    //  cout <<  structure.isMainchain(siteNumber)<<" Name "<<structure.site_info[siteNumber].atom_name<< endl;

    if( structure.isMainchain(siteNumber) == 101){ // For DNA, RNA structures

      alphaCarbon.push_back(siteNumber); // for per-residue RMSD, in case of DNA, RNA phosphorus is considered as backbone atom
      frodaAtom.at(siteNumber).isAlpha = true;
      if( parameters.verbose == 3 ){
	cerr << "  Noticing and labelling site " << siteNumber << " as Phosphorus." << endl;
      }
    }
    if( structure.isMainchain(siteNumber) == 2){      

      alphaCarbon.push_back(siteNumber); // for per-residue RMSD
      frodaAtom.at(siteNumber).isAlpha = true;
      if( parameters.verbose == 3 ){
	cerr << "  Noticing and labelling site " << siteNumber << " as C-alpha." << endl;
      }
      
    }
  }
   return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:  create the list of Ramachandran constraints.
////////////////////////////////////////////////////////////////////////////////
void Froda::makeRamaList( MolFramework &structure ){
 unsigned int anatom, anotheratom;
 Vector e1,e2,e3,e4;
 double pmet; // FIXME - unused - phi,psi
 double didot, diangle;
 bool posangle=0;

 int ntermini=0, ctermini=0, locked=0, calphas=0, mobiles=0, badstarts=0;
 int bin_phi, bin_psi; // will address the arrays
 // find C-alphas, put into dummy rama

 //Find N and C, put into dummy rama
 //Find Ci_1 and Ni1, put into dummy rama, or get terminus
 // check if phi or psi are locked
 // if both locked, don't bother
 //JULY 06 actually do include if locked, just don't recheck
 // if free, calculate phi0 and psi0, check status and set badStart flag
 // for now, don't bother with bad starts
 // if all the test are passed, push dummy rama into ramaList.

 for (unsigned int siteNumber =1; siteNumber<= structure.total_sites; siteNumber++){

   dummyRama.Calpha=dummyRama.Ni=dummyRama.Ci=dummyRama.Ci_1=dummyRama.Ni1=0;
   dummyRama.freePhi = dummyRama.freePsi=true;
   dummyRama.badStart=0;
   posangle=0;
   
   if( structure.isMainchain(siteNumber) == 2){
    calphas++;

    dummyRama.Calpha = siteNumber;

    for (unsigned int neighborNumber=0;neighborNumber<structure.site_info[siteNumber].neighbor_list.size(); neighborNumber++){
      anatom = structure.site_info[siteNumber].neighbor_list.at(neighborNumber);
      if  ( structure.isMainchain(anatom) == 1){
        dummyRama.Ni = anatom;
      }
      if  ( structure.isMainchain(anatom) == 3){
        dummyRama.Ci = anatom;
      }

    }

    if (!dummyRama.Ci && !dummyRama.Ni) continue; // not right bond

    anatom = dummyRama.Ni;
    for (unsigned int neighborNumber=0; neighborNumber< structure.site_info[anatom].neighbor_list.size();neighborNumber++){
      anotheratom = structure.site_info[anatom].neighbor_list.at(neighborNumber);
      if (anotheratom==siteNumber) continue; // don't loop back
      if  ( structure.isMainchain(anotheratom) == 3){
        dummyRama.Ci_1 = anotheratom;
      }

    }

    if (!dummyRama.Ci_1){ ntermini++; continue;} // no Ci-1, it's the N terminus

    anatom = dummyRama.Ci;
    for (unsigned int neighborNumber=0; neighborNumber< structure.site_info[anatom].neighbor_list.size();neighborNumber++){
      anotheratom = structure.site_info[anatom].neighbor_list.at(neighborNumber);
      if (anotheratom==siteNumber) continue; // don't loop back
      if  ( structure.isMainchain(anotheratom) == 1){
        dummyRama.Ni1 = anotheratom;
      }

    }

    if (!dummyRama.Ni1) {ctermini++; continue;} // no Ni+1, it's the C terminus

    // if we're here, we have all five atoms

    anatom = dummyRama.Ni;
    if (structure.site_info[siteNumber].rigid_label == structure.site_info[anatom].rigid_label) dummyRama.freePhi=0;
    anatom = dummyRama.Ci;
    if (structure.site_info[siteNumber].rigid_label == structure.site_info[anatom].rigid_label) dummyRama.freePsi=0;

    if (!dummyRama.freePhi && !dummyRama.freePsi) {
      locked++;
    } 
    else {
      mobiles++;
    }

    // evaluate phi0 psi0 and set badStart flag.
    // first phi:

    dummyvec.x = initialPos.at(dummyRama.Calpha ).x - initialPos.at(dummyRama.Ni ).x;
    dummyvec.y = initialPos.at(dummyRama.Calpha ).y - initialPos.at(dummyRama.Ni ).y;
    dummyvec.z = initialPos.at(dummyRama.Calpha ).z - initialPos.at(dummyRama.Ni ).z;
    dot = sqrt(dotProduct(dummyvec, dummyvec));
    pmet = 1.0/dot;
    e1.x=dummyvec.x*pmet;
    e1.y=dummyvec.y*pmet;
    e1.z=dummyvec.z*pmet; //e1 now unit vector along Ni to Calpha

    dummyvec.x = initialPos.at(dummyRama.Ci_1 ).x - initialPos.at(dummyRama.Ni ).x;
    dummyvec.y = initialPos.at(dummyRama.Ci_1 ).y - initialPos.at(dummyRama.Ni ).y;
    dummyvec.z = initialPos.at(dummyRama.Ci_1 ).z - initialPos.at(dummyRama.Ni ).z;
    //dummyvec is now the Ni to Ci-1 bond
    dot =dotProduct(e1,dummyvec);
    dummyvec.x -= dot*e1.x;
    dummyvec.y -= dot*e1.y;
    dummyvec.z -= dot*e1.z;
    dot = sqrt(dotProduct(dummyvec, dummyvec));
    pmet = 1.0/dot;
    e2.x=dummyvec.x*pmet;
    e2.y=dummyvec.y*pmet;
    e2.z=dummyvec.z*pmet; //e2 now unit projected vector for Ni to Ci-1 
    
    dummyvec.x = initialPos.at(dummyRama.Ci ).x - initialPos.at(dummyRama.Calpha ).x;
    dummyvec.y = initialPos.at(dummyRama.Ci ).y - initialPos.at(dummyRama.Calpha ).y;
    dummyvec.z = initialPos.at(dummyRama.Ci ).z - initialPos.at(dummyRama.Calpha ).z;
    //dummyvec is now the Calpha to Ci bond
    dot =dotProduct(e1,dummyvec);
    dummyvec.x -= dot*e1.x;
    dummyvec.y -= dot*e1.y;
    dummyvec.z -= dot*e1.z;
    dot = sqrt(dotProduct(dummyvec, dummyvec));
    pmet = 1.0/dot;
    e3.x=dummyvec.x*pmet;
    e3.y=dummyvec.y*pmet;
    e3.z=dummyvec.z*pmet; //e3 now unit projected vector for Calpha to Ci 
    
    e4 = crossProduct(e1,e2);
    didot = dotProduct(e3,e4);
    if(didot>0.0) posangle=1; // is in range positive
    didot = dotProduct(e2,e3);
    if(didot<-1.0) didot=-0.99999999;
    if(didot>1.0) didot=0.99999999;
    diangle = RAD_TO_DEG*acos(didot); // 0 to 180 degrees
    if (!posangle) diangle = -diangle;
     
    dummyRama.phi0 = diangle; // finally! 

    //now psi:
    posangle =0;

    dummyvec.x = initialPos.at(dummyRama.Ci ).x - initialPos.at(dummyRama.Calpha ).x;
    dummyvec.y = initialPos.at(dummyRama.Ci ).y - initialPos.at(dummyRama.Calpha ).y;
    dummyvec.z = initialPos.at(dummyRama.Ci ).z - initialPos.at(dummyRama.Calpha ).z;
    dot = sqrt(dotProduct(dummyvec, dummyvec));
    pmet = 1.0/dot;
    e1.x=dummyvec.x*pmet;
    e1.y=dummyvec.y*pmet;
    e1.z=dummyvec.z*pmet; //e1 now unit vector along Calpha to Ci
    
    dummyvec.x = initialPos.at(dummyRama.Ni ).x - initialPos.at(dummyRama.Calpha ).x;
    dummyvec.y = initialPos.at(dummyRama.Ni ).y - initialPos.at(dummyRama.Calpha ).y;
    dummyvec.z = initialPos.at(dummyRama.Ni ).z - initialPos.at(dummyRama.Calpha ).z;
    //dummyvec is now the Ni to Ci-1 bond
    dot =dotProduct(e1,dummyvec);
    dummyvec.x -= dot*e1.x;
    dummyvec.y -= dot*e1.y;
    dummyvec.z -= dot*e1.z;
    dot = sqrt(dotProduct(dummyvec, dummyvec));
    pmet = 1.0/dot;
    e2.x=dummyvec.x*pmet;
    e2.y=dummyvec.y*pmet;
    e2.z=dummyvec.z*pmet; //e2 now unit projected vector for Calpha to Ni 

    dummyvec.x = initialPos.at(dummyRama.Ni1 ).x - initialPos.at(dummyRama.Ci ).x;
    dummyvec.y = initialPos.at(dummyRama.Ni1 ).y - initialPos.at(dummyRama.Ci ).y;
    dummyvec.z = initialPos.at(dummyRama.Ni1 ).z - initialPos.at(dummyRama.Ci ).z;
    //dummyvec is now the  Ci to Ni+1 bond
    dot =dotProduct(e1,dummyvec);
    dummyvec.x -= dot*e1.x;
    dummyvec.y -= dot*e1.y;
    dummyvec.z -= dot*e1.z;
    dot = sqrt(dotProduct(dummyvec, dummyvec));
    pmet = 1.0/dot;
    e3.x=dummyvec.x*pmet;
    e3.y=dummyvec.y*pmet;
    e3.z=dummyvec.z*pmet; //e3 now unit projected vector for Ci to Ni+1 

    e4 = crossProduct(e1,e2);
    didot = dotProduct(e3,e4);
    if(didot>0.0) posangle=1; // is in range positive
    didot = dotProduct(e2,e3);
    if(didot<-1.0) didot=-0.99999999;
    if(didot>1.0) didot=0.99999999;
    diangle = RAD_TO_DEG*acos(didot); // 0 to 180 degrees
    if (!posangle) diangle = -diangle;
     
    dummyRama.psi0 = diangle; // finally! 

    //when we have arrays, do siteNumber check for bad start

    if (structure.site_info[dummyRama.Calpha].residue_name == "PRO"){
      dummyRama.isPro = true;
      bin_phi = (int) ((dummyRama.phi0+180.0)*0.125);
      bin_psi = (int) ((dummyRama.psi0+180.0)*0.125);     

      dummyRama.mapValue = prolineRamaPlot.at(bin_psi).at(bin_phi);
      if (dummyRama.mapValue > 2 ) dummyRama.badStart=true;
    } else if (structure.site_info[dummyRama.Calpha].residue_name == "GLY"){
      dummyRama.isGly = true;
      bin_phi = (int) ((dummyRama.phi0+180.0)*0.125);
      bin_psi = (int) ((dummyRama.psi0+180.0)*0.125);

      dummyRama.mapValue = glycineRamaPlot.at(bin_psi).at(bin_phi);
      if (dummyRama.mapValue > 2 ) dummyRama.badStart=true;
    } else{
      dummyRama.isGly = false;
      dummyRama.isPro = false;
      bin_phi = (int) ((dummyRama.phi0+180.0)*0.1);
      bin_psi = (int) ((dummyRama.psi0+180.0)*0.1);
           
      dummyRama.mapValue = generalRamaPlot.at(bin_psi).at(bin_phi);
      if (dummyRama.mapValue > 2 ) dummyRama.badStart=true;
    }

    if (dummyRama.badStart){
      badstarts++; 
    } 
    ramaList.push_back(dummyRama);
    
  } 


 }

 if( parameters.verbose >= 2 ){
   cout << "  Ramachandran data:" << endl;
   cout << "  " << calphas << " alpha carbons" << endl;
   cout << "  " << ntermini << " N termini" << endl;
   cout << "  " << ctermini << " C termini" << endl;
   cout << "  " << locked << " locked Phi-Psi pairs" << endl;
   cout << "  " << badstarts << " bad starting points" << endl;
   cout << "  " << mobiles << " variable pairs" << endl;
   cout << "  Final size of Ramachandran list: " << ramaList.size() << endl;
 }

 return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:  Rama plot for general residue.
////////////////////////////////////////////////////////////////////////////////
void Froda::getGeneralRamaPlot(){

 // Read the plot and the gradients from the text files into arrays
 // read from generalRamaPlot.txt, dphi_general.txt, dpsi_general.txt
 // Each line runs in 36 bins from PHI=-180/70 to PHI=+170/80
 // The first line is for PSI = -180/70
 // There are 36 lines
 // generalRamaPlot.at(0).at(0) is for PSi = -180/70, PHI = -180/70

vector<string> in_tokens;
string in_line, myfile;
vector<int> in_values;
int mysize=36;

generalRamaPlot.resize(mysize);

for (int a=0;a<mysize;a++){
 generalRamaPlot.at(a).resize(mysize);
}
myfile = parameters.path+ "/lib/generalRamaPlot.txt";
ifstream infile( myfile.c_str() );

  if( !infile ) {
    cout << "Error opening file." << endl;
    exit(1);
  }

  for (int a=0;a<mysize;a++){
    getline(infile, in_line);
    in_tokens=tokenize_string(in_line);


    for (int p=0;p<mysize;p++){
      in_values.push_back(atoi(in_tokens.at(p).c_str()));
    }

    for (int p=0;p<mysize;p++){
      generalRamaPlot.at(a).at(p) = in_values.at(p);
    }

    in_tokens.clear();
    in_values.clear();
  }

 infile.close();

 return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:  Rama plot for proline.
////////////////////////////////////////////////////////////////////////////////
void Froda::getProlineRamaPlot(){
 
 // Read the plot and the gradients from the text files into arrays
 // read from prolineRamaPlot.txt, dphi_proline.txt, dpsi_proline.txt
 // Each line runs in 45 bins from PHI=-180/70 to PHI=+170/80
 // The first line is for PSI = +180/70
 // There are 45 lines
 // prolineRamaPlot.at(0).at(0) is for PSI = -180/70, PHI = -180/70

vector<string> in_tokens;
string in_line, myfile;
vector<int> in_values;
int mysize=45;
int placer;

prolineRamaPlot.resize(mysize);

for (int a=0;a<mysize;a++){
 prolineRamaPlot.at(a).resize(mysize);
}

myfile = parameters.path+ "/lib/prolineRamaPlot.txt";
ifstream infile( myfile.c_str() );

  if( !infile ) {
    cout << "Error opening file." << endl;
    exit(1);
  }

  for (int a=0;a<mysize;a++){
    placer = mysize -a -1;
    getline(infile, in_line);
    in_tokens=tokenize_string(in_line);


    for (int p=0;p<mysize;p++){
      in_values.push_back(atoi(in_tokens.at(p).c_str()));
    }

    for (int p=0;p<mysize;p++){
      prolineRamaPlot.at(placer).at(p) = in_values.at(p);
    }

    in_tokens.clear();
    in_values.clear();
  }

 infile.close();

 return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:  Rama plot for glycine.
////////////////////////////////////////////////////////////////////////////////
void Froda::getGlycineRamaPlot(){
 // Read the plot and the gradients from the text files into arrays
 // read from glycineRamaPlot.txt, dphi_glycine.txt, dpsi_glycine.txt
 // Each line runs in 45 bins from PHI=-180/70 to PHI=+170/80
 // The first line is for PSI = +180/70
 // There are 45 lines
 // prolineRamaPlot.at(0).at(0) is for PSI = -180/70, PHI = -180/70

vector<string> in_tokens;
string in_line, myfile;
vector<int> in_values;
int mysize=45;
int placer;

glycineRamaPlot.resize(mysize);

for (int a=0;a<mysize;a++){
 glycineRamaPlot.at(a).resize(mysize);
}


myfile = parameters.path+ "/lib/glycineRamaPlot.txt";
ifstream infile( myfile.c_str() );

  if( !infile ) {
    cout << "Error opening file." << endl;
    exit(1);
  }

  for (int a=0;a<mysize;a++){
    placer = mysize -a -1;
    getline(infile, in_line);
    in_tokens=tokenize_string(in_line);


    for (int p=0;p<mysize;p++){
      in_values.push_back(atoi(in_tokens.at(p).c_str()));
    }

    for (int p=0;p<mysize;p++){
      glycineRamaPlot.at(placer).at(p) = in_values.at(p);
    }

    in_tokens.clear();
    in_values.clear();
  }


 infile.close();

 return;
}


////////////////////////////////////////////////////////////////////////////////
// Description:  find neighbouring atoms in target structure
////////////////////////////////////////////////////////////////////////////////
void Froda::setUpLocalTargeting(MolFramework &structure, MolFramework &target ){

  cerr << "Setting up local targetting: building friend lists." << endl;

  //the logic is:
  //atom I is targetted to atom J in target -- myTarget.at(I) = J
  //atoms Kn are near to atom J in target
  //atoms Ln are such that myTarget.at(Ln) = Kn
  //the friends of I are Ln

  Vector dr;
  double dist;

  double range = pow(4.5, 2); // radius squared for checker

  for (unsigned int firstSiteNumber = 1; firstSiteNumber <= structure.total_sites; firstSiteNumber++){

    unsigned int atomj = frodaAtom.at(firstSiteNumber).myTarget;
    if (atomj == 0 ) continue ; // skip absentees

    //cerr << "Structure atom " << firstSiteNumber << ", target atom " << atomj << endl;

    for (unsigned int secondSiteNumber = 1; secondSiteNumber <= target.total_sites; secondSiteNumber++){
      if ( secondSiteNumber == atomj) continue ; // skip self
      dr.x = target.site_info[secondSiteNumber].coords[0] - target.site_info[atomj].coords[0];
      dr.y = target.site_info[secondSiteNumber].coords[1] - target.site_info[atomj].coords[1];
      dr.z = target.site_info[secondSiteNumber].coords[2] - target.site_info[atomj].coords[2];

      dist = dotProduct(dr,dr);

      if (dist < range){
        //cerr << "Target friend " << secondSiteNumber << ", structure friend ";

        //find the structure atom corresponding
        if (!fancy_targetting){
          //the friend is secondSiteNumber
          //cerr << secondSiteNumber;
          frodaAtom.at(firstSiteNumber).friends.push_back(secondSiteNumber);
        } 
        else{ //more difficult
          for (unsigned int thirdSiteNumber = 1; thirdSiteNumber <= structure.total_sites; thirdSiteNumber++){
            if ( frodaAtom.at(thirdSiteNumber).myTarget == secondSiteNumber){
              //the friend is thirdSiteNumber
              //cerr << thirdSiteNumber;
              frodaAtom.at(firstSiteNumber).friends.push_back(thirdSiteNumber);
              break;
            }
          }
        }
        //cerr << endl;
      }

    }

  } //we now have all of the friends

  //position and offset calculation

  for ( unsigned int siteNumber =1; siteNumber <= structure.total_sites; siteNumber++){
    if ( frodaAtom.at(siteNumber).myTarget == 0 ) continue;
    if (frodaAtom.at(siteNumber).friends.size() == 0) continue;
    cerr << "Atom " << siteNumber << " has " << frodaAtom.at(siteNumber).friends.size() << " friends." << endl;

    Vector centralpoint = NULL_VEC;

    for (unsigned int friendNumber =0; friendNumber < frodaAtom.at(siteNumber).friends.size(); friendNumber++){
      int anatom =  frodaAtom.at( frodaAtom.at(siteNumber).friends.at(friendNumber)).myTarget;
      centralpoint.x += target.site_info[anatom].coords[0];
      centralpoint.y += target.site_info[anatom].coords[1];
      centralpoint.z += target.site_info[anatom].coords[2];
    }

    centralpoint.x /= frodaAtom.at(siteNumber).friends.size();
    centralpoint.y /= frodaAtom.at(siteNumber).friends.size();
    centralpoint.z /= frodaAtom.at(siteNumber).friends.size();

    dr.x = centralpoint.x - target.site_info[ frodaAtom.at(siteNumber).myTarget ].coords[0];
    dr.x = centralpoint.y - target.site_info[ frodaAtom.at(siteNumber).myTarget ].coords[1];
    dr.x = centralpoint.z - target.site_info[ frodaAtom.at(siteNumber).myTarget ].coords[2];

    dist = sqrt(dotProduct(dr,dr));

    frodaAtom.at(siteNumber).offset = dist;
    cerr << "Initial offset from friend sphere: " << frodaAtom.at(siteNumber).offset << endl;

    frodaAtom.at(siteNumber).offset += 0.5*lstep; //pad it for tolerance
    frodaAtom.at(siteNumber).offset *= frodaAtom.at(siteNumber).offset; // square it for later savings in compare

  }

  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: array initialisations for atom information
////////////////////////////////////////////////////////////////////////////////
void Froda::initialiseMainArrays( MolFramework &structure, MolFramework &restart ){

  frodaAtom.resize( structure.total_sites + 1 );

  currentPos.resize( structure.total_sites +1 );
  initialPos.resize( structure.total_sites +1 );
  cachedPos.resize( structure.total_sites +1 );

  overlapDegrees.resize(7);
 
  nMainchainAtoms=0; nSidechainAtoms=0; nMobileMainchainAtoms=0; nMobileSidechainAtoms=0;
  for( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ){

    frodaAtom.at(siteNumber).weight = 0;
    frodaAtom.at(siteNumber).isCore = false;
    frodaAtom.at(siteNumber).isMobile = true;
    tempvec = Vector(structure.site_info[siteNumber].coords[0],structure.site_info[siteNumber].coords[1], structure.site_info[siteNumber].coords[2]);
    initialPos.at(siteNumber) = tempvec;
    if (isRestart){
      tempvec = Vector(restart.site_info[siteNumber].coords[0],restart.site_info[siteNumber].coords[1], restart.site_info[siteNumber].coords[2]);
      currentPos.at(siteNumber) = tempvec;
    } else{
      currentPos.at(siteNumber) = tempvec;
    }

    frodaAtom.at(siteNumber).mainMismatch = NULL_VEC;
    //frodaAtom.at(siteNumber).addMismatch = NULL_VEC;
    frodaAtom.at(siteNumber).nMainMismatches = 0;
    //frodaAtom.at(siteNumber).nAddMismatches = 0;

    frodaAtom.at(siteNumber).myTarget = 0;

  }
  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: labels atoms as mainchain or sidechain
////////////////////////////////////////////////////////////////////////////////
void Froda::defineAtomsAsMainOrSide( MolFramework &structure){

  if( parameters.verbose >= 2 ){
    cout << "  Defining main and side chain atoms." << endl;
  }
  nMobileAtoms = 0;
  
  for (unsigned int siteNumber =1; siteNumber <= structure.total_sites; siteNumber++){
    if (structure.isMainchain(siteNumber)){
      frodaAtom.at(siteNumber).isMain = true; 
      nMainchainAtoms++;
      if( !frodaAtom.at(siteNumber).isCore ) {nMobileMainchainAtoms++; nMobileAtoms ++;} //core atoms not mobile
    } else{
      frodaAtom.at(siteNumber).isMain = false; 
      nSidechainAtoms++;
      if( !frodaAtom.at(siteNumber).isCore ) {nMobileSidechainAtoms++; nMobileAtoms++;} //core atoms not mobile
    }
  }

  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: labels clusters as mainchain or sidechain
////////////////////////////////////////////////////////////////////////////////
void Froda::defineClustersAsMainOrSide(){

  if( parameters.verbose ){
    cout << "Defining side and main clusters." << endl;
  }

  //define all clusters as main (includes backbone) or side 
  for (unsigned int clusterNumber=1; clusterNumber<=myNClusters; clusterNumber++){
    ghost[clusterNumber].isSidechain = 1; // true until proved otherwise
    
    for (unsigned int ghostAtomNumber=0; ghostAtomNumber< ghost[clusterNumber].atoms.size(); ghostAtomNumber++){
      if ( !ghost[clusterNumber].isSidechain ) break;
      if (  frodaAtom.at(ghost[clusterNumber].atoms.at(ghostAtomNumber)).isMain ) ghost[clusterNumber].isSidechain = 0;
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: array initialisation for RMSD by residue
////////////////////////////////////////////////////////////////////////////////
void Froda::setUpRMSDByResidueArrays(){

  if( parameters.verbose >= 2 ){
    cout << "  Setting up RMSD arrays." << endl;
  }

  //prepare resrmsd arrays for reporting
  //////////////////////////////////////////////////////////////////////
  currentResidueRMSD.resize( alphaCarbon.size() );
  runningResidueRMSD.resize( alphaCarbon.size() );
  for ( unsigned int carbon = 0; carbon < alphaCarbon.size(); carbon++){
    runningResidueRMSD.at(carbon) = 0.0;  //initialise to zero for later counting;
  }

  return;
  running_allmain_MSD = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// Description: labels polar atoms
// these are:
// N and O (negative)
// carbonyl C (positive)
// P and FE (positive)
// H on a polar atom
////////////////////////////////////////////////////////////////////////////////
void Froda::identifyPolarAtoms(MolFramework &structure){

  for ( unsigned int siteNumber = 0; siteNumber <= structure.total_sites; siteNumber++){
    frodaAtom.at(siteNumber).isPolar = false;
    frodaAtom.at(siteNumber).isHydrogen = false;
  }

    //find the charged atoms and populate the lists
    bool is_charged;
    int charge = 0;

    for (unsigned int siteNumber=1; siteNumber <= structure.total_sites; siteNumber++){
      is_charged = false;
      if (structure.site_info[siteNumber].element_name == "H "){
        frodaAtom.at(siteNumber).isHydrogen = true;
        for ( unsigned int neighborNumber = 0; neighborNumber < structure.site_info[siteNumber].neighbor_list.size(); neighborNumber++){
          if (structure.site_info[ structure.site_info[siteNumber].neighbor_list.at(neighborNumber) ].element_name == "N " ||
              structure.site_info[ structure.site_info[siteNumber].neighbor_list.at(neighborNumber) ].element_name == "O "){
            is_charged = true; charge = 1;
            frodaAtom.at(siteNumber).isPolar = true;
          }
        }
      }
      else if (structure.site_info[siteNumber].element_name == "N " || structure.site_info[siteNumber].element_name == "O "){
        is_charged = true; charge = -1;
        frodaAtom.at(siteNumber).isPolar = true;
      } else if (structure.site_info[siteNumber].element_name == "P " || structure.site_info[siteNumber].element_name == "FE"){
        is_charged = true; charge = 1;
        frodaAtom.at(siteNumber).isPolar = true;
      } else if ( structure.isMainchain(siteNumber) == 3){ //carbonyl backbone carbon
        is_charged = true; charge = 1;
        frodaAtom.at(siteNumber).isPolar = true;
      } else if ( structure.site_info[siteNumber].element_name == "C "){
        int n_charged_neighbors = 0;
        for ( unsigned int neighborNumber = 0; neighborNumber < structure.site_info[siteNumber].neighbor_list.size(); neighborNumber++){
          if (structure.site_info[ structure.site_info[siteNumber].neighbor_list.at(neighborNumber) ].element_name == "N " ||
              structure.site_info[ structure.site_info[siteNumber].neighbor_list.at(neighborNumber) ].element_name == "O "){
              n_charged_neighbors++;
          }
        }
        if (n_charged_neighbors > 1){ // siteNumber carbon bonded to two Os or Ns gets siteNumber positive charge to balance them
          is_charged = true; charge = 1;
          frodaAtom.at(siteNumber).isPolar = true;
        }
      }

      if (is_charged){

        chargedAtoms.push_back(siteNumber);
        frodaAtom.at(siteNumber).charge = charge;
      } else { frodaAtom.at(siteNumber).charge = 0;}

    }

    int n_neg =0, n_pos =0, n_null = 0;
    int n_polar_h = 0;

    for (unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++){
      if (frodaAtom.at(siteNumber).charge > 0 ) n_pos ++;
      if (frodaAtom.at(siteNumber).charge<0) n_neg ++;
      if (frodaAtom.at(siteNumber).charge==0) n_null ++;
      if ( frodaAtom.at(siteNumber).isPolar && frodaAtom.at(siteNumber).isHydrogen ) n_polar_h++;
    }
    
    if( parameters.verbose >= 2 ){
      cout << "  Found " << n_pos << " positively (+) charged atoms" << endl;
      cout << "  Found " << n_neg << " negatively (-) charged atoms" << endl;
      cout << "  Found " << n_null << " neutral atoms"  << endl;
      cerr << "  Found " << n_polar_h << " polar hydrogen atoms." << endl;
    }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: grab a list of all beta-carbon atoms from phobic groups
// For our purposes phobic groups are:
// ALA VAL LEU ILE MET PHE TRP TYR CYS
////////////////////////////////////////////////////////////////////////////////
void Froda::findPhobicBetaCarbons( MolFramework &structure ){
  phobicBetaCarbon.clear();
  for ( SiteID siteID = 1; siteID <= structure.total_sites; siteID ++ ) {
    if ( structure.site_info[siteID].atom_name == "CB" ) {
      //it's a beta carbon atom
      if ( structure.site_info[siteID].residue_name == "ALA" ||
           structure.site_info[siteID].residue_name == "VAL" ||
           structure.site_info[siteID].residue_name == "LEU" ||
           structure.site_info[siteID].residue_name == "ILE" ||
           structure.site_info[siteID].residue_name == "MET" ||
           structure.site_info[siteID].residue_name == "PHE" ||
           structure.site_info[siteID].residue_name == "TRP" ||
           structure.site_info[siteID].residue_name == "TYR" ||
           structure.site_info[siteID].residue_name == "CYS" ) {
        //it's a phobic residue
        phobicBetaCarbon.push_back(siteID);
      }
    }
  }

  cerr << "Found " << phobicBetaCarbon.size() << " phobic beta-carbons." << endl;
  cerr << "These are: " << endl;

  for ( unsigned int whichCB = 0; whichCB < phobicBetaCarbon.size(); whichCB++ ) {
    cerr << phobicBetaCarbon.at(whichCB);
    cerr << " " << structure.site_info[phobicBetaCarbon.at(whichCB)].residue_name << endl;
  }

  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: gets surface calculation parameters by atom type
//   this routine uses the isHydrogen label
//   it uses the neighbor counts to set atom type
//   parameters go into surfaceCk and surfaceRho
//   unrecognised atoms can get some defaults
////////////////////////////////////////////////////////////////////////////////
void Froda::obtainSurfaceParameters(MolFramework &structure){
  //atom type depends on element; n heavy neighbors and light
  string element;
  string element1;
  int nHeavy;
  int nLight;
  int nTotal; 

  struct AtomParamSet{
    string elem;
    string label;
    int nHeavy;
    int nLight;
    vector< double > surfaceCk;
    vector < double > HCk;
  };

  AtomParamSet C13, C22, C31, C21, C30;
  AtomParamSet N13, N22, N12, N21, N30, N20;
  AtomParamSet O11, O10, Ocarb10;
  AtomParamSet S11, S20;

  string filename = parameters.path + "/lib/surfaceParameters.txt"; 
  string linebuf;
  vector< string > words;
  ifstream surfaceParams;

  string elem;
  string label;
  vector< double > surfaceCk;
  vector< double > HCk;

  surfaceCk.resize(5);
  HCk.resize(5);

  //get all the params from a file in library

  surfaceParams.open( filename.c_str() );

  while ( !surfaceParams.eof() ) {

    getline( surfaceParams, linebuf );

    if (surfaceParams.eof() ) break;

    // Skip blank lines
    words = tokenize_string( linebuf );
    if ( words.size() == 0 ) continue;  

 
    if ( words[0] == "Atom" ) continue; //skip header line
    if ( words[0] == "C" || words[0] == "N" || words[0] == "O" || words[0] == "S" ) {
      elem = words[0];
      label = words[1];
      nHeavy = atoi( words[2].c_str() ); //heavy neighbors
      nLight = atoi( words[3].c_str() ); //hydrogens
      for ( int entry = 0; entry < 5; entry++ ) {
        surfaceCk[entry] = atof( words[entry + 5].c_str() );
      }
      //if there are hydrogens, read their params as well
      if ( nLight > 0 ) {
        getline( surfaceParams, linebuf );
        words = tokenize_string( linebuf );
        for ( int entry = 0; entry < 5; entry++ ) {
          HCk[entry] = atof( words[entry + 1].c_str() );
        }
        
      } 

      //we've now read a complete set; assign to the appropriate type
      AtomParamSet tempSet;
      tempSet.elem = elem;
      tempSet.label = label;
      tempSet.nHeavy = nHeavy;
      tempSet.nLight = nLight;
      tempSet.surfaceCk = surfaceCk;
      if ( nLight > 0 ) {
        tempSet.HCk = HCk;
      }

      if ( elem == "C" && nHeavy == 1 && nLight == 3 ) {
        C13 = tempSet;
      }
      else if ( elem == "C" && nHeavy == 2 && nLight == 2 ) {
        C22 = tempSet;
      }
      else if ( elem == "C" && nHeavy == 3 && nLight == 1 ) {
        C31 = tempSet;
      }
      else if ( elem == "C" && nHeavy == 2 && nLight == 1 ) {
        C21 = tempSet;
      }
      else if ( elem == "C" && nHeavy == 3 && nLight == 0 ) {
        C30 = tempSet;
      }
      else if ( elem == "N" && nHeavy == 1 && nLight == 3 ) {
        N13 = tempSet;
      }
      else if ( elem == "N" && nHeavy == 2 && nLight == 2 ) {
        N22 = tempSet;
      }
      else if ( elem == "N" && nHeavy == 1 && nLight == 2 ) {
        N12 = tempSet;
      }
      else if ( elem == "N" && nHeavy == 2 && nLight == 1 ) {
        N21 = tempSet;
      }
      else if ( elem == "N" && nHeavy == 3 && nLight == 0 ) {
        N30 = tempSet;
      }
      else if ( elem == "N" && nHeavy == 2 && nLight == 0 ) {
        N20 = tempSet;
      }
      else if ( elem == "O" && nHeavy == 1 && nLight == 1 ) {
        O11 = tempSet;
      }
      else if ( elem == "O" && nHeavy == 1 && nLight == 0 && label == "sp2" ) {
        O10 = tempSet;
      }
      else if ( elem == "O" && nHeavy == 1 && nLight == 0 && label == "carboxylate" ) {
        Ocarb10 = tempSet;
      }
      else if ( elem == "S" && nHeavy == 1 && nLight == 1 ) {
        S11 = tempSet;
      }
      else if ( elem == "S" && nHeavy == 2 && nLight == 0 ) {
        S20 = tempSet;
      }
      else {
        cout << "Reading unidentifiable atom type surface params!" << endl;
      }

    }

  }

  //now loop over the atoms assigning them their surface parameters

  for (unsigned int siteNumber=1; siteNumber <= structure.total_sites; siteNumber++){
    frodaAtom.at(siteNumber).surfaceCk.resize(5);
    frodaAtom.at(siteNumber).surfaceCk[0] = 0.0;
    frodaAtom.at(siteNumber).surfaceCk[1] = 0.0;
    frodaAtom.at(siteNumber).surfaceCk[2] = 0.0;
    frodaAtom.at(siteNumber).surfaceCk[3] = 0.0;
    frodaAtom.at(siteNumber).surfaceCk[4] = 0.0;

  }

  for (unsigned int siteNumber=1; siteNumber <= structure.total_sites; siteNumber++){
    element = structure.site_info[siteNumber].element_name;
    nTotal = structure.site_info[siteNumber].nCovalentNeighbors;
    nLight = structure.site_info[siteNumber].nHydrogens;
    nHeavy = nTotal - nLight;

    if ( nTotal == 1 ) {
      frodaAtom.at(siteNumber).isTerminal = true; //single bonded neighbour
    }

    //C, N O and S get their own values;
    //H gets a value depending on its parent;
    //others get a default
    if ( element == "H " ) {
      frodaAtom.at(siteNumber).surfaceRho = 1.20;
      continue; //wait for parent
    }
    else if ( element == "C " ) {
      frodaAtom.at(siteNumber).surfaceRho = 1.70;

      if ( nHeavy == 1 && nLight == 3 ) {
        frodaAtom.at(siteNumber).surfaceCk = C13.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = C13.HCk;
          }
        }        
      }
      else if ( nHeavy == 2 && nLight == 2 ) {
        frodaAtom.at(siteNumber).surfaceCk = C22.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = C22.HCk;
          }
        }        
      }
      else if ( nHeavy == 3 && nLight == 1 ) {
        frodaAtom.at(siteNumber).surfaceCk = C31.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = C31.HCk;
          }
        }        
      }
      else if ( nHeavy == 2 && nLight == 1 ) {
        frodaAtom.at(siteNumber).surfaceCk = C21.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = C21.HCk;
          }
        }        
      }
      else if ( nHeavy == 3 && nLight == 0 ) {
        frodaAtom.at(siteNumber).surfaceCk = C30.surfaceCk;       
      }
      else {
        cout << "Warning: atom " << siteNumber << " not a recognised class for area calcs." << endl;
      }
    }
    else if ( element == "N " ) {
      frodaAtom.at(siteNumber).surfaceRho = 1.55;
      if ( nHeavy == 1 && nLight == 3 ) {
        frodaAtom.at(siteNumber).surfaceCk = N13.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = N13.HCk;
          }
        }        
      }
      else if ( nHeavy == 2 && nLight == 2 ) {
        frodaAtom.at(siteNumber).surfaceCk = N22.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = N22.HCk;
          }
        }        
      }
      else if ( nHeavy == 1 && nLight == 2 ) {
        frodaAtom.at(siteNumber).surfaceCk = N12.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = N12.HCk;
          }
        }        
      }
      else if ( nHeavy == 2 && nLight == 1 ) {
        frodaAtom.at(siteNumber).surfaceCk = N21.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = N21.HCk;
          }
        }        
      }
      else if ( nHeavy == 3 && nLight == 0 ) {
        frodaAtom.at(siteNumber).surfaceCk = N30.surfaceCk;       
      }
      else if ( nHeavy == 2 && nLight == 0 ) {
        frodaAtom.at(siteNumber).surfaceCk = N20.surfaceCk;       
      }
      else {
        cout << "Warning: atom " << siteNumber << " not a recognised class for area calcs." << endl;
      }

    }
    else if ( element == "O " ) {
      frodaAtom.at(siteNumber).surfaceRho = 1.52;
      if ( nHeavy == 1 && nLight == 1 ) {
        frodaAtom.at(siteNumber).surfaceCk = O11.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = O11.HCk;
          }
        }        
      }
      else if ( nHeavy == 1 && nLight == 0 && structure.isMainchain( siteNumber ) ) {
        //this is the mainchain carbonyl oxygen, which uses data type O10
        frodaAtom.at(siteNumber).surfaceCk = O10.surfaceCk;       
      }
      else if ( nHeavy == 1 && nLight == 0 && structure.site_info[siteNumber].residue_name == "GLN" ) {
        //it's a glutamine amide carbonyl oxygen
        frodaAtom.at(siteNumber).surfaceCk = O10.surfaceCk;       
      }
      else if ( nHeavy == 1 && nLight == 0 && structure.site_info[siteNumber].residue_name == "ASN" ) {
        //it's an asparagine amide carbonyl oxygen
        frodaAtom.at(siteNumber).surfaceCk = O10.surfaceCk;       
      }
      else if ( nHeavy == 1 && nLight == 0 ) {
        //it's a carboxylate oxygen
        frodaAtom.at(siteNumber).surfaceCk = Ocarb10.surfaceCk;       
      }
      else {
        cout << "Warning: atom " << siteNumber << " not a recognised class for area calcs." << endl;
      }

    }
    else if ( element == "S " ) {
      frodaAtom.at(siteNumber).surfaceRho = 1.80;
      if ( nHeavy == 1 && nLight == 1 ) {
        frodaAtom.at(siteNumber).surfaceCk = S11.surfaceCk;       
        //label hydrogens too
        for ( int neighbor = 0; neighbor < nTotal; neighbor++ ) {
          unsigned int nay = structure.site_info[siteNumber].neighbor_list.at(neighbor);
          if ( frodaAtom.at(nay).isHydrogen ) {
            frodaAtom.at(nay).surfaceCk = S11.HCk;
          }
        }        
      }
      else if ( nHeavy == 2 && nLight == 0 ) {
        frodaAtom.at(siteNumber).surfaceCk = S20.surfaceCk;       
      }
      else {
        cout << "Warning: atom " << siteNumber << " not a recognised class for area calcs." << endl;
      }

    }
    else {
      frodaAtom.at(siteNumber).surfaceRho = 1.50;
      //placeholder
      cout << "Warning: atom " << siteNumber << " getting placeholder surface values." << endl;
    }

  }

}



////////////////////////////////////////////////////////////////////////////////
// Description: labels hydrophobic atoms
////////////////////////////////////////////////////////////////////////////////
void Froda::identifyPhobicAtoms(MolFramework &structure){

    //find the phobic atoms and populate the lists
    bool is_phobic;

    for (unsigned int siteNumber=1; siteNumber <= structure.total_sites; siteNumber++){
      is_phobic = false;

      if (structure.isHydrophobicAtom( siteNumber ) && structure.isConnectedToOnlyHydrophobicAtoms( siteNumber) ) is_phobic = true; 

      if (is_phobic){
        //cerr << "Atom " << siteNumber << ", species " << structure.site_info[siteNumber].element_name << " is phobic." << endl;
        phobicAtoms.push_back(siteNumber);
        frodaAtom.at(siteNumber).isPhobic = true;
      } else frodaAtom.at(siteNumber).isPhobic = false;

    }

    if( parameters.verbose >= 2 ){
      cerr << "  Found " << phobicAtoms.size() << " hydrophobic atoms" << endl;
    }

  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: finds matching atoms in the target structure
////////////////////////////////////////////////////////////////////////////////
void Froda::setUpTargetList(MolFramework &structure, MolFramework &target){

  initialToTargetL = 0;
  nTargetedMainchainAtoms = 0;
  Vector delta;

  if( fancy_targetting ){
    nAtomsWithTargets = 0;
    if( parameters.verbose ){
      cout << "Setting up fancy-target list: " << endl;
    }

    string temp1;
    const char* errmsg;
    int int_atom_number;
    int slen;
    map<int, int> structure_FIRSTnumber_from_orignumber;
    for (unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++){
      temp1= structure.site_info[siteNumber].orig_atom_number;
      slen = strlen ( temp1.c_str());
      errmsg = hy36decode(5, temp1.c_str(), slen , &int_atom_number);
      if ( errmsg )
	{
	  if (isNumber( temp1 )){
	    int_atom_number = atoi( temp1.c_str() );
	    structure_FIRSTnumber_from_orignumber.insert( 
	    pair<int, int> ( int_atom_number, structure.site_info[siteNumber].FIRST_number));
	  }
	  else{
	    cout << "ERROR: Atom [" << temp1 << "] has invalid number" << endl;
	    exit(1);
	  }
	}
      else{
	pair<int, int> ( int_atom_number, structure.site_info[siteNumber].FIRST_number) ;
      }
    }
    
    // for each atom in the target structure, find the matching
    // atom in the initial structure, store the pair in myTarget
    for (SiteID siteNumber = 1; siteNumber <= target.total_sites; siteNumber++){
      int target_FIRSTnumber = target.site_info[siteNumber].FIRST_number;
      int target_orignumber = atoi( target.site_info[siteNumber].orig_atom_number.c_str() );
      if (structure_FIRSTnumber_from_orignumber.count(target_orignumber) > 0) {
        frodaAtom.at(structure_FIRSTnumber_from_orignumber[target_orignumber]).myTarget = 
          target_FIRSTnumber;
        nAtomsWithTargets++;
        if ( frodaAtom.at(siteNumber).isMain ) {
          nTargetedMainchainAtoms++;
          delta.x = structure.site_info[siteNumber].coords[0] - target.site_info[siteNumber].coords[0]; 
          delta.y = structure.site_info[siteNumber].coords[1] - target.site_info[siteNumber].coords[1]; 
          delta.z = structure.site_info[siteNumber].coords[2] - target.site_info[siteNumber].coords[2]; 
          initialToTargetL += dotProduct( delta, delta );
        }
      }      
    }
    
    
  } else  {
    nAtomsWithTargets = structure.total_sites;

    for (unsigned int siteNumber=1; siteNumber <= structure.total_sites; siteNumber++){
      frodaAtom.at(siteNumber).myTarget = siteNumber;
      if ( frodaAtom.at(siteNumber).isMain ) {
        nTargetedMainchainAtoms++;
        delta.x = structure.site_info[siteNumber].coords[0] - target.site_info[siteNumber].coords[0]; 
        delta.y = structure.site_info[siteNumber].coords[1] - target.site_info[siteNumber].coords[1]; 
        delta.z = structure.site_info[siteNumber].coords[2] - target.site_info[siteNumber].coords[2]; 
        initialToTargetL += dotProduct( delta, delta );
      }
    }
  }

  if ( nTargetedMainchainAtoms > 0 ) {
    initialToTargetL = sqrt( initialToTargetL / nTargetedMainchainAtoms );
  }

  initialToTargetL = MolFramework::computePairRMSD(target, 
                                                    initialPos,
                                                    *preferredSiteSelector);
    
  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: label possible donor H atoms
//              not currently used, but potentially useful for scoring
////////////////////////////////////////////////////////////////////////////////
void Froda::findDonorHydrogens( MolFramework &structure){
  int dummydonor;

  for (unsigned int siteNumber =1; siteNumber <= structure.total_sites ; siteNumber++){
    if ( structure.isDonorHydrogen( siteNumber, &dummydonor) ) donorHydrogenList.push_back(siteNumber);
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: label the Hbond acceptor candidates
//              not currently used; potentially useful for scoring 
////////////////////////////////////////////////////////////////////////////////
void Froda::findHBondingAtoms( MolFramework &structure){

  for (unsigned int siteNumber =1; siteNumber <= structure.total_sites ; siteNumber++){
    if ( structure.isHydrogenBondAcceptor(siteNumber) ) {
      frodaAtom.at(siteNumber).isAcceptor = true;}
    else frodaAtom.at(siteNumber).isAcceptor = false;
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: identify the neighbors in the structure
////////////////////////////////////////////////////////////////////////////////
void Froda::findNeighbors( MolFramework &structure ){

  for (unsigned int atomA=1; atomA <= structure.total_sites; atomA++){
    for (unsigned int atomANeighborNumber=0; atomANeighborNumber < structure.site_info[atomA].neighbor_list.size(); atomANeighborNumber++){
      unsigned int atomB = structure.site_info[atomA].neighbor_list.at(atomANeighborNumber);
      frodaAtom.at(atomA).neighbors.push_back(atomB);
    }
  }

  return;
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Read in a list of user-defined pair associations. The format should be three white-space 
//   delimited columns listing the two atom numbers
//   and the desired closeness of associations.
////////////////////////////////////////////////////////////////////////////////
void Froda::readPairFile( MolFramework &structure){

  string filename = structure.path_name + "pair.in";
  string linebuf;
  ifstream pairList;
  myPairs.clear();

  // Try to open a file named "pair.in". If the file can't be found, 
  // prompt the user for the name of the file listing their associations. 
  //////////////////////////////////////////////////////////////////////
  if( parameters.verbose ){
    cout << "Reading user-defined atom pairs file:" << filename << endl;
  }
  pairList.open( filename.c_str() );
 
  if( !pairList ){

    clear_screen;
    cout << endl;
    cout << "  The file [" << filename << "] was not found in this directory." << endl << endl;
    do{
      pairList.clear();
      cout << "  Please enter the name of the file that contains the list of pair associations." << endl;
      cout << "  Filename = ";
      getline(cin, filename);
      
      pairList.open( filename.c_str(), ios::in );
      if( !pairList || did_press_enter(filename) ){
        clear_screen;
        cout << endl << "File not found." << endl << endl;
      }
      
    } while( !pairList || did_press_enter(filename) );
  }

  // Read each line of the file. The first and the second numbers are atoms.
  // The third number is how close to bring these two atoms 
  // if there exists a fourth tag, a or b, this affects which atom is perturbed
  //////////////////////////////////////////////////////////////////////
  unsigned int atom1 = 0;
  unsigned int atom2 = 0;
  float distance = 0;
  bool tagA = true;
  bool tagB = true;

  size_t field_start = 0;
  size_t field_end = 0;
  string number;

  vector< string > choppedLine;

  while( !pairList.eof() ){
    getline(pairList, linebuf);
    field_start = field_end = 0;

    // Skip blank lines and lines that begin with '#'.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos ||
        isComment(linebuf) )
      break;

    choppedLine = tokenize_string( linebuf );

    // Read in the first atom number
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", 0);
    field_end   = linebuf.find_first_of(" \t\n", field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    atom1 = atoi( number.c_str() );
    
    // Read in the second atom number
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    atom2 = atoi( number.c_str() );

    // Read the distance of the given association. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    distance = atof( number.c_str() );

    //is there a tag after the distance?
    if ( choppedLine.size() > 3 ) {
      if ( toLower( choppedLine.at(3) ) == "a" ) {
        tagB = false;
      }
      else if ( toLower( choppedLine.at(3) ) == "b" ) {
        tagA = false;
      }
    }

    // Map the original numbers into the numbers used internally by FIRST.
    //////////////////////////////////////////////////////////////////////
    if ( !parameters.use_first_numbering ) {
      atom1 = structure.orig_2_FIRST_atom_number[atom1];
      atom2 = structure.orig_2_FIRST_atom_number[atom2];
    }

    if( atom1 <= 0 || atom1 > structure.total_sites ||
        atom2 <= 0 || atom2 > structure.total_sites ){
      cout << "ERROR: An atom in the file [" << filename << "] does not exist." << endl;
      cout << "       The error occurred with the following line:" << endl;
      cout << "       " << linebuf << endl << endl;
      exit(1);
    }

    myPairs.push_back( myPair( atom1, atom2, distance, tagA, tagB ) );
 
    cerr << "Associating atoms " << atom1 << " and " << atom2 << " with r= " << distance << endl;
    if ( !tagA ) {
      cerr << "Only atom " << atom2 << " will be perturbed." << endl;
    }
    else if ( !tagB ) {
      cerr << "Only atom " << atom1 << " will be perturbed." << endl;
    }

  }

  nMyPairs = myPairs.size();

  if( parameters.verbose )
    cout << "Read " << nMyPairs << " pair associations from an external file." << endl;
  pairList.close();
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description: read a list of alpha carbons and elastic network vectors
// the format is atomNumber vectorX vectorY vectorZ
// atomNumber should be an alpha carbon
//scaling of the eigenvector will be done in the formBiasVectors routine
////////////////////////////////////////////////////////////////////////////////
void Froda::readElasticNetworkMode( MolFramework &structure ) {
  string filename = structure.path_name + "mode.in";
  string linebuf;
  ifstream modeList;
  nodalAtoms.clear();
  eigenvector.clear();
  

  // Try to open a file named "mode.in". If the file can't be found, 
  // prompt the user for the name of the file listing their associations. 
  //////////////////////////////////////////////////////////////////////
  cout << "Reading " << filename << endl;
  modeList.open( filename.c_str() );
 
  if( !modeList ){

    clear_screen;
    cout << endl;
    cout << " The file [" << filename << "] was not found in this directory." << endl << endl;
    do{
      modeList.clear();
      cout << " Please enter the name of the file that contains the list of node eigenvectors." << endl;
      cout << " Filename = ";
      getline(cin, filename);
      
      modeList.open( filename.c_str(), ios::in );
      if( !modeList || did_press_enter(filename) ){
        clear_screen;
        cout << endl << "File not found." << endl << endl;
      }
      
    } while( !modeList || did_press_enter(filename) );
  }

  // Read each line of the file. The first number is an atom.
  // The following vector is the eigenvector for that atoms residue. 
  //////////////////////////////////////////////////////////////////////
  unsigned int atom1 = 0;
  double inX, inY, inZ;
  size_t field_start = 0;
  size_t field_end = 0;
  string number;

  while( !modeList.eof() ){
    getline(modeList, linebuf);
    field_start = field_end = 0;

    // Skip blank lines and lines that begin with '#'.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos ||
        isComment(linebuf) )
      break;

    // Read in the first atom number
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", 0);
    field_end   = linebuf.find_first_of(" \t\n", field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    atom1 = atoi( number.c_str() );

    // Read the x component of the eigenvector. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    inX = atof( number.c_str() );
    // Read the y component of the eigenvector. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    inY = atof( number.c_str() );
    // Read the z component of the eigenvector. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    inZ = atof( number.c_str() );

    // Map the original numbers into the numbers used internally by FIRST.
    //////////////////////////////////////////////////////////////////////
    if ( !parameters.use_first_numbering ) {
      atom1 = structure.orig_2_FIRST_atom_number[atom1];
    }

    if( atom1 <= 0 || atom1 > structure.total_sites ) {
      cout << "ERROR: An atom in the file [" << filename << "] does not exist." << endl;
      cout << "       The error occurred with the following line:" << endl;
      cout << "       " << linebuf << endl << endl;
      exit(1);
    }

    nodalAtoms.push_back( atom1 );
    eigenvector.push_back( Vector( inX, inY, inZ ) );

    //cerr << "Bias read: " << atom1 << " " << inX << " " << inY << " " << inZ << endl; 
  }


  if( parameters.verbose )
    cout << "Read " << nodalAtoms.size() << " node eigenvectors from an external file." << endl;
  modeList.close();

}



void Froda::makeResidueList( MolFramework &structure) {
     
  unsigned int catom, crigidlabel;
  unsigned int oatom, rigidlabel;
    
  vector< unsigned int > atomsToCheck;
  vector< unsigned int > checkedAtoms;
          
  nodeResidues.resize( alphaCarbon.size() ); //one nodeResidue per alpha carbon
  for ( unsigned int carbon = 0; carbon < alphaCarbon.size(); carbon++){
    catom = alphaCarbon.at(carbon);
    nodeResidues.at(carbon).cAlpha = catom;
    for ( unsigned int k = 0; k < structure.site_info[catom].neighbor_list.size(); k++ ) {
      oatom = structure.site_info[catom].neighbor_list.at(k);
      if ( structure.isBackbone( oatom ) == 1 ) {
        nodeResidues.at(carbon).n1 = oatom;
      }
      else if ( structure.isBackbone( oatom ) == 3 ) {
        nodeResidues.at(carbon).n2 = oatom;
      }
    }
    if ( nodeResidues.at(carbon).n1 == 0 || nodeResidues.at(carbon).n2 == 0 ) {
      cerr << "DYING as c-alpha lacks neighbors C,N; ";
      cerr << nodeResidues.at(carbon).cAlpha << " ";
      cerr << nodeResidues.at(carbon).n1 << " ";
      cerr << nodeResidues.at(carbon).n2 << " ";
      cerr << endl;
      cerr.flush();
      exit(1);
    }
  //call old talksto logic to build nodeResidue.sameResidue array

    //cerr << "Checking carbon " << carbon << " / catom " << catom <<endl;
    crigidlabel = structure.site_info[catom].rigid_label;

    atomsToCheck.clear();
    checkedAtoms.clear();

    for ( unsigned int j=0; j < structure.site_info[catom].neighbor_list.size(); j++ ){
      atomsToCheck.push_back( structure.site_info[catom].neighbor_list.at(j) );
      //cerr << "Adding atom " << structure.site_info[catom].neighbor_list.at(j) << endl;
    }

    bool looking = true;

    while(looking) {
      //cerr << "Looping" << endl;
      looking = false;
      int nToCheck = atomsToCheck.size();
      //cerr << "Checking " << nToCheck << " atoms in list " << endl;
      for ( int j = 0; j < nToCheck; j++ ) {
        oatom = atomsToCheck.at(j);
        //cerr << "Checking " << oatom << endl;
        bool wasChecked = false;
        for ( unsigned int k=0; k < checkedAtoms.size(); k++ ) {
          if ( checkedAtoms.at(k) == oatom ) {
            wasChecked = true;
            break;
          }
        }
        if (wasChecked) continue;
        looking = true; //if we get here, at least one atom wasn't checked yet
        checkedAtoms.push_back( oatom );
 
        rigidlabel = structure.site_info[oatom].rigid_label;
        
        if ( structure.isDifferentResidue( catom, oatom ) ) continue;

        nodeResidues.at(carbon).sameResidue.push_back( oatom );

        for ( unsigned int  k =0; k < structure.site_info[oatom].neighbor_list.size(); k ++ ) {
          unsigned int patom = structure.site_info[oatom].neighbor_list.at(k);
          if ( patom == catom ) continue;
          bool pWasChecked = false;
          for ( unsigned int k=0; k < checkedAtoms.size(); k++ ) {
            if ( checkedAtoms.at(k) == patom ) {
              pWasChecked = true;
              break;
            }
          }
          if (pWasChecked) continue;
          atomsToCheck.push_back( patom );
        }        


      }

    }


  }
  //nodeResidues is now indexed as is alphaCarbons
  
}


////////////////////////////////////////////////////////////////////////////////
// Description: turns the ENM vectors into biases on the residues.
////////////////////////////////////////////////////////////////////////////////
void Froda::formBiasVectors() {
  unsigned int nNodes = nodalAtoms.size();
  
  //scale all the eigenvectors so that the longest is targetStep;
  double l, maxl;

  maxl = 0;
  for ( unsigned int which=0; which < nNodes; which++) {
    l = dotProduct( eigenvector.at(which), eigenvector.at(which) );
    if ( l > maxl ) maxl = l;
  }

  maxl = sqrt(maxl); //this is length of longest atom eigenvector
  //cerr << "Longest eigenvector found: " << maxl << endl;
  //cerr << "Scaling to match directed step " << targetStep << endl;

  double multiplier = targetStep/maxl;
  for ( unsigned int which=0; which < nNodes; which++) {
    eigenvector.at(which).x *= multiplier;
    eigenvector.at(which).y *= multiplier;
    eigenvector.at(which).z *= multiplier;
  }
  
  //now the eigenvectors are scales to dstep
  //and can be given to the residues

  bool foundCA;
  for ( unsigned int whichNode = 0; whichNode < nNodes; whichNode++ ) {
    //cerr << "Checking nodal atom: " << nodalAtoms.at(whichNode) << endl;
    foundCA = false;
    for ( unsigned int whichCA = 0; whichCA < nodeResidues.size(); whichCA++ ) {
      //cerr << "Comparing to residue CA: " << nodeResidues.at(whichCA).cAlpha << endl;
      if ( nodalAtoms.at(whichNode) == nodeResidues.at(whichCA).cAlpha ) {
        nodeResidues.at(whichCA).isBiased = true;
        nodeResidues.at(whichCA).originalBias = eigenvector.at(whichNode);
        nodeResidues.at(whichCA).makeInternalBias();
        foundCA = true;
      }
    }
    if ( !foundCA ) {
      cerr << "ERROR in mode bias routine." << endl;
      cerr << "Did not find nodal atom " << nodalAtoms.at(whichNode) << " among alpha carbons in residue list." << endl;
      exit(1);
    }

  }

}

////////////////////////////////////////////////////////////////////////////////
// Description: Load matrices used for symmetry
////////////////////////////////////////////////////////////////////////////////
void Froda::setMatrices(vector <Matrix> *matrices_) {
  symMat = new SymmetryMatrices(matrices_);
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: pull the variables from the parameter list into FRODA
////////////////////////////////////////////////////////////////////////////////
void Froda::obtainParameters( MolFramework &structure){
  if (parameters.using_restart) isRestart=1;
  if (parameters.using_target ) isTargeted =1;
  if (parameters.centri) centrifugal =1;
  centri_maxr = parameters.centri_maxr;

  localCentrifuge = false;
  if ( parameters.localCenter ) {

    cerr << "Using local centrifuge in setup." << endl;

    localCentrifuge = true;
    centerID.resize( parameters.localCenters.size() );
    centerDirection.resize( parameters.localCenters.size() );

    for ( unsigned int whichLC = 0; whichLC < centerID.size(); whichLC++ ) {
      if ( parameters.localCenters[whichLC] < 0 ) {
        centerDirection[whichLC] = -1;
        parameters.localCenters[whichLC] *= -1;
        centerID[whichLC] = (unsigned int) parameters.localCenters[whichLC];
      }
      else {
        centerDirection[whichLC] = 1;
        centerID[whichLC] = (unsigned int) parameters.localCenters[whichLC];
      }

      if ( !parameters.use_first_numbering ) {
        centerID[whichLC] = structure.orig_2_FIRST_atom_number[ centerID[whichLC] ];
      }

    }

  }
  

  bounce = parameters.step_size;
  outputEvery = parameters.output_freq;
  nConfigs = parameters.total_conformations; 
  if (parameters.body) { 
    moveByCluster =true;
  } else {
    moveByCluster =false;
  }

  moveDihedral = parameters.dihedral;
  makeSingleMove = parameters.makeSingleMove;
  targetStep = parameters.dstep;
  tolTarget = parameters.dtol;
  moveProb = parameters.prob;
  fancy_targetting = parameters.fancy_target;
  ghostTol = parameters.ghostTol;
  phobeTol = parameters.ph_radius;
  vdwTol = parameters.vdw_overlap;
  maxFitCycles = parameters.maxFitCycles;
  bodyResponseEvery = parameters.bodyResponseEvery;
  if (bodyResponseEvery ==0 ) use_bodyResponse = false;
  atLeastCycles = parameters.atleast;
  cylindrical = parameters.cylindrical;
  lowZRange = parameters.cyl_low_z;
  highZRange = parameters.cyl_high_z;

  localTargeting = parameters.local;
  lstep = parameters.lstep;

  vdwNotice = parameters.vdwNotice;
  doAnneal = parameters.doAnneal;
  annealScale = parameters.annealScale;
  constantRate = parameters.constantRate;
  acceptRate = parameters.acceptRate;

  propTotarget = parameters.propto;
  if (propTotarget) cerr << "Proportional" << endl;

  checkRama = parameters.rama;

  useSeed= parameters.useSeed;
  seed= parameters.seed;
  mobileRC1= parameters.mobileRC1;
  verbose= parameters.verbose;
  interactive= parameters.interactive;

  runMorph= parameters.runMorph;
  morphFrames= parameters.morphFrames;
  RMSDSpacing= parameters.RMSDSpacing;
  reportByRSMD= parameters.reportByRSMD;
  polarGeometry= parameters.polarGeometry;
  polarHRadius= parameters.polarHRadius;
  use_group_id= parameters.use_group_id;
  energy_cutoff= parameters.energy_cutoff;

  doScoringMC = parameters.scoreMC;
  frodaMCScale = parameters.frodaMCScale;

  mcRg = parameters.mcRg;
  minimumRg = parameters.minimumRg;
  RgWeight = parameters.RgWeight; 

  mcPhobicRg = parameters.mcPhobicRg;
  minimumPhobicRg = parameters.minimumPhobicRg;
  phobicRgWeight = parameters.phobicRgWeight; 

  mcSASA = parameters.mcSASA;
  minimumSASA = parameters.minimumSASA;
  SASAWeight = parameters.SASAWeight; 

  mcPhobicSASA = parameters.mcPhobicSASA;
  minimumPhobicSASA = parameters.minimumPhobicSASA;
  phobicSASAWeight = parameters.phobicSASAWeight; 

  mcPolarSASA = parameters.mcPolarSASA;
  polarSASAWeight = parameters.polarSASAWeight; 

  doHBScoring = parameters.doHBScoring;
  HBWellDepth = parameters.HBWellDepth; 
  SBWellDepth = parameters.SBWellDepth;

  doGetArea = parameters.doGetArea;
  doGetPhobicRg = parameters.doGetPhobicRg;
 
  nosterics = parameters.nosterics;

  chatty_fitting = parameters.chatty;

  useCM = parameters.useCM;
  useCMforward = parameters.useCMforward;
  useCMreverse = parameters.useCMreverse;
  cmforwardfactor = parameters.CMforward;
  cmreversefactor = parameters.CMreverse;
  cmlowerFitcycles = parameters.CMlower;
  cmupperFitcycles = parameters.CMupper;
  cmmaxstep = parameters.CMmax;

  persist = parameters.persist;
  nPersist = parameters.nPersist;

  usePairs = parameters.usePairs;
  stopWhenAllPairConstraintsSatisfied = usePairs && parameters.stopWhenAllPairConstraintsSatisfied;
  
  useElasticVector = parameters.useElasticVector;
  useInternalBias = parameters.useInternalBias;

  //useTalksto = parameters.useTalksto;

// for electron-density options
  useED = parameters.useED;
  useTheoMap = parameters.useTheoMap;
  useEZDMap = parameters.useEZDMap;
  trimMap = parameters.trimMap;
  edResFac = parameters.edResFac;
  theoMapFile = parameters.theoMapFile;
  ezdMapFile = parameters.ezdMapFile;
  tolED = parameters.edTol;
  edNoise = parameters.edNoise;
  trimMapFactor = parameters.trimMapFactor;
// for SAXS profiles
  useSAXS = parameters.useSAXS;
  saxsTol = parameters.saxsTol;
  saxsFile = parameters.saxsFile;

  //override mobRC1 setting if we're using the froda2Hybrid
  if ( froda2Hybrid ) {
    mobileRC1=true;
  }
  return;
}






