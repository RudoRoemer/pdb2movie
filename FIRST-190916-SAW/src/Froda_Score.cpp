#include "global_defs.h"
#include "interact.h"
#include "Froda.h"
#include "Sterics.h"
#include "mt19937ar.h"
#include "flexweb.h"
#include "EDStructure.h"
#include "GenericMap.h"
#include "AreaProximity.h"

#include <signal.h>

extern Parameters parameters;


////////////////////////////////////////////////////////////////////////////////
// Description: evaluate phi, psi for ramachandran constraints
//////////////////////////////////////////////////////////////////////////////
int Froda::checkRamaPlot( MolFramework &structure){
  int return_val=0;

  Vector bond1, bond2, bond3, bond4; // actual bonds
  Vector e1,e2,e3,e4; // unit vectors
  double phi,psi, pmet;
  double didot, diangle;
  bool negangle=false;
  int bin_phi, bin_psi; // will address the arrays
  int mapvalue;
  int dubious=0;
  int allowed=0;
  int core =0;
  int Ci_11, Ni1, Calpha1, Ci1, Ni11;


  for (unsigned int ramaNumber=0; ramaNumber<ramaList.size();ramaNumber++){
    dummyRama = ramaList.at(ramaNumber);
    mapvalue=0;

    if ( !dummyRama.freePhi && !dummyRama.freePsi) {
      continue; //no need to recalculate anything for a locked residue
    }

    // get the four bonds:
    
    Ci_11 = dummyRama.Ci_1 ;
    Ni1 = dummyRama.Ni ;
    Calpha1 = dummyRama.Calpha ;
    Ci1 = dummyRama.Ci ;
    Ni11 = dummyRama.Ni1 ;
    
    bond1.x = currentPos.at(Ni1).x - currentPos.at(Ci_11).x;
    bond1.y = currentPos.at(Ni1).y - currentPos.at(Ci_11).y;
    bond1.z = currentPos.at(Ni1).z - currentPos.at(Ci_11).z;
    bond2.x = currentPos.at(Calpha1).x - currentPos.at(Ni1).x;
    bond2.y = currentPos.at(Calpha1).y - currentPos.at(Ni1).y;
    bond2.z = currentPos.at(Calpha1).z - currentPos.at(Ni1).z;
    bond3.x = currentPos.at(Ci1).x - currentPos.at(Calpha1).x;
    bond3.y = currentPos.at(Ci1).y - currentPos.at(Calpha1).y;
    bond3.z = currentPos.at(Ci1).z - currentPos.at(Calpha1).z;
    bond4.x = currentPos.at(Ni11).x - currentPos.at(Ci1).x;
    bond4.y = currentPos.at(Ni11).y - currentPos.at(Ci1).y;
    bond4.z = currentPos.at(Ni11).z - currentPos.at(Ci1).z;


    // calculate phi:

    negangle=false;

    dot = sqrt(dotProduct(bond2, bond2));
    pmet = 1.0/dot;
    e1.x=bond2.x*pmet;
    e1.y=bond2.y*pmet;
    e1.z=bond2.z*pmet;

    dummyvec.x = -bond1.x;
    dummyvec.y = -bond1.y;
    dummyvec.z = -bond1.z;
 
    dot = dotProduct(dummyvec, e1);
    dummyvec.x -= dot*e1.x; 
    dummyvec.y -= dot*e1.y; 
    dummyvec.z -= dot*e1.z; 
    dot =sqrt(dotProduct(dummyvec,dummyvec));
    pmet = 1.0/dot;
    e2.x = dummyvec.x*pmet;
    e2.y = dummyvec.y*pmet;
    e2.z = dummyvec.z*pmet;

    dummyvec.x = bond3.x;
    dummyvec.y = bond3.y;
    dummyvec.z = bond3.z;
    
    dot = dotProduct(dummyvec, e1);
    dummyvec.x -= dot*e1.x; 
    dummyvec.y -= dot*e1.y; 
    dummyvec.z -= dot*e1.z; 
    dot =sqrt(dotProduct(dummyvec,dummyvec));
    pmet = 1.0/dot;
    e3.x = dummyvec.x*pmet;
    e3.y = dummyvec.y*pmet;
    e3.z = dummyvec.z*pmet;

    e4 = crossProduct(e1,e2);
    didot = dotProduct(e3,e4);
    if(didot<0.0) negangle=true; // is in range negative
    didot = dotProduct(e2,e3);

    if(didot<-1.0) didot=-0.99999999;
    if(didot>1.0) didot=0.99999999;

    diangle = RAD_TO_DEG*acos(didot); // 0 to 180 degrees
    if (negangle) diangle = -diangle;

    phi = diangle;    


    // calculate psi:

    negangle=false;
    dot = sqrt(dotProduct(bond3, bond3));
    pmet = 1.0/dot;
    e1.x=bond3.x*pmet;
    e1.y=bond3.y*pmet;
    e1.z=bond3.z*pmet;
    dummyvec.x = -bond2.x;
    dummyvec.y = -bond2.y;
    dummyvec.z = -bond2.z;
    dot = dotProduct(dummyvec, e1);
    dummyvec.x -= dot*e1.x; 
    dummyvec.y -= dot*e1.y; 
    dummyvec.z -= dot*e1.z; 
    dot =sqrt(dotProduct(dummyvec,dummyvec));
    pmet = 1.0/dot;
    e2.x = dummyvec.x*pmet;
    e2.y = dummyvec.y*pmet;
    e2.z = dummyvec.z*pmet;
    dummyvec.x = bond4.x;
    dummyvec.y = bond4.y;
    dummyvec.z = bond4.z;
    dot = dotProduct(dummyvec, e1);
    dummyvec.x -= dot*e1.x; 
    dummyvec.y -= dot*e1.y; 
    dummyvec.z -= dot*e1.z; 
    dot =sqrt(dotProduct(dummyvec,dummyvec));
    pmet = 1.0/dot;
    e3.x = dummyvec.x*pmet;
    e3.y = dummyvec.y*pmet;
    e3.z = dummyvec.z*pmet;
    e4 = crossProduct(e1,e2);
    didot = dotProduct(e3,e4);
    if(didot<0.0) negangle=true; // is in range negative
    didot = dotProduct(e2,e3);
    if(didot<-1.0) didot=-0.99999999;
    if(didot>1.0) didot=0.99999999;
    diangle = RAD_TO_DEG*acos(didot); // 0 to 180 degrees
    if (negangle) diangle = -diangle;

    psi = diangle;

    ramaList.at(ramaNumber).phi0 = phi;
    ramaList.at(ramaNumber).psi0 = psi;

   //get region:
    //if (structure.site_info[dummyRama.Calpha].residue_name == "PRO"){
    if (dummyRama.isPro){
      bin_phi = (int) ((phi+180.0)*0.125);
      bin_psi = (int) ((psi+180.0)*0.125);

      if (bin_phi==45) bin_phi=44;
      if (bin_psi==45) bin_psi=44;
      if (bin_phi==-1) bin_phi=0;
      if (bin_psi==-1) bin_psi=0;

      mapvalue = prolineRamaPlot.at(bin_psi).at(bin_phi);
    } 
    //else if (structure.site_info[dummyRama.Calpha].residue_name == "GLY"){
    else if (dummyRama.isGly){
      bin_phi = (int) ((phi+180.0)*0.125);
      bin_psi = (int) ((psi+180.0)*0.125);
      if (bin_phi==45) bin_phi=44;
      if (bin_psi==45) bin_psi=44;
      if (bin_phi==-1) bin_phi=0;
      if (bin_psi==-1) bin_psi=0;
      mapvalue = glycineRamaPlot.at(bin_psi).at(bin_phi);
    }
    else{
      bin_phi = (int) ((phi+180.0)*0.1);
      bin_psi = (int) ((psi+180.0)*0.1);
      if (bin_phi==36) bin_phi=35;
      if (bin_psi==36) bin_psi=35;
      if (bin_phi==-1) bin_phi=0;
      if (bin_psi==-1) bin_psi=0;
      mapvalue = generalRamaPlot.at(bin_psi).at(bin_phi);
    }

    ramaList.at(ramaNumber).mapValue = mapvalue;

    if (mapvalue==9) return_val++; // bad position
    if (mapvalue==2) dubious++; // edgy
    if (mapvalue==1) allowed++; // allowed region
    if (mapvalue==0) core++; // core region

 }

 nBadRama = return_val;
 nGenerousRama = dubious;
 nAllowedRama = allowed;
 nCoreRama = core;

 if (verbose && chatty_fitting) cerr << "RAMA " << setw(3) << return_val << "/"  << setw(3) << dubious << "/"  << setw(3) << allowed << "/" << setw(3) << core << "#";
 return 0; //removed reject-rama option
}


////////////////////////////////////////////////////////////////////////////////
// Description:  check distance to a particular structure.
////////////////////////////////////////////////////////////////////////////////
void Froda::distanceToTarget (MolFramework &restart, MolFramework &target){

 bool is_found;
 double tol;

 tol = tolTarget;

 is_found =1;
 targetMSD =0;
 restartMSD =0;

 double target_mainchain_MSD =0;

// check for distance to target 

 for (unsigned int siteNumber=1; siteNumber<=nTotalSites;siteNumber++){

   if ( frodaAtom.at(siteNumber).isCore ) continue;

   if (isRestart){
     dummyvec.x = currentPos.at(siteNumber).x - restart.site_info[siteNumber].coords[0];
     dummyvec.y = currentPos.at(siteNumber).y - restart.site_info[siteNumber].coords[1];
     dummyvec.z = currentPos.at(siteNumber).z - restart.site_info[siteNumber].coords[2];
     //vector from target to current position
   dot = (dotProduct(dummyvec, dummyvec));
   restartMSD += dot;
   }
   if (isTargeted) {  // targeted- check myTarget

     unsigned int targatom = frodaAtom.at(siteNumber).myTarget;
     if (targatom == 0 ) continue;

     dummyvec.x = currentPos.at(siteNumber).x - target.site_info[targatom].coords[0];
     dummyvec.y = currentPos.at(siteNumber).y - target.site_info[targatom].coords[1];
     dummyvec.z = currentPos.at(siteNumber).z - target.site_info[targatom].coords[2];
     dot = (dotProduct(dummyvec, dummyvec));
     targetMSD += dot;
     if (frodaAtom.at(siteNumber).isMain ) target_mainchain_MSD += dot;
   }

 }
 if (isRestart){
 restartRMSD = sqrt(restartMSD/nTotalSites);
 }
 if (isTargeted){
 targetRMSD = sqrt(targetMSD/nAtomsWithTargets);
   
 if ( nTargetedMainchainAtoms > 0 ){
   targmain_RMSD = sqrt(target_mainchain_MSD/nTargetedMainchainAtoms);
      
 }
 else {
   targmain_RMSD = 0;
 }
 if (targetRMSD > tol) is_found =0; // not fitted yet
 foundTarget = is_found;
 }
 
 currentRMSDfromInitial = MolFramework::computePairRMSD(initialPos, 
                                                         currentPos,
                                                         *preferredSiteSelector);

 currentRMSDfromTarget = MolFramework::computePairRMSD(target, 
                                                        currentPos,
                                                        *preferredSiteSelector);
 
 return;
}


////////////////////////////////////////////////////////////////////////////////
// Description: RMSD from starting point, by residue
////////////////////////////////////////////////////////////////////////////////
void Froda::getRMSDByResidue(){

  double r2;

  double allr2 = 0.0;

  for ( unsigned int alpha = 0; alpha < alphaCarbon.size(); alpha++){

    int carbon = alphaCarbon[alpha];

    dummyvec.x = currentPos.at(carbon).x - initialPos.at(carbon).x; 
    dummyvec.y = currentPos.at(carbon).y - initialPos.at(carbon).y; 
    dummyvec.z = currentPos.at(carbon).z - initialPos.at(carbon).z; 

    r2 = dotProduct(dummyvec, dummyvec); 

    allr2 += r2;

    runningResidueRMSD.at(alpha) += r2;
    currentResidueRMSD.at(alpha) = sqrt(r2);

  }

  //pass into running_allmain calculation

  //divide allr2 by the number of residues
  allr2 /= alphaCarbon.size();

  //add it onto the running MSD
  running_allmain_MSD += allr2;

}


////////////////////////////////////////////////////////////////////////////////
// Description: calculate RMSD from initial structure
////////////////////////////////////////////////////////////////////////////////
void Froda::getRMSD(){

      MSD = RMSD = 0.0;
      globalMSD = globalRMSD = 0.0;
      allmain_MSD = 0.0;
      allside_MSD = 0.0;
      mobmain_MSD = 0.0;
      mobside_MSD = 0.0;
      allTargetedMainchainMSD = 0.0;
      
      for( unsigned int siteNumber = 1; siteNumber <= nTotalSites; siteNumber++ ){
        dummyvec.x = currentPos.at(siteNumber).x - initialPos.at(siteNumber).x;
        dummyvec.y = currentPos.at(siteNumber).y - initialPos.at(siteNumber).y;
        dummyvec.z = currentPos.at(siteNumber).z - initialPos.at(siteNumber).z;
        dot = dotProduct(dummyvec,dummyvec);
        globalMSD += dot;
        if ( frodaAtom.at(siteNumber).isMain ){
          allmain_MSD += dot;
          if ( isTargeted ) {
            if ( frodaAtom.at(siteNumber).myTarget > 0 ) {
              allTargetedMainchainMSD += dot;
            }
          }
        } else{
          allside_MSD += dot;
        }

        if( !frodaAtom.at(siteNumber).isCore ){  // mobile atoms only
          MSD += dot;
          if ( frodaAtom.at(siteNumber).isMain ){
            mobmain_MSD += dot;
          } else{
            mobside_MSD += dot;
          }
        } 
      }
      
      globalRMSD = sqrt((globalMSD/nTotalSites));
      RMSD = sqrt(MSD/nMobileAtoms);

      if (nMainchainAtoms > 0 ){allmain_RMSD = sqrt(allmain_MSD/nMainchainAtoms);} else allmain_RMSD = 0.0;
      if (nSidechainAtoms > 0 ){allside_RMSD = sqrt(allside_MSD/nSidechainAtoms);} else allside_RMSD = 0.0;
      if (nMobileMainchainAtoms > 0){mobmain_RMSD = sqrt(mobmain_MSD/nMobileMainchainAtoms);} else mobmain_RMSD = 0.0;
      if (nMobileSidechainAtoms > 0){mobside_RMSD = sqrt(mobside_MSD/nMobileSidechainAtoms);} else mobside_RMSD = 0.0;
      if (nTargetedMainchainAtoms > 0 ) { 
        allTargetedMainchainRMSD = sqrt( allTargetedMainchainMSD/nTargetedMainchainAtoms);
      }
      else {
        allTargetedMainchainRMSD = 0.0;
      }
          
}


////////////////////////////////////////////////////////////////////////////////
// Description: decide to accept or reject a conformer
//              based on progress of RMSD in Monte-Carlo targeting
////////////////////////////////////////////////////////////////////////////////
void Froda::checkMonteCarloAnneal( MolFramework &structure){

        double deltarmsd = cachedRMSD - targetRMSD; // positive if we're improving
        if (deltarmsd < 0.0){
          double metropolis = genrand_real1();
          double chance = exp(deltarmsd/annealScale);
          if (metropolis > chance){ //reject conformer
			mcRejected++;
            currentPos = cachedPos;
            for( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ){

              structure.site_info[siteNumber].coords[0] = currentPos.at(siteNumber).x;
              structure.site_info[siteNumber].coords[1] = currentPos.at(siteNumber).y;
              structure.site_info[siteNumber].coords[2] = currentPos.at(siteNumber).z;
            }
          } 
        }
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Decide to accept or reject a conformer based on progress of correlation 
//   with an electron density map in Monte-Carlo targeting
// Parameters:  
//   EDStructure &pdb -- Structure containing the current FRODA conformer 
//                        to be evaluated
//   GenericMap &edMap -- Electron density map to which the current FRODA 
//                          conformer should be compared.
// Return Value List:
//   This function returns a value of type double containing the correlation
//   of the current FRODA conformer with the electron density map *after* the
//   accept/reject decision has been made.
////////////////////////////////////////////////////////////////////////////////
double Froda::checkMonteCarloED(EDStructure &currentStructure, GenericMap &edMap) {
  double oldCorrelation = currentStructure.cachedCorrelation;
  currentStructure.cachedCorrelation = edMap.correlate(currentStructure);
   // calculate correlation for current conformer
  double deltaCorrelation = currentStructure.cachedCorrelation - oldCorrelation;
   // will be > 0 if correlation improves
  if (deltaCorrelation < 0.0) {
    double metropolis = genrand_real1();
    double chance = exp(1.0*deltaCorrelation/annealScale);
    if (metropolis > chance) { //reject conformer
      mcRejected++;
      currentStructure.cachedCorrelation = oldCorrelation;
      currentPos = cachedPos;
      for(unsigned int siteNum = 1; siteNum <= currentStructure.structure->total_sites; siteNum++) {
        currentStructure.structure->site_info[siteNum].coords[0] = currentPos.at(siteNum).x;
        currentStructure.structure->site_info[siteNum].coords[1] = currentPos.at(siteNum).y;
        currentStructure.structure->site_info[siteNum].coords[2] = currentPos.at(siteNum).z;
      }
    } 
  }
  // check for convergence and return
  if (currentStructure.cachedCorrelation >= tolED) {
    foundTarget = 1;
  }
  return currentStructure.cachedCorrelation;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Decide to accept or reject a conformer based on progress of correlation 
//   with a SAXS profile in Monte-Carlo targeting
// Parameters:  
//   SAXS &saxs -- Structure containing the current FRODA conformer 
//                        to be evaluated
// Return Value List:
//   This function returns a value of type double containing the correlation
//   of the current FRODA conformer with the SAXS profile *after* the
//   accept/reject decision has been made.
////////////////////////////////////////////////////////////////////////////////
double Froda::checkMonteCarloSAXS(SAXS &saxs) {
  double oldCorrelation = saxs.cachedCorrelation;
  saxs.cachedCorrelation = saxs.correlate();
   // calculate correlation for current conformer
  double deltaCorrelation = saxs.cachedCorrelation - oldCorrelation;
  //cout << "Change in correlation: " << deltaCorrelation;
   // will be > 0 if correlation improves
  if (deltaCorrelation < 0.0) {
    double metropolis = genrand_real1();
    double chance = exp(deltaCorrelation/annealScale);
    //cout << ", metro = " << metropolis  << ", chance = " << chance;
    if (metropolis > chance) { //reject conformer
      //cout << ", rejected.\n";
      mcRejected++;
      saxs.cachedCorrelation = oldCorrelation;
      currentPos = cachedPos;
      for(unsigned int siteNum = 1; siteNum <= saxs.structure->total_sites; siteNum++) {
        saxs.structure->site_info[siteNum].coords[0] = currentPos.at(siteNum).x;
        saxs.structure->site_info[siteNum].coords[1] = currentPos.at(siteNum).y;
        saxs.structure->site_info[siteNum].coords[2] = currentPos.at(siteNum).z;
      }
    } else {
      //cout << ", accepted.\n";
    }     
  } else {
    //cout << ", automatic accept.\n";
  }
  // check for convergence and return
  if (saxs.cachedCorrelation >= saxsTol) {
    foundTarget = 1;
  }
  return saxs.cachedCorrelation;
}


////////////////////////////////////////////////////////////////////////////////
// Description: find new potential hydrogen bond interactions
////////////////////////////////////////////////////////////////////////////////
void Froda::findChargedInteractions( MolFramework &structure){
  //loop over the donor hydrogens

  HBScore = 0.0;
  vector<SiteID> nearbyList;

  for( unsigned int whichH = 0; whichH < donorHydrogenList.size(); whichH++ ) {
    SiteID hydrogen = donorHydrogenList[whichH];

    //cerr << "Working on hydrogen ";
    //cerr << structure.site_info[hydrogen].orig_atom_number << endl;

    //look for nearby acceptors
    nearbyList = sterics->getProximityList( hydrogen );

    for ( unsigned int whichA = 0; whichA < nearbyList.size(); whichA++ ) {
      SiteID acceptor = nearbyList[whichA];
      if ( !frodaAtom.at(acceptor).isAcceptor ) continue; //not an Acceptor
      //cerr << "  Working on acceptor ";
      //cerr << structure.site_info[acceptor].orig_atom_number;

      //cerr << "  at distance " << sqrt( getDist2( hydrogen, acceptor) ) << endl;

      //check if it's got any hbond status
      double energy = scoreHBond( hydrogen, acceptor, structure );

      if ( energy < 0.0 ) HBScore += energy;
    }

  }

  return;
} 



////////////////////////////////////////////////////////////////////////////////
// Description: accept/reject check on new interactions
////////////////////////////////////////////////////////////////////////////////
void Froda::checkMonteCarloBonding( MolFramework &structure){

  cerr << "Folding MC: ";
  double deltaScore = frodaMCScore - oldFrodaMCScore; // positive is bad

  cerr << "Old score " << oldFrodaMCScore << " ";
  cerr << "New score " << frodaMCScore << " ";

  cerr << " Delta: " << deltaScore;
  if (deltaScore > 0.0){

    double metropolis = genrand_real1();
    cerr << " random: " << metropolis << " ";
    double chance = exp( -deltaScore/frodaMCScale);
    cerr << " P(accept): " << chance << " ";
    if (metropolis > chance){ //reject conformer
      currentPos = cachedPos;
      cerr << "(Reject)";
      for( unsigned int a = 1; a <= structure.total_sites; a++ ){
        structure.site_info[a].coords[0] = currentPos.at(a).x;
        structure.site_info[a].coords[1] = currentPos.at(a).y;
        structure.site_info[a].coords[2] = currentPos.at(a).z;
      }
      restoreMCFoldingVariables();
      n_rejected++;
    } else { 
      cerr << "(Accept)";
      n_accepted++; 
    }
  }
  else { 
    cerr << "(Accept)";
    n_accepted++;
  }

  cerr << endl;
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: stored values of variables for monte-carlo folding
////////////////////////////////////////////////////////////////////////////////
void Froda::storeMCFoldingVariables(){
  oldFrodaMCScore = frodaMCScore;

  oldTotalSASA = totalSASA;
  oldPhobicSASA = phobicSASA;
  oldPolarSASA = polarSASA;
  oldPolarPairSA = polarPairSA;
  oldRg = Rg;
  oldPhobicRg = phobicRg;
 
  oldRgScore = RgScore;
  oldPhobicRgScore = phobicRgScore;
  oldSASAScore = SASAScore;
  oldPhobicSASAScore = phobicSASAScore;
  oldPolarSASAScore = polarSASAScore;
  oldHBScore = HBScore;

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: restore old values of variables for monte-carlo folding
////////////////////////////////////////////////////////////////////////////////
void Froda::restoreMCFoldingVariables(){
  frodaMCScore = oldFrodaMCScore;

  totalSASA = oldTotalSASA;
  phobicSASA = oldPhobicSASA;
  polarSASA = oldPolarSASA;
  polarPairSA = oldPolarPairSA;
  Rg = oldRg;
  phobicRg = oldPhobicRg;

  RgScore = oldRgScore;
  phobicRgScore = oldPhobicRgScore;
  SASAScore = oldSASAScore;
  phobicSASAScore = oldPhobicSASAScore;
  polarSASAScore = oldPolarSASAScore;
  HBScore = oldHBScore;

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: zero values of variables for monte-carlo folding
////////////////////////////////////////////////////////////////////////////////
void Froda::clearMCFoldingVariables(){

  frodaMCScore = 0.0;
  RgScore = 0.0;
  phobicRgScore = 0.0;
  SASAScore = 0.0;
  phobicSASAScore = 0.0;
  polarSASAScore = 0.0;
  HBScore = 0.0;

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: calculate radius of gyration for the protein
////////////////////////////////////////////////////////////////////////////////
double Froda::calculateRg() {
  double Rg; //radius of gyration of all atoms
  Vector centralPoint;
  Vector delta;

  Rg = 0;

  for ( unsigned int thisAtom = 1; thisAtom <= nTotalSites; thisAtom++ ) {
    centralPoint.x += currentPos.at(thisAtom).x;
    centralPoint.y += currentPos.at(thisAtom).y;
    centralPoint.z += currentPos.at(thisAtom).z;
  }

  centralPoint.x /= nTotalSites;
  centralPoint.y /= nTotalSites;
  centralPoint.z /= nTotalSites;

  for ( unsigned int thisAtom = 1; thisAtom <= nTotalSites; thisAtom++ ) {
    delta.x = currentPos.at(thisAtom).x - centralPoint.x;
    delta.y = currentPos.at(thisAtom).y - centralPoint.y;
    delta.z = currentPos.at(thisAtom).z - centralPoint.z;

    Rg += dotProduct( delta, delta );

  }

  Rg /= nTotalSites;
  Rg = sqrt(Rg);

  return Rg;
}

////////////////////////////////////////////////////////////////////////////////
// Description: calculate phobic radius of gyration for the protein
////////////////////////////////////////////////////////////////////////////////
double Froda::calculatePhobicRg() {
  double Rg; //radius of gyration of atoms
  Vector centralPoint;
  Vector delta;

  Rg = 0;

  if ( phobicBetaCarbon.size() == 0 ) return Rg;

  for ( unsigned int whichCB = 0; whichCB < phobicBetaCarbon.size(); whichCB ++ ) {
    unsigned int thisAtom = phobicBetaCarbon[whichCB];
    centralPoint.x += currentPos.at(thisAtom).x;
    centralPoint.y += currentPos.at(thisAtom).y;
    centralPoint.z += currentPos.at(thisAtom).z;
  }

  centralPoint.x /= phobicBetaCarbon.size();
  centralPoint.y /= phobicBetaCarbon.size();
  centralPoint.z /= phobicBetaCarbon.size();

  for ( unsigned int whichCB = 0; whichCB < phobicBetaCarbon.size(); whichCB ++ ) {
    unsigned int thisAtom = phobicBetaCarbon[whichCB];
    delta.x = currentPos.at(thisAtom).x - centralPoint.x;
    delta.y = currentPos.at(thisAtom).y - centralPoint.y;
    delta.z = currentPos.at(thisAtom).z - centralPoint.z;

    Rg += dotProduct( delta, delta );

  }

  Rg /= phobicBetaCarbon.size();
  Rg = sqrt(Rg);

  return Rg;
}

////////////////////////////////////////////////////////////////////////////////
// Description: calculate SASA of an atom, using Brooks/Guvench03 method
////////////////////////////////////////////////////////////////////////////////
void Froda::getAreaOfAtom( unsigned int atomID ) {
  double Ai;
  double kappa = 2.4;
  double Area;
  double myRadius = frodaAtom.at(atomID).surfaceRho + kappa;
  
  //first, get all the possibly relevant atoms from the monitor
  vector<SiteID> possibleInteractingAtoms;  
  possibleInteractingAtoms = areaProximity->getProximityList( atomID );

  Ai = 0.0;

  //get the rho, r, a parameters from each atom
  for ( unsigned int which = 0; which < possibleInteractingAtoms.size(); which++ ) {
    unsigned int other = possibleInteractingAtoms[which];
    double otherRadius = frodaAtom.at(other).surfaceRho;
    double rho2 = myRadius + otherRadius;
    rho2 *= rho2;

    Vector bond = currentPos.at(other);
    bond.x -= currentPos.at(atomID).x;
    bond.y -= currentPos.at(atomID).y;
    bond.z -= currentPos.at(atomID).z;
    double dist2 = dotProduct( bond, bond );

    double avalue = rho2 - dist2;
    if ( avalue > 0.0 ) { 
      Ai += avalue;
    }

  }

  Ai = pow ( Ai, 0.25 ); 

  Area = returnAreaCalculation( frodaAtom.at(atomID).surfaceCk[0], 
                                frodaAtom.at(atomID).surfaceCk[1], 
                                frodaAtom.at(atomID).surfaceCk[2], 
                                frodaAtom.at(atomID).surfaceCk[3], 
                                frodaAtom.at(atomID).surfaceCk[4], Ai );

  frodaAtom.at(atomID).SASA = Area;

}


////////////////////////////////////////////////////////////////////////////////
// Description: calculate SASA of the conformer, using Brooks/Guvench03 method
////////////////////////////////////////////////////////////////////////////////
void Froda::getAreaOfConformer() {
  areaProximity->update();

  totalSASA = 0;
  phobicSASA = 0;
  polarSASA = 0;
  polarPairSA = 0;  

  for ( unsigned int atom = 1; atom<=nTotalSites; atom++ ) {
    getAreaOfAtom( atom );
    totalSASA += frodaAtom.at(atom).SASA;
    if ( frodaAtom.at(atom).isPolar && frodaAtom.at(atom).isTerminal ) {
      if ( frodaAtom.at(atom).SASA > 0 ) {
        polarSASA += frodaAtom.at(atom).SASA;
      }
    }
    else {
      if ( frodaAtom.at(atom).SASA > 0 ) {
        phobicSASA += frodaAtom.at(atom).SASA;
      }
    }

  }

  //cerr << "Total area " << totalSASA;
  //cerr << "  Phobic area " << phobicSASA << "  Polar area " << polarSASA;
  //cerr << "  Paired polar area " << polarPairSA << endl;

}

double Froda::SBEnergy( SiteID donor, SiteID acceptor ) {
  double R_s = 3.2;
  double a = 0.375;
  double DA_distance;
  double ratioSquared;

  DA_distance = sqrt(getDist2( donor, acceptor ) );
  DA_distance += a;
  ratioSquared = R_s/DA_distance; 
  ratioSquared *= ratioSquared;
  double energy = SBWellDepth * ( 5 * pow( ratioSquared,6 ) - 6*pow( ratioSquared, 5 ) );

  return(energy);
}



double Froda::scoreHBond( SiteID hydrogen, SiteID acceptor, MolFramework &structure ) {
  SiteID donor = structure.site_info[hydrogen].neighbor_list[0];

  double energy = 0;

  //check if it's a salt bridge first
  if ( isSaltBridge( hydrogen, acceptor, structure ) ) {
    //rate it
    energy = SBEnergy( donor, acceptor );
    //cerr << "    Is salt bridge with energy " << energy << endl;
    return energy;
  }
  //now check if it's an hbond
  if ( !maybeHBond( hydrogen, acceptor, structure ) ) {
    //cerr << "    Can't be an hbond." << endl;
    return(0.0); //returning 0, no bond
  }

  double theta = computeAngle( donor, hydrogen, acceptor );

  //now we know it's OK for distance, theta
  int whatKind = HBType( hydrogen, acceptor, structure );
  //cerr << "    Potential hbond of kind " << whatKind;

  double phi;
  bool phiOK;

  scanPhi( hydrogen, acceptor, structure, whatKind, phi, phiOK );

  if ( !phiOK ) {
    //cerr << "    Unacceptable Phi value." << endl;
    return(0.0); //returning 0, no bond
  }

  //we might need gamma if it's type 1

  if ( whatKind == 1 ) {
    double gamma;
    bool itWorked = true;
    Vector donorNorm = getNormalVector( donor, itWorked, structure );
    if (!itWorked) {
      //cerr << "  PROBLEM with donor." << endl;
      return(0.0); //problem
    }
    Vector acceptorNorm;
    acceptorNorm = getNormalVector( acceptor, itWorked, structure );
    //acceptor doesn't give us a norm vector, use the base instead
    if (!itWorked){
      SiteID base = structure.site_info[acceptor].neighbor_list[0];
      acceptorNorm = getNormalVector( base, itWorked, structure );
    }
    if (!itWorked) {
      //cerr << "  PROBLEM with acceptor/base." << endl;
      return(0.0); //problem
    }
    double cosGamma = dotProduct( donorNorm, acceptorNorm );
    //remember that gamma is in the range pi to pi/2, not 0
    if ( cosGamma > 0 ) cosGamma *= -1;
    gamma = acos( cosGamma );
    if ( gamma > phi ) phi = gamma;// use gamma not phi
  }

  //we're now in position to rate the energy. First, the distance component 
  double R_s = 2.8;
  double DA_distance2;
  double ratioSquared;

  DA_distance2 = getDist2( donor, acceptor );
  ratioSquared = R_s*R_s/DA_distance2; 
  energy = HBWellDepth * ( 5 * pow( ratioSquared,6 ) - 6*pow( ratioSquared, 5 ) );
  //cerr << "  Distance function " << energy;  
 
  //thought: should I trap positive energies?

  //now the angular component prefactor
  double prefactor = pow( cos(theta), 2 ) * exp( -pow((PI - theta),6) ) ;

  //now the other angular terms, if any
  if (whatKind == 1) {
    //use best of phi or gamma; this best is already in phi
    prefactor *= pow( cos(phi), 2);
  }
  else if ( whatKind == 2) {
    prefactor *= pow( cos(phi), 2);
  }
  else if ( whatKind == 3) {
    prefactor *= pow( cos(theta), 2);
  }
  else if( whatKind == 4 ) {
    //uses the tetrahedral bonding
    phi -= (109.5/ RAD_TO_DEG );
    prefactor *= pow( cos(phi), 2);
  }
  energy *= prefactor; 
  //cerr << "  Final energy " << energy << endl;
  return ( energy );
}


////////////////////////////////////////////////////////////////////////////////
// Description: calculate and add up the MC scoring contributions
////////////////////////////////////////////////////////////////////////////////
void Froda::getMCScore( MolFramework &structure) {
    if (doHBScoring ) {
      findChargedInteractions( structure );
    }

    if ( doGetArea ) {
      getAreaOfConformer();
    }
    Rg = calculateRg();
    if ( doGetPhobicRg ) phobicRg = calculatePhobicRg();
    
    if ( mcRg ) {
      if ( minimumRg < 0 ) {
        RgScore = 0.5 * RgWeight * pow( ( Rg + minimumRg ), 2 ); //penalises large Rg
      } 
      else if ( Rg > minimumRg ) {
        RgScore = 0.5 * RgWeight * pow( ( Rg - minimumRg ), 2 ); //penalises large Rg
      }
    }
    if ( mcPhobicRg ) {
      if ( minimumPhobicRg < 0 ) {
        phobicRgScore = 0.5 * phobicRgWeight * pow( ( phobicRg + minimumPhobicRg ), 2 );
      } 
      else if ( phobicRg > minimumPhobicRg ) {
        phobicRgScore = 0.5 * phobicRgWeight * pow( ( phobicRg - minimumPhobicRg ), 2 );
      }
    }
    if ( mcSASA ) {
      if ( totalSASA > minimumSASA ) {
        SASAScore = SASAWeight * ( totalSASA - minimumSASA ); //penalises large exposed area 
      }
    }
    if ( mcPhobicSASA ) {
      if ( phobicSASA > minimumPhobicSASA ) {
        phobicSASAScore = phobicSASAWeight * ( phobicSASA - minimumPhobicSASA ); //penalises large exposed phobic area 
      }
    }
    if ( mcPolarSASA ) {
      //if ( > minimum ) {
        polarSASAScore = polarSASAWeight * polarSASA ; //penalises large exposed polar area 
        //probably use a negative weight on this one?
      //}
    }

    frodaMCScore = RgScore + phobicRgScore + SASAScore + phobicSASAScore + polarSASAScore + HBScore;

  return;
}

