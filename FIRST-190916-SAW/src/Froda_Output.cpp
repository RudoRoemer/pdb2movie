#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "Froda.h"
#include "mt19937ar.h"
#include "flexweb.h"
#include "hybrid_36_c.h"

#include <signal.h>

extern Parameters parameters;


////////////////////////////////////////////////////////////////////////////////
// Description: Outputting a conformer as pdb
////////////////////////////////////////////////////////////////////////////////
void Froda::outputConformer( MolFramework &structure, int i ){
  double runTime = static_cast<double>(timer.getutime())/1000.0;
  
  updateStatus();
  
  if ( parameters.report_resrmsd){
    stringstream residuefilename;
    residuefilename << structure.path_name << structure.base_name <<"_resrmsd_" << setfill('0') << setw(8) << i+1 << ".txt";
    ofstream resrmsd((residuefilename.str()).c_str() );

    double inconfs; //1 over nConfsFound, to normalise running array
    inconfs = 1.0/nConfsFound;

    int iterationWidth = (int) floor( log(i+1.0)/log(10.0) );
  
    resrmsd << "Residue # \t Iteration " << fixed << setw(iterationWidth) << (i+1) << endl;
  
    for ( unsigned int alpha = 0; alpha < alphaCarbon.size(); alpha++){
      int carbon = alphaCarbon[alpha];
      int resid = structure.site_info[carbon].seq_number;  

      resrmsd << fixed << setw(9) << resid << "\t ";
      resrmsd << fixed << setw(11 + iterationWidth) << setprecision(2) << sqrt ( inconfs* runningResidueRMSD.at(alpha) ) << endl;
    
    }

    resrmsd.close();
  }
 
  stringstream myfilepath;
  stringstream myfilename;
  
  myfilename << structure.base_name <<"_froda_" << setfill('0') << setw(8) << i+1 << ".pdb";
  myfilepath << structure.path_name << myfilename.str();
  
  // for flexweb - store the name of the current snapshot for the overlay 
  // TODO - refactor out to flexweb.cpp
  if (parameters.flexweb) { 
    int conformerNumber = i+1;
    
    if (parameters.numberOfFRODAwebOverlaySnapshots != 0) {
      int numberOfConformersPerOverlaySnapshot = nConfigs / parameters.numberOfFRODAwebOverlaySnapshots;
      
      int remainder = numberOfConformersPerOverlaySnapshot % outputEvery;
      
      // if numberOfConformersPerOverlaySnapshot isn't divisible by outputEvery,
      // choose the next 
      if (numberOfConformersPerOverlaySnapshot % outputEvery != 0) {
        numberOfConformersPerOverlaySnapshot -= remainder;
        
      }
      
      // ensure that numberOfConformersPerOverlaySnapshot != 0
      if (numberOfConformersPerOverlaySnapshot == 0) {
        numberOfConformersPerOverlaySnapshot = outputEvery;
        
      }
      
      if (conformerNumber % numberOfConformersPerOverlaySnapshot == 0) {
        string file_name = structure.path_name;
        file_name += ".pdbSnapshotList";
        ofstream currentSnapshotFile( file_name.c_str(), ios::app);
        currentSnapshotFile << endl << myfilename.str();
        currentSnapshotFile.close();
      }
    }
  }
  
  if( parameters.verbose ){
    cout <<  "Writing conformation to file " << myfilepath.str() << endl ;
  }
  
  stringstream new_conformer;
  new_conformer << "REMARK FRODA output" << endl;
  new_conformer << "REMARK initial structure from " << structure.base_name << endl;
  new_conformer << "REMARK FRODA conformer " << i+1 << endl;
  new_conformer << "REMARK GLOBAL RMS_D " << globalRMSD  ;
  new_conformer << "; MAIN " << allmain_RMSD;
  new_conformer << "; SIDE " << allside_RMSD;
  new_conformer << "; RUNTIME " << showpoint << setiosflags(ios::fixed) << setprecision(3) << runTime << endl;

  new_conformer << "REMARK MOBILE RMS_D " << RMSD ;
  new_conformer << "; MAIN " << mobmain_RMSD;
  new_conformer << "; SIDE " << mobside_RMSD;
  new_conformer << "; RUNTIME " << runTime << endl;
  new_conformer << "REMARK Rg= " << Rg << endl;
  if ( doGetPhobicRg ) {
    new_conformer << "REMARK Phobic Rg= " << phobicRg << endl;
  }
  if ( doGetArea ) {
    new_conformer << "REMARK PHOBIC AREA " << phobicSASA << " POLAR AREA " << polarSASA << " TOTAL AREA " << totalSASA << endl;
  }

  if (isRestart){
    new_conformer << "REMARK RMS-D FROM START: " << globalRMSD << "; FROM RESTART: " << restartRMSD << endl;
  }
  if (isTargeted){
    new_conformer << "REMARK RMS-D FROM START: " << globalRMSD << "; TO TARGET: " << targetRMSD << endl;
    new_conformer << "REMARK TARGETED MAINCHAIN RMS-D FROM START: " << allTargetedMainchainRMSD;
    new_conformer << "; TO TARGET: " << targmain_RMSD ;
    new_conformer << "; RUNTIME " << runTime << endl;

    double pathX;
    double pathH;
    double sigmaI = allTargetedMainchainRMSD;
    double sigmaF = targmain_RMSD; 
    
    sigmaI = currentRMSDfromInitial;
    sigmaF = currentRMSDfromTarget;
            
    pathX = ( 0.5/ initialToTargetL ) * ( pow( initialToTargetL, 2 ) + pow( sigmaI, 2 ) - pow( sigmaF, 2 ) );
    pathH = sqrt( pow(sigmaI, 2 ) - pow( pathX, 2 ) );

    new_conformer << "REMARK TARGETED PATHWAY x= " << pathX << "; h= " << pathH << endl;
  }
  
  if (checkRama){
    new_conformer << "REMARK RAMACHANDRAN " << nBadRama << " " << nGenerousRama << " " << nAllowedRama << " " << nCoreRama << endl;
    //provide detailed breakdown of rama by ramaList
    new_conformer << "REMARK RAMA CA PHI PSI MAP F/R TYPE" << endl;
    for ( unsigned int whichRama = 0; whichRama < ramaList.size(); whichRama++ ) {
      new_conformer << "REMARK RAMA " << ramaList.at(whichRama).Calpha << " ";
      new_conformer << ramaList.at(whichRama).phi0  << " ";
      new_conformer << ramaList.at(whichRama).psi0  << " ";

      //new JULY 06: more detail for rama data
      new_conformer << ramaList.at(whichRama).mapValue << " ";
      //label locked or variable
      if ( ramaList.at(whichRama).freePhi || 
           ramaList.at(whichRama).freePsi ) {
        new_conformer << "F "; //variable angle, "Free"
      }
      else {
        new_conformer << "R "; //locked angle, "Rigid"
      }
      //label glycine and proline
      if ( ramaList.at(whichRama).isPro ) {
        new_conformer << "PRO"; //proline residue
      }
      else if ( ramaList.at(whichRama).isGly ) {
        new_conformer << "GLY"; //proline residue
      }
      new_conformer << endl;
   
    }

    
  }

  if ( doScoringMC ){
    new_conformer << "REMARK MC SCORES: Rg " << RgScore << ", ";
    new_conformer << " Phobic Rg " << phobicRgScore << ", ";
    new_conformer << " SASA " << SASAScore << ", ";
    new_conformer << " Phobic SASA " << phobicSASAScore << ", ";
    new_conformer << " Polar SASA " << polarSASAScore << ": ";
    new_conformer << " HBonds " << HBScore << ": ";
    new_conformer << " TOTAL " << frodaMCScore << endl;
  }


  if (!nosterics){
    new_conformer << "REMARK CLASH BINS 0-0.1 0.1-0.2 0.2-0.3 0.3-0.4 0.4-0.5 0.5-0.6 0.6+ ANGSTROMS" << endl;
    new_conformer << "REMARK CLASH COUNT";
    for ( int bin=0; bin < 7; bin++){
      new_conformer << " " << overlapDegrees.at(bin);
    }
    new_conformer << endl;
  }
  

  if ( !froda2Hybrid ) {

    new_conformer << "REMARK WORST GHOST MISMATCH: " << worstGhostMismatch;
    new_conformer << " ; ATOM ";
    new_conformer << structure.site_info[worstGhostAtom].orig_atom_number;
    new_conformer << " ; GHOST " << worstGhost << endl;
    if ( Phobes.size() > 0 ) {
      new_conformer << "REMARK WORST TETHER STRETCH: " << worstTetherStretch;
      new_conformer << " ; ATOM1 ";
      new_conformer << structure.site_info[worstStretchAtom1].orig_atom_number;
      new_conformer << " ; ATOM2 ";
      new_conformer << structure.site_info[worstStretchAtom2].orig_atom_number;
      new_conformer << endl;
    }
    if ( !nosterics ) {
      new_conformer << "REMARK WORST CLASH: " << worstClash;
      new_conformer << " ; ATOM1 ";
      new_conformer << structure.site_info[worstClashAtom1].orig_atom_number;
      new_conformer << " ; ATOM2 ";
      new_conformer << structure.site_info[worstClashAtom2].orig_atom_number;
      new_conformer << endl;
    }
  }

  if ( usePairs ) {
    for ( unsigned int whichP = 0; whichP < myPairs.size(); whichP ++ ) {
      new_conformer << "REMARK PAIR ";
      new_conformer << structure.site_info[ myPairs[whichP].atomA ].orig_atom_number << " ";
      new_conformer << structure.site_info[ myPairs[whichP].atomB ].orig_atom_number << " ";
      new_conformer << myPairs[whichP].ABcurrent << " / ";
      new_conformer << myPairs[whichP].AB << " ";
      if ( myPairs[whichP].tagGood ) {
        new_conformer << "YES" << endl;
      }
      else{
        new_conformer << "NO" << endl;
      }
    }
  }

  float rigid_label  = 0.0;
  float coll_mode_label = 0.0;

  ostringstream FIRST_group_ID, atom_number, seg_id;
  string output_element_name;

  string temp1;
  const char* errmsg;
  int int_atom_number;
  char write_atom_number[6];

  for( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ){
    structure.site_info[siteNumber].rigid_label < 999 ? rigid_label = structure.site_info[siteNumber].rigid_label : rigid_label = 999.0;
    structure.site_info[siteNumber].coll_mode_label < 999 ? coll_mode_label = structure.site_info[siteNumber].coll_mode_label : coll_mode_label = 999.0;

    structure.site_info[siteNumber].coords[0] = currentPos.at(siteNumber).x;
    structure.site_info[siteNumber].coords[1] = currentPos.at(siteNumber).y;
    structure.site_info[siteNumber].coords[2] = currentPos.at(siteNumber).z;

    output_element_name = structure.site_info[siteNumber].element_name;
    if(output_element_name[1] == ' ') swap(output_element_name[0], output_element_name[1]);
    
    FIRST_group_ID.str("");
    if (use_group_id) {
      FIRST_group_ID << setw(6) << structure.site_info[siteNumber].FIRST_group_ID;
    }
        
    atom_number.str("");
    temp1 =   structure.site_info[siteNumber].orig_atom_number;
    if( isNumber( temp1))
      {
	int_atom_number = atoi( temp1.c_str());
	if(  int_atom_number <= 99999){
	  atom_number << setw(5) << structure.site_info[siteNumber].orig_atom_number;
	}
	else{
	  hy36encode(5, int_atom_number, write_atom_number);
	  atom_number << setw(5) <<  write_atom_number;
	}
      }
    else{
	atom_number << setw(5) << structure.site_info[siteNumber].orig_atom_number;
    }
    

    seg_id.str("");
    seg_id << setiosflags(ios_base::left)
    << setw(4) << structure.site_info[siteNumber].seg_id;
      
    new_conformer  << showpoint << setiosflags(ios::fixed) 
    << setw(6) << structure.site_info[siteNumber].record_name
    << atom_number.str()
    << setw(5) << structure.atomNamePDBFormat(structure.site_info[siteNumber].atom_name, structure.site_info[siteNumber].element_name)
    << setw(4) << structure.site_info[siteNumber].residue_name
    << setw(2) << char(structure.site_info[siteNumber].chain_ID)
    << setw(4) << structure.site_info[siteNumber].seq_number
    << setw(1) << structure.site_info[siteNumber].insert_code
    << "   "
    << setw(8) << setprecision(3) << structure.site_info[siteNumber].coords[X]
    << setw(8) << setprecision(3) << structure.site_info[siteNumber].coords[Y]
    << setw(8) << setprecision(3) << structure.site_info[siteNumber].coords[Z]
    << setw(6) << setprecision(2) << coll_mode_label
    << setw(6) << setprecision(2) << rigid_label // cols 61-66
    << "      " // cols 67-72
    << seg_id.str() // cols 73-76
    << setw(2) << output_element_name // cols 77-78
    << "   " // cols 79-81
    << FIRST_group_ID.str()
    << endl;
    if( find(structure.chainTermini.begin(), structure.chainTermini.end(), siteNumber) != structure.chainTermini.end() ){
      new_conformer << "TER" << endl;
    }

  }
  structure.outputCONECTRecords( new_conformer );
  new_conformer << "END" << endl;
  
  ofstream new_conformer_file((myfilepath.str()).c_str() );
  new_conformer_file << new_conformer.str();
  new_conformer_file.close();

}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description: Outputting the ghost member lists
////////////////////////////////////////////////////////////////////////////////
void Froda::outputGhosts(MolFramework &structure ){
  cerr << "Outputting ghost member lists." << endl;

  unsigned int atomsLeft;
  unsigned int atomsThisLine;
  unsigned int currentAtom;
  unsigned int atomsPerLine = 10;

  ofstream gfile;
  string filename = structure.path_name + "g.out";
  gfile.open( filename.c_str() );

  for ( unsigned int whichGhost = 1; whichGhost <= myNClusters; whichGhost++ ) {
    atomsLeft = ghost[whichGhost].atoms.size();
    bool stillWriting = true;
    currentAtom = 0;

    if (atomsLeft == 0 ) {
      gfile << "GHOST " << whichGhost << " GEN " << endl;
      //empty ghost
      stillWriting = false;
    }

    while (stillWriting) {
      if ( atomsLeft > atomsPerLine ) {
        atomsThisLine = atomsPerLine;
        gfile << "GHOST " << whichGhost << " GEN ";
        for ( unsigned int mycount = 0; mycount < atomsThisLine; mycount++ ) {
          //atom_number << setw(5) << structure.site_info[siteNumber].orig_atom_number;
          string thisAtom = structure.site_info[ ghost[whichGhost].atoms.at(currentAtom) ].orig_atom_number;
          gfile << thisAtom << " ";
          atomsLeft--;
          currentAtom++;
        }
        gfile << endl;
      }
      else if ( atomsLeft == 0 ) {
        stillWriting = false;
      }
      else {
        atomsThisLine = atomsLeft;
        gfile << "GHOST " << whichGhost << " GEN ";
        for ( unsigned int mycount = 0; mycount < atomsThisLine; mycount++ ) {
          string thisAtom = structure.site_info[ ghost[whichGhost].atoms.at(currentAtom) ].orig_atom_number;
          gfile << thisAtom << " ";
          atomsLeft--;
          currentAtom++;
        }
        gfile << endl;
        stillWriting = false;
      }
    }
  }

  gfile.close();
  return;
}
