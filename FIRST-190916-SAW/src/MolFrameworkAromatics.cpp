#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "MolFramework.h"
#include "VectorAlgebra.h"

extern const Parameters parameters;
inline float rad2deg( float a ){
  return(a*(180.0/PI));
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//     define the coordinates of the center of a given ring. 
// Parameters: 
//     ring_1 - aromatic ring (structure) under consideration.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::define_ring_center( aromatic_ring &ring_1) {

  vector<unsigned int>::iterator aromatic_atom = ( ring_1.ring_atoms ).begin();

  int rsize = (ring_1.ring_atoms).size();
  int uatom=0;
  int vatom =0;
  string aname = site_info[*aromatic_atom].atom_name;
  //  if(parameters.verbose) {
  // cout<<" "<<ring_1.residue_name<<" "<<ring_1.residue_number<<"\t"
  //  <<*aromatic_atom<<"   "<<aname<<"      "<< rsize<<endl;
  //}
  //////////////////////////////////////////////////////////////////////////////////////////

  while( aromatic_atom != (ring_1.ring_atoms ).end() ){
    string aname = site_info[*aromatic_atom].atom_name;
    // find geometric center
    // Consider only heavy atoms (nitrogens, oxygens & carbons), but consider them 
    // equivalently if aromatic has multiple rings, split into proper sub rings:
    //  i.e if 9-fold (TRP,G,A) = 5+6; if 12-fold = 5+6+5 ( YG).
    if(rsize==9) {
      if(site_info[*aromatic_atom].residue_name=="TRP") {
	if(aname=="CG" || aname=="CD1" || aname=="CE1") 
	  for(int a=0; a<3; a++) 
	    ring_1.second_cntr[a] += site_info[*aromatic_atom].coords[a];      
	else if(aname =="CD2" || aname =="CE2")
	  for(int a=0; a<3; a++) {
	    ring_1.geom_center[a] += site_info[*aromatic_atom].coords[a];      
	    ring_1.second_cntr[a] += site_info[*aromatic_atom].coords[a];      
	  }
	else 
	  for(int a=0; a<3; a++) 
	    ring_1.geom_center[a] += site_info[*aromatic_atom].coords[a];     
      }
      else {
	if(aname =="N9" || aname =="C8" || aname =="N7") 
	  for(int a=0; a<3; a++) 
	    ring_1.second_cntr[a] += site_info[*aromatic_atom].coords[a];
	
	else if(aname =="C4" || aname =="C5") 
	  for(int a=0; a<3; a++) {
	    ring_1.geom_center[a] += site_info[*aromatic_atom].coords[a];      
	    ring_1.second_cntr[a] += site_info[*aromatic_atom].coords[a];      
	  }
	
	else 
	  for(int a=0; a<3; a++) 
	    ring_1.geom_center[a] += site_info[*aromatic_atom].coords[a];     
      }
    }
    else { // set up to take the geometric center (average) of entire ring complex for larger residues
      for(int a=0; a<3; a++) 
	ring_1.geom_center[a] += site_info[*aromatic_atom].coords[a];
      /////////////////////////////////////////////////////////////////////////////////////////
      // Find the Local reference frame for each ring.
      //
      // Based on the work of Gabb, et.al. J. Mol. Graph. 1996 and utilized by 
      //        Gendron, et.al. JMB 2001.  
      // I have implemented the scheme such that only if it is a multi-ringed aromatic are the
      //   separate centers of each ring searched. The local reference frame involves the unit 
      //   vectors u and v which are between the anchor atoms: N1 in pyrimidines (CTU); N9 in 
      //   purines (AG); CG in amino acid sidechains; CA in amino acid mainchain and their 
      //   neighbors. These neighbors are identified as neighbor_atom1 & neighbor_atom2: 
      //   C2 & C6 in (CTU); C4 & C8 in (AG); CD2 & CD1 in sidechains; CB & N in mainchains
      //   in order to obtain u & v. 
      //
      //   Note that in the cases where there are multiple joined rings, .... 
      /////////////////////////////////////////////////////////////////////////////////////////
    }
    // identify the plane of the ring and the neighboring atoms (neighbor_atom1 & neighbor_atom2)
    if( ring_1.anchor_atom == *aromatic_atom ) {
      for(unsigned int neighborNumber=0; neighborNumber< (site_info[*aromatic_atom].neighbor_list).size(); neighborNumber++) {
	unsigned int atom_2 = site_info[*aromatic_atom].neighbor_list[neighborNumber];
	if(atom_2 != ring_1.external_atom && !is_ribose_atom(atom_2)) {
	  if(uatom && !vatom) vatom = atom_2;
	  if(!uatom) uatom = atom_2;
	}
      //cout<<"Anchor: "<<*aromatic_atom<<"-----"<<site_info[*aromatic_atom].neighbor_list[neighborNumber]<<" "
      //<<site_info[site_info[*aromatic_atom].neighbor_list[neighborNumber]].atom_name<<"||\t";
      } 
      //      cout<<endl;
//   for( int a = 0; a < (site_info[atom_1].neighbor_list).size(); a++ ) {
      //cout<<"---------------"<<*aromatic_atom<<"----------------"<<rsize<<endl;
      // for the case of known nucleic and amino acids.
      //////////////////////////////////////////////////////////	
//       if(rsize == 9) {
// 	if(ring_1.residue_name =="TRP") {
// 	  uatom = NthNearestNeighbor(*aromatic_atom,"CD2",1);
// 	  vatom = NthNearestNeighbor(*aromatic_atom,"CD1",1);
// 	}
// 	else {
// 	  uatom = NthNearestNeighbor(*aromatic_atom,"C4",1);
// 	  vatom = NthNearestNeighbor(*aromatic_atom,"C8",1);
// 	}
//       }
//       else if( rsize ==6) {
// 	if(ring_1.residue_name == "TYR" || ring_1.residue_name == "PHE" ) {
// 	  uatom = NthNearestNeighbor(*aromatic_atom,"CD1",1);
// 	  vatom = NthNearestNeighbor(*aromatic_atom,"CD2",1);
// 	}
// 	else {
// 	  uatom = NthNearestNeighbor(*aromatic_atom,"C2",1);
// 	  vatom = NthNearestNeighbor(*aromatic_atom,"C6",1);
// 	}
//       }
//       else if(rsize==5) {
// 	if(ring_1.residue_name =="PRO") {
// 	  uatom = NthNearestNeighbor(*aromatic_atom,"CB",1);
// 	  vatom = NthNearestNeighbor(*aromatic_atom,"N",1);	  
// 	}
// 	else if(ring_1.residue_name =="HIS") {
// 	  uatom = NthNearestNeighbor(*aromatic_atom,"CD2",1);
// 	  vatom = NthNearestNeighbor(*aromatic_atom,"ND1",1);
// 	}
//       }
//       else if(rsize==12) { // assuming a 5-6-5 configuration 
// 	cout<<"ring size of  "<< rsize<<" encountered. Check that this base "
// 	    <<ring_1.residue_name<<" is correctly connected."<<endl;
// 	uatom = NthNearestNeighbor(*aromatic_atom,"C4",1);
// 	vatom = NthNearestNeighbor(*aromatic_atom,"C8",1);
// 	if(uatom==0 ||vatom==0) {
// 	  uatom = NthNearestNeighbor(*aromatic_atom,"CD2",1);
// 	  vatom = NthNearestNeighbor(*aromatic_atom,"CD1",1);
// 	}
//       }
//       else {
// 	cout<<"ring size of  "<< rsize<<" is not yet defined. Add the proper base: "
// 	    <<ring_1.residue_name<<"."<<endl;
// 	while( aromatic_atom != (ring_1.ring_atoms ).end() ){
// 	  cout<<"\t Contains ring atoms  "<<*aromatic_atom<< "  "
// 	      <<site_info[*aromatic_atom].atom_name<<endl;
// 	  aromatic_atom++;
// 	} 
// 	return; // skips larger size rings (for now) 12.03.05 AJR
//       }
//      if(uatom) ring_1.neighbor_atom1 = uatom;
//       else {
// 	cout<<"problem with ring_1: "<< *aromatic_atom <<" "<<site_info[*aromatic_atom].atom_name
// 	    << " "<<ring_1.residue_name<<" "<<rsize<<endl;
// 	int i=(site_info[*aromatic_atom].atom_name).find("N9");
// 	int j=(site_info[*aromatic_atom].atom_name).find("N1");
// 	cout<<"\t"<< i <<" "<<j<<endl;
// 	string aname2;
// 	if(i!=string::npos) {
// 	  aname2.assign(site_info[*aromatic_atom].atom_name,0,i);
// 	  aname2.append("C4");
// 	  uatom = NthNearestNeighbor(*aromatic_atom,aname2,1);
// 	  aname2.replace(i,aname2.size(),"C8");
// 	  vatom = NthNearestNeighbor(*aromatic_atom,aname2,1);
// 	}
// 	else if(j!=string::npos) {
// 	  aname2.assign(site_info[*aromatic_atom].atom_name,0,i);
// 	  aname2.append("C2");
// 	  uatom = NthNearestNeighbor(*aromatic_atom,aname2,1);
// 	  aname2.replace(i,aname2.size(),"C6");
// 	  vatom = NthNearestNeighbor(*aromatic_atom,aname2,1);
// 	}
//       }
      //cout<<ring_1.residue_name<<" "<<rsize<<"\t"<<*aromatic_atom<<"   "<<uatom<<" "<<vatom;//<<endl;
      if(uatom&&vatom) {
	ring_1.neighbor_atom2 = vatom;
	ring_1.neighbor_atom1 = uatom;
      }
      else {
	cout<<"problem with ring_1: "<< *aromatic_atom<<" "<<site_info[*aromatic_atom].atom_name
 	    << " "<<ring_1.residue_name<<" "<<rsize<<endl;
	return;
      }
    }
    aromatic_atom++;   
  }
  for(int a=0; a<3; a++) {
    if(rsize!=9) {
      ring_1.geom_center[a] = ring_1.geom_center[a]/rsize;
      //cout<<"\t"<<a<<"======="<<ring_1.geom_center[a];
    }
    else if(rsize==9) {
      ring_1.geom_center[a] = ring_1.geom_center[a]/6.0;
      ring_1.second_cntr[a] = ring_1.second_cntr[a]/5.0;
      //cout<<"\t"<<a<<"======="<<ring_1.geom_center[a];
      //cout<<"\t"<<a<<"======="<<ring_1.second_cntr[a];
    }
  }
  if(uatom && vatom)
    getUnitNormalToPlane( ring_1.ring_norm, ring_1.anchor_atom, uatom, vatom);

  //  cout<<endl;
  //cout<<"\t"<<ring_1.anchor_atom<<": "<<ring_1.external_atom<<"\t["<<uatom<<","
  //<<vatom<<"]"<<endl;
  // cout<<"+++++"<<ring_1.ring_norm[0]<<"+++++"<<ring_1.ring_norm[1]<<"+++++"
  //<<ring_1.ring_norm[2]<<"+++++"<<endl;
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//     Return the coordinates of the center of a given ring. This function is for 
//     use once the ring_centers have been defined.
// Parameters: 
//     atom_1 - FIRST_number of a given atom
////////////////////////////////////////////////////////////////////////////////
void MolFramework::get_ring_center( float *center_coord, int ringsize, int atom_1, int atom_2, int atom_3) {

  //  int ringsize = 0;
  //  string neighbor_element;
  string aname = site_info[atom_1].atom_name;

// Find geometric center of the ring or ring complex
// Consider only heavy atoms (nitrogens, oxygens & carbons, etc?) but consider them equivalently
// if aromatic has multiple rings, split into proper ring. (i.e if 9-fold = 5+6, TRP, G, A.
//-----------------------------------------------------------------------------------------------
  if(ringsize==9) {
    if(site_info[atom_1].residue_name=="TRP") {
      if(aname=="CG" || aname=="CD1" || aname=="CE1") {
      }
    cout << "ringsize = " << ringsize << " isn't working yet. "<<endl;
    }
  }
  if(ringsize==6) {
    cout << "ringsize = " << ringsize << " isn't working yet. "<<endl;
  }
  

  return;
}

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Stacked rings can be identified in several ways. These rings include the 
//   base stacking in nucleotides and the interaction of other aromatic rings
//   such as TRP side chains. Using the parameter ringstack_fxn (not yet functional), 
//   one can easily adjust the way these are modeled.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::identify_stacked_rings() {

  if( parameters.verbose ) 
    cout << "Identifying stacked ring interactions... " << endl; 

  // 1. Specify the type of user-selected ring_stacking interation parameter to use.
  ////////////////////////////////////////////////////////////////////////////////
  //    a.  If the user selected "-noaromatics", then they don't want to include
  //        stacked ring interactions. Return to running the file. 
  ////////////////////////////////////////////////////////////////////////////////
  if( parameters.skip_identify_aromatics ) 
    return;

  ////////////////////////////////////////////////////////////////////////////////
  //    b.  If the user selected "-srin", then they want to include only the 
  //        specified stacked ring interactions. Read and return to running the file. 
  ////////////////////////////////////////////////////////////////////////////////
  if( parameters.srin ){
    read_ringstack_file();
    if( parameters.verbose )
      cout << " Reading stacked ring interactions... " << endl; 
    return;
  }
  
  // 2. Define a vector of struct aromatic_ring
  ////////////////////////////////////////////////////////////////////////////////
  //   a. identify aromatic_ring residues & place atoms into the proper rings.
  ////////////////////////////////////////////////////////////////////////////////
  int ring_count=0;
  int hold_resid=-99;
  int nres = total_residues+5; // TODO: 11.11.05 improve this allocation.
  // note this creates more than enough -- must check size for non-ring residues.
  struct aromatic_ring ring_info[nres]; 

  // reset_checked_label
  //  structure.reset_checked_label; // insures that all such values are tagged false.
  //////////////////////////////////////////////////////////////////////////
  // search over all atoms, but only place those that are in aromatic rings
  // into the aromatic_ring structure: ring_info, defined for each residue.
  ////////////////////////////////////////////////////////////////////////////
  for( SiteID current_atom = 1; current_atom <= total_sites; current_atom++ ) {

    if( is_ring_atom(current_atom) ){

      stringstream residue_id_string;
      residue_id_string << site_info[current_atom].seq_number << ";"
			<< site_info[current_atom].insert_code << ";"
			<< site_info[current_atom].chain_ID; 

      int current_resid = unique_res_id[residue_id_string.str()];
      int current_residue  = site_info[current_atom].seq_number;
      int current_chain_ID = site_info[current_atom].FIRST_chain_ID;

      if( current_residue != hold_resid ){ 
	ring_count++;

	// self-check to make sure no overfill errors
	//////////////////////////////////////////////////////////////////////
	if( ring_count > nres) {
	  cout<<"Error: Number of ringed residues exceeds the total number of residues (" << nres << ")" << endl;
	  exit(1);
	}

	//cout << " RRA: " << current_atom << " " << hold_resid << endl;
	//ring_info.push_back( aromatic_ring );
	//ring_info.push_back( aromatic_ring(ring_count, current_residue, current_chain_ID) ); 
	ring_info[ring_count].residue_number = current_residue;
	ring_info[ring_count].chain_ID = current_chain_ID;
	ring_info[ring_count].residue_name = site_info[current_atom].residue_name;
	ring_info[ring_count].resid = current_resid;
	hold_resid = current_residue;
      }

      ring_info[ring_count].ring_atoms.push_back(current_atom);
  
      int extra_atom = is_ring_anchor(current_atom);
      if( extra_atom ) {
	ring_info[ring_count].anchor_atom = current_atom;
	ring_info[ring_count].external_atom = extra_atom;
	//cout<<" RRR "<< ring_residues<< " "<<ring_info[ring_residues].anchor_atom<<":: "
	//<<extra_atom<< "  "<<site_info[extra_atom].atom_name<<endl;
      }
    }
  }

  if( parameters.verbose >= 2 ) 
    cout << ring_count << " aromatic_ring_residues found" << endl;

  // b. Find Ring geometric centers and defining vectors for each aromatic ring
  ////////////////////////////////////////////////////////////////////////////////
  for( int iring=1; iring<=ring_count; iring++) {
    //vector<int>::iterator aromatic_atom = ( ring_info[iring].ring_atoms ).begin();
    // int rsize = (ring_info[iring].ring_atoms).size(); // FIXME - warning: unused variable 'rsize' - does the call to size() have side-effects that we're depending on?
    //cout<<iring<<" "<<ring_info[iring].residue_name<<" "<<ring_info[iring].residue_number<<"  "
    //<<rsize<<"\t"<<ring_info[iring].anchor_atom<<"  "<<site_info[ring_info[iring].anchor_atom].atom_name<<" "<<*aromatic_atom<<endl;
    //<<*aromatic_atom<<"  "<<"       "<< rsize<<endl;
     define_ring_center( ring_info[iring]);
     /* cout<<iring<<"\tx\ty\tz"<<endl;
	cout<<"  "<<ring_info[iring].geom_center[0]<<","<<ring_info[iring].geom_center[1]<<","<<ring_info[iring].geom_center[2]<<endl;
	cout<<"\t"<<ring_info[iring].second_cntr[0]<<","<<ring_info[iring].second_cntr[1]<<","<<ring_info[iring].second_cntr[2]<<endl;
	cout<<"\t"<<ring_info[iring].ring_norm[0]<<","<<ring_info[iring].ring_norm[1]<<","<<ring_info[iring].ring_norm[2]<<endl;
     */
  } 
  // 3. Search for stacked aromatic rings (only over known aromatic ring atoms)
  //    To improve speed: only check anchor atoms from each ring.
  ////////////////////////////////////////////////////////////////////////////////
   vector<int> grids_to_check;
   for( int current_ring = 1; current_ring <= ring_count; current_ring++) {
     //
     if(ring_info[current_ring].neighbor_atom1) { // insures neighbor_atoms have been recognized
     vector<unsigned int>::iterator current_atom = (ring_info[current_ring].ring_atoms).begin();
    //  for( int current_atom = 1; current_atom <= total_sites; current_atom++ ){
    while( current_atom != (ring_info[current_ring].ring_atoms).end()) {    
      //if( *current_atom == ring_info[current_ring].anchor_atom) { // only searches over anchor atoms.
      //    if( is_ring_atom(current_atom) && !is_ribose_atom(current_atom)) { 
      // a. Create a list of coordinate grids to check. Check each grid. 
      //////////////////////////////////////////////////////////////////////
      grids_to_check = getGrids( *current_atom );
      current_grid = grids_to_check.begin();
      while( current_grid != grids_to_check.end() ){
	// b. Check each atom in the current grid. --> faster: just check anchor atoms.
	//////////////////////////////////////////////////////////////////////
	neighbor_atom = coordinate_grid[*current_grid].begin();
	while( neighbor_atom != coordinate_grid[*current_grid].end() ){
	  // c. Check the distance between current_atom and neighbor_atom. Only 
	  //    check atoms with larger site numbers, to prevent double counting. 
	  //////////////////////////////////////////////////////////////////////
	  if( *neighbor_atom > *current_atom &&
	      is_ring_atom(*neighbor_atom) && !is_ribose_atom(*neighbor_atom) &&
	      isDifferentResidue( *current_atom, *neighbor_atom) ) {
	    //	      && (getDistance(*current_atom, *neighbor_atom ) < 6.0) )  { // ajr TODO: check distance
	    // compare the centers to see if they are stacked.
	    // criteria: A) distances between ring centers is < 5.5 Ang. 
	    //           B) angle between 2 normals to base planes < 30 deg.
	    //           C) angle between normal of one base plane & vector to ring center < 40 deg.
	    // CHECK if already a neighbor.
	    int chk_ring = get_ringid(*neighbor_atom, ring_count, ring_info);
	    //	    int chk_ring = get_ringid(*neighbor_atom, site_info[*neighbor_atom].seq_number, ring_count, ring_info);
	    //cout<<"  RINGS:   "<<current_ring<<" "<<chk_ring<<endl;
	    if( current_ring == chk_ring ){
	      cout << " Warning: i to i ring interaction is being skipped. Probable bug in unique residue ID." << endl;
	      break;
	    }
	    if( *neighbor_atom == ring_info[chk_ring].anchor_atom) {
	      if( !is_stacked_ring_neighbor(chk_ring, ring_info[current_ring]) &&
		  !is_stacked_ring_neighbor(current_ring, ring_info[chk_ring]) 
		  && ring_info[chk_ring].neighbor_atom1 ) {
	   //cout<<"CHECKING ring_neighbors: "<< chk_ring << " "<<ring_info[chk_ring].residue_name
	    //<<" && "<< current_ring<<" "<<ring_info[current_ring].residue_name<< "\t "
	   //<<ring_info[chk_ring].neighbor_atom1 <<">> "<<*current_atom<<" "<<*neighbor_atom<<endl;
	   //cout<<"Arpchk::: "<<*current_atom<<" "<<*neighbor_atom<<"\t"<<ring_info[current_ring].residue_name<<ring_info[current_ring].residue_number;
		//cout<<"\t"<<ring_info[chk_ring].residue_name<<ring_info[chk_ring].residue_number<<endl;
		if(is_stacked_ring_pair(ring_info[current_ring],ring_info[chk_ring]) ) {
		  (ring_info[current_ring]).stacked_neighbor_list.push_back(chk_ring);
		  (ring_info[chk_ring]).stacked_neighbor_list.push_back(current_ring);
		}     
	      }
	    }
	  }
	  neighbor_atom++;
	}
	current_grid++;
      }
      current_atom++;
    }
    //}
     }
   
   }
}
////////////////////////////////////////////////////////////////////////////////
    

// ////////////////////////////////////////////////////////////////////////////////
// // Description: 
// //   Create a list of atoms in rings
// ////////////////////////////////////////////////////////////////////////////////
void MolFramework::list_ring_atoms2() {
   
  for( unsigned int current_atom =1; current_atom <= total_sites; current_atom++ ) { 
    string current_rname = site_info[current_atom].residue_name;
    if(!site_info[current_atom].checked)
      // Deal with known cases first
      // omit Hydrogens && Water  
      if( site_info[current_atom].element_name=="H " ||current_rname =="HOH" ) 
	site_info[current_atom].checked =true;
    // Handle Proline  
      else if( current_rname == "PRO" ) {
	if(site_info[current_atom].atom_name != "C" && site_info[current_atom].element_name != "O " )
	  ring_atom_list.push_back( current_atom );
	site_info[current_atom].checked=true;
      }
    // omit mainchain atoms in standard AA residues
      else if(isMainchain(current_atom)) site_info[current_atom].checked=true;
    
    // handle rings in standard amino acids:
      else if(current_rname == "HIS" || current_rname == "PHE" ||
	      current_rname == "TRP" || current_rname == "TYR" ) {
	if( site_info[current_atom].atom_name != "CB" && site_info[current_atom].atom_name == "OH" 
	    && site_info[current_atom].atom_name != "OXT" ) 
	  ring_atom_list.push_back( current_atom );
	site_info[current_atom].checked=true;
      }
      else if(current_rname == "ALA" || current_rname == "ARG" || 
	      current_rname == "ASN" || current_rname == "ASP" ||
	      current_rname == "CYS" || current_rname == "GLN" ||
	      current_rname == "GLU" || current_rname == "GLY" ||
	      current_rname == "ILE" || current_rname == "LEU" ||
	      current_rname == "LYS" || current_rname == "MET" ||
	      current_rname == "SER" || current_rname == "THR" ||
	      current_rname == "VAL" ) 
	site_info[current_atom].checked=true;
    // handle rings in standard nucleic acids:
      else if(current_rname == "  A" || current_rname == "  C" ||
	      current_rname == "  T" || current_rname == "  U" ||
	      current_rname == "  G" ) {
	if(!is_ribose_atom(current_atom) && !is_nucleic_backbone(current_atom)){
	  string element = site_info[current_atom].element_name;
	  int nbonds=0;
	  for( unsigned int neighborNumber=0; neighborNumber<(site_info[current_atom].neighbor_list).size(); neighborNumber++) 
	    if( site_info[(site_info[current_atom].neighbor_list[neighborNumber])].element_name!="H ") nbonds++;
	  if(nbonds>=2)  ring_atom_list.push_back( current_atom );
	}
	site_info[current_atom].checked = true;
      }
    // now deal with unknown residue types:
      else {
	int current_residue  = site_info[current_atom].seq_number;
	int current_chain_ID = site_info[current_atom].FIRST_chain_ID;
	
	cout<<"Unknown ring type: "<<current_atom<<" "<<site_info[current_atom].atom_name<<"("<<current_rname<<") "<<current_residue<<" "<<current_chain_ID<<endl;
	
// 	if( !is_ribose_atom(current_atom) && !is_nucleic_backbone(current_atom) && !is_) return(false);
// 	int atom_coord = (site_info[atom_1].neighbor_list).size();
// 	if(atom_coord<2) return(false); // if fewer than 2 bonds skip
	
// 	string element = site_info[atom_1].element_name;
// 	// if element is not N, O, C, skip    
// 	if( element != "N " && element != "C " && element != "O ") return(false);
// 	// count number of hydrogen neighbors
// 	int nhydr=0;
// 	//      cout << element<<" " <<atom_1<< " "<<atom_coord<<endl;
// 	for( int b=0; b<atom_coord; b++) {
// 	  if( site_info[(site_info[atom_1].neighbor_list[b])].element_name=="H ") nhydr++;
// 	}
// 	if( atom_coord - nhydr <2) return(false); // if all (or all but one) bond are to hydrogens not in a ring.
// 	else {
      
// 	  // if(site_info[atom_1].seq_number == 509) 
// 	  //    cout<<"\t"<<atom_1<<" "<<element<<" "<<site_info[atom_1].atom_name<<" "<<atom_coord<<" "<<nhydr<<"\t"
// 	  // 	    <<NthNearestNeighbor(atom_1,atom_1,5,atom_1)<<" "<<NthNearestNeighbor(atom_1,atom_1,6,atom_1)<<endl;
// 	  if( NthNearestNeighbor(atom_1,atom_1,5)) return(true);
// 	  else if(NthNearestNeighbor(atom_1,atom_1,6)) return(true);
// 	}
      }
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Look at atoms connected to "atom_1" to see if it belongs in a 5 or 6 fold 
//   ring. Otherwise return(0). 
// Parameters:
//   atom_1 - FIRST_number of a given atom.
//////////////////////////////////////////////////////////////////////////////// 
bool MolFramework::is_ring_atom(unsigned int atom_1){
  
  // omit Hydrogens && Water
  if(site_info[atom_1].element_name == "H ") return(false);
  else if(site_info[atom_1].residue_name == "HOH" ) return(false);

  // Handle Proline  
  //else if( site_info[atom_1].residue_name == "PRO" ) {
  //if(site_info[atom_1].element_name == "O " || site_info[atom_1].atom_name == "C") return(false);
  //else return(true);
  //}
  // omit standard AA residues
  else if(site_info[atom_1].residue_name == "ALA" || site_info[atom_1].residue_name == "ARG" || 
	  site_info[atom_1].residue_name == "ASN" || site_info[atom_1].residue_name == "ASP" ||
 	  site_info[atom_1].residue_name == "CYS" || site_info[atom_1].residue_name == "GLN" ||
 	  site_info[atom_1].residue_name == "GLU" || site_info[atom_1].residue_name == "GLY" ||
 	  site_info[atom_1].residue_name == "ILE" || site_info[atom_1].residue_name == "LEU" ||
 	  site_info[atom_1].residue_name == "LYS" || site_info[atom_1].residue_name == "MET" ||
 	  site_info[atom_1].residue_name == "SER" || site_info[atom_1].residue_name == "THR" ||
	  site_info[atom_1].residue_name == "VAL" || site_info[atom_1].residue_name == "PRO" ) 
    return(false);
  // omit mainchain atoms 
  else if(isMainchain(atom_1)) return(false);

  // handle rings in standard amino acids:
  else if(site_info[atom_1].residue_name == "HIS" || 
	  site_info[atom_1].residue_name == "PHE" ||
 	  site_info[atom_1].residue_name == "TRP" || 
	  site_info[atom_1].residue_name == "TYR" ) {

    if( site_info[atom_1].atom_name == "CB" || 
	site_info[atom_1].atom_name == "OH" ||
	site_info[atom_1].atom_name == "OXT" ) 
      return( false );

    else return(true);
  }

  // handle rings in standard nucleic acids:
  else if(site_info[atom_1].residue_name == "  A" || site_info[atom_1].residue_name == "  C" ||
  	  site_info[atom_1].residue_name == "  T" || site_info[atom_1].residue_name == "  U" ||
  	  site_info[atom_1].residue_name == "  G" ) {
    if(is_ribose_atom(atom_1)) return(false);
    if(is_nucleic_backbone(atom_1)) return(false);
    int atom_coord = (site_info[atom_1].neighbor_list).size();

    string element = site_info[atom_1].element_name;
    if(atom_coord<2) return(false); // if fewer than 2 bonds skip   
    // count number of hydrogen neighbors
    int nhydr=0;
    for( int b=0; b<atom_coord; b++) 
      if( site_info[(site_info[atom_1].neighbor_list[b])].element_name=="H ") nhydr++;
    // if all (or all but one) bond are to hydrogens not in a ring.
    // cout<<"ringatom:: "<<atom_1<<" "<<site_info[atom_1].atom_name<<"\t "<<atom_coord
    //<<"\t "<<nhydr<<endl;

    if( atom_coord - nhydr <2) return(false); 
    else
      return(true);
  }
  // now I have to deal with rings from non-standard groups
  else {
    
    //cout<<"OtherRing: "<<atom_1<<" "<<site_info[atom_1].atom_name<<"("<<site_info[atom_1].residue_name<<")"<<endl;
    //<<"\t "<<is_ribose_atom(atom_1)<<" "<<is_nucleic_backbone(atom_1)<<endl;
    if( is_ribose_atom(atom_1) ) return(false);
    if((site_info[atom_1].neighbor_list).size()<2) return(false); // if fewer than 2 bonds skip
    string element = site_info[atom_1].element_name;
    // if element is not N, O, C, skip    
    if( element != "N " && element != "C " && element != "O ") return(false);
    // count number of non-hydrogen neighbors
    int nbonds=0;
    int atom_coord = (site_info[atom_1].neighbor_list).size();
    for( int b=0; b <atom_coord; b++) 
      if( site_info[(site_info[atom_1].neighbor_list[b])].element_name!="H ") nbonds++;
    if( nbonds < 2) return(false);
    else {
      int rdepth=0;
      int ringflag = get_ring_atom(atom_1, rdepth);
      
      if(ringflag) return(true);
    }
  }
  return(false);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Look at "atom_1" to see if its atomname explicitly matches those used for 
//   dna/rna ribose sugar atoms. Hydrogen atoms are omitted.
// Parameters:
//   atom_1 - FIRST_number of a given atom.
//////////////////////////////////////////////////////////////////////////////// 
bool MolFramework::is_ribose_atom(unsigned int atom_1){

  // omit Hydrogens
  if(site_info[atom_1].element_name == "H ") return(false);
  else if(isBackbone(atom_1)) return(false);
  else {
    string aname = site_info[atom_1].atom_name;
    // note: this construction asumes that the sugars have been renamed with '''. 
    // AJR 04.03.08 to conform with PDBV3.0 * replaced by '.
    if(aname == "C1'" || aname =="C2'" || aname == "C3'" || aname == "C4'" ||
       aname == "C5'" || aname =="O2'" || aname == "O3'" || aname == "O4'" ||
       aname == "O5'") {
      if(NthNearestNeighbor(atom_1,"C1'",4)) {
	unsigned int c1atom = NthNearestNeighbor(atom_1,"C1'",4);
	// check that this atom && neighbor atom are in a ring of 5 atoms
	if( NthNearestNeighbor(atom_1,atom_1,5) == atom_1 &&
	    NthNearestNeighbor(c1atom,c1atom,5) == c1atom ){

	  if(NthNearestNeighbor(c1atom,"N1",1) || NthNearestNeighbor(c1atom,"N9",1) )
	    return(true); 	  // only accept those that are bonded to nucleobases
	  // TODO AJR 11.11.05 make this more generic.
	  else if(NthNearestNeighbor(c1atom, "C5",1) )
	    //else if(site_info[atom_1].residue_name =="PSU" && NthNearestNeighbor(c1atom,"C5",1))
	    return(true);
	}
      }
    }
    //Old PDB format
    if(aname == "C1*" || aname =="C2*" || aname == "C3*" || aname == "C4*" ||
       aname == "C5*" || aname =="O2*" || aname == "O3*" || aname == "O4*" ||
       aname == "O5*") {
      if(NthNearestNeighbor(atom_1,"C1*",4)) {
	unsigned int c1atom = NthNearestNeighbor(atom_1,"C1*",4);
	// check that this atom && neighbor atom are in a ring of 5 atoms
	if( NthNearestNeighbor(atom_1,atom_1,5) == atom_1 &&
	    NthNearestNeighbor(c1atom,c1atom,5) == c1atom ){

	  if(NthNearestNeighbor(c1atom,"N1",1) || NthNearestNeighbor(c1atom,"N9",1) )
	    return(true); 	  // only accept those that are bonded to nucleobases
	  else if(NthNearestNeighbor(c1atom, "C5",1) )
	    return(true);
	}
      }
    }

    else if(site_info[atom_1].record_name=="HETATM") {
     size_t namePosition=0; // FIXME - what is namePosition?
      string element = site_info[atom_1].element_name;
      if(element=="C " ) {
	//|| element=="O") {
	namePosition=aname.find("C1'");
	if(namePosition!=string::npos) return(true);
	namePosition=aname.find("C2'");
	if(namePosition!=string::npos) return(true);
	namePosition=aname.find("C3'");
	if(namePosition!=string::npos) return(true);
	namePosition=aname.find("C4'");
	if(namePosition!=string::npos) return(true);
	namePosition=aname.find("C5'");
	if(namePosition!=string::npos) return(true);
      }
      else if(element=="O ") {
	namePosition=aname.find("O2'");
	if(namePosition!=string::npos) return(true);
	namePosition=aname.find("O3'");
	if(namePosition!=string::npos) return(true);
	namePosition=aname.find("O4'");
	if(namePosition!=string::npos) return(true);
	namePosition=aname.find("O5'");
	if(namePosition!=string::npos) return(true);
	//	  cout<<namePosition<<" "<<string::npos<<"\t "<<aname<<" "<<element<<atom_1<<endl;
      }
    }
    else
      return(false);
  }

  return false;
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////////////////////////////////////////////////////////////////////////////////
// Description:
//   Look at "atom_1" to see if its atomname explicitly matches those used for 
//   dna/rna nucleobase atoms. Hydrogen atoms are omitted.
// Parameters:
//   atom_1 - FIRST_number of a given atom.
//////////////////////////////////////////////////////////////////////////////// 
bool MolFramework::is_nucleobase_atom(unsigned int atom_1){
  
  unsigned int c1atom;
  // omit Hydrogens
  if(site_info[atom_1].element_name == "H ") return(false);
  else if(isBackbone(atom_1)) return(false);
  else {
    string aname = site_info[atom_1].atom_name;
    if(aname == "N1" || aname == "N2" || aname =="C2" || aname == "N3" ||  aname == "C4" || aname == "N4" ||
       aname == "C5" || aname =="C6" || aname == "N6" || aname == "N7" || aname == "C8" ||
       aname == "N9" || aname == "C8" || aname == "O2" || aname == "O4" || aname == "O6"  ) {
      c1atom = NthNearestNeighbor(atom_1,"C1'",5);


      if(c1atom) {
	if(NthNearestNeighbor(c1atom,"N1",1) || NthNearestNeighbor(c1atom,"N9",1) )
	  return(true); 	  // only accept those that are bonded to nucleobases
	
	else if(NthNearestNeighbor(c1atom, "C5",1) )
	  return(true);
      }
      c1atom = NthNearestNeighbor(atom_1,"C1*",5); //Old PDB format
      if(c1atom) {
	if(NthNearestNeighbor(c1atom,"N1",1) || NthNearestNeighbor(c1atom,"N9",1) )
	  return(true); 	  // only accept those that are bonded to nucleobases
	
	else if(NthNearestNeighbor(c1atom, "C5",1) )
	  return(true);
      }
 
    }
  }

  return false;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Look at "atom_1" to see if its atomname explicitly matches those used for 
//   dna/rna MODIFIED nucleobase atoms. Hydrogen atoms are omitted.
// Parameters:
//   atom_1 - FIRST_number of a given atom.
//////////////////////////////////////////////////////////////////////////////// 
bool MolFramework::is_modifiedNucleobase_atom(unsigned int atom_1){
  
  unsigned int c1atom;
  // omit Hydrogens
  if(site_info[atom_1].element_name == "H ") return(false);
  else if(isBackbone(atom_1)) return(false);
  else {
    string aname = site_info[atom_1].atom_name;
    //For modified nucleotides 
    if(aname == "CM1" || aname == "CM2" || aname =="CM5" || aname == "C5M" ||  aname == "CM7" || 
       aname == "C10" || aname =="S4" ) {
      c1atom = NthNearestNeighbor(atom_1,"C1'",5);
      if(c1atom) {
	if(NthNearestNeighbor(c1atom,"N1",1) || NthNearestNeighbor(c1atom,"N9",1) )
	  return(true); 	  // only accept those that are bonded to nucleobases
	
	else if(NthNearestNeighbor(c1atom, "C5",1) )
	  return(true);
      }
      c1atom = NthNearestNeighbor(atom_1,"C1*",5);//Old PDB format
      if(c1atom) {
	if(NthNearestNeighbor(c1atom,"N1",1) || NthNearestNeighbor(c1atom,"N9",1) )
	  return(true); 	  
	
	else if(NthNearestNeighbor(c1atom, "C5",1) )
	  return(true);
      }
    }
  }

  return false;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Look at "atom_1" to determine if it is the ring atom that anchors the 
//   aromatic ring to the rest of the molecule. 
// Parameters:
//   atom_1 - FIRST_number of a given atom.
// Returns: 
//    0 if it is not an anchor atom
//    non-ring neighobor of anchor atom otherwise.
//////////////////////////////////////////////////////////////////////////////// 
int MolFramework::is_ring_anchor(unsigned int atom_1) {
  string aname = site_info[atom_1].atom_name;
  //if(site_info[atom_1].residue_name=="PRO") {  // case for PRO 
  //if(aname=="CA")   return( NthNearestNeighbor(atom_1,"C",1));
  //else return(0);
  //}
  if(aname == "N1" || aname == "N9")  // case for nucleotides
    return(NthNearestNeighbor(atom_1,"C1'",1));
  else if(aname == "CG")  // case for TRP,PHE,TYR, HIS,& non-std aromatic AA
    return(NthNearestNeighbor(atom_1,"CB",1));
  else if(int test = NthNearestNeighbor(atom_1,"C1'",1)) //case for non-std NT
    return(test);
  else if(int test = NthNearestNeighbor(atom_1,"C1*",1)) //case for non-std NT //Old PDB format
    return(test);
  else {	
    size_t n1location=aname.find("N1");
    size_t n9location=aname.find("N9");
    string aname2;
    if(n9location!=string::npos) {
      aname2.assign(aname,0,n9location);
      aname2.append("C1'");
      if(int test = NthNearestNeighbor(atom_1,aname2,1)) {
	//cout<<atom_1<<" "<<aname<<" "<<aname2<<endl;
	return(test);
      }
    }
    else if(n1location!=string::npos) {
      aname2.assign(aname,0,n1location);
      aname2.append("C1'");
      if(int test = NthNearestNeighbor(atom_1,aname2,1)) {
	cout<<atom_1<<" "<<aname<<" "<<aname2<<endl;
	return(test);
      }
    }
    //    else if
    //else
    return(0); 
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//    function to return the ring id number corresponding to a given atom
// Input: 
//    atom_1 --> query atom
//    number_of_rings --> 
//    *new_rings --> the aromatic_ring struct 
////////////////////////////////////////////////////////////////////////////////
int MolFramework::get_ringid(unsigned int atom_1, int number_of_rings, aromatic_ring *new_rings) {

////
  stringstream resnumid_string;
  resnumid_string << site_info[atom_1].seq_number << ";" << site_info[atom_1].insert_code << ";"
		  << site_info[atom_1].chain_ID;
  int chk_unique_id = unique_res_id[resnumid_string.str()];
  
  // loop over the atoms in each aromatic ring to find the query atom.
  for( int a = 1; a <= number_of_rings; a++ ){
    if(chk_unique_id == new_rings[a].resid) {
      vector<unsigned int>::iterator new_atoms=(new_rings[a].ring_atoms).begin();
      while(new_atoms != (new_rings[a].ring_atoms).end()) {
	if( *new_atoms == atom_1) return(a);
 	new_atoms++;
      }
    }
  }
  return(0);
}
////////////////////////////////////////////////////////////////////////////////
// Description: 
//    function to check if a given aromatic ring is already stacked to another 
//    aromatic ring.
//////////////////////////////////////////////////////////////////////////////// 
bool MolFramework::is_stacked_ring_neighbor(int ring_1, aromatic_ring check_aromatic) {

  // exit if check_ring has no ring neighbors
   if(check_aromatic.stacked_neighbor_list.size() == 0) return(false);
  vector<int>::iterator check_ring = (check_aromatic.stacked_neighbor_list).begin();
  while( check_ring != (check_aromatic.stacked_neighbor_list).end() ) {
    if( *check_ring == ring_1) 
      return(true);
    else
      check_ring++;
  }
  return(false);
}


////////////////////////////////////////////////////////////////////////////////  
bool MolFramework::is_stacked_ring_pair(aromatic_ring ring_1, aromatic_ring ring_2) {
  float ctr_ctr_distance;
  float snd_snd_distance;
  float snd_ctr_distance=0;
  float ctr_snd_distance=0;
  float center_vect[3];
  float eta,alpha1,alpha2;
  int size1 = (ring_1.ring_atoms).size();
  int size2 = (ring_2.ring_atoms).size();
  int typeFlag = 0;
  string type = "RS"; 
  bool nb_atom1 = is_nucleobase_atom(ring_1.anchor_atom);
  bool nb_atom2 = is_nucleobase_atom(ring_2.anchor_atom);
  bool ribose_atom1 = is_ribose_atom(ring_1.anchor_atom);
  bool ribose_atom2 = is_ribose_atom(ring_2.anchor_atom);
  float cutoff_stack_distance = 0;
 
  // Modification to take into account parameters for RNA, DNA. Source: paper of S.Fulle and H.Gohlke 
  // "Analyzing the Flexibility of RNA Structures by Constraint Counting" Biophys. J., 2008, 94,4202-4219
  // The stacked rings in RNAs, DNAs are taken into account through hydrophobic tethers. See  MolFramework::identifyHydrophobicTethers()
  if(nb_atom1 && nb_atom2 ) return(false);
  else if (nb_atom1 && ribose_atom2) return(false);
  else if (ribose_atom1 && nb_atom2 ) return(false);
  else if (ribose_atom1 && ribose_atom2 ) return(false);
 
  cutoff_stack_distance = parameters.cutoff_SR_distance;
  //cout <<ring_1.anchor_atom<<" "<<ring_2.anchor_atom <<" cutoff_stack_distance "<< cutoff_stack_distance <<endl;

  // 
  // define typeFlag for the ring sizes forming the RS pair interaction:
  // typeFlag = 4 for six-six, 3 for five-six, 2 for six-five and 1 for five-five 
  // a float with value of -1.0*typeFlag is entered into the stacked_ring new_bond structure
  //---------------------------------------------------------------------------------------

  ctr_ctr_distance = VectorAlgebra::distance( ring_1.geom_center, ring_2.geom_center);
  //
  // A. check the distance between the ring centers (six-fold ring centers) < 5.5 Ang. for protein-protein  or 
  // protein-nucleic acid, the distance < 3.55  for nucleic acid-nucleic acid rings
  //    distance set by parameters.cutoff_SR_distance and parameters.cutoff_SRnucleic_distance
  //---------------------------------------------------------------------------------
  if(ctr_ctr_distance <= cutoff_stack_distance ) { 
    for(int b=0; b<3; b++) 
      center_vect[b] = (ring_1.geom_center[b] - ring_2.geom_center[b])/ctr_ctr_distance;
     // B. check the angle between the two ring plane normals (eta) & the angles between 
    //    the center_to_center vector and the ring plane normals (alpha1 & alpha2).
    //    It doesn't matter the orientation so use absolute value of the dot product.
    //    normal_to_normal angle is called eta.
    //---------------------------------------------------------------------------------
    eta = acos( abs(VectorAlgebra::dot(ring_1.ring_norm,ring_2.ring_norm)) );

    if(rad2deg(eta) <= parameters.cutoff_SR_normal_plane_angle ){
      alpha1= acos(abs(VectorAlgebra::dot(center_vect,ring_1.ring_norm)));
      alpha2= acos(abs(VectorAlgebra::dot(center_vect,ring_2.ring_norm)));

      if(rad2deg(alpha1) <= parameters.cutoff_SR_center_normal_angle && 
	 rad2deg(alpha2) <= parameters.cutoff_SR_center_normal_angle ) {

	// C.  add to the stacked_rings bondlist: connecting the anchor_atoms 
	// of each base and classified by the typeFlag variable.
	// HIS typeFlag interaction fixed AJR 06.01.06
	//------------------------------------------------------------------
	typeFlag=4;
	if(size1==5 && size2==5) typeFlag=1;
	else if(size1==5) typeFlag=3;
	else if(size2==5) typeFlag=2;
	stacked_rings.push_back( new_bonds(ring_1.anchor_atom,ring_2.anchor_atom,
					   parameters.rs_bars, -1.0*typeFlag, type)); 
	//cout<<ring_1.residue_name<<ring_1.residue_number
	//<<">:<"<<ring_2.residue_name<<ring_2.residue_number
	//<<" CCd: "<<ctr_ctr_distance<<" ..norm2norm_angle: "
	// 	    <<rad2deg(eta)<<" , ctr2norm_angles: "<<rad2deg(alpha1)<<" , "
	//<<rad2deg(alpha2)<<endl;

	if(parameters.verbose == 3) {
	  cout<<ring_1.residue_name<<ring_1.residue_number<<" >"<<typeFlag<<"< "
	      <<ring_2.residue_name<<ring_2.residue_number
	      <<" CCd: "<<ctr_ctr_distance<<" ..norm2norm_angle: "
	      <<rad2deg(eta)<<" ,ctr2norm_angles: "<<rad2deg(alpha1)<<" "
	      <<rad2deg(alpha2)<<endl;
	}
	
	if(size1<=6 && size2<=6) // exit if both are single ring residues
	  return(true);
      }
    }
  }
  // --------- now check for cases where the residue of ring_1 has multiple rings
  if( size1 > 6) {
    snd_ctr_distance = VectorAlgebra::distance( ring_1.second_cntr, ring_2.geom_center);
  
    // A. check the distance between the ring centers (six-fold ring centers) < 5.5 Ang.
    //---------------------------------------------------------------------------------
    if(snd_ctr_distance <= cutoff_stack_distance) {
      // checks snd_ctr_distance only if shorter than ctr_ctr_distance or no stacked_ring
      // yet identified.
      if(!typeFlag || snd_ctr_distance < ctr_ctr_distance) {

	for(int b=0; b<3; b++) 
	  center_vect[b] = (ring_1.second_cntr[b] - ring_2.geom_center[b])/snd_ctr_distance;
	// B. check the angle between the two ring plane normals (eta) & the angles between 
	//    the center_to_center vector and the ring plane normals (alpha1 & alpha2).
	//    It doesn't matter the orientation so use absolute value of the dot product.
	//---------------------------------------------------------------------------------
	eta = acos(abs(VectorAlgebra::dot(ring_1.ring_norm,ring_2.ring_norm)));
	if(rad2deg(eta) <= parameters.cutoff_SR_normal_plane_angle ){
	  alpha1= acos(abs(VectorAlgebra::dot(center_vect,ring_1.ring_norm)));
	  alpha2= acos(abs(VectorAlgebra::dot(center_vect,ring_2.ring_norm)));
	  if(rad2deg(alpha1) <= parameters.cutoff_SR_center_normal_angle && 
	     rad2deg(alpha2) <= parameters.cutoff_SR_center_normal_angle &&
	     eta != 0.0 ) {
	    //
	    // C. add to the stacked_rings bondlist: connecting the anchor_atoms of each base. 
	    //    if this is a better interaction than before, remove the previously assigned
	    //    stacked_rings interaction and replace with this one.
	    //-----------------------------------------------------------------------------------
	    if(typeFlag) 
	      stacked_rings.pop_back(); 
	    typeFlag =3; // 
	    if(size2==5) typeFlag=1; //adjust for case where ring2 is HIS 
	    stacked_rings.push_back( new_bonds(ring_1.anchor_atom,ring_2.anchor_atom,
					       parameters.rs_bars, -1.0*typeFlag, type)); 

	    if(parameters.verbose) {
	      cout<<ring_1.residue_name<<ring_1.residue_number<<" >"<<typeFlag<<"< "
		  <<ring_2.residue_name<<ring_2.residue_number
		  <<" CCd: "<<snd_ctr_distance<<" ..norm2norm_angle: "
		  <<rad2deg(eta)<<" ,ctr2norm_angles: "<<rad2deg(alpha1)<<" "
		  <<rad2deg(alpha2)<<endl;
	    }  
	    // 	    cout<<ring_1.residue_name<<ring_1.residue_number
	    //<<">1<"<<ring_2.residue_name<<ring_2.residue_number
	    //<<" CCd: "<<snd_ctr_distance<<" ..norm2norm_angle: "
	    //<<rad2deg(eta)<<" , ctr2norm_angles: "<<rad2deg(alpha1)<<" , "
	    //<<rad2deg(alpha2)<<endl;

	    if(size2 <= 6) //exit if second ring residue has only 1 ring.
	      return(true);
	  }
	}
      }
    }
  }
  if(size2 >6){ 
    ctr_snd_distance = VectorAlgebra::distance(ring_1.geom_center,ring_2.second_cntr);
  // A. check the distance between the ring centers (six-fold ring centers) < 5.5 Ang.
  //---------------------------------------------------------------------------------
    if(ctr_snd_distance <= cutoff_stack_distance) {
      if(!typeFlag || ctr_snd_distance < ctr_ctr_distance) {// checks ctr_snd only for shorter values
	// for case where both residues are multiring
	if( (size1>6 && ctr_snd_distance < snd_ctr_distance) || size1<=6 ) { 

	  for(int b=0; b<3; b++) 
	    center_vect[b] = (ring_2.second_cntr[b] - ring_1.geom_center[b])/ctr_snd_distance;
      
	  // B. check the angle between the two ring plane normals (eta) & the angles between 
	  //    the center_to_center vector and the ring plane normals (alpha1 & alpha2).
	  //    It doesn't matter the orientation so use absolute value of the dot product.
	  //---------------------------------------------------------------------------------
	  eta = acos(abs(VectorAlgebra::dot(ring_1.ring_norm,ring_2.ring_norm)));
	  if(rad2deg(eta) <= parameters.cutoff_SR_normal_plane_angle ){
	    alpha1= acos(abs(VectorAlgebra::dot(center_vect,ring_1.ring_norm)));
	    alpha2= acos(abs(VectorAlgebra::dot(center_vect,ring_2.ring_norm)));
	    if(rad2deg(alpha1) <= parameters.cutoff_SR_center_normal_angle && 
	       rad2deg(alpha2) <= parameters.cutoff_SR_center_normal_angle &&
	       eta != 0.0 ) {
	      if(typeFlag) stacked_rings.pop_back();
	      typeFlag=2;
	      if(size1==5) typeFlag=1; //adjust for case where ring1 is HIS 
	      // add to the stacked_rings bondlist: connecting the anchor_atoms of each base. 
	      stacked_rings.push_back( new_bonds(ring_1.anchor_atom,ring_2.anchor_atom,
						 parameters.rs_bars, -1.0*typeFlag, type)); 
	      if(parameters.verbose) {
		cout<<ring_1.residue_name<<ring_1.residue_number<<" >"<<typeFlag<<"< "
		    <<ring_2.residue_name<<ring_2.residue_number
		    <<" CCd: "<<ctr_snd_distance<<" ..norm2norm_angle: "
		    <<rad2deg(eta)<<" ,ctr2norm_angles: "<<rad2deg(alpha1)<<" "
		    <<rad2deg(alpha2)<<endl;
	      }  
	      // 	      cout<<ring_1.residue_name<<ring_1.residue_number
	      //  <<">2<"<<ring_2.residue_name<<ring_2.residue_number
	      //  <<" CCd: "<<ctr_snd_distance<<" ..norm2norm_angle: "
	      //  <<rad2deg(eta)<<" , ctr2norm_angles: "<<rad2deg(alpha1)<<" , "
	      //  <<rad2deg(alpha2)<<endl;
	      
	      if(size1<=6) //exit if first residue has only 1 ring.
		return(true);
	    }
	  }
	}
      }
    }
  }
  if( size1 >6 &&  size2 > 6) {
    snd_snd_distance = VectorAlgebra::distance( ring_1.second_cntr, ring_2.second_cntr);
    // A. check the distance between the ring centers (five-fold ring centers) < 5.5 Ang.
    //---------------------------------------------------------------------------------
    if(snd_snd_distance <= cutoff_stack_distance) {
      switch(typeFlag) {
      case 1:
	if( snd_snd_distance > ctr_ctr_distance ) break;
      case 2: 
	if( snd_snd_distance > ctr_snd_distance ) break;
      case 3:
	if( snd_snd_distance > snd_ctr_distance ) break;
      default: 
	for(int b=0; b<3; b++) 
	  center_vect[b] = (ring_1.second_cntr[b] - ring_2.second_cntr[b])/snd_snd_distance;
	
	// B. check the angle between the two ring plane normals (eta) & the angles between 
	//    the center_to_center vector and the ring plane normals (alpha1 & alpha2).
	//    It doesn't matter the orientation so use absolute value of the dot product.
	//---------------------------------------------------------------------------------
	eta = acos(abs(VectorAlgebra::dot(ring_1.ring_norm,ring_2.ring_norm)));
	if(rad2deg(eta) <= parameters.cutoff_SR_normal_plane_angle ){
	  alpha1= acos(abs(VectorAlgebra::dot(center_vect,ring_1.ring_norm)));
	  alpha2= acos(abs(VectorAlgebra::dot(center_vect,ring_2.ring_norm)));
	  if(rad2deg(alpha1) <= parameters.cutoff_SR_center_normal_angle && 
	     rad2deg(alpha2) <= parameters.cutoff_SR_center_normal_angle ) {
    
	    if( typeFlag ) 
	      stacked_rings.pop_back();
	    typeFlag=1;
	    // add to the stacked_rings bondlist: connecting the anchor_atoms of each base. 
	    stacked_rings.push_back( new_bonds(ring_1.anchor_atom,ring_2.anchor_atom,
					       parameters.rs_bars, -1.0*typeFlag, type)); 
	    if(parameters.verbose) {
	      cout<<ring_1.residue_name<<ring_1.residue_number<<" >"<<typeFlag<<"< "
		  <<ring_2.residue_name<<ring_2.residue_number
		  <<" CCd: "<<snd_snd_distance<<" ..norm2norm_angle: "
		  <<rad2deg(eta)<<" ,ctr2norm_angles: "<<rad2deg(alpha1)<<" "
		  <<rad2deg(alpha2)<<endl;
	    }
 	    //cout<<ring_1.residue_name<<ring_1.residue_number
	    //<<">3<"<<ring_2.residue_name<<ring_2.residue_number
	    //<<" CCd: "<<snd_snd_distance<<" ..norm2norm_angle: "
	    //<<rad2deg(eta)<<" , ctr2norm_angles: "<<rad2deg(alpha1)<<" , "
	    //<<rad2deg(alpha2)<<endl;
	    //	     return(true);
	   }
	 }
       }
    }
  }
  if(typeFlag) 
    return(true);
  else
    return(false);
}


/////////////////////////////////////////////////////////////

bool MolFramework::is_nucleic_backbone(unsigned int atom_1) {
  // exclude protein mainchains
  if(isMainchain(atom_1) ) return(false);
  if(is_ribose_atom(atom_1) ) return(false);
  // cope with known NA residues: 
  if(site_info[atom_1].residue_name == "  A" || site_info[atom_1].residue_name == "  C" ||
     site_info[atom_1].residue_name == "  T" || site_info[atom_1].residue_name == "  U" ||
     site_info[atom_1].residue_name == "  G" || site_info[atom_1].residue_name == " DA" ||
     site_info[atom_1].residue_name == " DC" || site_info[atom_1].residue_name == " DT" ||
     site_info[atom_1].residue_name == " DG" ) {
    if( site_info[atom_1].atom_name == "P") return(true);
    else if(site_info[atom_1].element_name=="O ") {
      if(NthNearestNeighbor(atom_1,"P",1) && !is_ribose_atom(atom_1) ) return(true);
    }
    else return(false);
  }
  else if(site_info[atom_1].atom_name=="P" &&
	  NthNearestNeighbor( atom_1, "O1P", 1) &&
	  NthNearestNeighbor( atom_1, "O2P", 1) && 
	  NthNearestNeighbor( atom_1, "O5'", 1)) return(true);
  else if(site_info[atom_1].element_name=="O " && 
	  NthNearestNeighbor( atom_1, "P", 1) && NthNearestNeighbor( atom_1, "O5'",2) &&
	  (NthNearestNeighbor( atom_1, "O1P",2) || NthNearestNeighbor( atom_1, "O2P",2) ||
	   NthNearestNeighbor( atom_1, "O3P",2) )) return(true);
  //Possible different names of atoms
  else if(site_info[atom_1].atom_name=="P" &&
	  NthNearestNeighbor( atom_1, "OP1", 1) &&
	  NthNearestNeighbor( atom_1, "OP1", 1) && 
	  NthNearestNeighbor( atom_1, "O5'", 1)) return(true);
  else if(site_info[atom_1].element_name=="O " && 
	  NthNearestNeighbor( atom_1, "P", 1) && NthNearestNeighbor( atom_1, "O5'",2) &&
	  (NthNearestNeighbor( atom_1, "OP1",2) || NthNearestNeighbor( atom_1, "OP2",2) ||
	   NthNearestNeighbor( atom_1, "OP3",2) )) return(true);
  //Old PDB format
  else if(site_info[atom_1].atom_name=="P" &&
	  NthNearestNeighbor( atom_1, "OP1", 1) &&
	  NthNearestNeighbor( atom_1, "OP1", 1) && 
	  NthNearestNeighbor( atom_1, "O5*", 1)) return(true);
  else if(site_info[atom_1].element_name=="O " && 
	  NthNearestNeighbor( atom_1, "P", 1) && NthNearestNeighbor( atom_1, "O5*",2) &&
	  (NthNearestNeighbor( atom_1, "OP1",2) || NthNearestNeighbor( atom_1, "OP2",2) ||
	   NthNearestNeighbor( atom_1, "OP3",2) )) return(true);

  else
    return(false);

  return false;
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   A depth-first search from atom_1 through the covalent-bond network to find 
//   atom_2. The funtion searches recursively through the number of bonds set by
//   the variable degree. A record of the sites that have been checked is kept
//   in order to establish the existence of a loop. 
////////////////////////////////////////////////////////////////////////////////
int MolFramework::get_ring_atom( unsigned int atom_1, unsigned int ring_depth) {
  
  vector<unsigned int> ringatom_list;
  
  vector<unsigned int>::iterator atom_2 = (site_info[atom_1].neighbor_list).begin();
  while( atom_2 != (site_info[atom_1].neighbor_list).end() ) {
    if( site_info[*atom_2].element_name != "H ") {
      //      int atom_coord = (site_info[*atom_2].neighbor_list).size();
      for( unsigned int neighborNumber=0; neighborNumber<(site_info[*atom_2].neighbor_list).size(); neighborNumber++) {
	unsigned int atom_3 = site_info[*atom_2].neighbor_list[neighborNumber];
       	if(site_info[atom_3].element_name!="H " && atom_3 != atom_1 &&
	   !is_ribose_atom(atom_3) && !isMainchain(atom_3) ) {
	  ringatom_list.push_back(atom_3);
	  //cout<<"\t"<<atom_1<<"-- "<<*atom_2<<" bond "<<neighborNumber<<" is "
	  //    << atom_3 <<" "<< site_info[atom_3].atom_name<<endl;
	  if(ringatom_list.size() >1) {
	    vector<unsigned int>::iterator next_atom = ringatom_list.begin();
	    while( next_atom != ringatom_list.end()) {
	      for( unsigned int secondNeighborAtom=0; secondNeighborAtom<(site_info[atom_3].neighbor_list).size(); secondNeighborAtom++) {
		unsigned int atom_4 = site_info[atom_3].neighbor_list[secondNeighborAtom];
		if(atom_4 == *next_atom) {
		  ring_depth = 5;
		  return(*next_atom);
		}
	      } 
	      next_atom++;
	    }
	  }
		}
      }
    }
    atom_2++;
  }
  // now try to tackle 6-fold rings:

  vector<unsigned int> ringatom_list2;
  vector<unsigned int>::iterator next_atom = ringatom_list.begin();
  while( next_atom != ringatom_list.end()) {
    for( unsigned int neighborNumber=0; neighborNumber<(site_info[*next_atom].neighbor_list).size(); neighborNumber++) {
      unsigned int atom_4 = site_info[*next_atom].neighbor_list[neighborNumber];
      if(site_info[atom_4].element_name!="H " && !is_ribose_atom(atom_4) 
	 && !isMainchain(atom_4) && NthNearestNeighbor(atom_4,atom_1,1)!=atom_1 ) {
	  ringatom_list2.push_back(atom_4);
      }
    } 
    next_atom++;
  }
  unsigned int r2size = ringatom_list2.size();
  //  cout<<"  ringatom_list_size is : "<<ringatom_list2.size()<<endl;
  if(ringatom_list2.size() > 1){
    vector<unsigned int> atom = ringatom_list2;
    for(unsigned int i=0; i<r2size-1; i++) {
      for(unsigned int j=i+1;j<r2size; j++) {
	if(atom[j]==atom[i]) {
	  //    vector<int>::iterator next_atom = ringatom_list2.begin();
	  //    while( next_atom != ringatom_list2.end()) {
	  //      int listiterator = rinatom_list2.at;
	  //for( int c=0; c<(site_info[atom_3].neighbor_list).size(); c++) {
	  ring_depth = 6;
	 //cout<<"found match<6>: "<<site_info[atom_1].atom_name<<"\t"<<atom[i]<<" "<<ring_depth<<endl;
	  return(atom[i]);
	}
      }
    }
  }
  return(0);
}

////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is a depth-first recursive search through the site_info topology 
//   starting on site "current_atom". Only connections to sites with the same 
//   residue number and chain ID are followed. If the site is in a ring, it is 
//   added to "ring_atom_list" which is of type vector<int>. A list of the
//   sites that have been checked is kept to avoid cycles. 
// Parameters:
//   current_atom - 
////////////////////////////////////////////////////////////////////////////////////
 void MolFramework::list_ring_atoms( int current_atom, int current_residue, 
 						    int current_chain_ID ){

   search_list.push_back( current_atom );
   site_info[current_atom].checked = true;
   
   if( !is_ribose_atom(current_atom) && site_info[current_atom].element_name !="H " 
       && !is_nucleic_backbone(current_atom)) {
     // check the number of nonhydrogen neighbors:
     int neighbor = 0;
     int nbonds =0;
     for( unsigned int neighborNumber = 0; neighborNumber < (site_info[current_atom].neighbor_list).size(); neighborNumber++) {
       neighbor = site_info[current_atom].neighbor_list[neighborNumber];
       if( site_info[neighbor].element_name !="H ") nbonds++;
     }
     if(nbonds >= 2) {
     cout<<"RING ATOM " << current_atom<<" "<< current_residue<<"("<<site_info[current_atom].residue_name
	 <<") "<<site_info[current_atom].atom_name<<"  "<<nbonds<<endl;

     ring_atom_list.push_back( current_atom );
     }
   } 

  
   int neighbor = 0;

   for( unsigned int neighborNumber = 0; neighborNumber < (site_info[current_atom].neighbor_list).size(); neighborNumber++ ){
     neighbor = site_info[current_atom].neighbor_list[neighborNumber];
    
     if( site_info[neighbor].seq_number  == current_residue &&
 	site_info[neighbor].FIRST_chain_ID == current_chain_ID &&
 	!site_info[neighbor].checked ){
       list_ring_atoms( neighbor, current_residue, current_chain_ID );
     }
   } 
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description: check over the new_bond vector of stacked_rings to see if the
//   residues a pair of atoms belong to are in a stacked ring.
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::in_ring_stacked_pair(int atom_1, int atom_2) {
  
  if(parameters.skip_identify_aromatics) return(false);
  vector<new_bonds>::iterator current_stack = stacked_rings.begin();
  while( current_stack != stacked_rings.end() ){

    int residue_1 = site_info[current_stack->site_1].seq_number;
    int residue_2 = site_info[current_stack->site_2].seq_number;
    int chain_1 = site_info[current_stack->site_1].FIRST_chain_ID;
    int chain_2 = site_info[current_stack->site_2].FIRST_chain_ID;
    if( site_info[atom_1].seq_number == residue_1 && site_info[atom_1].FIRST_chain_ID == chain_1 && 
	site_info[atom_2].seq_number == residue_2 && site_info[atom_2].FIRST_chain_ID == chain_2 ) {
      return(true);
    }
    else if( site_info[atom_1].seq_number == residue_2 && site_info[atom_1].FIRST_chain_ID == chain_2 && 
	site_info[atom_2].seq_number == residue_1 && site_info[atom_2].FIRST_chain_ID == chain_1 ) {
      return(true);
    }
    current_stack++;
  }
  return(false);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::read_ringstack_file(){

  string filename = path_name + "stacked.in";
  string linebuf;
  ifstream stacked_list;
  //  ifstream old_file_name;

  // If we are in interactive mode, prompt the user for the name of 
  // the file listing their stacked aromatic rings. If we're in non-interactive
  // mode, look for a file named stacked.in. If this file does not 
  // exist, exit. 
  //////////////////////////////////////////////////////////////////////
  stacked_list.open( filename.c_str() );


  // 
  //////////////////////////////////////////////////////////////////////
  if( !stacked_list ){
    clear_screen;
    cout << endl;
    cout << " The file " << filename << " was not found in this directory." << endl << endl;
    do{
      stacked_list.clear();
      cout << " Please enter the name of the file that contains the list of stacked rings." << endl;
      cout << " Filename = ";
      getline(cin, filename);
      
      stacked_list.open( filename.c_str(), ios::in );
      if( !stacked_list || did_press_enter(filename) ){
	clear_screen;
	cout << endl << "File not found." << endl << endl;
      }
      
    } while( !stacked_list || did_press_enter(filename) );
  }
  
  // Read each line of the stacked_list file. The first number should
  // correspond to a ring anchor atom, and the second atom to the other
  // ring anchor atom. The third number allows the user to specify the
  // number of bars for the given interaction (the default is 3). 
  //////////////////////////////////////////////////////////////////////
  unsigned int atom_1 = 0;
  unsigned int atom_2 = 0;
  int bars   = parameters.rs_bars;
  size_t field_start = 0;
  size_t field_end = 0;
  string number;

  while( !stacked_list.eof() ){
    getline(stacked_list, linebuf);

    // Skip blank lines.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos ||
	isComment(linebuf) )
      break;

    // Read in site 1 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", 0);
    field_end   = linebuf.find_first_of(" \t\n", field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    atom_1 = atoi( number.c_str() );
    
    // Read in site 2 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    atom_2 = atoi( number.c_str() );
    
    // OPTIONAL Read the number of bars to use in the constraint, if it is listed. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    if( field_start != string::npos ){
      number = linebuf.substr( field_start, field_end-field_start );
      if( number.find("123456789") == string::npos)
	bars = atoi( number.c_str() );
      else
	bars = parameters.rs_bars;

      if( bars < 0 ||
	  bars > 6 ){
	cout << " ERROR: Found a stacked ring interaction with an invalid number of bars while reading" << endl;
	cout << "        the file " << filename << "." << endl;
	cout << "        INVALID LINE: " << linebuf << endl;
	exit(1);
      }
    }

    // Map the original numbers into the numbers used internally by FIRST.
    //////////////////////////////////////////////////////////////////////    
    atom_1 = orig_2_FIRST_atom_number[atom_1];
    atom_2 = orig_2_FIRST_atom_number[atom_2];
    
    if( atom_1 <= 0 || atom_1 > total_sites ||
	atom_2 <= 0 || atom_2 > total_sites ){
      cout << " ERROR: An atom in the file [" << filename << "] does not exist." << endl;
      cout << "        The error occurred with the following line:" << endl;
      cout << "        " << linebuf << endl << endl;
      exit(1);
    }
    cout << atom_1 << "  "<< atom_2 << " "<<bars<<endl;
    stacked_rings.push_back( new_bonds(atom_1, atom_2, bars) );
  }

  if( parameters.verbose )
    cout << " Read " << stacked_rings.size() << " stacked rings from external file." << endl;
  stacked_list.close();

  return;
}
////////////////////////////////////////////////////////////////////////////////
