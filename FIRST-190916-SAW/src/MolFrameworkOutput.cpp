#include "global_defs.h"
#include "Parameters.h"
#include "output.h"
#include "generalUtils.h"
#include "hybrid_36_c.h"

extern const Parameters parameters;
map<int,string> color_list;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Output the original PDB file with the occupancy and temperature factor fields
//   substituted by the collective motion label and rigid region labels, respectively. 
//   The PDB file format only allows for numbers <= 999.99 in these fields, any labels 
//   larger than 999 are reset to 999. An optional bool flag can be sent to this
//   routine, in which the columns the labels are printed are switched (ie. the 
//   rigid region label will now be in the occupancy field).
////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputRCD_PDBFormat( string file_extension, bool RCD_in_occu ){

  float rigid_label  = 0.0;
  float coll_mode_label = 0.0;

  string output_file = path_name + base_name;
  output_file += file_extension;

  ofstream R_decomp( output_file.c_str() );
  ostringstream FIRST_group_ID, atom_number, seg_id;
  float bvalue, occu;
  string output_element_name;

  string temp1;
  int int_atom_number;
  char write_atom_number[6];
  const char* errmsg;

  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    site_info[siteNumber].rigid_label < 999 ? rigid_label = site_info[siteNumber].rigid_label : rigid_label = 999.0;
    site_info[siteNumber].coll_mode_label < 999 ? coll_mode_label = site_info[siteNumber].coll_mode_label : coll_mode_label = 999.0;
    
    output_element_name = site_info[siteNumber].element_name;
    if(output_element_name[1] == ' ') swap(output_element_name[0], output_element_name[1]);
    
    FIRST_group_ID.str("");
    if (parameters.use_group_id) {
        FIRST_group_ID << setw(6) << site_info[siteNumber].FIRST_group_ID;
    }
      

    atom_number.str("");
    temp1 =   site_info[siteNumber].orig_atom_number;
    if( isNumber( temp1))
      {
	int_atom_number = atoi( temp1.c_str() );
	if(  int_atom_number <= 99999){
	  atom_number << setw(5) << site_info[siteNumber].orig_atom_number;
	}
	else{
	  hy36encode(5, int_atom_number, write_atom_number);
	  atom_number << setw(5) <<  write_atom_number;
	}
      }
    else{
      errmsg= hy36encode(5, int_atom_number, write_atom_number);
      if(!errmsg)
	atom_number << setw(5) << write_atom_number;
      else
	atom_number << setw(5) << site_info[siteNumber].orig_atom_number;
    }
    
    if (RCD_in_occu) {
      occu = rigid_label;
      bvalue = coll_mode_label;
    }
    else {
      bvalue = rigid_label;
      occu = coll_mode_label;
    }
    
    seg_id.str("");
    seg_id << setiosflags(ios_base::left)
    << setw(4) << site_info[siteNumber].seg_id;
      
    R_decomp  << showpoint << setiosflags(ios::fixed) 
    << setw(6) << site_info[siteNumber].record_name
    << atom_number.str()
    << setw(5) << atomNamePDBFormat(site_info[siteNumber].atom_name, site_info[siteNumber].element_name)
    << setw(4) << site_info[siteNumber].residue_name
    << setw(2) << char(site_info[siteNumber].chain_ID)
    << setw(4) << site_info[siteNumber].seq_number
    << setw(1) << site_info[siteNumber].insert_code
    << "   "
    << setw(8) << setprecision(3) << site_info[siteNumber].coords[X]
    << setw(8) << setprecision(3) << site_info[siteNumber].coords[Y]
    << setw(8) << setprecision(3) << site_info[siteNumber].coords[Z]
    << setw(6) << setprecision(2) << occu
    << setw(6) << setprecision(2) << bvalue // cols 61-66
    << "      " // cols 67-72
    << seg_id.str() // cols 73-76
    << setw(2) << output_element_name // cols 77-78
    << "   " // cols 79-81
    << FIRST_group_ID.str()
    << endl;

    if( find(chainTermini.begin(), chainTermini.end(), siteNumber) != chainTermini.end() ){
      R_decomp << "TER" << endl;
    }
  } 

  outputCONECTRecords(R_decomp);

  R_decomp << "END" << endl;
  R_decomp.close();
}    
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputPyMolScriptRCD_PDBFormat(){

  string output_file = path_name + base_name;
  output_file += "_RCD.pml";
  ofstream pymol_script( output_file.c_str() );

  //float stick_radius = 0.25; // FIXME - unused variable
  string render_style;

  ostringstream flexible; 
  flexible << "create Flexible, ";

  pymol_script << "# Rigid cluster decomposition coloring script for PyMol" << endl << "#" << endl;
  pymol_script << "# Created by Dan Farrell, Brandon Hespenheide." << endl;
  pymol_script << "# Department of Physics and Astronomy" << endl;
  pymol_script << "# Biophysics Theory Group" << endl;
  pymol_script << "# Arizona State University" << endl;
  pymol_script << "############################################################" << endl << endl;
 
  string pdb_file = base_name;
  pdb_file += parameters.run_number + "_RCD";

  pymol_script << "from pymol import cmd" << endl;
  pymol_script << "from pymol.cgo import *" << endl << endl;

  // Some final global attributes to set.
  //////////////////////////////////////////////////////////////////////
  pymol_script << "bg_color white" << endl;

  //Here, we either load all the froda conformers (if FRODA was being used)
  //Or we load just the single structure (if FRODA was not used)
  
  // Glob all the FRODA structures together for movie generation.
  //////////////////////////////////////////////////////////////////////
  if( parameters.runFRODA ||
      parameters.runMorph ){
    pymol_script << "# script to make a movie of the FRODA generated structures." << endl;
    pymol_script << endl << "from glob import glob" << endl;
    pymol_script << "filelist = glob (\"" << base_name << "_froda_*.pdb\") " << endl;
    pymol_script << "filelist.sort()" << endl;
    pymol_script << "cmd.load(\"" << pdb_file << ".pdb\", \"" << base_name << "_froda\")" << endl;
    pymol_script << "for file in filelist: cmd.load( file, \"" << base_name << "_froda\") " << endl;
    pymol_script << "show lines, " << base_name << "_froda" << endl;
    pymol_script << "color black" << endl;
  }
  else {
    // Color the whole structure black with thin lines. This represents the
    // underlying flexible network.
    //////////////////////////////////////////////////////////////////////
    pymol_script << "load " << pdb_file << ".pdb" << endl << endl;
    pymol_script << "set line_width = 1" << endl;
    pymol_script << "color black" << endl << endl;
  }

  // For each rigid cluster, calculate a unique color, and have pymol
  // assign that color to those atoms.
  // Also, we save each cluster's color in a map for use again later
  map< unsigned int, string > mapClusterIDtoColor;
  next_rcd_color_as_name(true);    
  for( unsigned int clusterLargerThanOneNumber = 0;
        clusterLargerThanOneNumber < total_clusters_larger_than_one;
        clusterLargerThanOneNumber++ ){
    if( RC_atom_list[clusterLargerThanOneNumber+1].size() >= parameters.min_output_cluster_size ){
      string color = next_rcd_color_as_name();
      mapClusterIDtoColor[clusterLargerThanOneNumber+1] = color;
      pymol_script << "color " << color << ", ( b > " << float(clusterLargerThanOneNumber+1-0.01) 
       << " and b < " << float(clusterLargerThanOneNumber+1+0.01) << ")" << endl;
    }
  }

  // Python commands to draw hbonds and hydrophobic tethers
  // Do not create the object if there are more
  // than 1000 hydrogen bonds in the structure.
  // If multiple FRODA conformers are loaded, then these pymol distance
  // objects will move with the movie.
  //////////////////////////////////////////////////////////////////////
  if( hydrogen_bonds.size() < 1000 ){
    pymol_script << "# Draw hbonds, hydrophoic tethers and stacked rings as distance objects" << endl;
    pymol_script << "set dash_gap, 0.1" << endl;

    int site_1, site_2;
    string origID1, origID2;
    
    for (unsigned int hydrogenBondNumber=0; hydrogenBondNumber < hydrogen_bonds.size(); hydrogenBondNumber++) {
      if (hydrogen_bonds[hydrogenBondNumber].bond_ranking <= parameters.energy_cutoff) {
        site_1 = hydrogen_bonds[hydrogenBondNumber].site_1;
        site_2 = hydrogen_bonds[hydrogenBondNumber].site_2;
        origID1 = site_info[site_1].orig_atom_number;
        origID2 = site_info[site_2].orig_atom_number;
        pymol_script << "distance hbonds = id " << origID1 << " , id " << origID2 << endl;
      }
    }
    pymol_script << "color red, hbonds" << endl;
    pymol_script << "hide labels, hbonds" << endl;
    pymol_script << "disable hbonds" << endl; 


    for (unsigned int hydrophobicTetherNumber=0; hydrophobicTetherNumber < hydrophobic_tethers.size(); hydrophobicTetherNumber++) {
      site_1 = hydrophobic_tethers[hydrophobicTetherNumber].site_1;
      site_2 = hydrophobic_tethers[hydrophobicTetherNumber].site_2;
      origID1 = site_info[site_1].orig_atom_number;
      origID2 = site_info[site_2].orig_atom_number;
      pymol_script << "distance hydrophobics = id " << origID1 << " , id " << origID2 << endl;
    } 
    if( hydrophobic_tethers.size() ) {
    pymol_script << "color green, hydrophobics" << endl;
    pymol_script << "hide labels, hydrophobics" << endl;
    pymol_script << "disable hydrophobics" << endl; 
    }


    for (unsigned int stackedRingNumber=0; stackedRingNumber < stacked_rings.size(); stackedRingNumber++) {
      site_1 = stacked_rings[stackedRingNumber].site_1;
      site_2 = stacked_rings[stackedRingNumber].site_2;
      origID1 = site_info[site_1].orig_atom_number;
      origID2 = site_info[site_2].orig_atom_number;
      pymol_script << "distance stackedrings = id " << origID1 << " , id " << origID2 << endl;
    } 
    if( stacked_rings.size() ){
      pymol_script << "color blue, stackedrings" << endl;
      pymol_script << "hide labels, stackedrings" << endl;
      pymol_script << "disable stackedrings" << endl; 
    }
  }
  // Create the rigid cluster objects for pymol, and color them. Only those
  // clusters larger than the min_output_cluster_size will have objects 
  // created for them. The remaining clusters will be binned into bulk
  // objects. 
  //////////////////////////////////////////////////////////////////////
  if( total_sites < 30000 )
    render_style = "sticks";
  else
    render_style = "lines";

  // for each rigid cluster, create a new pymol object of the cluster.
  // This overlaid pymol object will be drawn thicker, using sticks.
  // If FRODA was used, then these cluster objects will move with the movie.
  unsigned int total_RC_objects = 0;
  for( unsigned int clusterLargerThanOneNumber = 0; clusterLargerThanOneNumber < total_clusters_larger_than_one; clusterLargerThanOneNumber++ ){
    if( RC_atom_list[clusterLargerThanOneNumber+1].size() >= parameters.min_output_cluster_size ){
      pymol_script << "# Rigid Cluster " << clusterLargerThanOneNumber+1 << " has " 
                   << RC_atom_list[clusterLargerThanOneNumber+1].size() << " atoms." << endl;
      pymol_script << "create RC" << clusterLargerThanOneNumber+1 << ", ( b > " << float(clusterLargerThanOneNumber+1-0.01) 
                   << " and b < " << float(clusterLargerThanOneNumber+1+0.01) << ")" << endl;
      pymol_script << "show " << render_style << ", RC" << clusterLargerThanOneNumber+1 << endl;
      pymol_script << "set line_width = 3, " << "RC" << clusterLargerThanOneNumber+1 << endl;
      pymol_script << "color " << mapClusterIDtoColor[clusterLargerThanOneNumber+1] << ", RC" << clusterLargerThanOneNumber+1 << endl << endl;
      total_RC_objects++;
    }
  }

  // Bin the remaining rigid clusters into objects. Each bin will contain 
  // all the rigid clusters of a certain size (or range of sizes).
  //////////////////////////////////////////////////////////////////////
  //int bin_size = 10;// FIXME - unused variable
  unsigned int lower_bound = 0;
  int lower_b_cutoff = 0;
  
  if( parameters.min_output_cluster_size > 20 ){
    lower_bound = 20;
  }
  else if( parameters.min_output_cluster_size > 10 ){
    lower_bound = 10;
  }
  else{
    lower_bound = 1;
  }

  lower_bound += 10;

  do {
    lower_bound -= 10;

    lower_b_cutoff = total_RC_objects;

    while( RC_atom_list[total_RC_objects].size() > (unsigned int)lower_bound && // FIXME - warning: comparison between signed and unsigned integer expressions (no easy fix due to lower_bound = -1;)
      total_RC_objects < total_clusters ){
      total_RC_objects++;
    }

    int bin_label = lower_bound/10 +1;
    pymol_script << "# Rigid Cluster BIN" << bin_label << endl;
    if( total_RC_objects > 999 ){
      pymol_script << "create BIN" << bin_label << ", ( b > " << float(lower_b_cutoff+1-0.01) 
                   << " and b < 999.00 )" << endl;
    }
    else
      pymol_script << "create BIN" << bin_label << ", ( b > " << float(lower_b_cutoff+1-0.01) 
                   << " and b < "  << float(total_RC_objects+0.01) << ")" << endl;

    pymol_script << "show " << render_style << ", BIN" << bin_label << endl;
    pymol_script << "set line_width = 3, " << "BIN" << bin_label << endl;
    pymol_script << "color gray, BIN" << bin_label << endl;
    pymol_script << "disable BIN" << bin_label << endl << endl;

  } while ((lower_bound >= 10) && (total_RC_objects <= 999 ));

  
  pymol_script.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputPyMolScriptCMD_PDBFormat(){

  //float stick_radius = 0.25; // FIXME - unused variable
  string render_style;

  ostringstream flexible; 
  flexible << "create Flexible, ";

  string file_name = path_name + base_name;
  file_name += "_CMD.pml";
  ofstream pymol_script( file_name.c_str() );

  pymol_script << "# Rigid cluster decomposition coloring script for PyMol" << endl << "#" << endl;
  pymol_script << "# Created by Dan Farrell, Brandon Hespenheide." << endl;
  pymol_script << "# Department of Physics and Astronomy" << endl;
  pymol_script << "# Biophysics Theory Group" << endl;
  pymol_script << "# Arizona State University" << endl;
  pymol_script << "############################################################" << endl << endl;
 
  string pdb_file = base_name;
  pdb_file += parameters.run_number + "_RCD";

  pymol_script << "from pymol import cmd" << endl;
  pymol_script << "from pymol.cgo import *" << endl << endl;

  // Color the whole structure black with thin lines. This represents the
  // underlying flexible network.
  //////////////////////////////////////////////////////////////////////
  pymol_script << "load " << pdb_file << ".pdb" << endl << endl;
  pymol_script << "set line_width = 1" << endl;
  pymol_script << "color black" << endl << endl;

  // Python commands to draw hbonds and hydrophobic tethers as Compiled
  // Graphics Objects (CGO). Do not create the object if there are more
  // than 1000 hydrogen bonds in the structure.
  //////////////////////////////////////////////////////////////////////
  if( hydrogen_bonds.size() < 1000 ){
    pymol_script << "# Draw hbonds and hydrophoic tethers as CGO objects" << endl;
    pymol_script << "# The color of the objects can be modified below," << endl;
    pymol_script << "# by entering rgb values between 0.0 and 1.0." << endl;
    pymol_script << "hbondobj = [BEGIN,LINES,COLOR,0.0,0.0,0.0,]" << endl;
    pymol_script << "hydrophobicobj = [BEGIN,LINES,COLOR,0.0,0.0,0.0,]" << endl;
    
    for (unsigned int hydrogenBondNumber=0; hydrogenBondNumber < hydrogen_bonds.size(); hydrogenBondNumber++) {
      int site_1, site_2; float x1, y1, z1, x2, y2, z2;
      if (hydrogen_bonds[hydrogenBondNumber].bond_ranking <= parameters.energy_cutoff) {
        site_1 = hydrogen_bonds[hydrogenBondNumber].site_1;
        site_2 = hydrogen_bonds[hydrogenBondNumber].site_2;
        x1=site_info[site_1].coords[X];
        y1=site_info[site_1].coords[Y];
        z1=site_info[site_1].coords[Z];
        x2=site_info[site_2].coords[X];
        y2=site_info[site_2].coords[Y];
        z2=site_info[site_2].coords[Z];
        pymol_script << "hbondobj.extend([VERTEX," << x1 << "," << y1 << "," << z1 << ",";
        pymol_script << "VERTEX," << x2 << "," << y2 << "," << z2 << ",])" << endl;
      }
    }
    pymol_script << "hbondobj.extend([END])" << endl;
    pymol_script << "cmd.load_cgo(hbondobj,'hbonds')" << endl;
    pymol_script << "disable hbonds" << endl; 
    
    for (unsigned int hydrophobicTetherNumber=0; hydrophobicTetherNumber < hydrophobic_tethers.size(); hydrophobicTetherNumber++) {
      int site_1, site_2; float x1, y1, z1, x2, y2, z2;
      site_1 = hydrophobic_tethers[hydrophobicTetherNumber].site_1;
      site_2 = hydrophobic_tethers[hydrophobicTetherNumber].site_2;
      x1=site_info[site_1].coords[X];
      y1=site_info[site_1].coords[Y];
      z1=site_info[site_1].coords[Z];
      x2=site_info[site_2].coords[X];
      y2=site_info[site_2].coords[Y];
      z2=site_info[site_2].coords[Z];
      pymol_script << "hydrophobicobj.extend([VERTEX," << x1 << "," << y1 << "," << z1 << ",";
      pymol_script << "VERTEX," << x2 << "," << y2 << "," << z2 << ",])" << endl;
    }
    pymol_script << "hydrophobicobj.extend([END])" << endl;
    pymol_script << "cmd.load_cgo(hydrophobicobj,'hydrophobic')" << endl;
    pymol_script << "disable hydrophobic" << endl << endl;
  }

  // Create the rigid cluster objects for pymol, and color them. Only those
  // clusters larger than the min_output_cluster_size will have objects 
  // created for them. The remaining clusters will be binned into bulk
  // objects. 
  //////////////////////////////////////////////////////////////////////
  if( total_sites < 30000 )
    render_style = "sticks";
  else
    render_style = "lines";

  unsigned int total_RC_objects = 0;
  for( unsigned int siteLargerThanOneNumber = 0; siteLargerThanOneNumber < total_clusters_larger_than_one; siteLargerThanOneNumber++ ){
    if( RC_atom_list[siteLargerThanOneNumber+1].size() >= parameters.min_output_cluster_size ){
      pymol_script << "# Rigid Cluster " << siteLargerThanOneNumber+1 << " has " << RC_atom_list[siteLargerThanOneNumber+1].size() << " atoms." << endl;
      pymol_script << "create RC" << siteLargerThanOneNumber+1 << ", ( b > " << float(siteLargerThanOneNumber+1-0.01) << " and b < " << float(siteLargerThanOneNumber+1+0.01) << ")" << endl;
      pymol_script << "show " << render_style << ", RC" << siteLargerThanOneNumber+1 << endl;
      pymol_script << "set line_width = 3, " << "RC" << siteLargerThanOneNumber+1 << endl;
      pymol_script << "color " << next_rcd_color_as_name() << ", RC" << siteLargerThanOneNumber+1 << endl << endl;
      total_RC_objects++;
    }
  }

  // Bin the remaining rigid clusters into objects. Each bin will contain all
  // the rigid clusters of a certain size (or range of sizes).
  //////////////////////////////////////////////////////////////////////
  unsigned int lower_bound = 0;
  int lower_b_cutoff = 0;
  
  if( parameters.min_output_cluster_size > 20 ){
    lower_bound = 20;
  }
  else if( parameters.min_output_cluster_size > 10 ){
    lower_bound = 10;
  }
  else{
    lower_bound = 1;
  }

  lower_bound += 10;
  
  do {
    lower_bound -= 10;

    lower_b_cutoff = total_RC_objects;

    while( RC_atom_list[total_RC_objects].size() > lower_bound && // FIXME - warning: comparison between signed and unsigned integer expressions (no easy fix due to lower_bound = -1;)
      total_RC_objects < total_clusters ){
      total_RC_objects++;
    }

    int bin_label = lower_bound/10 +1;
    pymol_script << "# Rigid Cluster BIN" << bin_label << endl;
    if( total_RC_objects > 999 ){
      pymol_script << "create BIN" << bin_label << ", ( b > " << float(lower_b_cutoff+1-0.01) 
                   << " and b < 999.00 )" << endl;
    }
    else
      pymol_script << "create BIN" << bin_label << ", ( b > " << float(lower_b_cutoff+1-0.01) 
                   << " and b < "  << float(total_RC_objects+0.01) << ")" << endl;

    pymol_script << "show " << render_style << ", BIN" << bin_label << endl;
    pymol_script << "set line_width = 3, " << "BIN" << bin_label << endl;
    pymol_script << "color gray, BIN" << bin_label << endl;
    pymol_script << "disable BIN" << bin_label << endl << endl;

  } while ((lower_bound >= 10) && (total_RC_objects <= 999));

  // Some final global attributes to set.
  //////////////////////////////////////////////////////////////////////
  pymol_script << "bg_color white" << endl;
  pymol_script << "clip slab, 200" << endl;
  pymol_script << "center all" << endl;
  pymol_script << "zoom" << endl;

  pymol_script.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputCovalentBondList(){

  string filename;
  vector<unsigned int>::iterator next_bond;

  if( parameters.interactive &&
      !parameters.covout ){
    cout << " Please enter a file name for the list of covalent bonds: ";
    cin  >> filename;
    filename = path_name + filename;
  }
  else{
    filename = path_name + "cov.out";
  }

  ofstream output( filename.c_str() );
  SiteID site_1;
  SiteID site_2;
  string site_1o;
  string site_2o;

  for( unsigned int siteNumber = 0; siteNumber < total_sites; siteNumber++ ){

    sort( (site_info[siteNumber].neighbor_list).begin(), (site_info[siteNumber].neighbor_list).end() );
    next_bond = (site_info[siteNumber].neighbor_list).begin();
    while( next_bond != (site_info[siteNumber].neighbor_list).end() ){
      
      if( *next_bond > siteNumber ){ // FIXME - warning: comparison between signed and unsigned integer expressions 
				
				if (parameters.use_first_numbering) {
					site_1 = site_info[siteNumber].FIRST_number;
					site_2 = site_info[*next_bond].FIRST_number;
					output << setw(8) << site_1    << setw(8) << site_2
					       << setw(8) << site_info[siteNumber].number_of_bars[*next_bond]
					       << endl;

					
				} else {
					site_1o = site_info[siteNumber].orig_atom_number;
					site_2o = site_info[*next_bond].orig_atom_number;
					output << setw(8) << site_1o    << setw(8) << site_2o
					       << setw(8) << site_info[siteNumber].number_of_bars[*next_bond]
					       << endl;
				}
				
      }
      next_bond++;
    }
    
  } 

  output.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Create a file listing the hydrogen bonds identified by FIRST. The file will
//   be named "hbonds.out", and it will overwrite previous data. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputHydrogenBondList( bool full_output ){

  string filename;

  if( parameters.interactive &&
      !parameters.covout ){
    cout << " Please enter a file name for the list of hydrogen bonds: ";
    cin  >> filename;
    filename = path_name + filename;
  }
  else{
    filename = path_name + "hbonds.out";
  }
  ofstream output( filename.c_str() );

  stable_sort( hydrogen_bonds.begin(), hydrogen_bonds.end(), greater<new_bonds>() );

  for( unsigned int hydrogenBondNumber = 0; hydrogenBondNumber < hydrogen_bonds.size(); hydrogenBondNumber++ ){
		
		SiteID donorID;
		SiteID acceptorID;
		string str_donorID;
		string str_acceptorID;


		if (parameters.use_first_numbering) {
			donorID = hydrogen_bonds[hydrogenBondNumber].site_1;
			acceptorID = hydrogen_bonds[hydrogenBondNumber].site_2;
			output << setiosflags(ios::fixed) 
			       << setw(8)  << donorID
			       << setw(8)  << acceptorID
			       << setw(16) << setprecision(8) << hydrogen_bonds[hydrogenBondNumber].energy 
			       << setw(5)  << hydrogen_bonds[hydrogenBondNumber].bars
			       << endl;			
		} else {
			str_donorID = get_orig_atom_number( hydrogen_bonds[hydrogenBondNumber].site_1 );
			str_acceptorID = get_orig_atom_number( hydrogen_bonds[hydrogenBondNumber].site_2 );
			output << setiosflags(ios::fixed) 
			       << setw(8)  << str_donorID
			       << setw(8)  << str_acceptorID
			       << setw(16) << setprecision(8) << hydrogen_bonds[hydrogenBondNumber].energy 
			       << setw(5)  << hydrogen_bonds[hydrogenBondNumber].bars
			       << endl;

		}
		
  }

  output.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputHydrophobicTetherList(){
  
  string filename;
  if( parameters.interactive &&
      !parameters.covout ){
    cout << " Please enter a file name for the list of hydrophobic tethers: ";
    cin  >> filename;
    filename = path_name + filename;
  }
  else{
    filename = path_name + "hphobes.out";
  }

  ofstream output( filename.c_str() );
		
  for( unsigned int hydrophobicTetherNumber = 0; hydrophobicTetherNumber < hydrophobic_tethers.size(); hydrophobicTetherNumber++ ){
		SiteID site_1;
		SiteID site_2;
		string str_site_1;
		string str_site_2;

		if (parameters.use_first_numbering) {
			site_1 = hydrophobic_tethers[hydrophobicTetherNumber].site_1;
			site_2 = hydrophobic_tethers[hydrophobicTetherNumber].site_2;
			output << setw(8) << site_1
			       << setw(8) << site_2 
			       << setw(8) << hydrophobic_tethers[hydrophobicTetherNumber].bars
			       << endl;			
		} else {
			str_site_1 = get_orig_atom_number( hydrophobic_tethers[hydrophobicTetherNumber].site_1 );
			str_site_2 = get_orig_atom_number( hydrophobic_tethers[hydrophobicTetherNumber].site_2 );
			output << setw(8) << str_site_1
			       << setw(8) << str_site_2 
			       << setw(8) << hydrophobic_tethers[hydrophobicTetherNumber].bars
			       << endl;
		}
		
  }
  output.close();

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputStackedRingsList(){
  
  string filename;
  if( parameters.interactive &&
      !parameters.covout ){
    cout << " Please enter a file name for the list of stacked aromatic rings: ";
    cin  >> filename;
    filename = path_name + filename;
  }
  else{
    filename = path_name + "stacked.out";
  }

  ofstream output( filename.c_str() );
	
  for( unsigned int stackedRingNumber = 0; stackedRingNumber < stacked_rings.size(); stackedRingNumber++ ){
		SiteID site_1;
		SiteID site_2;
		string str_site_1;
		string str_site_2;
		
		if (parameters.use_first_numbering) {
			site_1 = stacked_rings[stackedRingNumber].site_1;
			site_2 = stacked_rings[stackedRingNumber].site_2;
			output << setw(8) << site_1
			       << setw(8) << site_2 
			       << setw(8) << stacked_rings[stackedRingNumber].bars
			       << endl;
		} else {
			str_site_1 = get_orig_atom_number( stacked_rings[stackedRingNumber].site_1 );
			str_site_2 = get_orig_atom_number( stacked_rings[stackedRingNumber].site_2 );
			output << setw(8) << str_site_1
			       << setw(8) << str_site_2 
			       << setw(8) << stacked_rings[stackedRingNumber].bars
			       << endl;
		}
  }
  output.close();

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::MakeTmpFileForHBDilute(){

  string output_file;

  if( parameters.bond_dilution ){
    output_file = path_name + "hbdil_temp.pdb";
    ofstream hbdil_tmp_file( output_file.c_str() );
    for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){

      hbdil_tmp_file << showpoint << setiosflags(ios::fixed) 
                     << setw(6) << site_info[siteNumber].record_name 
                     << setw(5) << site_info[siteNumber].FIRST_number
                     << setw(5) << site_info[siteNumber].atom_name
                     << setw(4) << site_info[siteNumber].residue_name
                     << setw(2) << char(site_info[siteNumber].chain_ID)
                     << setw(4) << site_info[siteNumber].seq_number
                     << setw(1) << site_info[siteNumber].insert_code
                     << "   "
                     << setw(8) << setprecision(3) << site_info[siteNumber].coords[X]
                     << setw(8) << setprecision(3) << site_info[siteNumber].coords[Y]
                     << setw(8) << setprecision(3) << site_info[siteNumber].coords[Z]
                     << endl;
    }
    hbdil_tmp_file.close();
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The textual output from FIRST is split into four sections: File Info, Input 
//   Parameters, Results, and Warnings. Within each section there are four levels 
//   of output detail. The dafualt is level 1, which contains the minimal amount
//   of information. The highest level is level 4, which is the most verbose output
//   level. 
// Parameters:
//   output_level - Sets how much information is to be included in the output file.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputText(){

  if( parameters.output_level == 0 )
    return;

  string output_file = path_name + base_name + "_results.txt";
  ofstream output( output_file.c_str(),  ios::app );


  // FILE INFO SECTION
  ////////////////////////////////////////////////////////////////////////////////
  output << hr() << endl;
  output << "# FILE INFO:" << endl << endl;
  output << " Input file name: " << infile_name << endl;
  output << " Date: " << parameters.start_time << endl;
  
  output << hr() << endl;

  // INPUT PARAMETERS SECTION
  ////////////////////////////////////////////////////////////////////////////////
  output << "# INPUT PARAMETERS:" << endl << endl;
  output << " Command line: \"" << parameters.command_line << "\"" << endl;
  if( is_nmr )
    output  << " NMR MODEL number: " << model_number << endl;
  output << " Energy cutoff: " << showpoint << setprecision(3) << parameters.energy_cutoff << " kcal/mol" << endl;

  if( parameters.output_level >= 2 ){
    output << " Interactive mode: " << boolalpha << parameters.interactive << endl;
    output << " Hydrophobic tether parameters:" << endl;
    output << " Hydrophobic function: " << parameters.hphobe_fxn << endl;
    output << " Hydrophobic distance cutoff: " << parameters.PH_cutoff << endl;
    output << " Hydrogen bond parameters:" << endl;
    output << " Hydrogen - Acceptor distance cutoff: " << parameters.cutoff_HB_hyd_accpt_dist << endl;
    output << " Donor - Acceptor distance cutoff: " << parameters.cutoff_HB_donor_accpt_dist << endl;
    output << " Donor - Hydrogen - Acceptor angle cutoff: " << parameters.cutoff_HB_donor_hyd_accpt_angle << endl;

    output << endl << " Hydrogen bonds included in the network." << endl;
    output << included_hbonds.str() << endl;

    output << endl << " Hydrophobic tethers included in the network." << endl;
    output << included_hphobes.str() << endl;

    output << endl << " Stacked aromatic rings included in the network." << endl;
    output << included_aromats.str() << endl;
  }

  if( parameters.output_level >= 3 ){
    output << endl << " Details of all hydrogen bonds." << endl;
    output << hbond_data.str() << endl;
  }

  output << endl << hr() << endl;

  // RESULTS SECTION
  ////////////////////////////////////////////////////////////////////////////////
  output << "# RESULTS:" << endl << endl;

  output << setw(15) << total_sites << " Sites included in the calculation." << endl;
  output << setw(15) << total_excluded_sites << " Sites excluded in the calculation." << endl;
  output << setw(15) << total_residues << " Unique residues found." << endl;

  if( isolated_sites == 1 )
    output << setw(15) << isolated_sites << " Isolated site." << endl;
  else
    output << setw(15) << isolated_sites << " Isolated sites." << endl;
  if( isolated_dimers == 1 )
    output << setw(15) << isolated_dimers << " Isolated dimer." << endl;
  else
    output << setw(15) << isolated_dimers << " Isolated dimers." << endl;
  output << setw(15) << total_hbonds << " Hydrogen bonds included." << endl;
  output << setw(15) << total_hphobes << " Hydrophobic tethers included." << endl;
  output << setw(15) << total_stacked_rings << " Stacked ring interactions included." << endl;
  output << setw(15) << total_clusters << " Rigid clusters." << endl;
  output << setw(15) << total_clusters_larger_than_min_cluster_size << " Rigid clusters larger than size " << parameters.min_output_cluster_size<<"."<<endl;
  int atoms_in_largestRC = 1;
  int atoms_in_second_largestRC = 1;
  int rcNumber_largestRC = 0;
  for( unsigned int rcNumber = 1; rcNumber < RC_atom_list.size(); rcNumber++ ){
    if( RC_atom_list[rcNumber].size() > atoms_in_largestRC ){
      atoms_in_largestRC =  RC_atom_list[rcNumber].size();
      rcNumber_largestRC = rcNumber;
    }
  }

  for( unsigned int rcNumber = 1; rcNumber < RC_atom_list.size(); rcNumber++ ){
    if(RC_atom_list[rcNumber].size() > atoms_in_second_largestRC && RC_atom_list[rcNumber].size() <= atoms_in_largestRC && rcNumber != rcNumber_largestRC ){
      atoms_in_second_largestRC = RC_atom_list[rcNumber].size();
    }
  }
  output << setw(15) << atoms_in_largestRC << " Sites in largest rigid cluster."<<endl;
  output << setw(15) << atoms_in_second_largestRC << " Sites in second largest rigid cluster."<<endl;


  if( !parameters.bond_dilution && !parameters.inputGhosts ){
    output << setw(15) << total_stressed_regions << " Stressed Regions." << endl;
    output << setw(15) << floppy_modes << " Total independent degrees of freedom." << endl;
  }
  else{
    output << "            N/A Stressed Regions." << endl;
    output << "            N/A Total independent degrees of freedom." << endl;
  }
  
  if( parameters.output_level >= 2 ){
    output << endl;
    output << "  Cluster |      Size (in atoms)" << endl;
    output << " -------------------------------" << endl;
    for( unsigned int rcNumber = 1; rcNumber < RC_atom_list.size(); rcNumber++ ){
      if( RC_atom_list[rcNumber].size() > 1 )
        output << " " << setw(8) << rcNumber << " | " << setw(8) << RC_atom_list[rcNumber].size() << endl; 
    }
    output << endl;
    output << "(" << setw(8) << total_clusters -total_clusters_larger_than_one << " clusters of size one)" << endl;
  }  

  if( parameters.output_level >= 3 ){
  }

  output << endl << hr() << endl;

  // WARNINGS SECTION
  ////////////////////////////////////////////////////////////////////////////////
  output << "# WARNINGS: LEVEL 1 warnings are most severe." << endl << endl;
  
  output << " LEVEL 1 WARNINGS: " << warning_level_1.size() << endl;
  for( unsigned int a = 0; a < warning_level_1.size(); a++ )
    output << "* " << warning_level_1.at(a);
  output << endl;

  output << " LEVEL 2 WARNINGS: " << warning_level_2.size() << endl;
  for( unsigned int a = 0; a < warning_level_2.size(); a++ )
    output << "* " << warning_level_2.at(a);
  output << endl;
  
  output << " LEVEL 3 WARNINGS: " << warning_level_3.size() << endl;
  for( unsigned int a = 0; a < warning_level_3.size(); a++ )
    output << "* " << warning_level_3.at(a);
  output << endl;
  
  output << endl << hr() << endl;
  output.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
 void MolFramework::outputRawData(){

  if( !parameters.make_raw_data_file )
    return;

  string output_file = path_name + base_name + "_data.txt";
  ofstream output( output_file.c_str() );

  output << "#  Raw data from FIRST5" << endl;
  output << "#" 
         << setw(11) << "FIRST num"
         << setw(12) << "orig num"
         << setw(12) << "rigid"
         << setw(12) << "stressed"
         << setw(12) << "coll mode" 
         << endl;
  output << hr() << endl;

  for( unsigned int a = 1; a <= total_sites; a++ ){
    output << setw(12) << a
           << setw(12) << site_info[a].orig_atom_number 
           << setw(12) << site_info[a].rigid_label
           << setw(12) << site_info[a].stressed_label 
           << setw(12) << site_info[a].coll_mode_label;
    if( site_info[a].excluded )
      output << setw(12) << "EXCLUDED";

    output << endl;
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Output the flexibility results for each bond. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputBondData(){

  unsigned int neighbor_site = 0;

  string output_file = path_name + base_name + "_bond.txt";
  ofstream output( output_file.c_str() );

  output << "# FIRST flexibility index results. Data refer to bonds." << endl;
  output << "#" 
         << setw(11) << "original"
         << setw(12) << "original"
         << setw(12) << "flex."
         << setw(12) << "FIRST"
         << setw(12) << "FIRST"
         << endl;
  output << "#" 
         << setw(11) << "atom 1"
         << setw(12) << "atom 2"
         << setw(12) << "index"
         << setw(12) << "atom 1"
         << setw(12) << "atom 2"
         << endl;
  output << hr() << endl;

  for( unsigned int current_site = 1; current_site <= total_sites; current_site++ ){
    if( !site_info[current_site].excluded ){
      for( unsigned int a = 0; a < site_info[current_site].neighbor_list.size(); a++ ){
        neighbor_site = site_info[current_site].neighbor_list[a];
        if( neighbor_site > current_site &&
            !site_info[neighbor_site].excluded ){
 
          output << setw(12) << site_info[current_site].orig_atom_number 
                 << setw(12) << site_info[neighbor_site].orig_atom_number 
                 << setw(12) << fixed << showpoint << setprecision(4) << site_info[current_site].bond_info[neighbor_site].flexibility_index
                 << setw(12) << current_site
                 << setw(12) << neighbor_site 
                 << endl;
        }
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputRCD_MOL2FORMAT(){

  string output_file = path_name + base_name + "_RCD.mol2";
  ofstream output( output_file.c_str() );

  output << "#" << endl 
         << "# MOL2 FILE CREATED BY FIRST5.0" << endl
         << "# CREATED ON: " << parameters.start_time << endl
         << "#" << endl;
  
  output << "@<TRIPOS>MOLECULE" << endl;
  output << base_name << endl;
  output << total_sites << " " << total_bonds << endl;
  output << "SMALL" << endl;
  output << "NO_CHARGES" << endl << endl;
  
  // Output the atoms in the file.
  ////////////////////////////////////////////////////////////////////////////////
  short int width_1 = int( ceil( log(double(total_sites))) );

  output << "@<TRIPOS>ATOM" << endl;
  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    if( RC_atom_list[site_info[siteNumber].rigid_label].size() >= parameters.min_output_cluster_size ){

      if( site_info[siteNumber].surface_site )
        output << showpoint << setiosflags(ios::fixed)
               << setw(width_1 +2) << left << site_info[siteNumber].orig_atom_number << right
               << setw(3) << site_info[siteNumber].element_name
               << setw(10) << setprecision(3) << site_info[siteNumber].coords[X]
               << setw(10) << setprecision(3) << site_info[siteNumber].coords[Y]
               << setw(10) << setprecision(3) << site_info[siteNumber].coords[Z]
               << " Du " 
               << setw(width_1 +2) << site_info[siteNumber].rigid_label
               << setw(10) << "RCS" << site_info[siteNumber].rigid_label 
               << " 0.0" 
               << endl;
      else
        output << showpoint << setiosflags(ios::fixed)
               << setw(width_1 +2) << left << site_info[siteNumber].orig_atom_number << right
               << setw(3) << site_info[siteNumber].element_name
               << setw(10) << setprecision(3) << site_info[siteNumber].coords[X]
               << setw(10) << setprecision(3) << site_info[siteNumber].coords[Y]
               << setw(10) << setprecision(3) << site_info[siteNumber].coords[Z]
               << " Du " 
               << setw(width_1 +2) << site_info[siteNumber].rigid_label
               << setw(10) << "RCB" << site_info[siteNumber].rigid_label 
               << " 0.0" 
               << endl;
    }
    else{
      output << showpoint << setiosflags(ios::fixed)
             << setw(width_1 +2) << left << site_info[siteNumber].orig_atom_number << right
             << setw(3) << site_info[siteNumber].element_name
             << setw(10) << setprecision(3) << site_info[siteNumber].coords[X]
             << setw(10) << setprecision(3) << site_info[siteNumber].coords[Y]
             << setw(10) << setprecision(3) << site_info[siteNumber].coords[Z]
             << " Du " 
             << setw(width_1 +2) << site_info[siteNumber].rigid_label
             << setw(10) << "FLEX"
             << " 0.0" 
             << endl;
    }
  }

  // Output the bonds in the file.
  ////////////////////////////////////////////////////////////////////////////////
  unsigned int counter = 0;
  int bonds = 0;
  output << endl << "@<TRIPOS>BOND" << endl;
  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    for( unsigned int neighborNumber = 0; neighborNumber < (site_info[siteNumber].neighbor_list).size(); neighborNumber++ ){
      unsigned int neighborSiteNumber = site_info[siteNumber].neighbor_list[neighborNumber];
      if( neighborSiteNumber > siteNumber ){
        output << ++counter << " "
               << siteNumber << " "
               << neighborSiteNumber 
               << " 1" << endl;
        bonds++;
      }
    }
    bonds = 0;
  }
  
  if( counter != total_bonds ){
    stringstream warning;
    warning << "Warning: total_bonds found differs from pebble_game: " << total_bonds << " != " << counter;
    add_warning( 1, warning.str() );
  }
  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputPyMolScriptRCD_MOL2Format(){

  //float stick_radius = 0.25;// FIXME - unused variable

  string output_file = path_name + base_name + "_RCD_mol2.pml";
  ofstream pymol_script( output_file.c_str() );
  
  pymol_script << "# Rigid cluster decomposition coloring script for PyMol" << endl << "#" << endl;
  pymol_script << "# Brandon Hespenheide." << endl;
  pymol_script << "# Department of Physics and Astronomy" << endl;
  pymol_script << "# Biophysics Theory Group" << endl;
  pymol_script << "# Arizona State University" << endl;
  pymol_script << "############################################################" << endl << endl;

  string mol2_file = base_name;
  mol2_file += parameters.run_number + "_RCD";
  pymol_script << "load " << mol2_file << ".mol2" << endl << endl;
  pymol_script << "hide everything" << endl << endl;

  // Create the rigid cluster objects for pymol, and color them.
  ////////////////////////////////////////////////////////////////////////////////
  int cluster_label = 0;
  while( RC_atom_list[++cluster_label].size() >= parameters.min_output_cluster_size ){
    pymol_script << "create RC" << cluster_label << ", resn RCS" << cluster_label << " | resn RCB" << cluster_label << endl;
    pymol_script << "show lines, RC" << cluster_label << endl;
    pymol_script << "set line_width=2.5, " << "RC" << cluster_label << endl;
    pymol_script << "color " << next_rcd_color_as_name() << ", RC" << cluster_label << endl << endl;
  }  
  
  // Color everything that isn't a rigid cluster gray, with thin lines.
  ////////////////////////////////////////////////////////////////////////////////
  pymol_script << endl
               << "create flexible, resn FLEX | resn RCS*" << endl
               << "show lines, flexible" << endl
               << "set line_width=1, flexible" << endl
               << "color white, flexible" << endl << endl;

  // Set some visual defaults
  ////////////////////////////////////////////////////////////////////////////////
  //pymol_script << "bg_color white" << endl;
  pymol_script << "clip slab, 100" << endl;
  pymol_script << "center all" << endl;
  pymol_script << "zoom" << endl;

  pymol_script.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string outputPDBFormat( Site_Info &site_info ){

  stringstream new_line;

  new_line << showpoint << setiosflags(ios::fixed) 
           << setw(6) << site_info.record_name 
           << setw(5) << site_info.orig_atom_number
           << setw(5) << site_info.atom_name
           << setw(4) << site_info.residue_name
           << setw(2) << char(site_info.chain_ID)
           << setw(4) << site_info.seq_number
           << setw(1) << site_info.insert_code
           << "   "
           << setw(8) << setprecision(3) << site_info.coords[X]
           << setw(8) << setprecision(3) << site_info.coords[Y]
           << setw(8) << setprecision(3) << site_info.coords[Z]
           << setw(6) << setprecision(2) << site_info.occupancy
           << setw(6) << setprecision(2) << site_info.temp_factor
           << endl;

  return( new_line.str() );
};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
// 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputHARLEMFormat(){

  string linebuf;

  // Open the harlem file
  //////////////////////////////////////////////////////////////////////
  ifstream harlem_file( infile_name.c_str(), ios::in );
  if( !harlem_file ){
    cout << " ERROR: Could not open file." << endl;
    exit(1);
  }

  // Open the output file
  //////////////////////////////////////////////////////////////////////
  string output_name = path_name + base_name + "_F.hlm";
  ofstream output_file( output_name.c_str() );
  if( !output_file ){
    cout << " ERROR: Could not open file." << endl;
    exit(1);
  }

  int counter = 0;
  while( !harlem_file.eof() ){

    getline(harlem_file, linebuf);
    if( linebuf.find("#ATOM LIST") == string::npos )
      output_file << linebuf << endl;
    else{
      output_file << linebuf << endl;

      for( unsigned int a = 1; a <= total_clusters; a++ ){
        counter = 0;
        output_file << " RIGIDCLST" << a << "  " << RC_atom_list[a].size() << " ";
        for( unsigned int atom_num = 0; atom_num < RC_atom_list[a].size(); atom_num++ ){
          counter++;
          output_file << harlem_names[RC_atom_list[a].at(atom_num)] << " " ;
          if( counter == 30 ){
            counter = 0;
            output_file << "###" << endl << "         ";
          }
        }
        output_file << endl;
      }

    }
  }

  harlem_file.close();
  output_file.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Replace the white-space on atom names from PDB formatted files. This white
//   space was stripped out while reading to adhere to standard and heteroatom
//   atom reading functions. 
////////////////////////////////////////////////////////////////////////////////
string MolFramework::atomNamePDBFormat( const string &orig_name, const string &element ){

  string new_name = "    ";

  if( orig_name.size() == 2 ){
    if( element[1] != ' ' )
      new_name.replace( 0, 2, orig_name );
    else{
      if( element[0] == orig_name[0] ) // if the element is the first character of the orig_atom string
        new_name.replace( 1, 2, orig_name );
      else if( element[0] == orig_name[1] )
        new_name.replace( 0, 2, orig_name );
      else{
        new_name.replace( 1, 2, orig_name );
        add_warning( 3, "Atom may be improperly named in output file.");
      }
    }
  }

  else if( orig_name.size() == 3 ){
    if( element[1] != ' ' )
      new_name.replace( 0, 3, orig_name );
    else{
      if( element[0] == orig_name[0] ) // if the element is the first character of the orig_atom string
        new_name.replace( 1, 3, orig_name );
      else if( element[0] == orig_name[1] )
        new_name.replace( 0, 3, orig_name );
      else{
        new_name.replace( 1, 3, orig_name );
        add_warning( 3, "Atom may be improperly named in output file.");
      }
    }
  }

  else if( orig_name.size() == 1 ){
    new_name.replace( 1, 1, orig_name );
  }

  else if( orig_name.size() == 4 ){
    new_name.replace( 0, 4, orig_name );
  }

  else{
    cout << " Error: &orig_name variable empty in atomNamePDBFormat()" << endl;
    exit(-1);
  }

  //cout << "[" << orig_name << "] [" << new_name << "] [" << element << "]" <<endl;

  return( new_name );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputFlexibilityIndexPDBFormat(){

  //float label = 0.0;// FIXME - unused variable

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  Description:
//    Output all of the CONECT lines which were read in from the original
//    PDB file to the RCD.pdb file.  The CONECT lines are written to the end
//    of the RCD.pdb file, just before the END tag.  The original atom numbers
//    from the CONECT records are output, not the FIRST atom numbers.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::outputCONECTRecords( ostream &output ){

  map< SiteID, vector<SiteID> >::const_iterator record = conect_records.begin();

  while( record != conect_records.end() ){
  
    int atom1 = record->first;
    const vector<SiteID> *neighbors = &record->second;
    size_t nNeighbors = neighbors->size();
    output << "CONECT" << setw(5) << atom1;
    //cout << "CONECT" << setw(5) << atom1;

    for(size_t j = 0; j < nNeighbors; j++){
      output << setw(5) << (*neighbors)[j];
      //cout << setw(5) << (*neighbors)[j];
    }

    //cout << endl;
    output << endl;
    record++;
  }
}
/////////////////////////////////////////////////////////////////////////////////
