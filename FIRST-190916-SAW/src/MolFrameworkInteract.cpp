#include "global_defs.h"
#include "Parameters.h"
#include "MolFramework.h"

extern Parameters parameters;

int choice = 1;
string input;
string temp_string;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Allow the user to choose which alternative side-chain conformation to 
//   use. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::userInteraction2( vector<int> mult_side_chain_conf_list,
					vector<string> alt_loc_label_list, 
					Site_Info *site_info ){

  int label_number = 1;
  string label;

  sort(alt_loc_label_list.begin(), alt_loc_label_list.end() );
  vector<string>::iterator pos = unique(alt_loc_label_list.begin(), alt_loc_label_list.end() );
  vector<string>::iterator unique_label = alt_loc_label_list.begin(); 

  clear_screen;

  cout << " -- Option - Multiple Side-Chain Conformations --------------" << endl << endl;
  cout << " The following alternative location labels were found in your file." << endl;
  while( unique_label < pos ){
    cout << " " << label_number << ") " << *unique_label << endl;
    unique_label++;
    label_number++;
  }

  cout << endl << endl;
  do{
    cout << " 1. Choose the number of the alternative location label you would like to use." << endl;
    cout << "    Choice = ";
    getline(cin, input);
    
    choice = (int) strtol(input.c_str(), 0, 10);

    if( choice > 0 && choice < label_number )
      label = alt_loc_label_list[choice-1];
    else
      cout << " Your choice should be between 1 and " << label_number-1 << endl;
    
  } while( choice <= 0 || choice >= label_number );


  int current_atom = 0;
  for( unsigned int sideChainConfNumber = 0; sideChainConfNumber < mult_side_chain_conf_list.size(); sideChainConfNumber++ ){

    current_atom = mult_side_chain_conf_list[sideChainConfNumber];
    if( site_info[current_atom].alt_location != label )
      site_info[current_atom].excluded = true;
  }

}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
int MolFramework::userInteraction3( unsigned int atom_1, unsigned int atom_2, 
						Site_Info *site_info, float distance ){

  clear_screen;
  
  cout << " -- Option - Bond Not Found in Residue Template Library -----" << endl << endl;
 
  cout.setf( ios::left );

  cout << spacing
       << setw(10) << "atom"
       << setw(10) << "atom"
       << setw(10) << "residue"
       << setw(10) << "sequence"
       << setw(8)  << "chain" << endl;

  cout << spacing
       << setw(10) << "number"
       << setw(10) << "name"
       << setw(10) << "name"
       << setw(10) << "number"
       << setw(8)  << "ID" << endl << endl;

  cout << spacing 
       << setw(10) << site_info[atom_1].orig_atom_number
       << setw(10) << site_info[atom_1].atom_name 
       << setw(10) << site_info[atom_1].residue_name
       << setw(10) << site_info[atom_1].seq_number
       << setw(8)  << char(site_info[atom_1].chain_ID) << endl

       << spacing << "           |" << endl

       << spacing 
       << setw(10) << site_info[atom_2].orig_atom_number
       << setw(10) << site_info[atom_2].atom_name 
       << setw(10) << site_info[atom_2].residue_name
       << setw(10) << site_info[atom_2].seq_number
       << setw(8)  << char(site_info[atom_2].chain_ID)
       << endl << endl

       << spacing << "Distance = " << setprecision(4) << distance << " Angstroms" << endl << endl;
  
  do{
    cout << " 1. CONNECT the atoms (Press Enter to connect)." << endl;
    cout << " 2. DO NOT CONNECT the atoms." << endl << endl;
    cout << "    Choice = ";
    getline(cin, input);
    choice = (int) strtol(input.c_str(), 0, 10);

    if( choice == 1 ||
	did_press_enter(input) ){
      return(1);
    }
    else if( choice == 2 )
      return(0);
    else
      cout << " Please choose 1 or 2." << endl;
    

    cout << "choice = " << choice << endl;
  } while( choice <= 0 || choice > 2 );

  return(1);
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
int userInteraction4( vector<string> model_list ){

  bool is_valid = true;
  vector<int> valid_model_numbers;

  do{
    clear_screen;
    cout << " -- Option - Multiple NMR models present --------------------" << endl << endl;

    cout << " Please choose a model to analyze from the following list.   " << endl;
    cout << " ------------------------------------------------------------" << endl;
    for( unsigned int modelNumber = 0; modelNumber < model_list.size(); modelNumber++ ){
      temp_string.assign( model_list[modelNumber], 11, 14 );
      int model_num = atoi(temp_string.c_str() );
      valid_model_numbers.push_back( model_num );
      if( model_num <= 9 )
	cout << " " << model_num << ")  " << model_list[modelNumber] << endl;
      else
	cout << " " << model_num << ") " << model_list[modelNumber] << endl;
    } 

    cout << endl << "     Choice = ";
    getline(cin, input);
    choice = (int) strtol(input.c_str(), 0, 10);
    
    is_valid = true;
    for( unsigned int modelNumber = 0; modelNumber < valid_model_numbers.size(); modelNumber++){
      if( valid_model_numbers[modelNumber] == choice )
	is_valid = false;
    }
        
  } while( choice <= 0 || 
	   (unsigned int)choice > model_list.size() || // FIXME - conversion from signed to unsigned integer 
	   did_press_enter(input) ||
	   is_valid );

  temp_string.assign( model_list[choice-1], 11, 14 );
  return( atoi(temp_string.c_str() ) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::userInteractionEnergyCutoff(){

  string input;

  //clear_screen;
  
  do{
    cout << endl << endl;
    cout << " Please select an energy cutoff to use (in kcal/mol): ";
    getline( cin, input );

  } while( invalid_float(input) ||
	   did_press_enter(input) );

  parameters.energy_cutoff = atof( input.c_str() );
}
////////////////////////////////////////////////////////////////////////////////
