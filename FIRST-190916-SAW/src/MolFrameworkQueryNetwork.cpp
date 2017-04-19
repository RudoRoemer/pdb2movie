#include "global_defs.h"
#include "Parameters.h"
#include "MolFramework.h"

extern const Parameters parameters;
vector<unsigned int> empty_vector;

////////////////////////////////////////////////////////////////////////////////
string bond_type_label( int bond_type ){
  if( bond_type == 1 ) return("User defined 1 bar bond.");
  if( bond_type == 2 ) return("User defined 2 bar bond.");
  if( bond_type == 3 ) return("User defined 3 bar bond.");
  if( bond_type == 4 ) return("User defined 4 bar bond.");
  if( bond_type == 5 ) return("User defined 5 bar bond.");
  if( bond_type == 6 ) return("User defined 6 bar bond.");
  if( bond_type == 7 ) return("Covalent Bond.");
  if( bond_type == 8 ) return("Covalent Bond.");
  if( bond_type == 9 ) return("Hydrogen Bond.");
  if( bond_type == 10) return("Salt Bridge.");
  if( bond_type == 11) return("Hydrophobic Tether.");
  if( bond_type == 12) return("Stacked Ring.");

  return ""; // FIXME - is "" a suitable default?
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   A macro to create a string that describe a user-defined bond between two 
//   atoms. The label consists of a dash-enclosed letter "U" followed by the 
//   number of bars between the two atoms. For example, -U5- is a user-defined 
//   5-bar bond.
// Parameters:
//   bars - The number of bars to be printed in the label. 
// Return Value:
//   string - Label for the bond type to be printed in the query_network menu
//     system.
////////////////////////////////////////////////////////////////////////////////
inline string make_bar_label( int bars ){
  string temp;
  stringstream s;
  s << bars;
  s >> temp;
  return( "-U"+temp+"-" );

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
void MolFramework::queryNetwork(){

  if( !parameters.interactive &&  
      !parameters.run_query ){

    //outputCovalentBondList();
    //outputHydrogenBondList();
    //outputHydrophobicTetherList();
    //outputStackedRingsList();
    return;
  }

  int choice = 0;
  string input;

  do{
    clear_screen;
    cout << " -- Menu A --------------------------------------------------" << endl;
    cout << " 1. Explore the current bond network." << endl << endl;
    cout << " <- Press the Enter key to continue running FIRST " << FIRST_VERSION_ID << endl << endl;
    cout << "    Choice = ";
    getline(cin, input);
    choice = (int) strtol(input.c_str(), 0, 10);

    if( choice == 1 )
      queryOptionsScreen();

  } while( did_not_press_enter(input) );

  clear_screen;
  cout << " Returning to the program." << endl;
  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
void MolFramework::queryOptionsScreen(){

  int choice = 0;
  string input;

  do{
    clear_screen;
    cout << " -- Menu B --------------------------------------------------" << endl << endl;
    cout << " 1. List all bonds to a specific atom number." << endl;
    cout << " 2. Add a constraint between two atoms." << endl << endl;
    cout << "  * The following output routines will create a list of the indicated bonds in a" << endl;
    cout << "  * format that is compatible with FIRST " << FIRST_VERSION_ID << ". If you require a more detailed" << endl;
    cout << "  * listing, please use the command-line option \"-o2\"." << endl << endl;
    cout << " 3. Output a list of all covalent bonds in [Atom_1 Atom_2] format." << endl;
    cout << " 4. Output a list of all hydrogen bonds in [Hydrogen Donor Energy Bars] format." << endl;
    cout << " 5. Output a list of all hydrophobic tethers in [Atom_1 Atom_2] format." << endl;
    cout << " 6. Output a list of all stacked aromatic rings in [Atom_1 Atom_2 Bars] format." << endl << endl;
    cout << " 7. List all bonds to a specific residue number." << endl << endl;

    if( choice == 3 )
      cout << "  * List of covalent bonds successfully written to file." << endl << endl;
    else if( choice == 4 )
      cout << "  * List of hydrogen bonds successfully written to file." << endl << endl;
    else if( choice == 5 )
      cout << "  * List of hydrophic tethers successfully written to file." << endl << endl;
    else if( choice == 6 )
      cout << "  * List of stacked aromatic rings successfully written to file." <<endl<<endl;
    cout << " <- Press enter to return to the previous menu." << endl << endl;
    cout << "    Choice = ";
    getline(cin, input);
    choice = (int) strtol(input.c_str(), 0, 10);

    if( choice == 1 )
      querySiteNumber( empty_vector );

    else if( choice == 2 )
      addConstraint();

    else if( choice == 3 )
      outputCovalentBondList();

    else if( choice == 4 )
      outputHydrogenBondList();

    else if( choice == 5 )
      outputHydrophobicTetherList();

    else if( choice == 6 )
      outputStackedRingsList();

    else if( choice == 7 )
      querySequenceNumber();

    else if( did_press_enter(input) )
      return;
    
  } while( did_not_press_enter(input) );
  
  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::printInfoLine1( int atom_number, int &counter, int position, string bond_type, string demark ){

  if( position == 0 ){
    cout.setf( ios::right );
    cout << setw(2) << counter << ") ";
    cout.unsetf( ios::right );
    cout.setf( ios::left );
    cout << setw(8) << site_info[atom_number].orig_atom_number 
	 << setw(8) << site_info[atom_number].atom_name 
	 << setw(8) << site_info[atom_number].residue_name
	 << setw(8) << site_info[atom_number].seq_number 
	 << setw(4) << char(site_info[atom_number].chain_ID)
	 << "   ";
  }
  else if( position == 1 ){
    cout.setf( ios::left );
    cout << setw(6) << bond_type 
	 << setw(8) << site_info[atom_number].orig_atom_number
	 << setw(8) << site_info[atom_number].atom_name 
	 << setw(8) << site_info[atom_number].residue_name
	 << setw(8) << site_info[atom_number].seq_number 
	 << setw(4) << char(site_info[atom_number].chain_ID)
	 << demark
	 << endl;
  }
  else{
    cout.setf( ios::right );
    cout << setw(2) << counter << ") ";
    cout.unsetf( ios::right );
    cout.setf( ios::left );
    cout << "                                       "
	 << setw(6) << bond_type 
	 << setw(8) << site_info[atom_number].orig_atom_number
	 << setw(8) << site_info[atom_number].atom_name 
	 << setw(8) << site_info[atom_number].residue_name
	 << setw(8) << site_info[atom_number].seq_number 
	 << setw(4) << char(site_info[atom_number].chain_ID)
	 << demark
	 << endl;
  }

  if( position != 1 )
    counter++;
 
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This function will dislpay information about a specific atom in the molecule.
//   The user enters the atom number as it appears in the input file. The number
//   is mapped to the internally used FIRST_number, and the atom information, 
//   including all the bonds to the given atom, is displayed.
// Parameters:
//   orig_atom_number - Atom number as found in the input file.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::querySiteNumber( vector<unsigned int> atom_numbers, bool interresidue_bonds_only ){

  bool invalid = false;
  unsigned int atom_number = 0;
  unsigned int FIRST_number = 0;
  int counter = 1;
  int position = 0;
  unsigned int bond = 0;
  int choice = 0;
  string input;
  string demark = "";
  vector<new_bonds> bond_list;

  if( atom_numbers.size() == 0 ){
    clear_screen; 
    do{
      cout << " -- Menu C ----------------------------------------------------------------------" << endl << endl;
      cout << " Please enter the atom number ";
      getline(cin, input);
      atom_number = (int) strtol(input.c_str(), 0, 10);
      FIRST_number = getFIRSTNumber( atom_number );
      
      if( did_press_enter(input) )
	return;

      if( FIRST_number <= 0 || FIRST_number > total_sites ){
	cout << endl << "*That atom number was not found. Please re-enter the atom number." << endl;
      }

      
    } while( FIRST_number <= 0 || FIRST_number > total_sites );

    atom_numbers.push_back( FIRST_number );
  }

  do{
    counter = 1;
    bond_list.clear();
    clear_screen;
    cout << " -- Menu C ----------------------------------------------------------------------" << endl << endl;
        
    cout.setf( ios::left );
    cout << "    "
	 << setw(8) << "atom"
	 << setw(8) << "atom"
	 << setw(8) << "res."
	 << setw(8) << "seq."
	 << setw(6) << "chain" 
	 << setw(7) << "|bond|"
	 << setw(8) << "atom"
	 << setw(8) << "atom"
	 << setw(8) << "res."
	 << setw(8) << "seq."
	 << setw(6) << "chain" 
	 << endl;
    
    cout << "    "
	 << setw(8) << "num"
	 << setw(8) << "name"
	 << setw(8) << "name"
	 << setw(8) << "num"
	 << setw(6) << "ID" 
	 << setw(7) << "|type|"
	 << setw(8) << "num"
	 << setw(8) << "name"
	 << setw(8) << "name"
	 << setw(8) << "num"
	 << setw(6) << "ID" 
	 << endl;  
    
    cout << "-------------------------------------------" 
	 << "-------------------------------------------" << endl;
    
    for( unsigned int atomNumber = 0; atomNumber < atom_numbers.size(); atomNumber++ ){
      atom_number = atom_numbers[atomNumber];
      position = 0;
      printInfoLine1( atom_number, counter, position++ );

      // 1. Check the covalent bond topology
      ////////////////////////////////////////////////////////////////////////////////
      neighbor_atom = (site_info[atom_number].neighbor_list).begin();
      while( neighbor_atom != (site_info[atom_number].neighbor_list).end() ){
	if( !isDifferentResidue( *neighbor_atom, atom_number) ) {
	  if( !interresidue_bonds_only )
	    printInfoLine1( *neighbor_atom, counter, position++, "-CB-" );
	}
	else
	  printInfoLine1( *neighbor_atom, counter, position++, "-CB-", " <---" );
	bond_list.push_back( new_bonds(atom_number, *neighbor_atom, 7) );
	neighbor_atom++;
      } 
      
      // 2. Check the hydrogen bond lists
      ////////////////////////////////////////////////////////////////////////////////
      for( unsigned int hydrogenBondNumber = 0; hydrogenBondNumber < hydrogen_bonds.size(); hydrogenBondNumber++ ){
	
	if( hydrogen_bonds[hydrogenBondNumber].site_1 == atom_number ){
	  if( !isDifferentResidue( hydrogen_bonds[hydrogenBondNumber].site_2, atom_number) &&
	      !interresidue_bonds_only )
	    printInfoLine1( hydrogen_bonds[hydrogenBondNumber].site_2, counter, position++, "-HB-" );
	  else
	    printInfoLine1( hydrogen_bonds[hydrogenBondNumber].site_2, counter, position++, "-HB-", " <---" );
	  bond_list.push_back( new_bonds(atom_number, hydrogen_bonds[hydrogenBondNumber].site_2, 9, hydrogen_bonds[hydrogenBondNumber].energy) );
	}
	else if( hydrogen_bonds[hydrogenBondNumber].site_2 == atom_number ){
	  if( !isDifferentResidue( hydrogen_bonds[hydrogenBondNumber].site_1, atom_number)  &&
	      !interresidue_bonds_only )
	    printInfoLine1( hydrogen_bonds[hydrogenBondNumber].site_1, counter, position++, "-HB-" );
	  else
	    printInfoLine1( hydrogen_bonds[hydrogenBondNumber].site_1, counter, position++, "-HB-", " <---" );
	  bond_list.push_back( new_bonds(atom_number, hydrogen_bonds[hydrogenBondNumber].site_1, 9, hydrogen_bonds[hydrogenBondNumber].energy) );
	}
      }
      
      // 3. Check the hydrophobic bond lists.
      ////////////////////////////////////////////////////////////////////////////////
      for( unsigned int hydrophobicTetherNumber = 0; hydrophobicTetherNumber < hydrophobic_tethers.size(); hydrophobicTetherNumber++ ){
	
	if( hydrophobic_tethers[hydrophobicTetherNumber].site_1 == atom_number ){
	  if( isDifferentResidue( hydrophobic_tethers[hydrophobicTetherNumber].site_2, atom_number ) )
	    printInfoLine1( hydrophobic_tethers[hydrophobicTetherNumber].site_2, counter, position++, "-PH-", " <---" );
	  else
	    printInfoLine1( hydrophobic_tethers[hydrophobicTetherNumber].site_2, counter, position++, "-PH-" );
	  bond_list.push_back( new_bonds(atom_number, hydrophobic_tethers[hydrophobicTetherNumber].site_2, 11) );
	}
	else if( hydrophobic_tethers[hydrophobicTetherNumber].site_2 == atom_number ){
	  if( isDifferentResidue( hydrophobic_tethers[hydrophobicTetherNumber].site_1, atom_number ) )
	    printInfoLine1( hydrophobic_tethers[hydrophobicTetherNumber].site_1, counter, position++, "-PH-", " <---" );
	  else
	    printInfoLine1( hydrophobic_tethers[hydrophobicTetherNumber].site_1, counter, position++, "-PH-" );
	  bond_list.push_back( new_bonds(atom_number, hydrophobic_tethers[hydrophobicTetherNumber].site_1, 11) );
	}
      } 
      // 4. Check the stacked ring lists.
      ////////////////////////////////////////////////////////////////////////////////
      for( unsigned int stackedRingNumber = 0; stackedRingNumber < stacked_rings.size(); stackedRingNumber++ ){
	
	if( stacked_rings[stackedRingNumber].site_1 == atom_number ){
	  if( isDifferentResidue( stacked_rings[stackedRingNumber].site_2, atom_number ) )
	    printInfoLine1( stacked_rings[stackedRingNumber].site_2, counter, position++, "-RS-", " <---" );
	  else
	    printInfoLine1( stacked_rings[stackedRingNumber].site_2, counter, position++, "-RS-" );
	  bond_list.push_back( new_bonds(atom_number, stacked_rings[stackedRingNumber].site_2, 12, 
					 stacked_rings[stackedRingNumber].energy) );
	}
	else if( stacked_rings[stackedRingNumber].site_2 == atom_number ){
	  if( isDifferentResidue( stacked_rings[stackedRingNumber].site_1, atom_number ) )
	    printInfoLine1( stacked_rings[stackedRingNumber].site_1, counter, position++, "-RS-", " <---" );
	  else
	    printInfoLine1( stacked_rings[stackedRingNumber].site_1, counter, position++, "-RS-" );
	  bond_list.push_back( new_bonds(atom_number, stacked_rings[stackedRingNumber].site_1, 12,
					stacked_rings[stackedRingNumber].energy) );
	}
      } 
      // 5. Check the user defined bond list.
      ////////////////////////////////////////////////////////////////////////////////
      for( unsigned int userDefinedConstraintNumber = 0; userDefinedConstraintNumber < user_defined_constraints.size(); userDefinedConstraintNumber++ ){
	
	if( user_defined_constraints[userDefinedConstraintNumber].site_1 == atom_number ){
	  string bars = make_bar_label( user_defined_constraints[userDefinedConstraintNumber].bars );
	  if( isDifferentResidue( user_defined_constraints[userDefinedConstraintNumber].site_2, atom_number ) )
	    printInfoLine1( user_defined_constraints[userDefinedConstraintNumber].site_2, counter, position++, bars, " <---" );
	  else
	    printInfoLine1( user_defined_constraints[userDefinedConstraintNumber].site_2, counter, position++, bars );
	  bond_list.push_back( new_bonds(atom_number, 
					 user_defined_constraints[userDefinedConstraintNumber].site_2,
					 user_defined_constraints[userDefinedConstraintNumber].bars) );
	}
	else if( user_defined_constraints[userDefinedConstraintNumber].site_2 == atom_number ){
	  string bars = make_bar_label( user_defined_constraints[userDefinedConstraintNumber].bars );
	  if( isDifferentResidue( user_defined_constraints[userDefinedConstraintNumber].site_1, atom_number ) )
	    printInfoLine1( user_defined_constraints[userDefinedConstraintNumber].site_1, counter, position++, bars, " <---" );
	  else
	    printInfoLine1( user_defined_constraints[userDefinedConstraintNumber].site_1, counter, position++, bars );
	  bond_list.push_back( new_bonds(atom_number, 
					 user_defined_constraints[userDefinedConstraintNumber].site_1, 
					 user_defined_constraints[userDefinedConstraintNumber].bars) );
	}
      }
      // If no neighbors.
      //////////////////////////////////////////////////////////////////////
      if( bond_list.size() == 0 )
	cout << "      NO NEIGHBORS" << endl;
      
      cout << "-------------------------------------------" 
	   << "-------------------------------------------" << endl;
    }
    
    bond = 0;
    cout << endl;
    cout << " 1. View details of a given bond." << endl;
    
    if( atom_numbers.size() == 1 ){
      cout << " 2. Add a constraint to atom: " << site_info[atom_numbers[0] ].orig_atom_number << endl;
      
      if( bond_list.size() )
	cout << " 3. Remove a constraint from atom: " << site_info[atom_numbers[0] ].orig_atom_number << endl;
    }
      
    //cout << " 4. Only show interresidue bonds." << endl << endl;
    
    cout << " <- Press enter to return to the previous menu." << endl << endl;
    if( invalid ){
      cout << "   *[" << input << "] is an invalid option. Please reenter." << endl;
      invalid = false;
    }
    cout << "    Choice = ";
    getline(cin, input);
    choice = (int) strtol(input.c_str(), 0, 10);
    
    if( choice == 1 ){
      while( !bond ){
	cout << endl << "   Enter the number of the bond (left-hand column): ";
	getline(cin, input);
	bond = (int) strtol(input.c_str(), 0, 10);
	if( bond <= 0 || bond > bond_list.size() ){
	  cout << "     Invalid bond number, please reenter." << endl;
	  bond = 0;
	}
      }
      viewBondDetails(bond_list[bond-1].site_1, 
			bond_list[bond-1].site_2, 
			bond_list[bond-1].bars,
			bond_list[bond-1].energy);
    }
    else if( choice == 2 && atom_numbers.size() ==1 ){
      addConstraint( atom_numbers.at(0) );
    }
    else if( choice == 3 && atom_numbers.size() ==1 ){
      while( !bond ){
	cout << endl << "    Enter the number of the bond (left-hand column): ";
	getline(cin, input);
	bond = (int) strtol(input.c_str(), 0, 10);
	if( did_press_enter(input) ){
	  input = "continue";
	  break;
	}
	if( bond <= 0 || bond > bond_list.size() ){
	  cout << "   *Invalid bond number, please reenter." << endl;
	  bond = 0;
	}
      }
      removeConstraint( bond_list[bond-1].site_1, bond_list[bond-1].site_2, bond_list[bond-1].bars );
    }
    //else if( choice == 4 ){
    //querySiteNumber( atom_numbers, true );
    //return;
    //}
    else if( did_press_enter(input) )
      return;
    else{
      invalid = true;
    }
  } while( did_not_press_enter(input) );
  
  cout << endl << " --------------------------------------------------------------------------------" << endl << endl;
  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
// Parameters:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::querySequenceNumber( int seq_number ){

  // bool invalid = false; // FIXME - unused variable
  //int atom_number = 0;  // FIXME - unused variable
  unsigned int chain_ID = 0;
  //int counter = 0;  // FIXME - unused variable
  //int spacer  = 0; // FIXME - unused variable
  //int bond = 0;  // FIXME - unused variable
  //int choice = 0; // FIXME - unused variable
  string input;
  pair<int,string> residue_and_chain;
  vector<new_bonds> bond_list;

  if( !seq_number ){
    clear_screen; 
    do{
      cout << " -- Menu E ----------------------------------------------------------------------" << endl << endl;
      cout << " Please enter a reside number: ";
      getline(cin, input);
      seq_number = (int) strtol(input.c_str(), 0, 10);

      if( did_press_enter(input) )
	return;

      if( seq_number <= 0 ){
	cout << endl << "*That residue number was not found. Please re-enter the atom number." << endl;
      }
      
    } while( seq_number <= 0 );
  }
  
  // Search the list of atoms for the given residue. Check for multiple chains
  ////////////////////////////////////////////////////////////////////////////////
  vector<unsigned int> atoms_in_this_residue;
  vector<unsigned int> chain_IDs;
  for( unsigned int a = 1; a <= total_sites; a++ ){
    if( site_info[a].seq_number == seq_number ){
      atoms_in_this_residue.push_back( a );
      chain_IDs.push_back( site_info[a].chain_ID );
    }
  }

  vector<unsigned int>::iterator pos = unique( chain_IDs.begin(), chain_IDs.end() );
  chain_IDs.erase( pos, chain_IDs.end() );

  // If multiple chains are found, prompt the user for a chain ID.
  //////////////////////////////////////////////////////////////////////
  if( chain_IDs.size() != 1 ){
    bool found = false;
    clear_screen;
    do{
      cout << " -- Menu E ----------------------------------------------------------------------" << endl << endl;
      cout << " Please select from the following chains" << endl;
      cout << " ------------------------------------------------------------" << endl;
      for( unsigned int a = 0; a < chain_IDs.size(); a++ )
	cout << " " << char(chain_IDs[a]) << endl;

      cout << endl << "     Choice = ";
      getline(cin, input);
      chain_ID = int(input[0]);

      if( did_press_enter(input) )
	return;

      if( find(chain_IDs.begin(), chain_IDs.end(), chain_ID) == chain_IDs.end() ){
	cout << endl << "*That chain ID was not found. Please select a chain ID from the list above." << endl;
      }
      else
	found = true;
      
    } while( !found );

    vector<unsigned int>::iterator not_in_chain = atoms_in_this_residue.begin();
    while( not_in_chain != atoms_in_this_residue.end() ){
      if( site_info[*not_in_chain].chain_ID != chain_ID )
	atoms_in_this_residue.erase( not_in_chain );
      else
	not_in_chain++;
    }
  }

  querySiteNumber( atoms_in_this_residue );

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::addConstraint( unsigned int orig_atom_1 ){
  
  unsigned int atom_1 = 0;
  int orig_atom_2 = 0;
  unsigned int atom_2 = 0;
  int bars   = 0;
  int option = 0;
  string input;

  // 1. Choose how to model the new constraint (determines how many bars to place 
  //    between the two atoms. 
  ////////////////////////////////////////////////////////////////////////////////
  do{ 
    clear_screen;
    cout << " -- Menu D -----------------------------------------------------------------------" << endl;
    cout << " Please choose how to model this constraint:" << endl;
    cout << " 1. 1 Bar"  << endl;
    cout << " 2. 2 Bars" << endl;
    cout << " 3. 3 Bars" << endl;
    cout << " 4. 4 Bars" << endl;
    cout << " 5. 5 Bars" << endl;
    cout << " 6. 6 Bars" << endl;
    cout << " 7. Rotatable Covalent Bond" << endl;
    cout << " 8. Non-rotatable Covalent Bond (like a peptide bond)" << endl;
    cout << " 9. Hydrogen Bond (Must enter Hydrogen and Acceptor atom numbers)" << endl;
    cout << "10. Salt Bridge (Must enter Hydrogen and Acceptor atom numbers)" << endl;
    cout << "11. Hydrophobic Tether (currently using " << parameters.hp_bars << " bars to model tethers)" << endl << endl;
    cout << " <- Press enter to return to the previous menu." << endl << endl;
    cout << "    Model constraint using option = ";
    getline(cin, input);
    option = (int) strtol(input.c_str(), 0, 10);
    
    if( did_press_enter(input) ){
      clear_screen;
      return;
    }

  } while( option <= 0 || option > 10 );

  // 2. Enter the atom number of the first atom, unless a non-zero atom number was
  //    passed as a formal argument. 
  ////////////////////////////////////////////////////////////////////////////////
  if( !orig_atom_1 ){
    cout << " Please enter the atom numbers (as in the original PDB file) for the atoms in this bond: " << endl;
    do{
      cout << "    Make connection to atom number = ";
      getline(cin, input);
      orig_atom_1 = (int) strtol(input.c_str(), 0, 10);
      atom_1 = getFIRSTNumber( orig_atom_1 );
      
      if( atom_1 <= 0 || atom_1 > total_sites )
	cout << " Error: The atom number must be between 1 and " << total_sites << endl << endl;
      
    } while( atom_1 <= 0 || atom_1 > total_sites );
  }
  else
    atom_1 = orig_atom_1;

  // 3. Enter the atom number of the second atom.
  //////////////////////////////////////////////////////////////////////
  do{
    cout << "    Make connection to atom number = ";
    getline(cin, input);
    orig_atom_2 = (int) strtol(input.c_str(), 0, 10);
    atom_2 = getFIRSTNumber( orig_atom_2 );

    if( atom_2 <= 0 || atom_2 > total_sites )
      cout << " Error: That atom number was not found." << endl << endl;

  } while( atom_2 <= 0 || atom_2 > total_sites );

  cout << endl;

  // Set atom_1 to have the smaller atom number. 
  //////////////////////////////////////////////////////////////////////
  if( atom_1 > atom_2 )
    swap( atom_1, atom_2 );


  // 4. Check to see if these atoms are already bonded. 
  //////////////////////////////////////////////////////////////////////
  if( alreadyBonded(atom_1, atom_2) ){
    do{
      cout << " ** Warning: " << orig_atom_1 << " and " << orig_atom_2 << " are alread bonded." << endl;
      cout << " ** Press 'y' to replace the bond, or press enter to cancel: ";
      getline(cin, input);
      
      if( did_press_enter(input) )
	return;
      
    } while( input != "y" && input != "Y" );

    remove_from_site_info_array( atom_1, atom_2 );
  }

  if( option == 7 || option == 8 ){ // Covalent bond
    bars = option -2;
    add_to_site_info_array( atom_1, atom_2, bars );
  }
  else if( option == 9 ){ // Hydrogen bond
    bars = 5;
    float HB_energy = hydrogenBondEnergy(atom_1, atom_2);
    hydrogen_bonds.push_back( new_bonds(atom_1, atom_2, bars, HB_energy) );
  }
  else if( option == 10 ){ // Salt bridge
    bars = 5;
    float SB_energy = saltBridgeEnergy(atom_1, atom_2);
    hydrogen_bonds.push_back( new_bonds(atom_1, atom_2, bars, SB_energy) );
  }
  else if( option == 11 ){ // Hydrophobic bond
    bars = parameters.hp_bars;
    hydrophobic_tethers.push_back( new_bonds(atom_1, atom_2, bars, PH_DEFAULT_ENERGY) );
  } 
  else{ // User defined constraint
    bars = option;
    user_defined_constraints.push_back( new_bonds(atom_1, atom_2, bars) );
  }

  return;    
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::removeConstraint( int atom_1, int atom_2, int bond_type ){

  if( bond_type >= 1 &&
      bond_type <= 6 ){
    if( removeFromBondList(atom_1, atom_2, user_defined_constraints) ){
      cout << "Warning: Bond not found in list." << endl;
      exit(1);
    }
  }      

  if( bond_type == 7 ||
      bond_type == 8 )
    remove_from_site_info_array( atom_1, atom_2 );

  if( bond_type == 9 ||
      bond_type == 10 ){
    if( removeFromBondList(atom_1, atom_2, hydrogen_bonds) ){
      cout << "Warning: Bond not found in list." << endl;
      exit(1);
    }    
  }

  if( bond_type == 11 ){
    if( removeFromBondList(atom_1, atom_2, hydrophobic_tethers) ){
      cout << "Warning: Bond not found in list." << endl;
      exit(1);
    }    
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::viewBondDetails( int atom_1, int atom_2, int bond_type, float energy ){

  string input;

  clear_screen;
  cout << " ---------------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << setw(8) << site_info[atom_1].orig_atom_number
       << setw(9) << site_info[atom_1].atom_name 
       << setw(7) << site_info[atom_1].residue_name
       << setw(5) << site_info[atom_1].seq_number
       << setw(9) << site_info[atom_1].chain_ID << endl;
  cout << "        |" << endl;
  cout << "        | Bond Type: " << bond_type_label(bond_type) << endl;
  cout << "        | Bond Energy = "<< energy << " kcal/mol" << endl;
  cout << "        | Distance = " << getDistance( atom_1, atom_2 ) << " Angstroms " << endl;
  cout << "        |" << endl;
  cout << setw(8) << site_info[atom_2].orig_atom_number
       << setw(9) << site_info[atom_2].atom_name 
       << setw(7) << site_info[atom_2].residue_name
       << setw(5) << site_info[atom_2].seq_number
       << setw(9) << site_info[atom_2].chain_ID;
  cout << endl << endl;
  do{
    cout << " <- Press the Enter key to return to the previous menu: ";
    getline(cin, input);
  } while( did_not_press_enter(input) );
  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
int MolFramework::alreadyBonded( unsigned int atom_1, unsigned int atom_2 ){

  neighbor_atom = (site_info[atom_1].neighbor_list).begin();
  while( neighbor_atom != (site_info[atom_1].neighbor_list).end() ){
    if( *neighbor_atom == atom_2 )
      return(1);
    neighbor_atom++;
  }

  return(0);
}
////////////////////////////////////////////////////////////////////////////////

