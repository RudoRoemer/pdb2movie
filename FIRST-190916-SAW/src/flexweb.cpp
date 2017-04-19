#include <unistd.h>

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>

#include "global_defs.h"
#include "output.h"
#include "generalUtils.h"

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Save the PID of an instance of FIRST in a file named 'PID'
////////////////////////////////////////////////////////////////////////////////
void output_PID( string path ) {

  string fileName = path;
  fileName += "PID";
  ofstream pidfile( fileName.c_str() );
  
  pidfile << getpid() << endl;
  
  pidfile.close();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Populate the vector taskNames with those "tasks" that should appear in the
//   runStatus tile when a job is submitted to FlexWeb. The task names should
//   be added to the vector in the order in which they are to appear in the 
//   runStatus tile. This vector of task names is actually a vector of pairs,
//   each pair is <string,int>, where string is the task name, and int is a number 1
//   or 2, that affect the size of font used when presenting the information 
//   in HTML. 1 is for sections, 2 is for subsections (which will have smaller 
//   font).
////////////////////////////////////////////////////////////////////////////////
void setTaskNames(){

  vector< pair<string,int> > taskNames;

  if( parameters.run_first ){
    parameters.taskNames.push_back( pair<string,int>("Build molecular framework",1) );
    parameters.taskNames.push_back( std::make_pair("Analyze flexibility",1) );
  }

  if( parameters.runFRODA ){
    parameters.taskNames.push_back( make_pair("Initialize FRODA dynamics",1) );
    parameters.taskNames.push_back( make_pair("Generate conformations",1) );
    parameters.taskNames.push_back( make_pair("Conformations generated",2) );
    parameters.taskNames.push_back( make_pair("Current All-Atom RMSD",2) );
    parameters.taskNames.push_back( make_pair("Current Main Chain RMSD",2) );
 }

  if( parameters.run_timme ){
    parameters.taskNames.push_back( make_pair("Build molecular framework",1) );
	parameters.taskNames.push_back( make_pair("Compute rigid cluster hierarchy",1) );
	parameters.taskNames.push_back( make_pair("Compute pair deviations",2) );
	parameters.taskNames.push_back( make_pair("Decompose into rigid clusters",2) );
	parameters.taskNames.push_back( make_pair("Assign rigid labels",2) );
    parameters.taskNames.push_back( make_pair("Generate dilution plot",2) );
    parameters.taskNames.push_back( make_pair("Calculate flexibility",2) );
//    parameters.taskNames.push_back( make_pair("Calculate mobility",1) );
    
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void outputStatus() {
  if (parameters.flexweb) {
    
    // TODO - more this to flexweb.cpp (pass an ordered map from task(s) to status)
    // TODO - DRY - this belongs in a common utility 
    int last_slash = parameters.infile_name.find_last_of("/\\");
    string path_name;
    path_name.assign( parameters.infile_name, 0, last_slash+1 );
    // end TODO 
    
    // internally store status in a stringstream 
    stringstream status;
    
    string task;
    
    status << "<table class='status'>" << endl;
    
    for( unsigned int taskNumber = 0; taskNumber < parameters.taskNames.size(); taskNumber++ ){
      
      task = parameters.taskNames[taskNumber].first;
      
      status << "\t<tr>" << endl;
      status << "\t\t<td class=\"statusTaskName_"<< parameters.taskNames[taskNumber].second << "\" >" 
	       << task
	       << "</td> <td>:</td>" << endl;
      status << "\t\t<td class=\"statusTaskValue_"<< parameters.taskNames[taskNumber].second << "\" >"
	       << parameters.mapFromTaskNameToStatus[task]
	       << "</td>" << endl; 
      status << "\t</tr>" << endl;
    }
    
    status << "</table>" << endl;
    
    // actually save the file
    string statusFilename = path_name + ".jobStatus.html";
    ofstream statusFile( statusFilename.c_str() );
    statusFile << status.str();
    statusFile.close();
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void output_details_flexweb() {
  
  string file_name = parameters.path_to_working_dir + ".runDetails.html";
  ofstream output_details( file_name.c_str() );
  
  output_details << "<h3>";
    
  if (parameters.runFRODA) {
    output_details << "FRODA";
  } else if (parameters.run_timme) {
    output_details << "TIMME";
  } else {
    output_details << "FIRST";
  }
  
  output_details << " run</h3><br/>" << endl;
  
  time_t current_time;
  time(&current_time);
  struct tm *timeptr = localtime(&current_time);  

  const size_t currentTimeLength = 128;
  char currentTime[currentTimeLength];
  strftime( currentTime, currentTimeLength,
            "Date: %A, %d %B %Y", timeptr );
  
  output_details << currentTime << "<br />" << endl;
  int last_slash = parameters.infile_name.find_last_of("/\\");
  string path_name;
  path_name.assign( parameters.infile_name, 0, last_slash+1 );
  string basename = parameters.infile_name.substr( last_slash+1 );

  output_details << "Input file name: " << basename << "<br />" << endl;
      
  if (parameters.run_first || parameters.runFRODA) {
    output_details << "Input Parameters:<br />" << endl;

    output_details << "Energy cutoff: " << parameters.energy_cutoff << endl;

  }   
  
  if (parameters.run_timme) {
    
  } 
  
  output_details.close();
  
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void output_Jmol_script(ofstream &jmol_script, 
                        MolFramework &structure ){
  
  jmol_script << "viewProperties.colorByRigidClusterJmolScript = \"\";";
  // Create the rigid cluster objects for Jmol, and color them. Only those
  // clusters larger than the min_output_cluster_size will have objects 
  // created for them. The remaining clusters will be binned into bulk
  // objects. 
  //////////////////////////////////////////////////////////////////////
  int total_RC_objects = 0;
  for( unsigned int clusterNumber = 0; clusterNumber < structure.total_clusters; clusterNumber++ ){
    
    if( structure.RC_atom_list[clusterNumber+1].size() >= parameters.min_output_cluster_size ){
      if( total_RC_objects == 0 )
        jmol_script << "viewProperties.colorByRigidClusterJmolScript += \"select temperature>=" << clusterNumber+1 << " and temperature<" << clusterNumber+2 
          << "; color " << next_jmol_color_as_name(true) << "; \"" << endl;
      else
        jmol_script << "viewProperties.colorByRigidClusterJmolScript += \"select temperature>=" << clusterNumber+1 << " and temperature<" << clusterNumber+2 
          << "; color " << next_jmol_color_as_name() << "; \"" << endl;
      total_RC_objects++;
    }
  } 
  
  jmol_script << "viewProperties.showWireframeRigidClusterJmolScript = \"\";" << endl;
  // Create the rigid cluster objects for Jmol, and color them. Only those
  // clusters larger than the min_output_cluster_size will have objects 
  // created for them. The remaining clusters will be binned into bulk
  // objects. 
  //////////////////////////////////////////////////////////////////////
  total_RC_objects = 0;
  for( unsigned int clusterNumber = 0; clusterNumber < structure.total_clusters; clusterNumber++ ){
    
    if( structure.RC_atom_list[clusterNumber+1].size() >= parameters.min_output_cluster_size ){
      if( total_RC_objects == 0 )
        jmol_script << "viewProperties.showWireframeRigidClusterJmolScript += \"select temperature>=" << clusterNumber+1 << " and temperature<" << clusterNumber+2 
          << "; wireframe 0.2; \"" << endl;
      else
        jmol_script << "viewProperties.showWireframeRigidClusterJmolScript += \"select temperature>=" << clusterNumber+1 << " and temperature<" << clusterNumber+2 
          << "; wireframe 0.2; \"" << endl;
      total_RC_objects++;
    }
  } 
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void output_Jmol_script( MolFramework &structure ){
  string file_name = structure.path_name + structure.base_name;
  file_name += "_RCD.js";
  ofstream jmol_script( file_name.c_str() );
  
  string replacementFileName = structure.path_name + ".RCD.js"; // TODO - replace file_name with this when we migrate to the cluster
  ofstream replacementJmolScript(replacementFileName.c_str());
//  replacementJmolScript << "var rcdPDBfilename = "+structure.base_name+"_RCD.pdb;";
  
  // Print some initialization information for Jmol.
  //////////////////////////////////////////////////////////////////////
  jmol_script << "jmolInitialize(\"/jmol\");" << endl;
  jmol_script << "jmolCheckBrowser(\"popup\", \"../../browsercheck\", \"onClick\"); " << endl;
  jmol_script << "jmolSetAppletColor(\"#FFFFFF\"); " << endl;
  jmol_script << "jmolSetAppletCssClass(\"jmol\"); " << endl;
  
  output_Jmol_script(jmol_script, structure);
  output_Jmol_script(replacementJmolScript, structure);
  
  jmol_script << "viewProperties.showBackbone = false;" << endl;
  jmol_script << "viewProperties.showOverlay = false;" << endl;
  jmol_script << "viewProperties.showNetwork = true;" << endl;
  jmol_script << "viewProperties.showAnimation = false;" << endl;
  jmol_script << "viewProperties.modelURI = \"" << structure.base_name << "_RCD.pdb\";" << endl;
  
  // TODO - removeme
  jmol_script << "jmolApplet( 600, \"load " << structure.base_name << "_RCD.pdb; select all; cpk off; color black; wireframe 0.1; \\" << endl;
      
  // Create the rigid cluster objects for Jmol, and color them. Only those
  // clusters larger than the min_output_cluster_size will have objects 
  // created for them. The remaining clusters will be binned into bulk
  // objects. 
  //////////////////////////////////////////////////////////////////////
  int total_RC_objects = 0;
  for( unsigned int clusterNumber = 0; clusterNumber < structure.total_clusters; clusterNumber++ ){
    
    if( structure.RC_atom_list[clusterNumber+1].size() >= parameters.min_output_cluster_size ){
      if( total_RC_objects == 0 ) {
        jmol_script << "\t\tselect temperature>=" << clusterNumber+1 << " and temperature<" << clusterNumber+2 
        << "; color BONDS " << next_jmol_color_as_name(true) << "; wireframe 0.2; \\" << endl;
      }
        
      else {
        jmol_script << "\t\tselect temperature>=" << clusterNumber+1 << " and temperature<" << clusterNumber+2 
        << "; color BONDS " << next_jmol_color_as_name() << "; wireframe 0.2; \\" << endl;

      }
      
      total_RC_objects++;
    }
  }
  
  jmol_script << "\");" << endl;
  // END TODO 
  
  jmol_script.close();
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void output_FRODA_TIMME_Jmol_script( MolFramework &structure ){ // FIXME - Depreciate
  output_Jmol_script(structure);
  /*
//  string file_name = structure.base_name;
  string file_name = ".RCD.js";
  ofstream jmol_script( file_name.c_str() );
  // Print some initialization information for Jmol.
  //////////////////////////////////////////////////////////////////////
  jmol_script << "jmolInitialize(\"/jmol\");" << endl;
  jmol_script << "jmolCheckBrowser(\"popup\", \"../../browsercheck\", \"onClick\"); " << endl;
  jmol_script << "jmolSetAppletColor(\"#FFFFFF\"); " << endl;
  jmol_script << "jmolSetAppletCssClass(\"jmol\"); " << endl;
  
  output_Jmol_script(jmol_script, structure);
  
  jmol_script << "viewProperties.showBackbone = true;" << endl;
  jmol_script << "viewProperties.showOverlay = true;" << endl;
  jmol_script << "viewProperties.showNetwork = false;" << endl;
  jmol_script << "viewProperties.showAnimation = false;" << endl;
  jmol_script << "viewProperties.modelURI = \".ensemble.pdb\";" << endl;
  
  jmol_script << "jmolApplet(600);" << endl;

  // TODO - request an update of the Jmol view
  
  jmol_script.close();*/
} 
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void output_RCD_XML( MolFramework &structure ){
  
  string file_name = structure.path_name + structure.base_name;
  file_name += "_RCD_flexweb.xml";
  ofstream rcd_xml( file_name.c_str() );
  
  rcd_xml << "<html>" << endl;
  rcd_xml << "  <head>" << endl;
  rcd_xml << "    < title>Flexweb > FIRST online > Results > " << structure.base_name << "</title>" << endl;
  rcd_xml << "  <head>" << endl;
  rcd_xml << "    <color_scheme color=\"yellow\"/>" << endl;
  rcd_xml << "  </head>" << endl << endl;
  
  rcd_xml << "  <menus>" << endl;
  rcd_xml << "    <menu_flexweb/>" << endl;
  rcd_xml << "    <menu_first_online/>" << endl;
  rcd_xml << "  </menus>" << endl << endl;
  
  rcd_xml << "  <body>" << endl;
  rcd_xml << "  <div class=\"section\">" << endl;
  rcd_xml << "  <div class=\"section_title\">Results for " << structure.base_name << "</div>" << endl;
  rcd_xml << "  <form> <!-- form tags should always be used around UI controls -->" << endl;
  rcd_xml << "    <script src=\"test_jmol.js\">&#160;</script>" << endl;
  rcd_xml << "  </form>" << endl;
  rcd_xml << "  </div>" << endl << endl;
  
  rcd_xml << "  </body>" << endl;
  rcd_xml << "</html>" << endl;
  
  rcd_xml.close();
}
////////////////////////////////////////////////////////////////////////////////
