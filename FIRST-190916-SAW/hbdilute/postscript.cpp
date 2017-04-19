/************************************************************/
/* Various bits and pieces of the output postscript file    */
/* need to be included. The lines were hard coded in to eli-*/
/* minate the need for the buch of files in a specific      */
/* directory that need to be catted to the final file in a  */
/* certain order.                                           */
/************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "types.h"

#define TRUE  1
#define FALSE 0
/************************************************************/
/*            The postscript header info                    */
/* Includes standard stuff to make the file "legal" level 2 */
/* postscript, as well as definitions for the colors used   */
/* and the geometric primitives.                            */
/************************************************************/
void print_header(FILE *ps_file, float scale ){

  time_t current_time;

  current_time = time(NULL);
  
  fprintf(ps_file, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(ps_file, "%%%%Title: HBD.ps\n");
  fprintf(ps_file, "%%%%Creator: Brandon Hespenheide\n");
  fprintf(ps_file, "%%%%CreationDate: 1999\n");
  fprintf(ps_file, "%%%%BoundingBox:0 0 620 800\n");
  fprintf(ps_file, "%%%%DocumentFonts: Times-Roman\n");
  fprintf(ps_file, "%%%%EndComments\n");
  fprintf(ps_file, "%%%%BeginProlog\n");
  fprintf(ps_file, "/Col { sethsbcolor } bind def\n");
  fprintf(ps_file, "%% These are the colors for the flexibility scale and the\n");
  fprintf(ps_file, "%% lines that display the hydrogen bonds.\n\n");

  fprintf(ps_file, "/Blue {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/Red  {1.0000  0.0000  0.0000 setrgbcolor } def\n");

  fprintf(ps_file,"/Col00 {1.0000 0.0000 0.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col01 {0.0000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col02 {0.3000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col03 {0.6000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col04 {0.9000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col05 {0.2000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col06 {0.5000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col07 {0.8000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col08 {0.1000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col09 {0.4000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col10 {0.7000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col11 {0.0000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col12 {0.3000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col13 {0.6000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col14 {0.9000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col15 {0.2000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col16 {0.5000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col17 {0.8000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col18 {0.1000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col19 {0.4000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col20 {0.7000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col21 {0.0250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col22 {0.3250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col23 {0.6250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col24 {0.9250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col25 {0.2250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col26 {0.5250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col27 {0.8250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col28 {0.1250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col29 {0.4250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col30 {0.7250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col31 {0.0250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col32 {0.3250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col33 {0.6250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col34 {0.9250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col35 {0.2250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col36 {0.5250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col37 {0.8250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col38 {0.1250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col39 {0.4250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col40 {0.7250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col41 {0.0500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col42 {0.3500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col43 {0.6500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col44 {0.9500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col45 {0.2500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col46 {0.5500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col47 {0.8500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col48 {0.1500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col49 {0.4500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col50 {0.7500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col51 {0.0500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col52 {0.3500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col53 {0.6500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col54 {0.9500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col55 {0.2500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col56 {0.5500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col57 {0.8500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col58 {0.1500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col59 {0.4500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col60 {0.7500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col61 {0.0750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col62 {0.3750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col63 {0.6750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col64 {0.9750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col65 {0.2750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col66 {0.5750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col67 {0.8750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col68 {0.1750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col69 {0.4750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col70 {0.7750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col71 {0.0750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col72 {0.3750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col73 {0.6750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col74 {0.9750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col75 {0.2750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col76 {0.5750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col77 {0.8750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col78 {0.1750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col79 {0.4750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col80 {0.7750 0.8000 0.8000 sethsbcolor } def\n");
  
  fprintf(ps_file, "/FlexCol01 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol02 {0.1000  0.1000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol03 {0.2000  0.2000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol04 {0.3000  0.3000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol05 {0.4000  0.4000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol06 {0.4500  0.4500  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol07 {0.5000  0.5000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol08 {0.6000  0.6000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol09 {0.7000  0.7000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol10 {0.8000  0.8000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol11 {1.0000  0.8000  0.8000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol12 {1.0000  0.7000  0.7000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol13 {1.0000  0.6000  0.6000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol14 {1.0000  0.5000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol15 {1.0000  0.4500  0.4500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol16 {1.0000  0.4000  0.4000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol17 {1.0000  0.3000  0.3000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol18 {1.0000  0.2000  0.2000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol19 {1.0000  0.1000  0.1000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol20 {1.0000  0.0000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol21 {0.5000  0.5000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol30 {1.0000  0.0000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol31 {0.7500  0.7500  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol32 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol33 {1.0000  0.5000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol34 {0.0000  1.0000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol35 {0.9000  0.9000  0.9000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol36 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol37 {0.5000  1.0000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol38 {0.5000  0.5000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol39 {0.5000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol40 {1.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol41 {1.0000  0.5000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol42 {0.2500  1.0000  0.2500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol43 {0.0000  0.7500  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol44 {0.8000  1.0000  0.2500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol45 {0.8000  1.0000  0.2500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol46 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol47 {0.5000  1.0000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol48 {0.5000  0.5000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol49 {0.5000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol50 {1.0000  0.0000  1.0000 setrgbcolor } def\n");

  
  fprintf(ps_file, "/Poly4 { moveto lineto lineto lineto fill } bind def\n");
  fprintf(ps_file, "/Pl4 { 8 copy Poly4 moveto moveto moveto moveto closepath stroke } bind def\n");
  fprintf(ps_file, "/Poly3 { moveto lineto lineto fill } bind def\n");
  fprintf(ps_file, "/Pl3 { 6 copy Poly3 moveto moveto moveto closepath stroke } bind def\n");
  fprintf(ps_file, "/Print { /Times-Roman findfont exch scalefont setfont show } bind def  \n");
  fprintf(ps_file, "/sequence_line { moveto lineto stroke } bind def\n\n");
  fprintf(ps_file, "/CenterRot90 {\n");
  fprintf(ps_file, "  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file, "  exch stringwidth pop -2 div exch 3 div exch rmov\n");
  fprintf(ps_file, " } bind def\n");
  fprintf(ps_file, "/UncenterRot90 {\n");
  fprintf(ps_file, "  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file, "  exch stringwidth } bind def\n");
  fprintf(ps_file, "/Rot90 { gsave currentpoint translate 90 rotate } bind def\n");
  fprintf(ps_file, "%%%%EndProlog\n\n");
  fprintf(ps_file, "%%%%Page:    1   1\n");

  fprintf( ps_file, "420 770 moveto\n");
  fprintf( ps_file, "(%s) 10.0 Print\n", ctime(&current_time) ); 

  fprintf(ps_file, "%3.2f %3.2f scale\n", scale, scale );
  fprintf(ps_file, "save\n\n");
}
/************************************************************/


/************************************************************/
/* print a line of information at the bottom of the last    */
/* page of postscript output. for PORTRAIT style output.    */
/************************************************************/
void print_portrait_footer( FILE *ps_file, char protein_name[128], int ps_line_number ){

  ps_line_number -= 10;
  fprintf(ps_file, "Blue\n");
  fprintf(ps_file, "40 %d moveto\n", ps_line_number );
  fprintf(ps_file, "(Blue:donor) 10 Print\n");
  fprintf(ps_file, "Red\n");
  fprintf(ps_file, "90 %d moveto\n", ps_line_number );
  fprintf(ps_file, "(Red:acceptor) 10 Print\n");
  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "180 %d moveto\n", ps_line_number);
  fprintf(ps_file, "(M:main-chain   S:side-chain   W:water   H:hetero-atom) 10 Print\n");
  fprintf(ps_file, "500 %d moveto\n", ps_line_number);
  fprintf(ps_file, "(%s) 12 Print\n", protein_name );
}
/************************************************************/


/************************************************************/
/* print a line of information at the bottom of the last    */
/* page of postscript output. for LANDSCAPE style output.   */
/************************************************************/
void print_landscape_footer( FILE *ps_file, char protein_name[128], int ps_line_number ){

  ps_line_number += 10;
  fprintf(ps_file, "Col07\n");
  fprintf(ps_file, "%d 110 moveto\n", ps_line_number);
  fprintf(ps_file, "(Blue:donor) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(Blue:donor) 10 Print\n");
  fprintf(ps_file, "grestore\n");
  fprintf(ps_file, "Col01\n");
  fprintf(ps_file, "%d 170 moveto\n", ps_line_number);
  fprintf(ps_file, "(Red:acceptor) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(Red:acceptor) 10 Print\n");
  fprintf(ps_file, "grestore\n");
  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "%d 260 moveto\n", ps_line_number);
  fprintf(ps_file, "(M:main-chain   S:side-chain   W:water   H:hetero-atom) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(M:main-chain   S:side-chain   W:water   H:hetero-atom) 10 Print\n");
  fprintf(ps_file, "grestore\n");
  fprintf(ps_file, "%d 675 moveto\n", ps_line_number);
  fprintf(ps_file, "(%s) 12 UncenterRot90 Rot90\n", protein_name );
  fprintf(ps_file, "(%s) 12 Print\n", protein_name );
  fprintf(ps_file, "grestore\n");

}
/************************************************************/


/************************************************************/
/************************************************************/
void print_data_headings( FILE *ps_file, int long_output, int number_of_residues,
			  int number_of_chains, int chain_size[10] ){
  
  /* Print portrait */
  if( !number_of_chains && number_of_residues <= 125 ){ 
	
    fprintf(ps_file, "20 740 moveto\n");
    fprintf(ps_file, "(H-bond) 10.0 Print\n");
    fprintf(ps_file, "20 730 moveto\n");
    fprintf(ps_file, "(number) 10.0 Print\n");
    fprintf(ps_file, "58 740 moveto\n");
    fprintf(ps_file, "(E) 10.0 Print\n");
    fprintf(ps_file, "78 740 moveto\n");
    fprintf(ps_file, "(<r>) 10.0 Print\n");
  }    

  /* Print landscape */
  else{ 
    fprintf(ps_file, "33 20 moveto\n");
    fprintf(ps_file, "(H-bond) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(H-bond) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
    fprintf(ps_file, "43 20 moveto\n");
    fprintf(ps_file, "(number) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(number) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
    fprintf(ps_file, "43 57 moveto\n");
    fprintf(ps_file, "(E) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(E) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
    fprintf(ps_file, "43 79 moveto\n");
    fprintf(ps_file, "(<r>) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(<r>) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
  }
}
/************************************************************/


/************************************************************/
/************************************************************/
void print_current_Hbond_info_portrait( FILE *ps_file, int long_output, int short_output, 
					int first_line, int *ps_line_number, int hb_number, 
					float hb_energy, int cluster_counter, int colors[30], 
					int HB_start, int HB_end, int NUMBER_OF_RESIDUES, 
					float scale_factor, int *ps_page_number, int missing_one, 
					int HBtwo_start, int HBtwo_end, int y_translate[10],
					char donor_atom_type[4], char accept_atom_type[4], 
					int number_of_chains, char donor_chain_ID, 
					char accept_chain_ID, char chain_IDs[10],
					char protein_name[128], int side_chain_test, float mean_coordination,
					int fit_it_all_on_one_page ){

  int   
    a = 0,
    this_chain = 0;

  //float // FIXME - unused variables
    //triangle_x = 0, 
    //triangle_y = 0;

  //float // FIXME - unused variables
    //y_start = 0.0,
    //y_end = 0.0;

  fprintf(ps_file, "Col00\n");

  if( !side_chain_test ){
    if( !first_line ){
      fprintf( ps_file, "21 %3d moveto\n", *ps_line_number-4 );
      fprintf( ps_file, "(%4d) 9.0 Print\n", hb_number );
      fprintf( ps_file, "44 %3d moveto\n", *ps_line_number-4 );
      fprintf( ps_file, "(%7.3f) 9.0 Print\n\n", hb_energy );
      fprintf( ps_file, "72 %3d moveto\n", *ps_line_number-4 );
      fprintf( ps_file, "(%7.3f) 9.0 Print\n\n", mean_coordination );
    }
    else{
      fprintf( ps_file, "21 %3d moveto\n", *ps_line_number-3 );
      fprintf( ps_file, "(%4d) 9.0 Print\n\n", hb_number);
      fprintf( ps_file, "72 %3d moveto\n", *ps_line_number-3 );
      fprintf( ps_file, "(%7.3f) 9.0 Print\n\n", mean_coordination);
    }
  }
  
  for( a = 0; a <= number_of_chains; a++ ){
    
    fprintf(ps_file,"103 %3d moveto\n", *ps_line_number-((10*a)+3) );
    fprintf(ps_file,"(%c) 9.0 Print\n", chain_IDs[a] );
    
    /************************************************************/
    /* Print the red and blue triangles that denote where the   */
    /* broken Hbond occured, to the postscript file.            */
    /************************************************************/      
    /*if( !first_line ) {
      if( donor_chain_ID == chain_IDs[a] ) {
	this_chain = a;
	triangle_x = ( HB_start * 3 ) + 111.5;
	triangle_y = *ps_line_number-(3+(10*a));
	
	fprintf(ps_file,"Blue\n");
	fprintf(ps_file,"%3d %3d %3d %3d %3d %3d Pl3\n", triangle_x, triangle_y, (triangle_x - 2),
		(triangle_y - 3), (triangle_x + 2), (triangle_y - 3) );
      }
      if( accept_chain_ID == chain_IDs[a] ){
	
	triangle_x = ( HB_end * 3 ) + 111.5;
	triangle_y = *ps_line_number-(3+(10*a));
	
	fprintf(ps_file,"Red\n");
	fprintf(ps_file,"%3d %3d %3d %3d %3d %3d Pl3\n", triangle_x, triangle_y, (triangle_x - 2),
		(triangle_y - 3), (triangle_x + 2), (triangle_y - 3) );
      }
      }*/
    /************************************************************/      
  }
  
  HB_end   += (y_translate[this_chain] - 1);
  HB_start += (y_translate[this_chain] - 1);
  
  if( !first_line && ( short_output || missing_one || fit_it_all_on_one_page ) ){
    fprintf(ps_file,"Blue\n");
    
    if( donor_chain_ID == 'W' ) {      
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
      fprintf(ps_file,"(W) 9.0 Print\n");
    }
    else{
      if( donor_chain_ID == 'H' ) {      
	fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
	fprintf(ps_file,"(H) 9.0 Print\n");
      }
      else{
	if( !strncmp(donor_atom_type, "N", 2 ) ||
	    !strncmp(donor_atom_type, "O", 2 ) ) {
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
	  fprintf(ps_file,"(M) 9.0 Print\n");
	}
	else{
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
	  fprintf(ps_file,"(S) 9.0 Print\n");
	}
      }
    }
    
    fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+120, *ps_line_number-4 );
    fprintf(ps_file,"(%3d) 9.0 Print\n", HB_start);
    
    if( donor_chain_ID != 'W' && donor_chain_ID != 'H' ){
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+131, *ps_line_number-4 );
      fprintf(ps_file,"(%c) 8.0 Print\n", donor_chain_ID );
    }
    
    fprintf(ps_file,"Red\n");
    
    if( accept_chain_ID == 'W' ) {      
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );
      fprintf(ps_file,"(W) 9.0 Print\n");
    }
    else{
      if( accept_chain_ID == 'H' ) {      
	fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );	
	fprintf(ps_file,"(H) 9.0 Print\n");
      }
      else{
	if( !strcmp(accept_atom_type, "N" ) ||
	    !strcmp(accept_atom_type, "O" ) ) {
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );
	  fprintf(ps_file,"(M) 9.0 Print\n");
	}
	else{
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );
	  fprintf(ps_file,"(S) 9.0 Print\n");
	}
      }
    }
    
    fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+148, *ps_line_number-4 );
    fprintf(ps_file,"(%3d) 9.0 Print\n", HB_end );
    
    if( accept_chain_ID != 'W' && accept_chain_ID != 'H' ){
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+162, *ps_line_number-4 );
      fprintf(ps_file,"(%c) 8.0 Print\n", accept_chain_ID );
    }
    
  }
  
  if( long_output )
    *ps_line_number -= 60;
  else
    *ps_line_number -= 10;

  if( short_output || missing_one ) {
    if(*ps_line_number < 60 || ( *ps_line_number - ( number_of_chains*10 ) ) < 60){
      *ps_line_number = 710;
      (*ps_page_number)++;
      fprintf(ps_file, "\nshowpage\n" );
      fprintf(ps_file,"%%%%Page:    %d   %d\n\n", *ps_page_number, 
	      *ps_page_number );
    }
  }
  
}

/*************************************************************************/
/* If the protein has more than 125 amino acids, the output is displayed */
/* in landscape format. 250 amino-acids will fit on a page, multiple land*/
/* scape pages are used if neccesary.                                    */
/*************************************************************************/
void print_current_Hbond_info_landscape( FILE *ps_file, int long_output, int short_output, 
					 int first_line, int *ps_line_number, int hb_number, 
					 float hb_energy, int cluster_counter, int colors[30], 
					 int HB_start, int HB_end, int NUMBER_OF_RESIDUES, 
					 float scale_factor, int *ps_page_number, int missing_one, 
					 int HBtwo_start, int HBtwo_end, int y_translate[10],
					 char donor_atom_type[4], char accept_atom_type[4],
					 int number_of_chains, char donor_chain_ID, 
					 char accept_chain_ID, char chain_IDs[10], int chain_size[10],
					 char protein_name[128], int side_chain_test,
					 float mean_coordination, int fit_it_all_on_one_page ){

  int   
    a = 0, 
    b = 0,
	// triangle_x = 0, // FIXME - warning: unused variable 'triangle_x'
    summed_length = 0, 
    end_pt = 0, 
    // start_pt = 0, // FIXME - warning: unused variable 'start_pt'
    beggining_of_current_chain = 0,
    this_chain = 0;

  //float 
    // y_start=0.0, // FIXME - warning: unused variable 'y_start'
    //y_end=0.0, // FIXME - warning: unused variable 'y_start'
    // triangle_y = 0.0; // FIXME - warning: unused variable 'triangle_y'

  fprintf(ps_file, "Col00\n");

  if( first_line ){
    fprintf(ps_file, "%3d 36 moveto\n", *ps_line_number+3 );
    fprintf(ps_file, "(All %3d Hbonds) 9.0 UncenterRot90 Rot90\n", hb_number);
    fprintf(ps_file, "(All %3d Hbonds) 9.0 Print\n", hb_number);
    fprintf(ps_file, "grestore\n" );
  }
  else{
    fprintf(ps_file, "%3d 21 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%4d) 9.0 UncenterRot90 Rot90\n", hb_number );
    fprintf(ps_file, "(%4d) 9.0 Print\n", hb_number );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 44 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", hb_energy );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", hb_energy );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 72 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", mean_coordination );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", mean_coordination );
    fprintf(ps_file, "grestore\n" );
    
    for( a = 0; a <= number_of_chains; a++ ){

      summed_length += ( (chain_size[a]*3) + 18 ); /* end point for the current chain on the output file */

      for( b = 0; b < a; b++ )
	beggining_of_current_chain += (chain_size[b]*3) + 18;

      beggining_of_current_chain += 110;
      
      fprintf(ps_file, "Col00\n");
      fprintf(ps_file, "%3d %3d moveto\n", *ps_line_number+3, beggining_of_current_chain - 8 );
      fprintf(ps_file, "(%c) 9.0 UncenterRot90 Rot90\n", chain_IDs[a] );
      fprintf(ps_file, "(%c) 9.0 Print\n", chain_IDs[a] );
      fprintf(ps_file, "grestore\n" );
      
      /************************************************************/
      /* Print the red and blue triangles that denote where the   */
      /* broken Hbond occured, to the postscript file.            */
      /************************************************************/
      /*
      if( donor_chain_ID == chain_IDs[a] ) {
	
	triangle_x = *ps_line_number + 3;
	triangle_y = ( (HB_start-1) * 3 ) + 1.5 + beggining_of_current_chain;
	
	fprintf(ps_file, "Blue\n");
	fprintf(ps_file, "%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x, triangle_y, (triangle_x + 3),
		(triangle_y - 2), (triangle_x + 3), (triangle_y + 2) );
      }
      
      if( accept_chain_ID == chain_IDs[a] ){
	
	if( HB_start == HB_end){
	  triangle_x = *ps_line_number - 3;
	  triangle_y = ( (HB_end-1) * 3 ) + 1.5 + beggining_of_current_chain;
	  
	  fprintf(ps_file,"Red\n");
	  fprintf(ps_file,"%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x, triangle_y, (triangle_x - 3),
		  (triangle_y - 2), (triangle_x - 3), (triangle_y + 2) );
	}
	else{
	  triangle_x = *ps_line_number + 3;
	  triangle_y = ( (HB_end-1) * 3 ) + 1.5 + beggining_of_current_chain;
	  
	  fprintf(ps_file,"Red\n");
	  fprintf(ps_file,"%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x, triangle_y, (triangle_x + 3),
		  (triangle_y - 2), (triangle_x + 3), (triangle_y + 2) );
	}
	this_chain = a;
      }
      */
      /************************************************************/
      beggining_of_current_chain = 0;
    }

    HB_start += (y_translate[this_chain] - 1);
    HB_end   += (y_translate[this_chain] - 1);
    
    end_pt = summed_length-15;
    
    if( !first_line && ( short_output || fit_it_all_on_one_page ) ){
      fprintf(ps_file,"Blue\n");
      
      if( donor_chain_ID == 'W' ) {      
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+120 );
	fprintf(ps_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(ps_file,"(W) 9.0 Print\n");
	fprintf(ps_file,"grestore\n" );
      }
      else{
	if( donor_chain_ID == 'H' ) {      
	  fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+120 );
	  fprintf(ps_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(ps_file,"(H) 9.0 Print\n");
	  fprintf(ps_file,"grestore\n" );
	}
	else{
	  if( !strcmp(donor_atom_type, "N" ) ||
	      !strcmp(donor_atom_type, "O" ) ) {
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+120 );
	    fprintf(ps_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(M) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
	  else{
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+120 );
	    fprintf(ps_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(S) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
	}
      }
      
      fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+129 );
      fprintf(ps_file,"(%3d) 9.0 UncenterRot90 Rot90\n", HB_start );
      fprintf(ps_file,"(%3d) 9.0 Print\n", HB_start);
      fprintf(ps_file,"grestore\n" );
      if( donor_chain_ID != 'W' ){
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+144 );
	fprintf(ps_file,"(%c) 8.0 UncenterRot90 Rot90\n", donor_chain_ID );
	fprintf(ps_file,"(%c) 8.0 Print\n", donor_chain_ID );
	fprintf(ps_file,"grestore\n" );
      }
      
      fprintf(ps_file,"Red\n");
      
      if( accept_chain_ID == 'W' ) {      
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+155 );
	fprintf(ps_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(ps_file,"(W) 9.0 Print\n");
	fprintf(ps_file,"grestore\n" );
      }
      else{
	if( accept_chain_ID == 'H' ) {      
	  fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+155 );
	  fprintf(ps_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(ps_file,"(H) 9.0 Print\n");
	  fprintf(ps_file,"grestore\n" );
	}
	else{
	  if( !strcmp(accept_atom_type, "N" ) ||
	      !strcmp(accept_atom_type, "O" ) ) {
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+155 );
	    fprintf(ps_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(M) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
	  else{
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+155 );
	    fprintf(ps_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(S) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
	}
      }
      
      fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+164 );
      fprintf(ps_file,"(%3d) 9.0 UncenterRot90 Rot90\n", HB_end );
      fprintf(ps_file,"(%3d) 9.0 Print\n", HB_end );
      fprintf(ps_file,"grestore\n" );
      if( accept_chain_ID != 'W' ){
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+177 );
	fprintf(ps_file,"(%c) 8.0 UncenterRot90 Rot90\n", accept_chain_ID );
	fprintf(ps_file,"(%c) 8.0 Print\n", accept_chain_ID );
	fprintf(ps_file,"grestore\n" );
      }
      
      if( side_chain_test ){
	fprintf(ps_file, "%3d 90 moveto\n", *ps_line_number+4 );
	fprintf(ps_file, "(%4d) 9.0 UncenterRot90 Rot90\n", hb_number );
	fprintf(ps_file, "(%4d) 9.0 Print\n", hb_number );
	fprintf(ps_file, "grestore\n" );
	fprintf(ps_file, "%3d 60 moveto\n", *ps_line_number+4 );
	fprintf(ps_file, "((%4d)) 9.0 UncenterRot90 Rot90\n", HB_start );
	fprintf(ps_file, "((%4d)) 9.0 Print\n", HB_start );
	fprintf(ps_file, "grestore\n" );
	fprintf(ps_file, "%3d 30 moveto\n", *ps_line_number+4 );
	fprintf(ps_file, "(%4s) 9.0 UncenterRot90 Rot90\n", donor_atom_type );
	fprintf(ps_file, "(%4s) 9.0 Print\n", donor_atom_type );
	fprintf(ps_file, "grestore\n" );
      }
      
    }
  }
  
  if( long_output )
    *ps_line_number += 60;
  if( side_chain_test )
    *ps_line_number += 8  + ( ( 10 * number_of_chains ) );
  else
    *ps_line_number += 12;
  
      
  /* comment out the following code if you want all the output to be on a single page */
  /*
    if( short_output || missing_one || side_chain_test ) {
    if(*ps_line_number > 560 || ( ( *ps_line_number + ( number_of_chains * 10 ) ) > 560 ) ) {
    *ps_line_number = 60;
    (*ps_page_number)++;
    fprintf(ps_file, "\nshowpage\n" );
    fprintf(ps_file,"%%%%Page:    %d   %d\n\n", *ps_page_number, 
    *ps_page_number );
    }
    }
  */
      
}
/**********************************************************************/

/**********************************************************************/
/* prints the info in landscape format, but on the right page. Allows */
/* for multipage output for proteins with more than 200 amino-acids.  */
/**********************************************************************/
void print_current_Hbond_info_landscape_multipage( FILE *ps_file, int long_output, int short_output, 
						   int first_line, int *ps_line_number, int hb_number, 
						   float hb_energy, int cluster_counter, int colors[30], 
						   int donor_res_num, int accept_res_num,  int NUMBER_OF_RESIDUES, 
						   float scale_factor, int *ps_page_number, int missing_one, 
						   int HBtwo_start, int HBtwo_end, int y_translate[10],
						   char donor_atom_type[4], char accept_atom_type[4],
						   FILE *file_list[5], int number_of_pages, 
						   int number_of_chains, char donor_chain_ID, 
						   char accept_chain_ID, char chain_IDs[10],
						   char protein_name[128], int chain_size[10],
						   float mean_coordination ) {

  int 
    a = 0, 
    b = 0, 
    last_page = 0, 
    end_point = 0, 
    output_page = 0,
    triangle_x = 0, 
    //end_pt=0, // FIXME - unused variable
    start_pt = 0, 
    current_pg = 0, 
    summed_length = 0,
    previous_length = 0, 
    donor_pointer = 0, 
    accept_pointer = 0;
  
  float 
    //y_start = 0.0,  // FIXME - unused variable
    //y_end = 0.0,  // FIXME - unused variable
    y_position = 0.0, 
    triangle_y = 0.0;

  char  
    linebuf[300], 
    file_name[6];
  
  FILE 
    *current_file;

  /**********************************************************************/
  /* Print first 3 columns of output (hb number, energy, mean coord).   */
  /* or just the number of HBonds if this is the first line of results. */
  fprintf(ps_file, "Col00\n");
  if( !first_line ){
    fprintf(ps_file, "%3d 21 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%4d) 9.0 UncenterRot90 Rot90\n", hb_number );
    fprintf(ps_file, "(%4d) 9.0 Print\n", hb_number );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 44 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", hb_energy );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", hb_energy );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 72 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", mean_coordination );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", mean_coordination );
    fprintf(ps_file, "grestore\n" );
  }
  else{
    fprintf(ps_file, "%3d 36 moveto\n", *ps_line_number+3 );
    fprintf(ps_file, "(All %3d Hbonds) 9.0 UncenterRot90 Rot90\n", hb_number);
    fprintf(ps_file, "(All %3d Hbonds) 9.0 Print\n", hb_number);
    fprintf(ps_file, "grestore\n" );
  }
  /**********************************************************************/
  
  /**********************************************************************/
  /* Print the chain ID before each chain.                              */
  for( b = 0; b <= number_of_chains; b++ ){
    
    summed_length += ( chain_size[b] + ( b*6 ) );
    start_pt = ( summed_length -chain_size[b] );
    current_pg = start_pt/200;
    current_file = file_list[current_pg];

    printf("summed %d start_pt %d chain_size %d chain %d\n", summed_length, start_pt, chain_size[b], b );
    while( start_pt > 200 )
      start_pt -= 200;

    fprintf(current_file, "Col00\n");
    fprintf(current_file, "%3d %3d moveto\n", *ps_line_number+3, (start_pt*3) +110 -9 );
    fprintf(current_file, "(%c) 9.0 UncenterRot90 Rot90\n", chain_IDs[b] );
    fprintf(current_file, "(%c) 9.0 Print\n", chain_IDs[b] );
    fprintf(current_file, "grestore\n" );
    summed_length -= ( b*6 );
      
  }
  /**********************************************************************/

  summed_length = 0;
  /**********************************************************************/
  /* FIRST outputs the donor and acceptor residue numbers as starting   */
  /* from 1, regardless of the res number the chains actually start from*/
  /* Modify these numbers to coincide with the real residue number.     */
  for( b = 0; b <= number_of_chains; b++ ){
    if( donor_chain_ID == chain_IDs[b] )
      donor_res_num  += y_translate[b] -1;
    if( accept_chain_ID == chain_IDs[b] )
      accept_res_num += y_translate[b] -1;
  }
  /**********************************************************************/

  /**********************************************************************/
  /* Print the red and blue triangles under each line showing the loc-  */
  /* ation of the hydrogen bond along the primary structure.            */
  if( !first_line ) {
    
    for( b = 0; b <= number_of_chains; b++ ){

      summed_length += (chain_size[b] + ( 6*b ));

      /**********************************************************************/
      /* The second condition will stop hetatm's from having triangles.     */
      /**********************************************************************/
      if( donor_chain_ID == chain_IDs[b] &&
	  donor_res_num <= chain_size[b] ) {
	
	output_page = int (( summed_length -chain_size[b] +donor_res_num -y_translate[b]) / 200.0); 
	current_file = file_list[output_page];
	
	previous_length = summed_length - chain_size[b];
	
	while( previous_length > 200 ) {
	  previous_length -= 200;
	}

	donor_pointer = donor_res_num -y_translate[b];
	while( donor_pointer >= 200 ) {
	  donor_pointer -= 200;
	}
	
	if( ( 200 - previous_length ) > ( chain_size[b] ) )
	  y_position = ( (donor_pointer +previous_length)*3) + 111.5;
	else {
	  if( donor_pointer >= ( 200-previous_length) )
	    y_position = ( (donor_pointer -(200 -previous_length))*3 ) + 111.5;
	  else
	    y_position = ( (donor_pointer +previous_length)*3) + 111.5;
	}
	
	triangle_x = *ps_line_number + 3;
	triangle_y = y_position;
	
	fprintf(current_file, "Blue\n");
	fprintf(current_file, "%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x, triangle_y, (triangle_x + 3),
		(triangle_y - 2), (triangle_x + 3), (triangle_y + 2) );
	
      }
      if( accept_chain_ID == chain_IDs[b] &&
	  accept_res_num <= chain_size[b] ) {      
	
	output_page = int( ((summed_length -chain_size[b] +accept_res_num -y_translate[b])*1.0) / 200.0 ); 
	current_file = file_list[output_page];

	printf("res num %3d output_page %d, summed %d chain_size %d y_trans %d\n", accept_res_num, 
	       output_page, summed_length, chain_size[b], y_translate[b] );
	previous_length = summed_length - chain_size[b];
	
	while( previous_length > 200 ) {
	  previous_length -= 200;
	}
	
	accept_pointer = accept_res_num -y_translate[b];
	while( accept_pointer >= 200 ) {
	  accept_pointer -= 200;
	}
	
	if( ( 200 - previous_length ) > ( chain_size[b] ) )
	  y_position = ( (accept_pointer +previous_length) *3) + 111.5;
	else {
	  if( accept_pointer >= ( 200-previous_length) )
	    y_position = ( (accept_pointer -(200-previous_length)) *3) + 111.5;
	  else
	    y_position = ( (accept_pointer +previous_length) *3) + 111.5;
	}
	triangle_x = *ps_line_number + 3;
	triangle_y = y_position;
	
	fprintf(current_file, "Red\n");
	fprintf(current_file, "%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x, triangle_y, (triangle_x + 3),
		(triangle_y - 2), (triangle_x + 3), (triangle_y + 2) );
	
      }
    }
  }
  /**********************************************************************/

  last_page = summed_length / 200;
  current_file = file_list[last_page];

  if( number_of_chains ){
    
    while( summed_length > 200 ) {
      summed_length -= 200;
    }
    end_point = summed_length;

  }
  else
    end_point = NUMBER_OF_RESIDUES - (last_page * 200);

  if( !first_line && short_output ){
    
    fprintf(current_file, "Blue\n");

      if( donor_chain_ID == 'W' ) {      
	fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+120 );
	fprintf(current_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(current_file,"(W) 9.0 Print\n");
	fprintf(current_file,"grestore\n" );
      }
      else{
	if( donor_chain_ID == 'H' ) {      
	  fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+120 );
	  fprintf(current_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(current_file,"(H) 9.0 Print\n");
	  fprintf(current_file,"grestore\n" );
	}
	else{
	  if( !strcmp(donor_atom_type, "N" ) ||
	      !strcmp(donor_atom_type, "O" ) ) {
	    fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+120 );
	    fprintf(current_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(M) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
	  else{
	    fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+120 );
	    fprintf(current_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(S) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
	}
      }
      
      fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+129 );
      fprintf(current_file,"(%3d) 9.0 UncenterRot90 Rot90\n", donor_res_num );
      fprintf(current_file,"(%3d) 9.0 Print\n", donor_res_num );
      fprintf(current_file,"grestore\n" );
      if( donor_chain_ID != 'W' ){
	fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+144 );
	fprintf(current_file,"(%c) 8.0 UncenterRot90 Rot90\n", donor_chain_ID );
	fprintf(current_file,"(%c) 8.0 Print\n", donor_chain_ID );
	fprintf(current_file,"grestore\n" );
      }
    
    fprintf(current_file, "Red\n");

      if( accept_chain_ID == 'W' ) {      
	fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+155 );
	fprintf(current_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(current_file,"(W) 9.0 Print\n");
	fprintf(current_file,"grestore\n" );
      }
      else{
	if( accept_chain_ID == 'H' ) {      
	  fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+155 );
	  fprintf(current_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(current_file,"(H) 9.0 Print\n");
	  fprintf(current_file,"grestore\n" );
	}
	else{
	  if( !strcmp(accept_atom_type, "N" ) ||
	      !strcmp(accept_atom_type, "O" ) ) {
	    fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+155 );
	    fprintf(current_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(M) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
	  else{
	    fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+155 );
	    fprintf(current_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(S) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
	}
      }
      
      fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+164 );
      fprintf(current_file,"(%3d) 9.0 UncenterRot90 Rot90\n", accept_res_num );
      fprintf(current_file,"(%3d) 9.0 Print\n", accept_res_num );
      fprintf(current_file,"grestore\n" );
      if( accept_chain_ID != 'W' ){
	fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+177 );
	fprintf(current_file,"(%c) 8.0 UncenterRot90 Rot90\n", accept_chain_ID );
	fprintf(current_file,"(%c) 8.0 Print\n", accept_chain_ID );
	fprintf(current_file,"grestore\n" );
      }
    
    fprintf(current_file, "Col00\n");
  }
  /**********************************************************************/

  /**********************************************************************/
  /* Incremenet the line number. Also, if we just printed the last line */
  /* that will fit on this page, reset all variables to print at the top*/
  /* of the page, and print appropriate page markers to all postscript  */
  /* files.                                                             */
  *ps_line_number += 12;  
  
  if(*ps_line_number > 560) {
    *ps_line_number = 60;
    (*ps_page_number)++;
    
    for( a = 1; a <= number_of_pages; a++ ){
      current_file = file_list[a];
      fprintf(current_file, "\nshowpage\n" );
      fflush( current_file );
      fclose( current_file );
    }
    fprintf(ps_file, "\nshowpage\n" );
    fflush( ps_file );
    
    for( a = 1; a <= number_of_pages; a++ ) {
      sprintf( file_name, "%d", a );
      current_file = fopen( file_name, "r" );
      while( fgets( linebuf, sizeof(linebuf), current_file ) != NULL ){
	fprintf(ps_file, "%s", linebuf );
      }
    }
    
    fprintf(ps_file, "%%%%Page:    %d   %d\n", *ps_page_number,
	    *ps_page_number );
  }
  /**********************************************************************/
  
}
/**********************************************************************/

/**********************************************************************/
/* Prints the actual colorbars corresponding to the the rigid cluster */
/* decomposition. The information for each cluster in included as a   */
/* linked list. The lists for each cluster are indexed from the array */
/* nc_count[]. This routine is for single page output only (full size */
/* portrait, full size landscape, scaled "fit-it-all-on-one-page").   */
/**********************************************************************/
void print_decomp( clusters *nc_count[120], int cluster_counter,
		   int ps_line_number, FILE *ps_file, int colors[120],
		   int y_translate[10], int number_of_residues, 
		   int number_of_chains, int chain_size[10], char chain_IDs[10],
		   int total_insertions, int insertion_res_num[100],
		   int number_of_insertions[100], char insertion_chain_ID[100], int lrc[1000] ){

  int 
    a = 0, 
    b = 0,
    c = 0,
    d = 0,
    e = 0,
    this_chain = 0,
    top = 0, 
    bottom = 0, 
    y_position = 0,
    start_pt = 0, 
    end_pt = 0,
    previous_chains_size = 0,
    next_residue = 0,
    //print_this_number = 0,  // FIXME - unused variable
    print_number = 1,
    rr,
    rrhold=-1,
    file_flag=0,
    move_up = 0;

  FILE
    *map_file = NULL;

  clusters  
    *nc_current = NULL;

  /*************************************************************/
  /* First need to compute how the output will fit on the page */
  /* If there are <= 125 residues, use the portrait format. If */
  /* it's 125 < residues >= 250, use landscape. Greater than   */
  /* 250 residues will require multiple page printouts, there  */
  /* is an option to print it all on one page.                 */
  /* In the output, the height of the color bars is always the */
  /* same, only the width should change.                       */
  /*************************************************************/
  // checks if the print_decomp has been run before without using an additional parameter
  if(!lrc[0]) { 
    map_file   = fopen("residue_map.txt" , "w" );
    //    map_file   = fopen(filename, "w" );
    fprintf(map_file,"#res_index    res_number  chain_id\n");
    file_flag =1;
    //    printf("Initial decomp.  %d  %d\n",number_of_chains,number_of_residues);
  }
  for(a=0;a<1000;a++) {
    lrc[a]=0;
  }
  if( !number_of_chains && (number_of_residues <= 125) ) {
    
    start_pt = 110;

    /**********************************************************************/
    /* Print residue numbers, sequence lines, asterisks, etc, ..          */
    /**********************************************************************/
    for( a = 0; a <= number_of_chains; a++ ){
      
      end_pt = ( chain_size[a] * 3 ) + start_pt;

      /************************************************************/
      /* Draw the black line representing the primary structure.  */
      fprintf(ps_file, "\nCol00\n");
      fprintf(ps_file, "%3d %3d %3d %3d sequence_line\n", start_pt,
	      ps_line_number, end_pt, ps_line_number );
      fprintf(ps_file, "grestore\n");
      /************************************************************/

      /**********************************************************************/
      /* Print the sequence numbers at the top of the page. Also, if the    */
      /* proteins contained labeled insertions, print an asterisk over the  */
      /* residue.                                                           */
      /**********************************************************************/
      if( ps_line_number == 710 ) {
	fprintf(ps_file, "%d 722 moveto\n", 113 );
	fprintf(ps_file, "(%d) 8.0 UncenterRot90 Rot90\n", y_translate[a] );
	fprintf(ps_file, "(%d) 8.0 Print\n", y_translate[a] );
	fprintf(ps_file, "grestore\n");
	fprintf(ps_file, "%5.2f 715 %5.2f 720 sequence_line\n", start_pt + 1.5, start_pt + 1.5 );
	fprintf(ps_file, "grestore\n\n");

	/**********************************************************************/
	/* Display an asterik over each inserted residue in the protein.      */
	/**********************************************************************/
	for( d = 0; d < total_insertions; d++ ){
	  if( insertion_chain_ID[d] == chain_IDs[a] ){

	    for( e = 0; e < number_of_insertions[d]; e++ ){
	      
	      if( insertion_res_num[d] != 1 ){
		fprintf( ps_file, "%5.1f 711.5 moveto\n",  start_pt + 7.5 + (3* (insertion_res_num[d] + e - y_translate[a] ))  );
		printf("%d %d %d %d\n", insertion_res_num[d], e, start_pt + 5 + (3* (insertion_res_num[d] + e - y_translate[a] )), 
		       number_of_insertions[d] );
		fprintf(ps_file, "(*) 10.0 Print\n" );
	      }
	      else{
		fprintf( ps_file, "%5.1f 711.5 moveto\n",  start_pt - 1.0 + (3 * e)  );
		printf("%d %d %d %d\n", insertion_res_num[d], e, start_pt + 5 + (3* (insertion_res_num[d] + e - y_translate[a] )), 
		       number_of_insertions[d] );
		fprintf(ps_file, "(*) 10.0 Print\n" );
	      }
		
	    }
	  }
	}
	/**********************************************************************/

	next_residue = (y_translate[a] + 10) / 10; 
	next_residue *= 10; 
	for( b = next_residue; b < (chain_size[a]+ y_translate[a]); b+=10 ){
	  fprintf(ps_file, "%d 722 moveto\n", start_pt + 3 + (3*(b - y_translate[a] )));
	  fprintf(ps_file, "(%d) 8.0 UncenterRot90 Rot90\n", b );
	  fprintf(ps_file, "(%d) 8.0 Print\n", b );
	  fprintf(ps_file, "grestore\n\n");
	  fprintf(ps_file, "%5.2f 715 %5.2f 720 sequence_line\n", start_pt + 1.5 + (3*(b - y_translate[a] )), 
		  111.5 + (3*(b - y_translate[a] )) );
	  fprintf(ps_file, "grestore\n\n");
	}
      }
      
      start_pt = end_pt + 18;    
    }
    /**********************************************************************/
    /**********************************************************************/
    
    for( a = 0; a < cluster_counter; a++ ){

      nc_current = nc_count[a];
      while( nc_current ) {
	
	if( colors[a] < 10 )
	  fprintf(ps_file,"Col0%d\n", colors[a] );
	else if( colors[a] <=80 )
	  fprintf(ps_file,"Col%d\n", colors[a] );
	else
	  fprintf(ps_file,"Col%d\n", int( (colors[a]-40)%40 +40 ) );

	if( number_of_chains ){
	  for( b = 0; b <= number_of_chains; b++ ){

	    if( nc_current->chain_ID == chain_IDs[b] ) {
	      top    = ps_line_number + 3 - (10*b);
	      bottom = ps_line_number - 3 - (10*b);
	      y_position = nc_current->residue_number - y_translate[b];
	    }
	  }
	}
	
	else{
	  top    = ps_line_number + 3;
	  bottom = ps_line_number - 3;
	  y_position = nc_current->residue_number - y_translate[0];
	}

	/**********************************************************************/
	/* Modify spacing to include insertions that share the same residue   */
	/* numbering.                                                         */
	/**********************************************************************/
	for( d = 0; d < total_insertions; d++ ){
	  
	  if( nc_current->residue_number == insertion_res_num[d] &&
	      nc_current->chain_ID == insertion_chain_ID[d] )
	    y_position += nc_current->insertion_space;
	  
	  if( nc_current->residue_number > insertion_res_num[d] &&
	      nc_current->chain_ID == insertion_chain_ID[d] )
	    y_position += number_of_insertions[d];

	}
	/**********************************************************************/
	
	/*printf("y_position %3d  res_num %3d,  %d %d %d\n", y_position, nc_current->residue_number, 
	  total_insertions, insertion_res_num[0], number_of_insertions[0],
	  insertion_chain_ID[0]);*/
	rr = y_position;	
	if(file_flag && rr != rrhold) {
	  fprintf(map_file,"%5d %5d %c\n",rr,nc_current->residue_number,nc_current->chain_ID);
	  rrhold = rr;
	}

	if( !strcmp( nc_current->atom_type, "N" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  (y_position*3)+110, top,
		  (y_position*3)+111, top,
		  (y_position*3)+111, bottom,
		  (y_position*3)+110, bottom );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}
	
	if( !strcmp( nc_current->atom_type, "CA" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  (y_position*3)+111, top,    
		  (y_position*3)+112, top,    
		  (y_position*3)+112, bottom, 
		  (y_position*3)+111, bottom );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}
	
	if( !strcmp( nc_current->atom_type, "C" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  (y_position*3)+112, top,    
		  (y_position*3)+113, top,     
		  (y_position*3)+113, bottom, 
		  (y_position*3)+112, bottom );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}
	//For nucleic acids
	if( !strcmp( nc_current->atom_type, "P" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  (y_position*3)+110, top,
		  (y_position*3)+111, top,
		  (y_position*3)+111, bottom,
		  (y_position*3)+110, bottom );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}
	
	if( !strcmp( nc_current->atom_type, "O5'" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  (y_position*3)+111, top,
		  (y_position*3)+112, top,
		  (y_position*3)+112, bottom,
		  (y_position*3)+111, bottom );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}

	if( !strcmp( nc_current->atom_type, "C5'" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  (y_position*3)+112, top,
		  (y_position*3)+113, top,
		  (y_position*3)+113, bottom,
		  (y_position*3)+112, bottom );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}


	nc_current = nc_current->next_element;
      }
    }
    fflush( ps_file );
  }

  /**********************************************************************/
  /* print the decomp in single page landscape for proteins with        */
  /* 125 > amino-acids <= 175. The following code is also executed when */
  /* hbdilute.c is run with the "b" option for, fit_it_all_on_one_page. */
  /**********************************************************************/
  else{
     start_pt = 110;

    for( a = 0; a <= number_of_chains; a++ ){
      
      end_pt = ( chain_size[a] * 3 ) + start_pt;

      /************************************************************/
      /* Draw a black line representing the primary structure.    */
      fprintf(ps_file, "\nCol00\n");
      fprintf(ps_file, "%3d %3d %3d %3d sequence_line\n", ps_line_number,
	      start_pt, ps_line_number, end_pt );
      fprintf(ps_file, "grestore\n");
      /************************************************************/

      /************************************************************/
      /* For landscape output, the top line of the page x=60. If  */
      /* ps_line_number == 60, print the residue numbers for each */
      /* chain.                                                   */
      if( ps_line_number == 60 ) {
	fprintf(ps_file, "38 %d  moveto\n", start_pt -1 );
	fprintf(ps_file, "(%d) 8.0 Print\n", y_translate[a] );
	fprintf(ps_file, "50 111.5 55 111.5 sequence_line\n"); 

	/**********************************************************************/
	/* Display an asterik over each inserted residue in the protein.      */
	/**********************************************************************/
	for( d = 0; d < total_insertions; d++ ){

	  if( insertion_chain_ID[d] == chain_IDs[a] ){

	    for( e = 0; e < number_of_insertions[d]; e++ ){
	      if( insertion_res_num[d] != 1 ){
		fprintf( ps_file, "50 %5.1f  moveto\n",  start_pt + 4.5 + ( 3*(insertion_res_num[d] + e - y_translate[a] )) -4 );
		fprintf(ps_file, "(*) 10.0 Print\n" );
	      }
	      else{
		fprintf( ps_file, "50 %5.1f  moveto\n",  start_pt + 1.5 + (3 * e) -4  );
		fprintf(ps_file, "(*) 10.0 Print\n" );
	      }
		
	    }
	  }
	}
	/**********************************************************************/

	/************************************************************/
	/* Because the first residue in the protein need not be 1 or*/
	/* a multiple of 10, find the next residue number to print  */
	/* to the top of the page such that it is at least 10       */
	/* greater than the first residue. ie. if the first residue */
	/* in the protein is 34, the next res num to print wound be */
	/* 50.                                                      */
	next_residue = (y_translate[a] + 10) / 10; 
	next_residue *= 10; 

	for( b = next_residue; b < (chain_size[a] + y_translate[a]); b+=10 ){
	  
	  /*
	  print_this_number = b;
	  print_number = TRUE;
	  */
	  
	  move_up = 0;
	  for( d = 0; d < total_insertions; d++ ){

	    if( insertion_chain_ID[d] == chain_IDs[a] ){
	      if( insertion_res_num[d] < b ){
		move_up += number_of_insertions[d];
	      }
	      //if( insertion_res_num[d] == b ){
	      //b -= 9; 
	      //print_number = FALSE;
	      //}
	    }

	  }

	  if( print_number ){
	    if( b < 100 ){
	      fprintf(ps_file, "38 %d  moveto\n",  start_pt -1 +(3*(b - y_translate[a] +move_up)) );
	      fprintf(ps_file, "(%d) 8.0 Print\n", b );
	      fprintf(ps_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt + 1.5 + (3*(b - y_translate[a] +move_up)), 
		      start_pt + 1.5 + (3*(b - y_translate[a] +move_up)) );
	    }
	    else{
	      fprintf(ps_file, "34 %d  moveto\n",  start_pt -1 +(3*(b - y_translate[a] +move_up)) );
	      fprintf(ps_file, "(%d) 8.0 Print\n", b );
	      fprintf(ps_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt + 1.5 + (3*(b - y_translate[a] +move_up)), 
		    start_pt + 1.5 + (3*(b - y_translate[a] +move_up)) );
	    }
	  }

	}
      }
      /************************************************************/

      start_pt = end_pt + 18;    
    }
    
    for( a = 0; a < cluster_counter; a++ ){
      
      nc_current = nc_count[a];

      while( nc_current ) {
	
	//if( nc_current != NULL )
	//printf("how did i get here if nc_current is null? %d\n", nc_current->atom_number);
	
	if( colors[a] < 10 )
	  fprintf(ps_file,"Col0%d\n", colors[a] );
	else if( colors[a] <=80 )
	  fprintf(ps_file,"Col%d\n", colors[a] );
	else
	  fprintf(ps_file,"Col%d\n", int( (colors[a]-40)%40 +40 ) );
	
	if( number_of_chains ){


	  for( b = 0; b <= number_of_chains; b++ ){

	    if( nc_current->chain_ID == chain_IDs[b] ) {

	      this_chain = b;
	      top    = ps_line_number + 3;
	      bottom = ps_line_number - 3;
	      for( c = 0; c < b; c++ )
		previous_chains_size += chain_size[c];
	      y_position = (nc_current->residue_number - y_translate[b] + (6*b) + previous_chains_size) * 3;
	      previous_chains_size = 0;
	    }

	    /**********************************************************************/
	    /* Modify spacing to include insertions that share the same residue   */
	    /* numbering.                                                         */
	    /**********************************************************************/
	    for( d = 0; d < total_insertions; d++ ){
	      if( insertion_chain_ID[d] == nc_current->chain_ID &&
		  insertion_chain_ID[d] == chain_IDs[b] ){

		if( nc_current->residue_number > insertion_res_num[d] ){
		  y_position += ((number_of_insertions[d]) * 3);
		}

		if( nc_current->residue_number == insertion_res_num[d] ){

		  if( nc_current->residue_number == 1 ){
		    if( nc_current->insertion_space )
		      y_position += ((nc_current->insertion_space - 1) * 3);
		    else
		      y_position += (number_of_insertions[d]) * 3;
		  }
		  else
		    y_position += (nc_current->insertion_space) * 3;
		}
	      }

	    }
	    /**********************************************************************/

	  }
	}
	
	else{
	  top    = ps_line_number + 3;
	  bottom = ps_line_number - 3;
	  y_position = ((nc_current->residue_number)*3) - (y_translate[0]*3);
	}

	/*printf("y_position %3d  res_num %3d,  %d %d %d  %d\n", y_position, nc_current->residue_number, 
	  total_insertions, insertion_res_num[d], number_of_insertions[d], y_translate[this_chain] );*/
	rr = y_position/3;	
	if( number_of_chains ){
	  for( b = 1; b <= number_of_chains; b++ ){
	    if( nc_current->chain_ID == chain_IDs[b] ) {
	      rr -= 6*b;
	    }
	  }
	}
  //printf("resval %4d  res_num %4d  %c \n", rr, nc_current->residue_number, nc_current->chain_ID);
	if(file_flag && rr != rrhold) {
	  fprintf(map_file,"%5d %5d %c\n",rr,nc_current->residue_number,nc_current->chain_ID);
	  //	  printf("%5d %5d %c\n",rr,nc_current->residue_number,nc_current->chain_ID);
	  rrhold = rr;
	}
	if( !strcmp( nc_current->atom_type, "N" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+110,
		  top,    (y_position)+111, 
		  bottom, (y_position)+111, 
		  bottom, (y_position)+110 );
	  if(a==0) {
	    lrc[rr]++;
	  }
	}
	
	if( !strcmp( nc_current->atom_type, "CA" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+111,
		  top,    (y_position)+112, 
		  bottom, (y_position)+112, 
		  bottom, (y_position)+111 );
	  if(a==0) {
	    //	    rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}
	
	if( !strcmp( nc_current->atom_type, "C" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+112,
		  top,    (y_position)+113, 
		  bottom, (y_position)+113,
		  bottom, (y_position)+112 );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}

	if( !strcmp( nc_current->atom_type, "P" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+110,
		  top,    (y_position)+111, 
		  bottom, (y_position)+111,
		  bottom, (y_position)+110 );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}

	if( !strcmp( nc_current->atom_type, "O5'" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+111,
		  top,    (y_position)+112, 
		  bottom, (y_position)+112,
		  bottom, (y_position)+111 );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}

	if( !strcmp( nc_current->atom_type, "C5'" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+112,
		  top,    (y_position)+113, 
		  bottom, (y_position)+113,
		  bottom, (y_position)+112 );
	  if(a==0) {
	    //rr = nc_current->residue_number;
	    lrc[rr]++;
	  }
	}

	nc_current = nc_current->next_element;
      }
    }
    //    fclose(map_file);
  }

}  

/**********************************************************************/
/* Prints the actual colorbars corresponding to the the rigid cluster */
/* decomposition. The information for each cluster in included as a   */
/* linked list. The lists for each cluster are indexed from the array */
/* nc_count[].                                                        */
/**********************************************************************/
void print_multipage_decomp( clusters *nc_count[120], int cluster_counter,
			     int ps_line_number, FILE *ps_file, int colors[120],
			     int y_translate[10], int number_of_residues,
			     FILE *file_list[5], int number_of_pages,
			     int number_of_chains, int chain_size[10], 
			     char chain_IDs[10], int total_insertions, 
			     int insertion_res_num[100], int number_of_insertions[100], 
			     char insertion_chain_ID[100] ){

  int 
    a = 0, 
    b = 0, 
    c = 0, 
    d = 0,
    top = 0, 
    bottom = 0, 
    y_position = 0, 
    file_number = 0, 
    //end_point = 0,  // FIXME - unused variable
    start_pt=0, 
    summed_length = 0, /* length of all the chains so far */
    end_pt = 0, 
    next_residue = 0,
    last_res_num_on_this_page = 0,
    current_pg = 0, 
    previous_length = 0, 
    seq_page = 0,
    seq_num = 0, 
    //space = 0,  // FIXME - unused variable
    times = 0, 
    insertion_space = 0, 
    residue_pointer = 0,
    printed_first_residue_number_yet = 0; 

  clusters  
    *nc_current = NULL;

  FILE 
    *current_file = NULL; 
  
  /*****************************************************/
  /* print the sequence line (black line) on all pages */
  /* and a scale at the top of the page indicating the */
  /* residue number.                                   */
  /*****************************************************/
  start_pt = 110;
  current_pg = 0;

  for( a = 0; a <= number_of_chains; a++ ){
    
    printed_first_residue_number_yet = FALSE;
    seq_num  = 0;
    seq_page = 0;    
    current_file = file_list[current_pg];
    end_pt = ( chain_size[a] * 3 ) + start_pt;
    
    /**********************************************************************/
    /* If end_pt > 710, then the sequence line extends off the page. Make */
    /* accomodations to print up to y=710, then put the rest on the next  */
    /* page.                                                              */
    /**********************************************************************/
    if( end_pt > 710 ){
      
      /**********************************************************************/
      /* Print what you can on the current page. Keep printing the current  */
      /* chain until the end_pt will fit on the page, then code to the next */
      /* snipet of code.                                                    */
      /**********************************************************************/
      while( end_pt > 710 ){
	
	fprintf(current_file, "\nCol00\n");
	fprintf(current_file, "%3d %3d %3d 710 sequence_line\n", ps_line_number,
		start_pt, ps_line_number );
	fprintf(current_file, "grestore\n");
	
	/**********************************************************************/
	/* print the residue numbers on the top of the page.                  */
	/**********************************************************************/
	if( ps_line_number == 60 ){
	  
	  if( !printed_first_residue_number_yet ){
	    fprintf(current_file, "38 %d  moveto\n", start_pt -1 );
	    fprintf(current_file, "(%d) 8.0 Print\n", y_translate[a] );
	    fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt +1.5, start_pt +1.5); 
	    
	    /************************************************************/
	    /* Because the first residue in the protein need not be 1 or*/
	    /* a multiple of 10, find the next residue number to print  */
	    /* to the top of the page such that it is at least 10       */
	    /* greater than the first residue. ie. if the first residue */
	    /* in the protein is 34, the next res num to print wound be */
	    /* 50.                                                      */
	    /************************************************************/
	    next_residue = (y_translate[a] + 10) / 10; 
	    next_residue *= 10; 
	    printed_first_residue_number_yet = TRUE;
	    
	    /************************************************************/
	    /* set the upper bound for residue numbering on this page.  */
	    /************************************************************/
	    last_res_num_on_this_page = chain_size[a] - ((end_pt - 710)/3) + y_translate[a];

	    for( b = next_residue; b < last_res_num_on_this_page; b+=10 ){
	      if( b < 100 ){
		fprintf(current_file, "38 %d  moveto\n",  start_pt -1 +(3*(b - y_translate[a] )) );
		fprintf(current_file, "(%d) 8.0 Print\n", b );
		fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", 
			start_pt +1.5 + (3*(b - y_translate[a])), 
			start_pt + 1.5 + (3*(b - y_translate[a] )) );
	      }
	      else{
		fprintf(current_file, "34 %d  moveto\n",  start_pt -1 +(3*(b - y_translate[a] )) );
		fprintf(current_file, "(%d) 8.0 Print\n", b );
		fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", 
			start_pt +1.5 + (3*(b - y_translate[a])), 
			start_pt + 1.5 + (3*(b - y_translate[a] )) );
	      }
	    }
	  }
	  else{
	    last_res_num_on_this_page = chain_size[a] - ((end_pt - 710)/3) + y_translate[a];

	    for( b = seq_num; b < last_res_num_on_this_page; b+=10 ){
	      if( b < 100 ){
		fprintf(current_file, "38 %d  moveto\n",  
			start_pt -1 +(3*(b - (last_res_num_on_this_page -200 )) ) );
		fprintf(current_file, "(%d) 8.0 Print\n", b );
		fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", 
			start_pt +1.5 + (3*(b - (last_res_num_on_this_page -200) ) ), 
			start_pt + 1.5 + (3*(b -  (last_res_num_on_this_page -200 ))) );
	      }
	      else{
		fprintf(current_file, "34 %d  moveto\n",  
			start_pt -1 +(3*(b - (last_res_num_on_this_page -200 ) )) );
		fprintf(current_file, "(%d) 8.0 Print\n", b );
		fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", 
			start_pt +1.5 + (3*(b - (last_res_num_on_this_page - 200) ) ), 
			start_pt + 1.5 + (3*(b - (last_res_num_on_this_page -200 )) ) );
	      }
	    }
	    
	  }
	}
	
	current_pg++;
	start_pt = 110;
	end_pt   = ( end_pt - 600 );
	current_file = file_list[current_pg];
	seq_page++;
	seq_num = b;
      }
      
      /**********************************************************************/
      /* Print the rest on the next page. The variable "seq_num" will set a */
      /* marker for where to begin the sequence numbering.                  */
      /**********************************************************************/      
      fprintf(current_file, "\nCol00\n");
      fprintf(current_file, "%3d %3d %3d %3d sequence_line\n", ps_line_number,
	      start_pt, ps_line_number, end_pt );
      fprintf(current_file, "grestore\n");

      if( ps_line_number == 60 ){
	
	for( c = seq_num; c < (chain_size[a]+y_translate[a]); c+=10 ){
	  if( c < 100 ){
	    fprintf(current_file, "38 %d  moveto\n", start_pt -1 +(3*(c - last_res_num_on_this_page)) );
	    fprintf(current_file, "(%d) 8.0 Print\n", c );
	    fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt +1.5 + (3*(c - last_res_num_on_this_page)), 
		    start_pt + 1.5 + (3*(c - last_res_num_on_this_page)) );
	  }
	  else{
	    fprintf(current_file, "34 %d  moveto\n", start_pt -1 +(3*(c - last_res_num_on_this_page)) );
	    fprintf(current_file, "(%d) 8.0 Print\n", c );
	    fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt +1.5 + (3*(c - last_res_num_on_this_page)), 
		    start_pt + 1.5 + (3*(c - last_res_num_on_this_page)) );
	  }
	}
      }
      start_pt = end_pt + 18;
    }
    
    /**********************************************************************/
    /* use this code if the sequence line will fit on the current page.   */
    /* (That is, end_pt <= 710).                                          */
    /**********************************************************************/
    else{
      fprintf(current_file, "\nCol00\n");
      fprintf(current_file, "%3d %3d %3d %3d sequence_line\n", ps_line_number,
	      start_pt, ps_line_number, end_pt );
      fprintf(current_file, "grestore\n");

      if( ps_line_number == 60 ){
	fprintf(current_file, "38 %d  moveto\n", start_pt -1 );
	fprintf(current_file, "(%d) 8.0 Print\n", y_translate[a] );
	fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt +1.5, start_pt +1.5); 
	/************************************************************/
	/* Because the first residue in the protein need not be 1 or*/
	/* a multiple of 10, find the next residue number to print  */
	/* to the top of the page such that it is at least 10       */
	/* greater than the first residue. ie. if the first residue */
	/* in the protein is 34, the next res num to print wound be */
	/* 50.                                                      */
	next_residue = (y_translate[a] + 10) / 10; 
	next_residue *= 10; 
	for( b = next_residue; b < (chain_size[a] + y_translate[a]); b+=10 ){
	  if( b < 100 ){
	    fprintf(current_file, "38 %d  moveto\n",  start_pt -1 +(3*(b - y_translate[a] )) );
	    fprintf(current_file, "(%d) 8.0 Print\n", b );
	    fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt +1.5 + (3*(b - y_translate[a])), 
		    start_pt + 1.5 + (3*(b - y_translate[a] )) );
	  }
	  else{
	    fprintf(current_file, "34 %d  moveto\n",  start_pt -1 +(3*(b - y_translate[a] )) );
	    fprintf(current_file, "(%d) 8.0 Print\n", b );
	    fprintf(current_file, "50 %5.2f 55 %5.2f sequence_line\n", start_pt +1.5 + (3*(b - y_translate[a])), 
		    start_pt + 1.5 + (3*(b - y_translate[a] )) );
	  }
	}
	
      }
    }    

    start_pt = end_pt + 18;    
  }
  /**********************************************************************/
  /* End sequence line, numbering, tick marks, etc...                   */
  /**********************************************************************/
  /**********************************************************************/
  

  /**********************************************************************/
  /* The following code prints the actual colorbar text to the          */
  /* postscript file.                                                   */
  /**********************************************************************/
  for( a = 0; a < cluster_counter; a++ ){
    
    nc_current = nc_count[a];

    while( nc_current ) {
      
      summed_length = previous_length = times = 0;
      y_position = insertion_space = 0;

      /**********************************************************************/
      /* If there are 2 or more chains, find the right output file.         */
      /**********************************************************************/
      if( number_of_chains ){

	for( b = 0; b <= number_of_chains; b++ ){

	  summed_length += (chain_size[b] + ( 6*b ));

	  if( nc_current->chain_ID == chain_IDs[b] ) {

	    /**********************************************************************/
	    /* Modify spacing to include insertions that share the same residue   */
	    /* numbering.                                                         */
	    /**********************************************************************/
	    for( d = 0; d < total_insertions; d++ ){
	      if( insertion_chain_ID[d] == nc_current->chain_ID &&
		  insertion_chain_ID[d] == chain_IDs[b] ){

		if( nc_current->residue_number > insertion_res_num[d] ){
		  y_position += ((number_of_insertions[d]) * 3);
		  insertion_space += number_of_insertions[d];
		}

		if( nc_current->residue_number == insertion_res_num[d] ){
		  
		  if( nc_current->residue_number == 1 ){
		    if( nc_current->insertion_space ){
		      y_position += ((nc_current->insertion_space - 1) * 3);
		      insertion_space += nc_current->insertion_space -1;
		    }
		    else{
		      y_position += (number_of_insertions[d]) * 3;
		      insertion_space += number_of_insertions[d];
		    }
		  }
		  else{
		    y_position += (nc_current->insertion_space) * 3;
		    insertion_space += nc_current->insertion_space;
		  }
		}
	      }

	    }
	    /**********************************************************************/
	    
	    /**********************************************************************/
	    /* Select which temp file to print to.                                */
	    file_number = ( summed_length -chain_size[b] +nc_current->residue_number -y_translate[b] +insertion_space) / 200; 
	    printf("file number %d\n", file_number );
	    current_file = file_list[file_number];
	    /**********************************************************************/
	    
	    /**********************************************************************/
	    /* Tell the postscript file what color to print the current atom.     */
	    if( colors[a] < 10 )
	      fprintf(current_file,"Col0%d\n", colors[a] );
	    else
	      fprintf(current_file,"Col%d\n", colors[a] );
	    /**********************************************************************/

	    /**********************************************************************/
	    /* Calculate where on the current page the rigid cluster information  */
	    /* should be printed. This will depend on how many chains have been   */
	    /* printed previously (summed_length) and how much of the current     */
	    /* chain was printed on the last page.                                */
	    previous_length = summed_length - chain_size[b];
	    while( previous_length > 200 ) {
	      previous_length -= 200;
	    }

	    residue_pointer = (nc_current->residue_number) -y_translate[b];
	    while( residue_pointer >= 200 ) {
	      residue_pointer -= 200;
	    }
	    
	    if( ( 200 - previous_length ) > ( chain_size[b] ) )
	      y_position += (residue_pointer +1 +previous_length) *3;

	    else {
	      if( residue_pointer >= ( 200-previous_length) )
		y_position += (residue_pointer -( 200 -previous_length ) ) *3;
	      /*will fit on page with the previous length */
	      else 
		y_position += (residue_pointer +previous_length) *3;
	    }
	    /**********************************************************************/

	  }
	}
	top    = ps_line_number + 3;
	bottom = ps_line_number - 3;
	
      }
      else{
	file_number = ( nc_current->residue_number -y_translate[0] ) / 200; 
	current_file = file_list[file_number];
        top    = ps_line_number + 3;
	bottom = ps_line_number - 3;
	y_position = ( nc_current->residue_number -y_translate[0] -(file_number*200) ) *3;	

	if( colors[a] < 10 )
	  fprintf(current_file,"Col0%d\n", colors[a] );
	else
	  fprintf(current_file,"Col%d\n", colors[a] );
      }
	
      if( !strcmp( nc_current->atom_type, "N" ) ){
	fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		top,    y_position +110,
		top,    y_position +111, 
		bottom, y_position +111, 
		bottom, y_position +110 );
      }
	
	if( !strcmp( nc_current->atom_type, "CA" ) ){
	  fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    y_position +111,
		  top,    y_position +112, 
		  bottom, y_position +112, 
		  bottom, y_position +111 );
	}
	
	if( !strcmp( nc_current->atom_type, "C" ) ){
	  fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    y_position +112,
		  top,    y_position +113, 
		  bottom, y_position +113,
		  bottom, y_position +112 );
	}
	
	if( !strcmp( nc_current->atom_type, "P" ) ){
	  fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    y_position +110,
		  top,    y_position +111, 
		  bottom, y_position +111,
		  bottom, y_position +110 );
	}

	if( !strcmp( nc_current->atom_type, "O5'" ) ){
	  fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    y_position +111,
		  top,    y_position +112, 
		  bottom, y_position +112,
		  bottom, y_position +111 );
	}

	if( !strcmp( nc_current->atom_type, "C5'" ) ){
	  fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    y_position +112,
		  top,    y_position +113, 
		  bottom, y_position +113,
		  bottom, y_position +112 );
	}

	nc_current = nc_current->next_element;
    }
  }
}

/**********************************************************************/
