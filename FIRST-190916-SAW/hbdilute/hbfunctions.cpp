/**********************************************************************/
/* These are the main functions that read in the decomp data, sort it */
/* and produce the color consistency in the decomposition output.     */
/*                                                                    */
/* BMH 6.19.02 Added code in read_chem_file to read insertion codes.  */
/* The logic isn't perfect, it requires that all three backbone atoms */
/* be present in the file (N, CA, C).                                 */
/**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "types.h"
#include "postscript.h"

#define column(n)   (linebuf+(n)-1) /* manipulate linebuf by column */

/***********************************************************************/
/* This routine reads the data for the main-chain atoms in the protein */
/* from the *_FIRSTdataset file and stores them in a linked list named */
/* chem_file_head. A pointer to the beginning of the list is passed    */
/* back to hbdilute.c                                                  */
/***********************************************************************/
void read_chem_file( atom_list **chem_file_head, FILE *pdb_file, int y_translate[10],
		     int chain_size[10], char chain_IDs[10], int *number_of_chains, 
		     int *total_insertions, int insertion_res_num[100], 
		     int number_of_insertions[100], char insertion_chain_ID[100] ){
  
  int  
    a = 0, 
    is_atom = 0,
    known_insertion = 0;

  char 
    linebuf[150], 
    atom[4], 
    chain_label = '\0';

  atom_list *temp = NULL,
            *current = NULL;

  for( a = 0; a < 20; a++ )
    number_of_insertions[a] = 0;


  while( fgets( linebuf, sizeof(linebuf), pdb_file ) != NULL ) {
    
    is_atom   = !strncmp( linebuf, "ATOM", 4 );
    sscanf( linebuf+13, "%3s", atom );

    if( is_atom && 
	( !strcmp( atom, "N" ) ||
	  !strcmp( atom, "CA") ||
	  !strcmp( atom, "C" ) ||
	  !strcmp( atom, "O5'") ||
	  !strcmp( atom, "C5'") ||
	  !strcmp( atom, "P") ) ) {

      temp = ( atom_list *) malloc( sizeof( atom_list ));
      temp->atom_number = atoi( (linebuf+5)  );
      sscanf( linebuf+13, "%3s", temp->atom_type );
      sscanf( linebuf+22, "%4d", &temp->residue );
      temp->chain_ID = *column(22);
      temp->atom_or_hetatm = 1;

      /**********************************************************************/
      /* The following code was added to handle insertions that may be found*/
      /* in the PDB file. Assumes PDB file format standards.                */
      /**********************************************************************/
      if( *(linebuf+26) != ' ' ){
	
	if( *total_insertions == 0 ){
	  if( !strcmp( atom, "N" ) ){
	    insertion_res_num[*total_insertions] = temp->residue;
	    insertion_chain_ID[*total_insertions] = temp->chain_ID;
	    temp->insertion_space = 1; 
	  }
	  if( !strcmp( atom, "CA" ) )
	    temp->insertion_space = 1; 
	  if( !strcmp( atom, "C" ) ){
	    temp->insertion_space = 1; 
	    (number_of_insertions[*total_insertions])++;
	    (*total_insertions)++;
	  }
	}

      	else{
	  for( a = 0; a < *total_insertions; a++ ){
	    if( temp->residue == insertion_res_num[a] &&
		temp->chain_ID == insertion_chain_ID[a] ){
	      if( !strcmp( atom, "N" ) )
		temp->insertion_space = number_of_insertions[a] + 1; 
	      if( !strcmp( atom, "CA" ) )
		temp->insertion_space = number_of_insertions[a] + 1; 
	      if( !strcmp( atom, "C" ) ){
		temp->insertion_space = number_of_insertions[a] + 1; 
		(number_of_insertions[a])++;
	      }
	      known_insertion = -1; 
	    }
	  }
	  
	  if( !known_insertion ){
	    if( !strcmp( atom, "N" ) ){
	      insertion_res_num[*total_insertions] = temp->residue;
	      insertion_chain_ID[*total_insertions] = temp->chain_ID;
	      temp->insertion_space = 1; 
	    }
	    if( !strcmp( atom, "CA" ) )
	      temp->insertion_space = 1; 
	    if( !strcmp( atom, "C" ) ){
	      temp->insertion_space = 1; 
	      (number_of_insertions[*total_insertions])++;;
	      (*total_insertions)++;
	    }
	  }
	  known_insertion = 0; 
	}
      }
      /**********************************************************************/

      /*printf("residue: %d   total %d: number_of_insertions at residue %d: %3d\n", temp->residue, 
	     *total_insertions, insertion_res_num[(*total_insertions)-1],
	     number_of_insertions[(*total_insertions)-1] );*/

      if( *chem_file_head == NULL ){
	*chem_file_head = current = temp;
	y_translate[*number_of_chains] = current->residue;
	chain_label = temp->chain_ID;
	chain_IDs[*number_of_chains] = temp->chain_ID;
      }
      else{
	current->next_atom = temp;
	current = current->next_atom;
	current->next_atom = NULL;
      }
      
      if( !strcmp( atom, "CA") ) {

	if( current->chain_ID == chain_label ) {
	  (chain_size[*number_of_chains])++;
	}
	else{
	  chain_label = current->chain_ID;
	  (*number_of_chains)++;
	  y_translate[*number_of_chains] = current->residue;
	  (chain_size[*number_of_chains])++;
	  chain_IDs[*number_of_chains] = current->chain_ID;
	}

	if( current->residue < y_translate[*number_of_chains] )
	  y_translate[*number_of_chains] = current->residue;
	
      }

      if( !strcmp( atom, "P") ) {

	if( current->chain_ID == chain_label ) {
	  (chain_size[*number_of_chains])++;
	}
	else{
	  chain_label = current->chain_ID;
	  (*number_of_chains)++;
	  y_translate[*number_of_chains] = current->residue;
	  (chain_size[*number_of_chains])++;
	  chain_IDs[*number_of_chains] = current->chain_ID;
	}

	if( current->residue < y_translate[*number_of_chains] )
	  y_translate[*number_of_chains] = current->residue;
	
      }

      
    }
  }
}
/**********************************************************************/
/**********************************************************************/

/**********************************************************************/
/* This rountine will determine the number of clusters in the current */
/* decomposition. Since only main-chain atoms are counted, a lower    */
/* bound on cluster size must be set. The variable is minimum_cluster */
/* _size. I have been using 3, corresponding to a single residue. The */
/* date for this calculation comes from the decomp_list file.         */
/**********************************************************************/
  int compute_number_of_clusters( int atom_label[10000], atom_list *chem_file,
				int cluster_label[120], int number_of_atoms ){
  
  int 
    a = 0,
    tally = 0, 
    //current_number = 0, // FIXME - unused variable
    total_clusters = 0, 
    //current_label = 1,  // FIXME - unused variable
    //renew_list = 0,     // FIXME - unused variable
    //check = 0,    // FIXME - unused variable
    minimum_cluster_size = 3,
    upper_bound_on_number_of_clusters = 0;
  
  atom_list 
    *start_of_list = NULL;


  /**********************************************************************/
  /* The upper bound on the number of clusters is set to limit the      */
  /* number of times the decomp file needs to be read. Mainchain rigid  */
  /* clusters will not necessarily be labelled by consecutive integers. */
  /* side-chains and hetero groups can form rigid clusters whose label  */
  /* will (possibly) be smaller than a main-chain rigid cluster. It is  */
  /* therefore impossible to simply stop checking the decomp file the   */
  /* first time it fails to find a main-chain rigid cluster. For example*/
  /* rigid cluster 23 may consist of 10 main-chain atoms, so it will be */
  /* added to the cluster_label array. However, cluster 24 might be all */
  /* heteratoms, so we don't add it. cluster 25 could then be 9 main-   */
  /* chain atoms, which we want to pick up. Therefore, i set the upper  */
  /* bound to be 300, as a precaution to make sure we find all the main */
  /* chain clusters of size "minimum_cluster_size" or greater.          */
  /* For the dimeric protein 2cts, with 876 total residues, the largest */
  /* cluster label found was 167, so 300 should be approriate, for now. */
  /*     The size of nc_count and oc_count is related to the total      */
  /* number of clusters found ("total_clusters"), and does not depend on*/
  /* the actual value of the cluster label.                             */
  /**********************************************************************/
  upper_bound_on_number_of_clusters = 3000;

  start_of_list = chem_file;

  for( a = 1; a < upper_bound_on_number_of_clusters; a++ ){

    while( chem_file ){

      if( atom_label[chem_file->atom_number] == a )
	tally++;
      
      if( tally == minimum_cluster_size ) {
	//printf("cluster %d, label %d %d\n", total_clusters, a, chem_file->atom_number );
	cluster_label[total_clusters] = a;
	total_clusters++;
	tally = 0;
	chem_file = NULL;
      }
      
      if( chem_file != NULL )
	chem_file = chem_file->next_atom;
      
    }
    chem_file = start_of_list;
    tally = 0;
  }
  
  
  return( total_clusters );
}    

/**********************************************************************/
/* Using the chem_file_head list, this routine stores the data from   */
/* initial decompostion in the array of structures, nc_count. The ini */
/* decomposition refers to the protein with all the H-bonds present.  */
/**********************************************************************/
void first_decomp( int long_output, int total_clusters, int bulk_atom_label[10000],
		   atom_list *chem_file, clusters *new_clusters[120], int colors[120],
		   int cluster_labels[1000] ) {

  int current_cluster=0, cluster_number=0;

  atom_list    *start_of_list = NULL;

  clusters     *nc_temp    = NULL,
               *nc_current = NULL;
  
  start_of_list = chem_file;

  if( !long_output ){

    for( current_cluster = 0; current_cluster < total_clusters; current_cluster++ ) {

      colors[current_cluster] = current_cluster+1;
      chem_file = start_of_list;

      while( chem_file ){
	
	/************************************************************/
	/* the array bulk_atom_label contains the cluster label     */
	/* assigned by first. Clusters are number consecutively, the*/
	/* largest cluster is 1, the next largest 2, ... . Which    */
	/* cluster an atom in the protein belongs to is indexed by  */
	/* bulk_atom_label. We only look at clusters of a certain   */
	/* size. The variable total_clusters indicates how many     */
	/* clusters we have in the current decomp.                  */
	/************************************************************/	  
	cluster_number = bulk_atom_label[chem_file->atom_number];

	if(cluster_number == cluster_labels[current_cluster] ){

	  nc_temp = (clusters *) malloc( sizeof( clusters ));
	 
	  if( new_clusters[current_cluster] == NULL ){
	    new_clusters[current_cluster] = nc_temp;
	    nc_current = nc_temp;
	  }
	  
	  nc_temp->atom_number     = chem_file->atom_number;
	  nc_temp->residue_number  = chem_file->residue;
	  nc_temp->insertion_space = chem_file->insertion_space;
	  sscanf( chem_file->atom_type, "%3s", nc_temp->atom_type );
	  nc_temp->chain_ID = chem_file->chain_ID;
	  nc_temp->cluster_color   = current_cluster+1;
	  
	  nc_current->next_element = nc_temp;
	  nc_current = nc_current->next_element;
	  nc_current->next_element = NULL;
	  
	}
	
	
	chem_file = chem_file->next_atom;
	
      }
    }
  }
  
}

/**********************************************************************/
/* store the atom info from the chem_file corresponding to the new    */
/* clusters in an array of linked lists called nc_count.              */
/**********************************************************************/
void set_new_decomp_info( int total_clusters, int bulk_atom_label[10000], 
			  int atom_num,  clusters *new_decomp[120], 
			  atom_list *chem_file, int cluster_labels[1000] ) {
  
  int cluster_number = 0,
      current_cluster = 0;
      //found_cluster = 0, // FIXME - unused variable
      //this_cluster = 0;  // FIXME - unused variable
  
  atom_list *start_of_list = NULL;
  
  clusters *nc_current = NULL,
           *nc_temp    = NULL;

  start_of_list = chem_file;

  for( current_cluster = 0; current_cluster < total_clusters; current_cluster++ ) {
    
    chem_file = start_of_list;
    
    while( chem_file ){

      cluster_number = bulk_atom_label[chem_file->atom_number];

      if( cluster_number == cluster_labels[current_cluster] ){
	nc_temp = (clusters *) malloc( sizeof( clusters ));

	if( new_decomp[current_cluster] == NULL ){
	  new_decomp[current_cluster] = nc_temp;
	  nc_current = nc_temp;
	}

	nc_temp->atom_number      = chem_file->atom_number;
	nc_temp->residue_number   = chem_file->residue;
	nc_temp->insertion_space  = chem_file->insertion_space;
	sscanf( chem_file->atom_type, "%3s", nc_temp->atom_type );
	nc_temp->chain_ID = chem_file->chain_ID;
	nc_temp->next_element = NULL;
	
	nc_current->next_element = nc_temp;
	nc_current = nc_current->next_element;
	nc_current->next_element = NULL;
      }
      
      chem_file = chem_file->next_atom;
      
    }
  }

  /* error check. make sure there are no empty cluster lists. */
  for( current_cluster = 0; current_cluster < total_clusters; current_cluster++ ) {
    if( new_decomp[current_cluster] == NULL ){
      printf(" Error: Found empty cluster list.\n");
      printf("        Cluster = %d\n", current_cluster);
      exit(1);
    }
  }

  //if( total_clusters > 32 ){
  //nc_temp = new_decomp[32];
  //printf("atom list for cluster 32: ");
  //while( nc_temp ){
  //  printf(" %d", nc_temp->atom_number );
  //  nc_temp = nc_temp->next_element;
  //}
  //printf("\n");
  //}
  
}
/**************************************************/

/**************************************************/
void swap_cluster_info( clusters *array[120], int element_A, int element_B, 
			int size[120] ){
  
  int temp_size[120];
  
  clusters *temp_struct[2];
  
  temp_size[0] = size[element_A];
  temp_struct[0] = array[element_A];

  size[element_A] = size[element_B];
  array[element_A] = array[element_B];

  size[element_B] = temp_size[0];
  array[element_B] = temp_struct[0];
  
}
/**************************************************/

/**************************************************/
void bsort( clusters *array[120], int size[120], int start, int end ){
  
  int a=0, last=0;
  
  if( start >= end )
    return;  
  swap_cluster_info( array, start, (start+end)/2, size );
  last = start;
  
  for( a = start+1; a <= end; a++ ){
    if( size[a] > size[start] ){
      swap_cluster_info( array, ++start, a, size );
    }
  }
  swap_cluster_info( array, start, last, size );
  bsort( array, size, start, last-1 );
  bsort( array, size, last+1, end );
}
/*****************************************************/

/*****************************************************/
void swap_colors( int array[120], int A, int B ){
  int temp;
  temp = array[A];
  array[A] = array[B];
  array[B] = temp;
}
/*****************************************************/

/*****************************************************/
void csort( int colors[120], int start, int end ) {

  int 
    a = 0, 
    last = 0;

  if( start >= end )
    return;
  swap_colors( colors, start, (start+end)/2 );
  last = start;
  for( a = start+1; a <= end; a++ ){
    if( colors[a] > colors[start] )
      swap_colors( colors, ++start, a );
  }
  swap_colors( colors, start, last );
  csort( colors, start, last-1 );
  csort( colors, last+1, end   );
}
/****************************************************/

/****************************************************/
/* This array is designed to compute the cluster    */
/* coloring, so as to keep it consistent through-   */
/* out the H-bond dilution plots. It an old cluster */
/* has more than one subset in the new decomposition*/
/* the largest of the new subsets is given the color*/
/* of the previous color, and a new color is assig- */
/* ned to the remaining subcluster(s)               */
/****************************************************/
void compute_new_cluster_list(clusters *nc_list[120], int num_new,
			      clusters *oc_list[120], int num_old,
			      int old_colors[120], int colors[120]){
  
  int 
    a = 0, 
    b = 0, 
    c = 0, 
    d = 0, 
    //e = 0,  // FIXME - unused variable
    total_count = 0,
    color_count = 1,
    is_a_subset = 0, 
    subset_index = 0, 
    subset_element = 0, 
    total_atoms[120], 
    current_color_index = 0,
    atom_count = 0, 
    sort_color[120],
    t = 0, 
    new_colors[120];
  
  clusters 
    *new_current,
    *old_current,
    *subset_list[120],
    *new_list[120];

  /**********************************************************************/
  /* Initialize variables.                                              */
  /**********************************************************************/
  for( d = 0; d < num_old; d++ )
    sort_color[d] = old_colors[d];
  
  for( d = 0; d < 120; d++)
    new_colors[d] = 99999;
  /**********************************************************************/

  csort( sort_color, 0, num_old-1 );
  current_color_index = sort_color[0];

  /***********************************************/
  /* First, find all new clusters which are sub- */
  /* sets of the first old cluster. If there are */
  /* more than one, sort by size and recolor. The*/
  /* new color will be 1 + the number of old     */
  /* clusters.                                   */
  /***********************************************/
  for( a = 0; a < num_old; a++ ){
    
    subset_index   = 0;
    subset_element = 0;

    /* Find all the new clusters that are subsets of old cluster "a" */
    for( b = 0; b < num_new; b++ ){
      
      is_a_subset = 0;
      new_current = nc_list[b];
      old_current = oc_list[a];
      
      while(old_current){
	
	if(new_current->residue_number == old_current->residue_number &&
	   !strcmp( new_current->atom_type, old_current->atom_type )  &&
	   ( new_current->chain_ID == old_current->chain_ID  ) &&
	   new_current->atom_number == old_current->atom_number ) {
	  is_a_subset = -1;
	  old_current = NULL;
	}
	if( old_current )
	  old_current = old_current->next_element;
      }

      if( is_a_subset ){
	subset_list[subset_index] = nc_list[b];
	new_current = nc_list[b];
	while(new_current){
	  atom_count++;
	  new_current = new_current->next_element;
	}
	total_atoms[subset_index] = atom_count;
	subset_index++;	
	atom_count = 0;
      }
    }

    /* if there is at least one subset of old cluster "a" */
    if( subset_index ){     
      //printf(" Cluster %d has a subset. first atom is: %d\n", a, subset_list[subset_index-1]->atom_number ); 
      
      bsort( subset_list, total_atoms, 0, subset_index-1 );

      for( t = 0; t < subset_index; t++ ){
	if( t == 0){
	  new_list[total_count] = subset_list[t];
	  new_colors[total_count] = old_colors[a];
	  total_count++;
	}
	else{    /* need to sort color_nums, and take 1+ the largest */
	  new_list[total_count] = subset_list[t];
	  new_colors[total_count] = current_color_index + color_count;
	  total_count++;
	  color_count++;
	}
      }
    }
    
  }

  if( total_count != num_new ){
    printf(" ERROR: Number of subclusters incorrectly calculated.\n");
    //exit(1);
  }

  //printf("num_new %d, num_old %d, total_count %d\n", num_new, num_old, total_count );
  for( c = 0; c < num_new; c++ ) {
    nc_list[c] = new_list[c];
    colors[c]  = new_colors[c];
  }

}
/**********************************************************************/

/**********************************************************************/
/* Tallies the number of main-chain atoms in new clusters and old     */
/* clusters. If the main-chain count did not change, the output of the*/
/* current decomp is not produced.                                    */
/**********************************************************************/
int compare_old_and_new_decomps( clusters *new_clst[120], clusters *old_clst[120], 
				 int cluster_counter, int cluster_sizes[120], 
				 int old_cluster_count ) {

  int 
    a = 0, 
    size_new = 0, 
    size_old = 0, 
    this_cluster_size = 0;
  
  clusters  
    *nc_current = NULL,
    *oc_current = NULL;
  

  for( a = 0; a < cluster_counter; a++ ){
    
    nc_current = new_clst[a];
    
    while( nc_current ){
      
      size_new++;
      /*      printf("size: %4d (%5d %4d %c)\n",size_new, nc_current->atom_number, nc_current->residue_number, nc_current->chain_ID );*/
      fflush(stdout);
      this_cluster_size++;
      nc_current = nc_current->next_element;
    }
    cluster_sizes[a] = this_cluster_size - 1;
    this_cluster_size = 0;
    //printf("new cluster sizes: Cluster %d = %d  (%d)\n", a, cluster_sizes[a], cluster_counter );
  }
  for( a = 0; a < old_cluster_count; a++ ){
    
    oc_current = old_clst[a];
    
    while( oc_current != NULL ){
      size_old++;
      oc_current = oc_current->next_element;
    }
   
  }
  
  /*printf("size new %d (%3d clusters) size old %d (%3d clusters)\n", size_new, cluster_counter, size_old, old_cluster_count );*/
  
  if( size_new == size_old ) {
    return(1);
  }
  return(0);
}

/**********************************************************************/
/* If there are more than 220 residues in the protein, the output is  */
/* displayed on multiple pages, in landscape format. In order to keep */
/* a clean postscript file, the output for the extra pages is sent to */
/* temporary files, one for each extra page needed. When the last     */
/* decomp of the page has finished, the tempfiles will be catted to   */
/* the output file cluster.ps, and deleted. The array file_list holds */
/* the addresses of the tempfiles.                                    */
/**********************************************************************/
void set_temp_files( FILE *file_list[5], int *number_of_pages, int *ps_page_number) {

  int a = 0;

  char 
    file_name[15];

  FILE *fp;
  
  for( a = 1; a <= *number_of_pages; a++ ) {
    printf("setting up temp files...\n");
    sprintf( file_name, "%d", a );
    fp = fopen( file_name,"w+");
    (*ps_page_number)++;
    fprintf(fp, "%%%%Page:     %d   %d\n\n", *ps_page_number, *ps_page_number );
    fprintf(fp, "0 109 620 109 sequence_line\n");
    file_list[a] = fp;
  }
}

/**********************************************************************/
/* Clean up unnecessary files and stuff.                              */
/**********************************************************************/    
void clean_up( FILE *ps_file, int number_of_pages ){

  int  a=0;
  char file_name[4], linebuf[300];
  FILE *current_file;

  fprintf( ps_file, "showpage\n\n" );
  fflush( ps_file );

  for( a = 1; a <= number_of_pages; a++ ) {
    sprintf( file_name, "%d", a );
    current_file = fopen( file_name, "r" );
    while( fgets( linebuf, sizeof(linebuf), current_file ) != NULL ){
      fprintf(ps_file, "%s", linebuf );
    }
    fprintf(ps_file, "showpage\n");
  }
  
  fprintf(ps_file, "%%%%EOF\n\n");
  
  if( number_of_pages )  
    system("\\rm -f 1 2 3 4 5");
}
