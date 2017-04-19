#ifndef _POSTSCRIPT_
#define _POSTSCRIPT_

#include "types.h"

void print_header( FILE *, float );
void print_data_headings( FILE *, int, int, int, int[]  );
void print_current_Hbond_info_portrait( FILE *, int, int, int, int *, int, float, 
					int, int [], int, int, int, float, int *, 
					int, int, int, int[], char [], char [], int,
					char, char, char[], char[], int, float, int );

void print_current_Hbond_info_landscape( FILE *, int, int, int, int *, int, float, 
					 int, int [], int, int, int, float, int *, 
					 int, int, int, int[], char [], char [], int,
					 char, char, char[], int[], char[], int, float, int );

void print_current_Hbond_info_landscape_multipage( FILE *, int, int, int, int *, int, float, 
						   int, int [], int, int, int, float, int *, 
						   int, int, int, int[], char [], char [],
						   FILE *[], int, int, char, char, char[], 
						   char[], int[], float );
void print_decomp( clusters *[], int, int, FILE *, int [], int[], int, int, int[], char[], int, 
		   int[], int[], char[], int [] );
void print_multipage_decomp( clusters *[], int, int, FILE *, int [], int[], int, FILE *[], int, 
			     int, int[], char[], int, int[], int[], char[] );
void print_portrait_footer( FILE *, char[], int );
void print_landscape_footer( FILE *, char[], int );

#endif
