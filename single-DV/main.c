/*************************************************************************
main.c:  Change point detection of diffusion and velocity (non-parallelized)
         Originally written using OpenMPI v1.4.2
         Arguments are filename, timestep, alpha

filename:   File with Time Series Data
alpha:      Type-I error
*************************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <gsl/gsl_version.h>
#include "changepoint.h"
#include "critical_values.h"
#include <math.h>

#define TOL 1.0e-10

int main(int argc, char *argv[])
{
	FILE *fpin, *fpout;    // File pointers
	char in_name[255],  filename[255], *endptr[255], out_name[255];	// Store strings
	double alpha = 0.075;  // Type-I error, mis-specify transition
	double delta_t = 1;    // Time unit between measurements
	int L = 0;             // Total number of data points
	int NA_BASE = 10000;   // Maximum length that can be processed before recursion
	int NA_OVERLAP = 200;  // Overlap length
	int ui;    	       // Dummy index
	int cpl, cpr, cp1;     // FindCP left start, right end, and return value
	int Ncp = 0, Ncpdlt=0;   // Change points found, deleted
	double tmp;            // Temporary variable for input
	double *traj;          // Time series
	struct changepoint *cp_root=NULL; // Change point binary tree
	clock_t time0, time1;  // Variables for timing
	int h,i,j,k;           // Dummy indices
	int trace=0;           // Debug flag
	enum bool done = false; // Boolean for has a change point been added to tree

	int** cpArray = (int**) malloc(sizeof(int*)); // Change points found
	int** lbArray = (int**) malloc(sizeof(int*)); // Left bounds of confidence interval
	int** rbArray = (int**) malloc(sizeof(int*)); // Right bounds of confidence interval

	// Get command line arguments
	if (argc != 4) {
		// No file name entered in the command line
		printf( "\nchangepoint %s%s build %s (GSL version %s)\n", CHANGEPOINT_VERSION, PLATFORM, COMPILE_DATE, GSL_VERSION);
		printf( "Syntax : changepoint filename delta_t alpha beta\n");
		printf( "Example: changepoint myfile 0.001 0.05 0.95\n");
		printf( "         produces myfile.cp with type-I error (alpha) of 5%% and\n");
		printf( "         confidence interval of 0.95 \n" );
		printf( "BUG    : Please send emails to hawyang-at-princeton.edu\n\n");
		exit(1);
	}
	else if (argc == 4) {
		strcpy(in_name, argv[1]);
		strcpy(out_name, argv[1]);
		delta_t = strtod(argv[2], endptr);
		alpha = strtod(argv[3], endptr);
	}
	// Print title
	printf( "\n" );
	printf( "********************************************************\n");
	printf( "					Change Point Detection			   	 \n");
	printf( "********************************************************\n");
	printf( "Version: %s%s build %s (GSL version %s)\n\n", 
		   CHANGEPOINT_VERSION, PLATFORM, COMPILE_DATE, GSL_VERSION );					
	
	// Read in data
	if ( (fpin=fopen(in_name,"r")) == NULL ) {
			// File does not exist, exit with error
			printf("File [%s] does not exist. Exiting ...\n", in_name );
			exit(1);
	}
	else {
		// File exists, start reading in data
		printf("] Reading in data file ...\n");
		traj = (double*) malloc(sizeof(double*));
		traj[0] = 0.0;
		while (!feof(fpin) ) {
			fscanf( fpin, "%le", &tmp );
			L++;
			traj = (double *) realloc(traj, (L+1)*sizeof(double));
			traj[L] = tmp;
		}
		L--;
		fclose( fpin );
		printf("  Total of %lu data points read in.\n", (long unsigned int) L);
		// Find change points recursively
		printf("] Finding change points recursively ...\n"); 
		time0 = clock();
		if (NA_BASE > L )
			cp1 = FindCP(&cp_root, traj, delta_t, 0, L+1, alpha, &Ncp, L, 1);
		else {
			// If length of trajectory exceeds bounds of critical value
			cpl = 0;
			cpr = NA_BASE;
			done = false;
			while (!done ) {
				cp1 = FindCP(&cp_root, traj, delta_t, cpl, cpr, alpha, &Ncp, L, 1);
				if (cpr < L) {
					// No change points found
					if (cp1 == 0 ) cpl = cpr - NA_OVERLAP;
					else cpl = cp1;
					cpr = cpl + NA_BASE - 1;
					if (cpr > L) cpr = L+1;
					done = false;
				}
				else 
					done = true;
			}
		}
		time1 = clock();
		printf("  %i change points found. [%.0f s]\n", Ncp, (double)(time1-time0)/CLOCKS_PER_SEC);
		// Remove spurious change points
		printf("] Examining %i change points sequentially ...\n",Ncp);
		time0 = clock();
		if (Ncp > 0)
			CheckCP(&cp_root, traj, delta_t, alpha, &Ncpdlt, trace, L);
		time1 = clock();
		printf("  %i spurious change points removed. [%.0f s]\n", Ncpdlt, 
			(double)(time1-time0)/CLOCKS_PER_SEC);
		Ncp = 0;
		MakeCPArray(cp_root, &cpArray[0], &lbArray[0], &rbArray[0], &Ncp);
		// Save results...
		strcpy(filename, out_name);
		strcat(filename, ".cp" );
		printf("  saving file: %s \n", filename);
		fpout = fopen( filename, "w");
		for (i = 1; i <= Ncp; i++)
			fprintf(fpout, "%d %d %d\n", cpArray[0][i], lbArray[0][i], rbArray[0][i]);
		fclose(fpout);
		
		free(traj);
		free(cpArray[0]);
		free(lbArray[0]);
		free(rbArray[0]);
		free(cpArray);
		free(lbArray);
		free(rbArray);
	}
}
