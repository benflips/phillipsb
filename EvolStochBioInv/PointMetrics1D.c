/* Based on PointMetrics1D.c in the KBGrad functions, but modified to get pairwise matings...

Estimates density of individuals, mean trait values, and trait variances at each individual's location
Assumes a Gaussian kernel as the neighbourhood size.
Requires an array of 1D locations (X[i]), trait values (H[i], D[i])
bandwidth currently fixed at bw=1 (easy enough to make variable)
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
//#include <Rdefines.h>
#include <math.h>
#include "Rinterface.h"

#define Pi 3.141593


//--------------------------------------//
// Function declarations //

SEXP mate (SEXP R_X, SEXP R_dens, SEXP R_n, SEXP R_bw);
double norm (double x, double bw);
SEXP metrics (SEXP R_X, SEXP R_H, SEXP D_H, SEXP R_n, SEXP R_ev);
SEXP pmetrics (SEXP R_X, SEXP R_H, SEXP R_D, SEXP R_n, SEXP R_ev, SEXP R_ncores);
SEXP sum_metrics (SEXP R_X, SEXP R_H, SEXP R_D, SEXP R_n, SEXP R_bins, SEXP R_nbins, SEXP R_ev, SEXP R_bw);


//--------------------------------------//
//  Function definitions  //

// Given a popmatrix, finds a mate for each individual
SEXP mate (SEXP R_X, SEXP R_dens, SEXP R_n, SEXP R_bw){
	R_X=coerceVector(R_X, REALSXP);
	//R_dens=coerceVector(R_dens, REALSXP);
	int bw = INTEGER(coerceVector(R_bw, INTSXP))[0];
	int n = INTEGER(coerceVector(R_n, INTSXP))[0];
	
	SEXP mate; PROTECT(mate=allocVector(INTSXP, n));
	
	const double small_num = 1e-8;
	double *X, *D;
	int *M;
	X = REAL(R_X);
	M = INTEGER(mate);
	//D=REAL(R_dens);
	double w, rnum, Dtemp;
	
	GetRNGstate();
	for (int ii=0; ii<n; ii++){ // find a mate for each individual
		Dtemp = 0; //D[ii]-norm(0, bw);	
		for (int jj=0; jj<n; jj++){
			if (jj==ii) continue;
			Dtemp += norm(X[ii]-X[jj], bw);
		}
		if (Dtemp==0) {M[ii]=-99; continue;} // catch the lonely ones
	
		rnum = unif_rand();//Random number between 0-1
		w=0;
		int jj = 0;
		do
		{  //finds the R vector index of the individual's mate
			if (jj==ii) {jj++; continue;}
			w += norm(X[ii]-X[jj], bw)/Dtemp; //cumulant, should add to one eventually
			//printf("ind=%i rnum=%e cumulant=%e lessthan=%i catch=%i\n", ii, rnum, w, w<rnum, (rnum-w)>small_num);
			jj++;
		} while (w<rnum && (rnum-w)>small_num ); 

		M[ii] = jj;	
	}
	PutRNGstate();
	
	UNPROTECT(1);
	return(mate);
}

double norm (double x, double bw){
	return exp(-pow(x,2)/(2*pow(bw,2)))/(bw*sqrt(2*Pi));
}

// Calculates population metrics back to individual locations, 1D
SEXP metrics (SEXP R_X, SEXP R_H, SEXP R_D, SEXP R_n, SEXP R_ev){
	int n = INTEGER(coerceVector(R_n, INTSXP))[0]; //grab vector length
	int ev = INTEGER(coerceVector(R_ev, INTSXP))[0];
	R_X=coerceVector(R_X, REALSXP); //digest R_X
	R_H=coerceVector(R_H, REALSXP); //digest R_H
	R_D=coerceVector(R_D, REALSXP); //digest R_D
	SEXP outmat; PROTECT(outmat=allocMatrix(REALSXP, n, 5)); // a place to put all the output
	
	double *X, *RH, *RD, *out;  //pointer variables
	double w; //to take weights
	
	X = REAL(R_X); //pointers to real parts of R vectors
	RH = REAL(R_H);
	RD = REAL(R_D);
	out = REAL(outmat);
	
	//calculate density and mean trait values
	for (int ii=0; ii<n; ii++){
		out[ii] = 0; //column 1 for density
		out[ii+n] = 0; //column 2 for meanD
		out[ii+2*n] = 0; //column 3 for meanH
		for (int jj=0; jj<n; jj++){
			w = norm(X[ii]-X[jj], 1);
			out[ii] += w;
			out[ii+n] += w*(RD[jj]);
			out[ii+2*n] += w*(RH[jj]);
		}
		out[ii+n] = out[ii+n]/out[ii];
		out[ii+2*n] = out[ii+2*n]/out[ii];
	}
	
	//calculate trait variances
	if (ev==1){
		for (int ii=0; ii<n; ii++){
			out[ii+3*n] = 0;
			out[ii+4*n] = 0;
			for (int jj=0; jj<n; jj++){
				w = norm(X[ii]-X[jj], 1);
				out[ii+3*n] += w*pow(out[ii+n]-RD[jj], 2);
				out[ii+4*n] += w*pow(out[ii+2*n]-RH[jj], 2);
			}
			out[ii+3*n] = sqrt(out[ii+3*n]/out[ii]);
			out[ii+4*n] = sqrt(out[ii+4*n]/out[ii]);
		}
	}
	
	UNPROTECT(1);
	return(outmat);
}



// A serial version of the metrics function that calculates metrics back to fixed points rather than individuals
SEXP sum_metrics (SEXP R_X, SEXP R_H, SEXP R_D, SEXP R_n, SEXP R_bins, SEXP R_nbins, SEXP R_ev, SEXP R_bw){
	int n = INTEGER(coerceVector(R_n, INTSXP))[0]; //grab vector length
	int nbins = INTEGER(coerceVector(R_nbins, INTSXP))[0]; //grab bin length
	int ev = INTEGER(coerceVector(R_ev, INTSXP))[0];
	int bw = INTEGER(coerceVector(R_bw, INTSXP))[0];
	
	R_X=coerceVector(R_X, REALSXP); //digest R_X
	R_H=coerceVector(R_H, REALSXP); //digest R_H
	R_D=coerceVector(R_D, REALSXP); //digest R_D
	R_bins=coerceVector(R_bins, INTSXP);
	SEXP outmat; PROTECT(outmat=allocMatrix(REALSXP, nbins, 5)); // a place to put all the output
	
	double *X, *RH, *RD, *out;  //pointers variables
	double w; //to take weights
    int *Rb, ii, jj; // iterations variables
	
	X = REAL(R_X); //pointers to real parts of R vectors
	RH = REAL(R_H);
	RD = REAL(R_D);
	out = REAL(outmat);
	Rb = INTEGER(R_bins);
		
	//calculate density and mean trait values
	for (ii=0; ii<nbins; ii++){
		out[ii] = 0; //column 1 for density
		out[ii+nbins] = 0; //column 2 for meanD
		out[ii+2*nbins] = 0; //column 3 for meanH
		for (jj=0; jj<n; jj++){
			w = norm(Rb[ii]-X[jj], bw);
			out[ii] += w;
			out[ii+nbins] += w*(RD[jj]);
			out[ii+2*nbins] += w*(RH[jj]);
		}
		out[ii+nbins] = out[ii+nbins]/out[ii];
		out[ii+2*nbins] = out[ii+2*nbins]/out[ii];
	}
	
	//calculate trait variances
	if (ev==1){
			for (ii=0; ii<nbins; ii++){
				out[ii+3*nbins] = 0;
				out[ii+4*nbins] = 0;
				for (jj=0; jj<n; jj++){
					w = norm(Rb[ii]-X[jj], bw);
					out[ii+3*nbins] += w*pow(out[ii+nbins]-RD[jj], 2);
					out[ii+4*nbins] += w*pow(out[ii+2*nbins]-RH[jj], 2);
				}
				out[ii+3*nbins] = sqrt(out[ii+3*nbins]/out[ii]);
				out[ii+4*nbins] = sqrt(out[ii+4*nbins]/out[ii]);
			}
	}
	
	UNPROTECT(1);
	return(outmat);
}


// Calculates population metrics back to individual locations, 2D
SEXP metrics2D (SEXP R_X, SEXP R_Y, SEXP R_H, SEXP R_D, SEXP R_n, SEXP R_ev){
	int n = INTEGER(coerceVector(R_n, INTSXP))[0]; //grab vector length
	int ev = INTEGER(coerceVector(R_ev, INTSXP))[0];
	R_X=coerceVector(R_X, REALSXP); //digest R_X
	R_Y=coerceVector(R_Y, REALSXP); //digest R_Y
	R_H=coerceVector(R_H, REALSXP); //digest R_H
	R_D=coerceVector(R_D, REALSXP); //digest R_D
	SEXP outmat; PROTECT(outmat=allocMatrix(REALSXP, n, 5)); // a place to put all the output
	
	double *X, *Y, *RH, *RD, *out;  //pointer variables
	double w; //to take weights
	
	X = REAL(R_X); //pointers to real parts of R vectors
	Y = REAL(R_Y);
	RH = REAL(R_H);
	RD = REAL(R_D);
	out = REAL(outmat);
	
	//calculate density and mean trait values
	for (int ii=0; ii<n; ii++){
		out[ii] = 0; //column 1 for density
		out[ii+n] = 0; //column 2 for meanD
		out[ii+2*n] = 0; //column 3 for meanH
		for (int jj=0; jj<n; jj++){
			w = norm(pow((pow(X[ii]-X[jj],2)+pow(Y[ii]-Y[jj],2)),0.5), 1);
			out[ii] += w;
			out[ii+n] += w*(RD[jj]);
			out[ii+2*n] += w*(RH[jj]);
		}
		out[ii+n] = out[ii+n]/out[ii];
		out[ii+2*n] = out[ii+2*n]/out[ii];
	}
	
	//calculate trait variances
	if (ev==1){
		for (int ii=0; ii<n; ii++){
			out[ii+3*n] = 0;
			out[ii+4*n] = 0;
			for (int jj=0; jj<n; jj++){
				w = norm(pow((pow(X[ii]-X[jj],2)+pow(Y[ii]-Y[jj],2)),0.5), 1);
				out[ii+3*n] += w*pow(out[ii+n]-RD[jj], 2);
				out[ii+4*n] += w*pow(out[ii+2*n]-RH[jj], 2);
			}
			out[ii+3*n] = sqrt(out[ii+3*n]/out[ii]);
			out[ii+4*n] = sqrt(out[ii+4*n]/out[ii]);
		}
	}
	
	UNPROTECT(1);
	return(outmat);
}

// Given a popmatrix, finds a mate for each individual, 2D case
SEXP mate2D (SEXP R_X, SEXP R_Y, SEXP R_dens, SEXP R_n, SEXP R_bw){
	R_X=coerceVector(R_X, REALSXP);
	R_Y=coerceVector(R_Y, REALSXP);
	//R_dens=coerceVector(R_dens, REALSXP);
	int bw = INTEGER(coerceVector(R_bw, INTSXP))[0];
	int n = INTEGER(coerceVector(R_n, INTSXP))[0];
	
	SEXP mate; PROTECT(mate=allocVector(INTSXP, n));
	
	const double small_num = 1e-8;
	double *X, *Y, *D;
	int *M;
	X = REAL(R_X);
	Y = REAL(R_Y);
	M = INTEGER(mate);
	//D=REAL(R_dens);
	double w, rnum, Dtemp;
	
	GetRNGstate();
	for (int ii=0; ii<n; ii++){ // find a mate for each individual
		Dtemp = 0; //D[ii]-norm(0, bw);	
		for (int jj=0; jj<n; jj++){
			if (jj==ii) continue;
			Dtemp += norm(pow((pow(X[ii]-X[jj],2)+pow(Y[ii]-Y[jj],2)),0.5), bw);
		}
		if (Dtemp==0) {M[ii]=-99; continue;} // catch the lonely ones
	
		rnum = unif_rand();//Random number between 0-1
		w=0;
		int jj = 0;
		do
		{  //finds the R vector index of the individual's mate
			if (jj==ii) {jj++; continue;}
			w += norm(pow((pow(X[ii]-X[jj],2)+pow(Y[ii]-Y[jj],2)),0.5), bw)/Dtemp; //cumulant, should add to one eventually
			//printf("ind=%i rnum=%e cumulant=%e lessthan=%i catch=%i\n", ii, rnum, w, w<rnum, (rnum-w)>small_num);
			jj++;
		} while (w<rnum && (rnum-w)>small_num ); 

		M[ii] = jj;	
	}
	PutRNGstate();
	
	UNPROTECT(1);
	return(mate);
}


