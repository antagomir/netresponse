/*

  This file is a part of the NetResponse R package.

  Copyright (C) 2008-2010 Leo Lahti, Olli-Pekka Huovilainen and
  António Gusmão. Contact: Leo Lahti <leo.lahti@iki.fi>

  This file is based on the Agglomerative Independent Variable Group
  Analysis package, Copyright (C) 2001-2007 Esa Alhoniemi, Antti
  Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and Paul
  Wagner.
 
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
 
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

#define POW2(x) ((x) * (x))

void compute_nc(int ncentroids, long datalen, double *trueNc,
		double *qOFz, double *Nc)
{
  register int i, j;

  /* Go through the clusters */
  for (i=0; i<ncentroids; i++) {
    trueNc[i] = 0.0; /* initialize cluster size with 0 */
    for (j=0; j<datalen; j++) {
      trueNc[i] += qOFz[i*datalen + j]; /* sum up the individual memberships */
    }
    Nc[i] = trueNc[i];
  }
  /* treat separately the last cluster: ensure it's empty */
  Nc[ncentroids-1] = 0.0;
  for (j=0; j<datalen; j++) {
    qOFz[(ncentroids-1)*datalen+j] = 0;
  }

  return;
}


void
update_centroids(long datalen, int ncentroids, int dim1, int dim2,
		 double *data1, int **data2_int,
                 double *Nc, double *qOFz, double *Mumu, double *S2mu,
                 double *Mubar, double *Mutilde,
		 double *KsiAlpha, double *KsiBeta, double *AlphaKsi,
		 double *BetaKsi, double implicitnoisevar,
		 double *U_p, double ***U_hat_table, double *Ns) {
  register int i, k;
  register long ind, j, t;
  double term, term2, term3, s2x, s2x_new;
  
  
  if (dim1){
    for (k = 0; k < dim1; k++) {
      s2x = BetaKsi[k] / AlphaKsi[k];
      for (i = 0; i < ncentroids; i++) {
	term = 0.0;
	ind  = k * ncentroids + i;
	for (j = 0; j < datalen; j++)
	  term += qOFz[i*datalen + j] * data1[k * datalen + j];
	term2         = s2x + S2mu[k] * Nc[i];
	Mubar[ind]   = ((s2x * Mumu[k]) + (S2mu[k] * term))/term2;
	Mutilde[ind] = (s2x * S2mu[k])/term2;
	KsiAlpha[ind] = AlphaKsi[k] + 0.5 * Nc[i];
	term3 = 0.0;
	for (j = 0; j < datalen; j++)
	  term3 += qOFz[i*datalen + j] *
	    (Mutilde[ind] + POW2(data1[k * datalen + j] - Mubar[ind]) +
	     implicitnoisevar);
	KsiBeta[ind] = BetaKsi[k] + 0.5 * term3;

	s2x_new = KsiBeta[ind] / KsiAlpha[ind];
	term2         = s2x_new + S2mu[k] * Nc[i];
	Mubar[ind]   = ((s2x_new * Mumu[k]) + (S2mu[k] * term))/term2;
	Mutilde[ind] = (s2x_new * S2mu[k])/term2;
      }
    }
  }

  for(j=0;j<dim2;j++){
    for(i=0;i<ncentroids;i++){   
      for(k=0;k<(int)(Ns[j]);k++)
      	U_hat_table[j][i][k]=U_p[j];
      for(t=0;t<datalen;t++){	     /***************/
	U_hat_table[j][i][data2_int[j][t]] += qOFz[i*datalen+t];
      }
    }
  }  
  return;
}


void update_gamma(int ncentroids, double *trueNc, double *prior_alpha,
		  double *post_gamma)
{
  register int i;
  double ncsum, nccumsum;

  ncsum = 0.0;
  for (i=0; i<ncentroids; i++) {
    ncsum += trueNc[i];
  }
  nccumsum = 0.0;
  for (i=0; i<ncentroids; i++) {
    nccumsum += trueNc[i];
    post_gamma[2*i] = 1 + trueNc[i];
    post_gamma[2*i+1] = *prior_alpha + ncsum - nccumsum;
  }

  return;
}



void
allocate_memory_A(long datalen,int ncentroids,int dim2,
                double ****U_hat_table,int ***data2_int,double *Ns) {
  register int i,j;
  
  if (dim2) {
    *U_hat_table=(double ***)malloc(dim2 * sizeof(double*));
    *data2_int=(int **)malloc(dim2 * sizeof(int*));
  }
  for (j=0;j<dim2;j++){
    (*data2_int)[j]  = (int *)malloc(datalen * sizeof(int));
    (*U_hat_table)[j]=(double**)malloc(ncentroids * sizeof(double*));
    for (i=0;i<ncentroids;i++) {
      (*U_hat_table)[j][i] =(double *)malloc(((int)(Ns[j]))*sizeof(double));
    }
  }
    
  return;
}
 
/************************************************************/

void
free_memory_A(int ncentroids, int dim2, 
	    double ****U_hat_table, int ***data2_int) {
  register int i,j;
  for (j=0;j<dim2;j++){
    for (i=0;i<ncentroids;i++) {
      free((*U_hat_table)[j][i]);
    }
    free((*data2_int)[j]);
    free((*U_hat_table)[j]);
  }

  if (dim2) {
    free(*U_hat_table);
    free(*data2_int);
  }
  return;
}


void
vdp_mk_HPposterior(double *Mumu, double *S2mu, double *Mubar, double *Mutilde, 
		    double *AlphaKsi, double *BetaKsi, 
		    double *KsiAlpha, double *KsiBeta, 
		    double *post_gamma, double *prior_alpha,
		    double *U_p, SEXP *U_hat,
		    long datalen, int dim1, int dim2, double *data1, double *data2, 
		    double *Ns, int ncentroids, 
		    double implicitnoisevar, double *qOFz,
		    double *Nc, double *trueNc) {
  register long i, j, t;
  register int k;
  double  *U_hat_j;
  SEXP U_hat_j_SEXP;
  double       ***U_hat_table;
  int          **data2_int;

  allocate_memory_A(datalen, ncentroids, dim2,
		  &U_hat_table,&data2_int,Ns );

  for (j=0;j<dim2;j++){
    for(t=0;t<datalen;t++)
      data2_int[j][t]=((int)(data2[j*datalen+t]))-1;
  }

  compute_nc(ncentroids, datalen, trueNc, qOFz, Nc);

  update_centroids(datalen, ncentroids, dim1, dim2,
		   data1, data2_int,
		   Nc, qOFz, Mumu, S2mu,
		   Mubar, Mutilde, 
		   KsiAlpha, KsiBeta, AlphaKsi,
		   BetaKsi, implicitnoisevar,
		   U_p, U_hat_table, Ns);

  update_gamma(ncentroids, trueNc, prior_alpha, post_gamma);

  for (j=0;j<dim2;j++){
    PROTECT(U_hat_j_SEXP = NEW_NUMERIC(ncentroids*Ns[j]));
    U_hat_j=NUMERIC_POINTER(U_hat_j_SEXP);
    SET_ELEMENT(*U_hat, j, U_hat_j_SEXP); // U_hat has a list of size dim2
    for (i=0;i<ncentroids;i++)
      for (k=0;k<(int)(Ns[j]);k++)
	      U_hat_j[k*ncentroids+i]=U_hat_table[j][i][k];
  }

  free_memory_A(ncentroids, dim2, &U_hat_table, &data2_int);

  return;
}


/************************************************************/
/* bridge function                                          */

SEXP mHPpost(SEXP X1, SEXP X1_Columns, SEXP X1_Rows, SEXP X2, SEXP X2_Columns,
        SEXP realS, SEXP OPTSimplicitnoisevar,
        SEXP HP_PRIOR_Mumu, SEXP HP_PRIOR_S2mu,
	SEXP HP_PRIOR_AlphaKsi, SEXP HP_PRIOR_BetaKsi,
	SEXP HP_PRIOR_U_p, SEXP HP_PRIOR_prior_alpha,
	SEXP QOFZ, SEXP QOFZ_Columns)
{
  long datalen;
  int  i, dim1, dim2, ncentroids;
  double *Mumu, *S2mu, *Mubar, *Mutilde, 
    *AlphaKsi, *BetaKsi, *KsiAlpha, *KsiBeta, *U_p, *prior_alpha,
    *post_gamma;
  double *data1;
  double *data2;
  double *Ns;
  double *qOFz_in, *qOFz, *Nc, *trueNc;
  double implicitnoisevar;
  const char *posterior_fields[]={"Mubar","Mutilde",
				  "KsiAlpha","KsiBeta",
				  "gamma","Nc","trueNc","qOFz","Uhat"};
  const char *prior_fields[]={"Mumu","S2mu",
			      "AlphaKsi","BetaKsi",
			      "alpha","U_p"};
  SEXP list, list_names;
  SEXP oMubar, oMutilde,  oKsiAlpha, oKsiBeta, opost_gamma,
       oNc    , otrueNc ,  oqOFz   , oU_hat;
  SEXP* U_hat;

  /************ CONVERTED input variables ******************/
  PROTECT(X1 = AS_NUMERIC(X1));  
  data1   = NUMERIC_POINTER(X1);
  dim1    = INTEGER_VALUE(X1_Columns);
  datalen = INTEGER_VALUE(X1_Rows);

  PROTECT(X2 = AS_NUMERIC(X2));  
  data2   = NUMERIC_POINTER(X2);
  dim2    = INTEGER_VALUE(X2_Columns);

  Ns = NUMERIC_POINTER(realS);
  implicitnoisevar = NUMERIC_VALUE(OPTSimplicitnoisevar);

  /************ CONVERTED initial values of model parameters ******************/
  Mumu        = NUMERIC_POINTER(HP_PRIOR_Mumu);
  S2mu        = NUMERIC_POINTER(HP_PRIOR_S2mu);
  AlphaKsi    = NUMERIC_POINTER(HP_PRIOR_AlphaKsi);
  BetaKsi     = NUMERIC_POINTER(HP_PRIOR_BetaKsi);
  U_p         = NUMERIC_POINTER(HP_PRIOR_U_p);
  prior_alpha = NUMERIC_POINTER(HP_PRIOR_prior_alpha);
  qOFz_in     = NUMERIC_POINTER(QOFZ);
  ncentroids  = INTEGER_VALUE(QOFZ_Columns);

  /********* CONVERTED output variables ***********************/
  /* Necessary to allocate memory for the output variables :| */
  /*Allocation **/

  //printf("b0");
  PROTECT(oMubar      = NEW_NUMERIC(ncentroids*dim1));

  //printf("b1");
  PROTECT(oMutilde    = NEW_NUMERIC(ncentroids*dim1));
  PROTECT(oKsiAlpha   = NEW_NUMERIC(ncentroids*dim1));
  PROTECT(oKsiBeta    = NEW_NUMERIC(ncentroids*dim1));
  PROTECT(opost_gamma = NEW_NUMERIC(2*ncentroids));
  PROTECT(oNc         = NEW_NUMERIC(1*ncentroids));
  PROTECT(otrueNc     = NEW_NUMERIC(1*ncentroids));
  PROTECT(oqOFz       = NEW_NUMERIC(datalen*ncentroids));
  PROTECT(oU_hat      = NEW_LIST(dim2)); /* CHECK This should be a Cell Matrix??? */

  Mubar      = NUMERIC_POINTER(oMubar);
  Mutilde    = NUMERIC_POINTER(oMutilde);
  KsiAlpha   = NUMERIC_POINTER(oKsiAlpha);
  KsiBeta    = NUMERIC_POINTER(oKsiBeta);
  post_gamma = NUMERIC_POINTER(opost_gamma);
  Nc         = NUMERIC_POINTER(oNc);
  trueNc     = NUMERIC_POINTER(otrueNc);
  qOFz       = NUMERIC_POINTER(oqOFz);
  U_hat      = &oU_hat;

  for (i=0; i<datalen*ncentroids; i++) {
    qOFz[i] = qOFz_in[i];
  }

  vdp_mk_HPposterior(Mumu, S2mu, Mubar, Mutilde, 
		      AlphaKsi, BetaKsi, KsiAlpha, KsiBeta, 
		      post_gamma, prior_alpha,
		      U_p, U_hat,
		      datalen, dim1, dim2, data1, data2, 
		      Ns, ncentroids, implicitnoisevar, qOFz, Nc, trueNc);
  
  /****************** CREATE A LIST WITH THE OUTPUT *****************/
  // Creating a character string vector 
  // of the "names" attribute of the
  // objects in out list:

  PROTECT(list_names = NEW_CHARACTER(9));    

  for(i = 0; i < 9; i++)   
    SET_STRING_ELT(list_names,i,mkChar(posterior_fields[i])); 

  // Creating a list with 9 vector elements:
  PROTECT(list = NEW_LIST(9)); 

  // attaching elements to the list:
  SET_ELEMENT(list, 0, oMubar);
  SET_ELEMENT(list, 1, oMutilde);
  SET_ELEMENT(list, 2, oKsiAlpha);
  SET_ELEMENT(list, 3, oKsiBeta);
  SET_ELEMENT(list, 4, opost_gamma);
  SET_ELEMENT(list, 5, oNc);
  SET_ELEMENT(list, 6, otrueNc);
  SET_ELEMENT(list, 7, oqOFz);
  SET_ELEMENT(list, 8, oU_hat);

  // and attaching the vector names:
  SET_NAMES(list, list_names);

  UNPROTECT(2+9+2+dim2);
  //CHECK UPDATE NUMBER OF UNPROTECTS : UNPROTECT(4);
  return list;
}
 
/************************************************************/
 
#define DIGAMMA_S 1e-5
#define DIGAMMA_C 8.5
#define DIGAMMA_S3 1.0/12
#define DIGAMMA_S4 1.0/120
#define DIGAMMA_S5 1.0/252
#define DIGAMMA_D1 -0.5772156649

#define POW2(x) ((x) * (x))


/************************************************************/
/* Read Elements of a R list in C.                         */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

/************************************************************/
/* digamma function (by Antti Honkela)                      */

double 
digamma(double x) {
  double y = 0.0, r = 0.0, xn = x;

  if (xn <= 0) {
    return R_NaN;
  }
  
  if (xn <= DIGAMMA_S)
    y = DIGAMMA_D1 - 1.0 / xn;
  else {
    while (xn < DIGAMMA_C) {
      y  -= 1.0 / xn;
      xn += 1.0;
    }
    r = 1.0 / xn;
    y += log(xn) - .5 * r;
    r = POW2(r);
    y -= r * (DIGAMMA_S3 - r * (DIGAMMA_S4 - r*DIGAMMA_S5));
  }
 
  return y;
}                                                                              


/************************************************************/

void
compute_variance(int ncentroids, int dim1, double *KsiAlpha, 
                 double *KsiBeta, double **S2_x, double **Ksi_log) {
  register int i, ind, k;
  
  for (i = 0; i < ncentroids; i++) 
    for (k = 0; k < dim1; k++) {
      ind = k * ncentroids + i;
      S2_x[i][k]    = KsiBeta[ind]/KsiAlpha[ind];
      Ksi_log[i][k] = digamma(KsiAlpha[ind])-log(KsiBeta[ind]);
      
      if( S2_x[i][k] < 1e-100 ) S2_x[i][k] = 1e-100;
    }
  
  return;
}

/************************************************************/
void
compute_tempmat(long datalen, int dim1, int dim2, int ncentroids,
		double **Temp, double *data1, int **data2_int,
		double *Mubar, double *Mutilde, double **S2_x,
                double **Ksi_log, double ***U_hat_table, double *Ns,
		double implicitnoisevar, double *log_lambda) {
  register int i, k;
  long         ind, j,t;
  double term;
  
  for (i = 0; i < ncentroids; i++) {
    for (j = 0; j < datalen; j++) {
      Temp[i][j] = 0.0;
      for (k = 0; k < dim1; k++) {
	ind  = k * ncentroids + i;
	Temp[i][j] += ((Mutilde[ind]+POW2(data1[k*datalen + j]-Mubar[ind]) + implicitnoisevar)/
		       S2_x[i][k]) - Ksi_log[i][k];
      }
      Temp[i][j] /= 2.0;
    }
  }
  for(j=0;j<dim2;j++){
    for(i=0;i<ncentroids;i++){
      term=0.0;
      for(k=0;k<(int)(Ns[j]);k++){
	term += U_hat_table[j][i][k]; 
	U_hat_table[j][i][k]=digamma(U_hat_table[j][i][k]);
      }
      term=digamma(term);
      for (t=0;t<datalen;t++){
	Temp[i][t] += (term - U_hat_table[j][i][data2_int[j][t]]);
      }
    }
  }

  for (i = 0; i < ncentroids; i++) {
    for (j = 0; j < datalen; j++) {
      log_lambda[i * datalen + j] += -dim1*log(2*M_PI)/2 - Temp[i][j];
    }
  }
  return;
}


void log_p_of_z_given_other_z_c(int datalen, long ncentroids,
				double *post_gamma, double *log_lambda)
{
  register int c, i;
  double E_log_p;

  for (c=0; c<ncentroids; c++) {
    E_log_p = digamma(post_gamma[2*c]) - digamma(post_gamma[2*c] + post_gamma[2*c+1]);
    for (i=0; i<c; i++) {
      E_log_p += digamma(post_gamma[2*i+1]) - digamma(post_gamma[2*i] + post_gamma[2*i+1]);
    }
    for (i=0; i<datalen; i++) {
      log_lambda[c*datalen+i] = E_log_p;
    }
  }

  return;
}



void fix_lambda(int ncentroids, long datalen, double *prior_alpha, double *log_lambda)
{
  register int i;
  double correction;

  correction = log(1 - exp(digamma(*prior_alpha) - digamma(1 + *prior_alpha)));
  for (i=0; i<datalen; i++) {
    log_lambda[(ncentroids-1)*datalen + i] -= correction;
  }

  return;
}



void
allocate_memory_B(long datalen,int ncentroids,int dim1,int dim2,double ***S2_x,
                double ***Ksi_log, double ***Temp,
		double ****U_hat_table,int ***data2_int,double *Ns) {
  register int i,j;
  
  *Temp    = (double **)malloc(ncentroids * sizeof(double*));
  if (dim1) {
    *S2_x    = (double **)malloc(ncentroids * sizeof(double*));
    *Ksi_log = (double **)malloc(ncentroids * sizeof(double*));
  }
  if (dim2) {
    *U_hat_table=(double ***)malloc(dim2 * sizeof(double*));
    *data2_int=(int **)malloc(dim2 * sizeof(int*));
  }
  for (i = 0; i < ncentroids; i++) {
    (*Temp)[i]    = (double *)malloc(datalen * sizeof(double));
    if (dim1) {
      (*S2_x)[i]    = (double *)malloc(dim1 * sizeof(double));
      (*Ksi_log)[i] = (double *)malloc(dim1 * sizeof(double));
    }
  }
  for (j=0;j<dim2;j++){
    (*data2_int)[j]  = (int *)malloc(datalen * sizeof(int));
    (*U_hat_table)[j]=(double**)malloc(ncentroids * sizeof(double*));
    for (i=0;i<ncentroids;i++) {
      (*U_hat_table)[j][i] =(double *)malloc(((int)(Ns[j]))*sizeof(double));
    }
  }
    
  return;
}
 
/************************************************************/

void
free_memory_B(int ncentroids,int dim1, int dim2, double ***Temp, double ***W, 
            double ***S2_x, double ***Ksi_log,
		double ****U_hat_table, int ***data2_int) {
  register int i,j;
  for (j=0;j<dim2;j++){
    for (i=0;i<ncentroids;i++) {
      free((*U_hat_table)[j][i]);
    }
    free((*data2_int)[j]);
    free((*U_hat_table)[j]);
  }

  for (i = 0; i < ncentroids; i++) { 
    free((*Temp)[i]); 
    if (dim1){
      free((*S2_x)[i]); 
      free((*Ksi_log)[i]);
    }
  }
  
  free(*Temp);
  if (dim1) {
    free(*S2_x);
    free(*Ksi_log);
  }
  if (dim2) {
    free(*U_hat_table);
    free(*data2_int);
  }
  return;
}



void
vdp_mk_log_lambda(double *Mumu, double *S2mu, double *Mubar, double *Mutilde, 
		  double *AlphaKsi, double *BetaKsi, 
		  double *KsiAlpha, double *KsiBeta, 
		  double *post_gamma, double *log_lambda, double *prior_alpha,
		  double *U_p, SEXP *U_hat,
		  long datalen, int dim1, int dim2, double *data1, double *data2, 
		  double *Ns, int ncentroids, 
		  double implicitnoisevar) {
  register long i, j, t;
  register int k;
  double  *U_hat_j;
  SEXP U_hat_j_SEXP;
  double       **W, **Temp, **S2_x,**Ksi_log,***U_hat_table;
  int          **data2_int;

  allocate_memory_B(datalen, ncentroids, dim1,dim2, &S2_x, &Ksi_log, &Temp,
		  &U_hat_table,&data2_int,Ns );

  for (j=0;j<dim2;j++){
    for(t=0;t<datalen;t++)
      data2_int[j][t]=((int)(data2[j*datalen+t]))-1;

    U_hat_j_SEXP = VECTOR_ELT(*U_hat, (int)j);
    U_hat_j=NUMERIC_POINTER(U_hat_j_SEXP);

    for(i=0;i<ncentroids;i++)
      for(k=0;k<Ns[j];k++)
        U_hat_table[j][i][k]=U_hat_j[k*ncentroids+i];
  }

  
  if (dim1) 
    compute_variance(ncentroids, dim1, KsiAlpha, KsiBeta, S2_x, Ksi_log);

  log_p_of_z_given_other_z_c(datalen, ncentroids, post_gamma, log_lambda);

  compute_tempmat(datalen,dim1,dim2,ncentroids,Temp,data1,data2_int,
		  Mubar,Mutilde,S2_x,Ksi_log,U_hat_table,Ns,
		  implicitnoisevar, log_lambda);
    
  fix_lambda(ncentroids, datalen, prior_alpha, log_lambda);
    
  free_memory_B(ncentroids, dim1,dim2, &Temp, &W, &S2_x, &Ksi_log, &U_hat_table, &data2_int);
  return;
}





/************************************************************/
/* bridge function                                          */

SEXP
mLogLambda(SEXP X1, SEXP X1_Columns, SEXP X1_Rows, 
             SEXP X2, SEXP X2_Columns,
             SEXP realS, SEXP OPTSimplicitnoisevar,
             SEXP hp_prior, SEXP HPposterior) {
  long datalen;
  int  dim1, dim2, ncentroids;
  double *Mumu, *S2mu, *Mubar, *Mutilde, 
    *AlphaKsi, *BetaKsi, *KsiAlpha, *KsiBeta, *U_p, *prior_alpha,
    *post_gamma, *log_lambda;
  double *data1;
  double *data2;
  SEXP olog_lambda, oU_hat;
  SEXP* U_hat;

  double *Ns;
  double implicitnoisevar;
  
  /******************** input variables ********************/
  
  
  /************ CONVERTED input variables ******************/
  /* data */
  PROTECT(X1 = AS_NUMERIC(X1));  
  data1   = NUMERIC_POINTER(X1);
  dim1    = INTEGER_VALUE(X1_Columns);
  datalen = INTEGER_VALUE(X1_Rows);

  PROTECT(X2 = AS_NUMERIC(X2));  
  data2   = NUMERIC_POINTER(X2);
  dim2    = INTEGER_VALUE(X2_Columns);

  Ns = NUMERIC_POINTER(realS);
  implicitnoisevar = NUMERIC_VALUE(OPTSimplicitnoisevar);
  

  /* Converted Initial Values of Model Parameters */

  if(dim1) {
    Mumu       = NUMERIC_POINTER(getListElement(hp_prior,"Mumu"));
    S2mu       = NUMERIC_POINTER(getListElement(hp_prior,"S2mu"));
    AlphaKsi   = NUMERIC_POINTER(getListElement(hp_prior,"AlphaKsi"));
    BetaKsi    = NUMERIC_POINTER(getListElement(hp_prior,"BetaKsi"));
    Mubar      = NUMERIC_POINTER(getListElement(HPposterior,"Mubar"));
    Mutilde    = NUMERIC_POINTER(getListElement(HPposterior,"Mutilde"));
    KsiAlpha   = NUMERIC_POINTER(getListElement(HPposterior,"KsiAlpha"));
    KsiBeta    = NUMERIC_POINTER(getListElement(HPposterior,"KsiBeta"));
  }
  if(dim2) {
    U_p         = NUMERIC_POINTER(getListElement(hp_prior,"U_p"));
    oU_hat      = getListElement(HPposterior,"Uhat");
    U_hat       = &oU_hat;
  }
  
  prior_alpha = NUMERIC_POINTER(getListElement(hp_prior,"alpha"));
  post_gamma  = NUMERIC_POINTER(getListElement(HPposterior,"gamma"));

  ncentroids = INTEGER_POINTER( GET_DIM(getListElement(HPposterior,"Mubar")) )[0];

  /******************** output variables ********************/
  PROTECT(olog_lambda     = NEW_NUMERIC(datalen*ncentroids));
  log_lambda = NUMERIC_POINTER(olog_lambda);


  vdp_mk_log_lambda(Mumu, S2mu, Mubar, Mutilde, 
		    AlphaKsi, BetaKsi, KsiAlpha, KsiBeta, 
		    post_gamma, log_lambda, prior_alpha,
		    U_p, U_hat,
		    datalen, dim1, dim2, data1, data2, 
		    Ns, ncentroids, implicitnoisevar);

  UNPROTECT(3);

  return olog_lambda;
}
 
/************************************************************/


void softmax(int dim1, int dim2, double *in, double *out)
{
  register int i, j;
  double rowsum, rowmax;

  for (i=0; i<dim1; i++) {
    rowmax = DBL_MIN;
    for (j=0; j<dim2; j++) {
      if (in[j*dim1 + i] > rowmax)
	rowmax = in[j*dim1 + i];
    }
    rowsum = 0;
    for (j=0; j<dim2; j++) {
      out[j*dim1 + i] = exp(in[j*dim1 + i] - rowmax);
      rowsum += out[j*dim1 + i];
    }
    for (j=0; j<dim2; j++) {
      out[j*dim1 + i] /= rowsum;
    }
  }
}




/************************************************************/
/* bridge function                                          */

SEXP
vdpSoftmax(SEXP matrix_M) {
  int dim1, dim2;
  double *in, *out;
  SEXP output, dims;
  
  /******************** input variables ********************/
  in      = NUMERIC_POINTER(matrix_M);
  dim1    = INTEGER_POINTER(GET_DIM(matrix_M))[0];
  dim2    = INTEGER_POINTER(GET_DIM(matrix_M))[1];
  PROTECT(dims = allocVector(INTSXP, 2));
  INTEGER(dims)[0] = dim1; INTEGER(dims)[1] = dim2;

  /******************** output variables ********************/
  PROTECT(output = NEW_NUMERIC(dim1*dim2));
  SET_DIM(output, dims);
  out = NUMERIC_POINTER(output);

  softmax(dim1, dim2, in, out);
  
  UNPROTECT(2);
  return output;
}
 
/************************************************************/

void sumlogsumexp(int dim1, int dim2, double *in, double *out)
{
  register int i, j;
  double rowsum, rowmax;

  *out = 0.0;
  for (i=0; i<dim1; i++) {
    rowmax = DBL_MIN;
    for (j=0; j<dim2; j++) {
      if (in[j*dim1 + i] > rowmax)
	rowmax = in[j*dim1 + i];
    }
    rowsum = 0.0;
    for (j=0; j<dim2; j++) {
      rowsum += exp(in[j*dim1 + i] - rowmax);
    }
    *out += rowmax + log(rowsum);
  }
}


/************************************************************/
/* bridge function                                          */

SEXP
vdpSumlogsumexp(SEXP matrix_M) {
  int dim1, dim2;
  double *in, *out;
  SEXP output, dims;
  
  /******************** input variables ********************/
  in      = NUMERIC_POINTER(matrix_M);
  dim1    = INTEGER_POINTER(GET_DIM(matrix_M))[0];
  dim2    = INTEGER_POINTER(GET_DIM(matrix_M))[1];
  PROTECT(dims = allocVector(INTSXP, 2));
  INTEGER(dims)[0] = 1; INTEGER(dims)[1] = 1;

  /******************** output variables ********************/
  PROTECT(output = NEW_NUMERIC(1));
  SET_DIM(output, dims);
  out = NUMERIC_POINTER(output);

  sumlogsumexp(dim1, dim2, in, out);
  
  UNPROTECT(2);
  return output;
}
 
/************************************************************/




