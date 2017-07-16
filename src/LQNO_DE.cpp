// Illustrate map feature of TMB to perform likelihood ratio tests on a nRNAged array dataset.
#define TMB_LIB_INIT R_init_isomiRs
#include <TMB.hpp>
// LQNO regression family for the estimation of bias factors

template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(y);	// read counts
	DATA_MATRIX(X);   // library membership for each read count
	DATA_FACTOR(u_X); // RNA membership for each read count
	DATA_FACTOR(u_G); // Group (experimental condition) membership for each read count
	PARAMETER_VECTOR(b);  // overal level (group mean) for each group (mean submodel)
	PARAMETER_VECTOR(s_b); // overal level (group mean) for each group (sigma submodel)
	PARAMETER_MATRIX(u_m); // expression of microRNA relative to the mean (mean submodel)
	PARAMETER_MATRIX(u_s); // expression of microRNA relative to the mean (sigma submodel)
	PARAMETER_VECTOR(sigmu);  // standard deviation of the random effects (u_m) on the mean submodel
	PARAMETER_VECTOR(sigsig); // standard deviation of the random effects (u_s) on the sigma submodel

	Type res = 0;  // negative log-likelihood 

	vector<Type> eta(y.size());   // linear predictor of the mean submodel
	vector<Type> logsigma(y.size()); // linear predictor of the sigma submodel
	int nRNA = u_m.rows();	// number of RNAs in the libraries
	int nG = u_m.cols();		// number of treatment groups

	array<Type> Dmu(nRNA, nG - 1);	// log-fold change relative to baseline (mean submodel)
	array<Type> Dsig(nRNA, nG - 1);  // log-fold change relative to baseline (sigma submodel)

	logsigma = X*s_b;
	eta = X*b;
	// construction of log-fold expression changes as transformed variables
	for (int i = 0; i<nRNA; i++){
		for (int j = 0; j<nG - 1; j++) {
			Dmu(i, j) = b[j + 1] + u_m(i, j + 1) - u_m(i, 0);
			Dsig(i, j) = s_b[j + 1] + u_s(i, j + 1) - u_s(i, 0);
		}
	}
	ADREPORT(Dmu);
	ADREPORT(Dsig);

	// linear predictor contribution for each miRNA
	for (long int i = 0; i<y.size(); i++){
		eta[i] += u_m(u_X[i] - 1, u_G[i] - 1);
		logsigma[i] += u_s(u_X[i] - 1, u_G[i] - 1);
	}
	// construction of negative log-likelihood

  for(long int i=0;i<y.size();i++){
    res -= dnorm(y[i],exp(eta[i]),exp(0.5*(eta[i]+log(1+exp(logsigma[i]+eta[i])))),true);
  }
  // random effect contribution to the negative log-likelihood
  for(int i=0;i<nRNA;i++){
    for(int j=0;j<nG;j++) {
      res -= dnorm(u_m(i,j),Type(0.0),sigmu[j],true);
      res -= dnorm(u_s(i,j),Type(0.0),sigsig[j],true);
    }
  }
  return res;
}
