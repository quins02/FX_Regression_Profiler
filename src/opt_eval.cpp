#include <cmath>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <iostream>

using namespace std;

double opt_call(double S, double K, double r, double t){
	double ND = max(S-K,(double) (0));
	double D = ND * exp(-r*t);
	return D;
}

double opt_put(double S, double K, double r, double t){
	double ND = max(K-S,(double) (0));
	double D = ND * exp(-r*t);
	return D;
}

double opt_dig_call(double S, double K, double r, double t){
	double ND;
	if(S>K){ND = 1;}
	else{ND=0;}
	double D = ND * exp(-r*t);
	return D;
}

double NORMFUNC(double x, void * params){
	double alpha = *(double *) params;
	double f = (1/sqrt(2*M_PI))*exp(-(alpha*x*x)/2);
	return f;
}

double NORMDIST(double x){
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	
	double result, error;
	double alpha = 1.0;

	gsl_function F;
	F.function = &NORMFUNC;
	F.params = &alpha;

	gsl_integration_qag (&F, -100, x, 0, 1e-7, 1000, 3, w, &result, &error);

	return result;

	gsl_integration_workspace_free(w);

}

double call_CF(double S, double K, double R, double t, double vol){
	double d1 = (1/(vol * sqrt(t)))*(log(S/K) +(R+(vol*vol*0.5))*(t));
	//cout<<d1<<endl;
	double d2 = d1 - vol*sqrt(t);

	return NORMDIST(d1)*S - NORMDIST(d2)*K*exp(-R*t);
}

double dig_call_CF(double S, double K, double R, double t, double vol){
	double d1 = (1/(vol * sqrt(t)))*(log(S/K) +(R+(vol*vol*0.5))*(t));
	double d2 = d1 - vol*sqrt(t);

	return 10*NORMDIST(d2)*exp(-R*t);
}