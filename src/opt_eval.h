#ifndef OPT_EVAL_H
#define OPT_EVAL_H

double opt_call(double S, double K, double r, double t);
double opt_put(double S, double K, double r, double t);
double opt_dig_call(double S, double K, double r, double t);
double NORMFUNC(double x, void * params);
double NORMDIST(double x);
double call_CF(double S, double K, double R, double t, double vol);
double dig_call_CF(double S, double K, double R, double t, double vol);

#endif