#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <gsl/gsl_statistics.h>

using namespace std;

vector <vector <double> > rand_gsl(int seed, int path, int T, double dt){
	double DRIFT=0;
	const gsl_rng_type * Q = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Q);
	gsl_rng_set(r, seed);

	gsl_rng_env_setup();

	double t = 0;
	int i = 0, k = 1;
	// vector <double> A;
	// vector <vector <double> > B;
	// vector <vector <vector <double> > > C;
	vector <double> A1;
	vector <vector <double> > B1;
	B1.clear();

	double tmp;

	//output in order PATH, TIME
	while(i<path){
		A1.push_back(100);
		while(t<T-dt){
			tmp = DRIFT*dt+(gsl_rng_uniform_pos(r) - 0.5);
			A1.push_back(A1[k-1]+(log(A1[k-1])*tmp));
			// A.push_back(tmp);
			t+=dt;
			k++;
		}
		t=0;
		k=1;
		// C.push_back(B);
		// B.clear();
		B1.push_back(A1);
		A1.clear();
		i++;
	}

	gsl_rng_free(r);

	return B1;
}

vector <vector <double> > CorMat(vector <vector <double> > Data){
	size_t Row = Data.size();
	size_t Col = Data[0].size();

	vector < vector <double> > Mat;
	vector <double> tmpVec;

	vector <double> tmp1;
	vector <double> tmp2;
	double * a;
	double * b;
	for (size_t i = 0; i < Row; i++){
		for (size_t j = 0; j < Row; j++){
			for (size_t k = 0; k < Col; k++){
				tmp1.push_back(Data[i][k]);
				tmp2.push_back(Data[j][k]);
			}
			a = &tmp1[0];
			b = &tmp2[0];
			tmpVec.push_back(gsl_stats_correlation(a,1,b,1,Col));
			tmp1.clear();
			tmp2.clear();
		}
		Mat.push_back(tmpVec);
		tmpVec.clear();
		
	}

	return Mat;
}

vector< vector <double> > PathGen(double seed, int PATH, double T, double dt, string HistData){
	
	
	double k = 170;		//Mean reversion parameter
	double V = 0.01;	//Average interest rate
	double o = 0.05;	//Interest rate volatility

	//Stochastic Interest rate 1 parameters
	double Sk = 170;	//Mean reversion parameter
	double SV = 0.06;	//Average interest rate
	double So = 0.05;	//Interest rate volatility

	const gsl_rng_type * Q = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Q);
	gsl_rng_set(r, seed);

	const gsl_rng_type * q = gsl_rng_default;
	gsl_rng * W = gsl_rng_alloc(q);
	gsl_rng_set(W, seed*seed);

	const gsl_rng_type * O = gsl_rng_default;
	gsl_rng * S = gsl_rng_alloc(O);
	gsl_rng_set(S, seed*seed*seed);

	gsl_rng_env_setup();

	vector< vector <double> > MultOut;

	vector<double> R1_vec;
	vector<double> A;
	vector<double> R2_vec;

	vector<double> Hist;

	ifstream data (HistData.c_str());
	string line;
	double TEMP;
	
	while(getline(data,line)){
		stringstream Lin(line);
		Lin >> TEMP;
		Hist.push_back(TEMP);
	}	

	double R1, R2, tmp;
	//int m;

	for(int COUNT = 0 ; COUNT < PATH ; COUNT++){
		R1 = -0.005;	//Initial Interest Rate T=0
		R2 = 0.006;	//Initial Volatility T=0
		
		//m=1;	
		
		R2_vec.push_back(R2);		
		R1_vec.push_back(R1);

		//tmp = Hist[Hist.size()-1];	//Initial Asset Price T=0
		tmp =1;
		A.push_back(tmp);

		for( int i = 1; i<(T/dt) ; i++ ){
			//if(i==(int)((m/20)*(T/dt))){
			R1+=k*(V-R1)*dt + R1*(gsl_ran_gaussian(W,o));
			//	m++;
			//}	
			R2 += Sk*(SV-R2)*dt + R2*(gsl_ran_gaussian(W,So));
			R2=abs(R2);
			R2_vec.push_back(R2);
			R1_vec.push_back(R1);
			//tmp += tmp*(R*dt+(gsl_ran_gaussian(r,0.006)));
			tmp += tmp*((R1-R2)*dt+(gsl_ran_gaussian(r,0.01)));
			A.push_back(tmp);
		}
		MultOut.push_back(A);
		R1_vec.clear();
		A.clear();
		R2_vec.clear();
	}

	return MultOut;
}