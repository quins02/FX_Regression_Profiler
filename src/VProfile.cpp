#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "rng_func.h"
#include "opt_eval.h"	//Include optiom calculations functions
#include "mat_opp.h"	//Include matrix operations functions
#include <time.h>
#include <chrono>
#include <algorithm>

using namespace std;

vector < vector < vector < double > > >  BetaGen();
vector< vector<double> > Val();
//vector< vector <double> > PathGen(double seed, int PATH, double T, double dt, string HistData);
vector <double> Reg(vector< vector< vector <double> > > X, double DMONTH);
double parab(vector<double> Coef, double X);
double hyp(vector<double> Coef, double X);
double mean(vector <double> vec);
vector <double> X_percent(vector <double> vec, double Percentile);

//DECLARE GLOBAL VARIABLES
double STRIKE = 1.2;
double dt = 0.01 ;
int Path = 10;
int T = 1;
int Poly = 2;
bool OUT = 0;

int main(){
	//Timing Start
	clock_t start;
	double duration;
	start = clock();

	vector < vector < vector < double > > > X = BetaGen();
	vector< vector<double> > V = Val();

	//Declare Variables necessary for valutaion
	size_t Path_Length = V[0].size();
	size_t Num_Path = V.size();
	size_t Num_Buckets = X[0].size();

	// cout<<	Path_Length <<endl;
	// cout<< Num_Path <<endl;
	// cout<< Num_Buckets <<endl;
	// cout<<X.size() <<endl;
	// cout<<X[0].size() <<endl;
	// cout<<X[0][0].size() <<endl;
	// cout<<X[10][19][4]<<endl;


	int Tv = 0;
	int Interval = 12;

	double tmp;

	string CVANA = "tmp/VPath.dat";
	ofstream CVAA;
	CVAA.open(CVANA.c_str());

	for(size_t i = 0 ; i < V[0].size() ; i++){
		for(size_t j = 0 ; j < V.size() ; j++){
			CVAA<<V[j][i]<<"\t";
		}
		CVAA<<endl;
	}

	// vector < double > tmpBeta;
	vector < double > Val_Path;
	vector < vector < double > > Val_T;

	int i = (int) (Path_Length/Interval) ;

	for( Tv =  0 ; Tv < X.size()   ; Tv++ ){
	
		for( size_t j = 0 ; j < Num_Path ; j++ ){
			tmp = V[j][i];
			for ( size_t k = 0 ; k < Num_Buckets ; k++ ){
				vector < double > tmpBeta((X[Tv][k].begin() + 2), X[Tv][k].end());
				if( X[Tv][k][0] < tmp && tmp < X[Tv][k][1] ){
					Val_Path.push_back(parab(tmpBeta,tmp));
					// cout<<tmpBeta[0]<<endl;
					k = Num_Buckets;
					tmpBeta.clear();
				}
				else{}
			}
		}
		Val_T.push_back(Val_Path);
		// cout<<Val_T.size()<<"\t"<<Val_T[Val_T.size()-1].size()<<endl;
		Val_Path.clear();
		i += (int) (Path_Length/Interval);
	}

	string CVAN = "tmp/CVA.dat";
	ofstream CVA;
	CVA.open(CVAN.c_str());

	string CVAM = "tmp/CVAMean.dat";
	ofstream Mean;
	Mean.open(CVAM.c_str());

	string CVAPM = "tmp/CVAPMean.dat";
	ofstream PMean;
	PMean.open(CVAPM.c_str());

	for( size_t i = 0 ; i < Val_T.size() ; i++ ){
		for( size_t j = 0 ; j < Val_T[0].size() ; j++ ){
			CVA<<i+1<<"\t"<<Val_T[i][j]<<endl;
		}
		Mean<<i+1<<"\t"<<mean(Val_T[i])<<endl;
		PMean<<i+1<<"\t"<<mean(X_percent(Val_T[i],0.95))<<endl;
	}

	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to generate risk profile: "<<duration<<"s"<<endl;

	cout<<"Press any key to exit...";
	cin.ignore();
	cin.get();

	return 0;
}

vector < vector < vector < double > > >  BetaGen(){
	//Set regression options from user
	cout<<"Regression Pricing for FX Exotics\n----------------------------\n\n";
	cout<<"Strike price: ";
	cin>>STRIKE;
	// cout<<"Month to discount to: ";
	// cin>> DMONTH;
	cout<<"Number of Paths: ";
	cin>> Path;
	cout<<"Number of Buckets: ";
	int Buck;
	cin>> Buck;
	vector<double> Buck_vec;
	cout<<"Override bucket postions (equal number)? (Y=1, N=0): ";
	bool TF;
	cin>>TF;
	if(TF == 1){
		int TEMP;
		cout<<"Bucket positons: \n";
		for(int i = 0; i<Buck ; i++){
			cin>>TEMP;
			Buck_vec.push_back(TEMP);
		}
	}
	cout<<"Order of polynomial: ";
	cin>>Poly;
	cout<<"Output all paths to file (slow)(Y=1, N=0): ";
	cin>>OUT;


	string strikeS = "tmp/strike.dat";
	ofstream SK;
	SK.open(strikeS.c_str());

	SK<<STRIKE<<endl;

	//Timing Start
	clock_t start;
	double duration;
	start = clock();


	vector< vector <double> > X =  PathGen(time(NULL), Path, T, dt, "Hist.txt");
	vector< vector < vector <double> > > T = ExpPolyReg(X,Poly,1);

	//Declare variables for Bucketing
	int min, max;
	vector<double> v, Bucket;
	
	vector < vector < vector < double > > > RetVec;
	vector < vector < double > > BUCKETVEC;

	for(double DMONTH = 1; DMONTH<=12 ; DMONTH++){

		//Extract values for asset price from training paths at time t. 
		//Add to vector v
		for(size_t i=0; i<T[0].size(); i++){
			v.push_back(T[1][i][(int)((T[0][0].size()/12)*DMONTH)]);
		}

		//Setting auto buckets equal number in each
		vector<double> tmp_v= v;
		sort(tmp_v.begin(),tmp_v.end());
		if(TF == 0){
			for (int i = 1; i<Buck+1 ; i++){
				Buck_vec.push_back(tmp_v[(i*(int)(v.size()/(Buck+1)))]);
			}
		}

		//Find minimum and maximum values of the asset price and add buffer
		max = (int) *max_element(v.begin(),v.end()) + 1;
		min = (*min_element(v.begin(),v.end()) > 0.02)? 0.02: *min_element(v.begin(),v.end());

		//Add min, max & user elements to Bucket
		//sort in ascending order.
		Bucket.push_back(min);
		Bucket.push_back(max);
		for(int i = 0; i<Buck ; i++){
			Bucket.push_back(Buck_vec[i]);
		}

		sort(Bucket.begin(),Bucket.end());

		//SPLIT TO BUCKETS
		//And perform individual regressions, finding beta coefficents
		//and applying to function
		//Output to BucketFit.dat
		//**************NOT OPTIMISED******************
		string BUCK = "tmp/BucketFit" + static_cast<ostringstream*>( &(ostringstream() << DMONTH) )->str() + string(".dat");
		ofstream BF;
		BF.open(BUCK.c_str());

		//cout<<"Calculating bucketed regression coefficents...\n";

		bool YN;
		vector <double> tmp1;
		vector < vector <double> > tmp2;
		vector < vector < vector <double> > > tmp3;
		vector <double> Btmp;
		for(size_t j = 1 ; j<Bucket.size() ; j++){
			for(size_t a = 0 ; a<T.size(); a++){
				for( size_t i = 0 ; i<v.size() ; i++ ){
					YN = bool(Bucket[j-1]<v[i] && v[i]<=Bucket[j]);
					if(YN==1){
						for (size_t b = 0 ; b<T[0][0].size(); b++){
							tmp1.push_back(T[a][i][b]);
						}
						tmp2.push_back(tmp1);
						tmp1.clear();
					}
				}
				tmp3.push_back(tmp2);
				tmp2.clear();
			}
			//APPLY REGRESSION TO INDIVIDUAL BUCKETS
			Btmp = Reg(tmp3, DMONTH);
			// Generate function and output to BucketFit.dat
			for(double s = (Bucket[j-1]) ; s < ((Bucket[j])) ; s+=0.001){
				BF<<setprecision(4)<<s<<"\t"<<parab(Btmp,s)<<endl;
			}
			Btmp.insert(Btmp.begin(), Bucket[j]);
			Btmp.insert(Btmp.begin(), Bucket[j-1]);
			BUCKETVEC.push_back(Btmp);
			//Output beta coefficents to screen
			// cout<<"Between S= "<<Bucket[j-1]<<" and "<<Bucket[j]<<" ("<<tmp3[0].size()<<" paths): "<<endl;
			// for(size_t F = 0; F<Btmp.size(); F++){
			// 	cout<<"Beta coefficent "<<j<<"["<<F<<"] at "<<DMONTH<<": "<<Btmp[F]<<endl;
			// }
			Btmp.clear();
			tmp3.clear();
		}
		RetVec.push_back(BUCKETVEC);
		Bucket.clear();
		Buck_vec.clear();
		v.clear();
		tmp_v.clear();
		BUCKETVEC.clear();
		BF.close();

	}
		
		// string Full = "AllOut.dat";
		// ofstream AllOut;
		// AllOut.open(Full.c_str());

		// for (size_t i = 0 ; i < A[0].size(); i++){
		// 	for (size_t j = 0 ; j < A.size(); j++){
		// 		AllOut<<A[j][i]<<"\t";
		// 	}
		// 	AllOut<<endl;
		// }

	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to calculate Beta coefficents: "<<duration<<"s"<<endl;

	//Generate closed form solution for Call option
	//output to CFValue.dat
	// string name = "CFValue.dat";
	// ofstream CFS;
	// CFS.open(name.c_str());

	// for(double S  = 0.001 ; S<STRIKE*2 ; S+=0.001){
	// 	 CFS<<S<<"\t"<<call_CF(S, STRIKE, 0.01, ((12-DMONTH)/12), 0.06)<<endl;
	// }

	T.clear();

	return RetVec;
}

vector< vector<double> > Val(){
	clock_t start;
	double duration;
	start = clock();


	vector <vector <double> > X = PathGen((time(NULL)*time(NULL)),2000,T,dt, "Hist.txt");


	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to generate valuation paths: "<<duration<<"s"<<endl;

	return X;

}

vector <double> Reg(vector< vector< vector <double> > > X, double DMONTH){
	//Timing Start
	// clock_t start;
	// double duration;
	// start = clock();


	//Initilise counting variables
	size_t i=0, j=0;

	//Evaluate CF at time t 
	double tau = 0., R = 0, n=0, tmp=0;
	vector<double> Phi;
	while(i<X[0].size()){
			//Exotic to be priced (i.e. option/spread/barrier SEE opt_eval.cpp for possible options)
			// n = -2*opt_put(X[1][i][T/dt - 1],STRIKE,0.01,1)
			// 	+1*opt_put(X[1][i][T/dt - 1],STRIKE-0.1,0.01,1)
			// 	+1*opt_put(X[1][i][T/dt - 1],STRIKE+0.1,0.01,1);
			n = barrier_call(X[1][i],1.05,1,0,STRIKE,0.01,1-tau);
			Phi.push_back(n);
			tau+=dt;
		//Opt<<n<<endl;
		j=0;
		i++;
	}

	i=0;

	//Extract Variables at time t
	vector<double> tmp_vec;
	vector<double> tmp_phi;
	vector< vector<double> > tmp_mat;
	vector< vector<double> > phi_T;

	i=0;
	while(i<X.size()){
		j=0;
		tmp_phi.clear();
		while(j<X[0].size()){
			tmp_vec.push_back(X[i][j][(int)(((T/dt)/12)*DMONTH)-1]);
		 	j++;
		}
		tmp_mat.push_back(tmp_vec);
		tmp_vec.clear();
		i++;
	}

	phi_T.push_back(Phi);

	vector< vector<double> > MATRIX = transpose(tmp_mat);
	vector< vector<double> > MATRIXMultD = MatMult (tmp_mat,MATRIX);
	vector< vector<double> > MATRIXMultN = MatMult(tmp_mat ,transpose(phi_T));
	vector< vector<double> > Beta = MatMult(MatInv(MATRIXMultD), MATRIXMultN);

	vector <double> BetaVec;
	for(size_t i = 0; i < (Beta).size(); i++){
		for(size_t j = 0; j < (Beta)[0].size(); j++){
			tmp = (Beta)[i][j];
		}
		BetaVec.push_back(tmp);
	}

	// duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	// cout<<"Time for regression & training: "<<duration<<"s"<<endl;

	return BetaVec;
}

double parab(vector<double> Coef, double X){
	double Y = 0;
	for(size_t i =0 ; i<Coef.size(); i++){
		Y+=pow(X,i)*Coef[i];
	}
	// if(Y<0){Y=0;}
	// if(Y>1){Y=1;}
	return Y;
}

double hyp(vector<double> Coef, double X){
	double TOT = 0;
	for(size_t i =0 ; i<Coef.size(); i++){
		TOT+=pow(X,i)*Coef[i];
	}
	return TOT;
}

double mean(vector <double> vec){
	int S = vec.size();
	double SUM = 0;
	for(int i = 0 ; i < S ; i++){
		SUM+=vec[i];
	}
	return SUM/S;
}

vector <double> X_percent(vector <double> vec, double Percentile){
	sort(vec.begin(),vec.end());
	int S = vec.size() - 1;
	int Bracket = S-(int)(S*Percentile);
	vector <double> Percent (vec.end()-Bracket, vec.end());
	return Percent;
}
