#include "../src/rng_func.h"
#include "../src/mat_opp.h"
#include "../src/opt_eval.h"
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

vector < double > Extract(vector<double> Vector, vector<int> Index);
vector <vector <vector <double> > > Mat3Ext(vector <vector <vector <double> > > S);

int main(){
	//Timing Start
	clock_t start;
	double duration;
	start = clock();

	vector <vector <vector <double> > > X = PathGen(time(NULL),400000,1, 0.01, 1);
	vector <vector <vector <double> > > S = ExpPolyReg(X[0],2,1);
	X.clear();
	// vector <vector <vector <double> > > Y = Mat3Ext(S);

	// for(int i = 0 ; i<Y[0][0].size() ; i++){
	// 	cout<<Y[1][0][i]<<endl;
	// }

	//Record duration
	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout<<"Time to generate and extract paths: "<<duration<<"s"<<endl;

	return 0;
}

vector < double > Extract(vector<double> Vector, vector<int> Index){
	vector<double> v;
	for (size_t i = 0 ; i<Index.size() ; i++){
		 v.push_back(Vector[Index[i]]);
		 // cout<<v<<endl;
	}
	return v;
}

vector <vector <vector <double> > > Mat3Ext(vector <vector <vector <double> > > S){ 
	vector <int> tmp1;	

	int Time = S[0][0].size();
	int Paths = S[0].size();
	int Factor = S.size();

	for( int i = 0 ; i<Time ; i+=Time/12){
		for( int j = 0 ; j < Paths ; j++ ){
			if ( opt_call(S[1][j][i],1,0.0,1)>0 ){
				// cout<<i<<"\t"<<S[1][j][i]<<endl;
				tmp1.push_back(i);

			}
			// else { cout <<"OTM"<<endl;}
		}
	}

	vector <double> TEMP;
	vector < vector <double> > Xtra;
	vector < vector < vector <double> > > Xtract;
	for (int i = 0 ; i < Factor ; i++){
		for(int j = 0 ; j < Paths ; j++ ){
			TEMP = Extract(S[i][j],tmp1);
			Xtra.push_back(TEMP);
		}
		Xtract.push_back(Xtra);
		Xtra.clear();
	}
	tmp1.clear();
	S.clear();
	return Xtract;
}

