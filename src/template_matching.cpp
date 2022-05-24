
/**
 * @file template_matching.cpp
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#include "template_matching.h"

#include <algorithm>
#include <cmath> 
#include <fstream>
#include <iostream>
#include <limits>
#include <random>

#include "omp.h"

#include "BIC.h"
#include "FDR.h"
#include "load.h"
#include "model.h"

using namespace std;

double nINF	=-exp(1000);
double INF 	= exp(1000);

/**
 * @brief 
 * 
 * @param X 
 * @param window 
 * @return vector<vector<double>> 
 */
vector<vector<double>> merge(vector<vector<double>> X, double window){
  int N = X.size(), t=0;
  vector<vector<double> > nX;
  while (t < N){
    double start = X[t][0]-window, stop = X[t][1], n=0;
    vector<double> S = {0.0,0.0,0.0};
    while (t < N and stop+window > X[t][0]-window ){
      stop = X[t][1];
      S[0]+=X[t][2], S[1]+=X[t][3], S[2]+=X[t][4];
      t+=1,n+=1;
    }
    vector<double> row = {start , stop + window, S[0]/n, S[1]/n, S[2]/n };
    nX.push_back(row);
  }
  return nX;
}


vector<vector<double>> bubble_sort3(vector<vector<double>> X){ //sort vector of vectors by second
	bool changed=true;
	if (X.size()<2){
		return X;
	}
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i][2] < X[i+1][2]){
				vector<double> copy 	= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}

int sample_centers(vector<double> centers, double p){
	random_device rd;
	mt19937 mt(rd());
	// default_random_engine generator;    // Unused...
	uniform_int_distribution<int> distribution(0,centers.size()-1);
	int i 	= distribution(mt);
	return i;//required in model.o (ugh...)
}

/**
 * @brief 
 * 
 * @param data 
 * @param BIC_values 
 * @param densities 
 * @param densities_r 
 * @param window 
 * @param sigma 
 * @param lambda 
 * @param foot_print 
 * @param pi 
 * @param w 
 * @param thr 
 */
void BIC_template(segment * data,  double * BIC_values, double * densities, double * densities_r, double window, 
		  double sigma, double lambda, double foot_print, double pi, double w, int thr){
  double vl;
  int NN 	= int(data->XN);
  int threads  	= thr; // omp_get_max_threads();
  int counts 		= NN / threads;
  #pragma omp parallel num_threads(threads)
  {
    int tid 	= omp_get_thread_num();
    int start 	= tid*counts;
    int stop 	= (tid+1)*(counts);
    if (tid+1 == threads){
      stop 	= NN;
    }
    int j = start, k =start;
    double N_pos=0,N_neg=0;
    double total_density;
    for (int i = start; i < stop; i++) {
      while (j < data->XN and (data->Coordinate(j) - data->Coordinate(i)) < -window) {
        N_pos -= data->ForwardCoverage(j);
        N_neg -= data->ReverseCoverage(j);
        j++;
      }
      while (k < data->XN and (data->Coordinate(k) - data->Coordinate(i)) < window) {
        N_pos += data->ForwardCoverage(k);
        N_neg += data->ReverseCoverage(k);
        k++;
      }
      int aa=k < data->XN;
      int bb=j < data->XN;
      int cc=N_neg > 0;
      int dd=N_pos > 0;
      
      if (k < data->XN  and j < data->XN and k!=j and N_neg > 0 and N_pos > 0 and (data->Coordinate(k) - data->Coordinate(j)) > 1.75*window  ){
	total_density 	= (N_pos / (data->Coordinate(k) - data->Coordinate(j))) + (N_neg / (data->Coordinate(k) - data->Coordinate(j)));
	densities[i] 	= N_pos ;
	densities_r[i] 	= N_neg ;
	
	double mu 	= (data->Coordinate(k) + data->Coordinate(j)) /2.;
		
	BIC_values[i] 	= BIC3(data->X,  j,  k,  i, N_pos,  N_neg, 
			       sigma, lambda, foot_print, pi, w);
      }else{
	BIC_values[i] 	= 0;
	densities[i] 	= 0;
	densities_r[i] 	= 0;
      }
    }
  }
}

/**
 * @brief 
 * 
 * @param a 
 * @param b 
 * @param c 
 * @param x 
 * @param y 
 * @param z 
 * @return true 
 * @return false 
 */
bool check_hit(double a, double b, double c, double x, double y, double z){
  if (a > x and b > y and c > z){
    return true;
  }
  return false;
} 

//================================================================================================
/**
 * @brief The template matching algorithm from Azofeifa 2018
 * 
 * @param segments 
 * @param P     the ubiquitous parameters hash
 * @param SC    the slice ratio (some preset values)
 * @return double 
 */
double run_global_template_matching(vector<segment*> segments, 
            params * P, slice_ratio SC){
	
  double CTT                    = 5; //filters for low coverage regions, WHY hard coded?!!?

  double ns 			= stod(P->p["-ns"]);
  double window 		= stod(P->p["-pad"])/ns;
  double sigma, lambda, foot_print, pi, w;
  double ct 			= stod(P->p["-bct"]);
  
  sigma 	= stod(P->p["-sigma"])/ns , lambda= ns/stod(P->p["-lambda"]);
  foot_print= stod(P->p["-foot_print"])/ns , pi= stod(P->p["-pi"]), w= stod(P->p["-w"]);
  
  bool SCORES 		= not P->p["-scores"].empty();
  
  ofstream FHW_scores;
  
  if (SCORES){ FHW_scores.open(P->p["-scores"]); }
  string annotation;
  int prev, prev_start, stop;
  int N;
  int start, center;

  for (int i = 0; i < segments.size(); i++){
    double * BIC_values 	= new double[int(segments[i]->XN)];
    double * densities 		= new double[int(segments[i]->XN)];
    double * densities_r 	= new double[int(segments[i]->XN)];

    double l 		=  segments[i]->getXLength(); // maxX-segments[i]->minX;
    double ef 		= segments[i]->fN*( 2*(window*ns)*0.05  /(l*ns ));
    double er 		= segments[i]->rN*( 2*(window*ns)*0.05 /(l*ns ));
    double stdf 	= sqrt(ef*(1- (  2*(window*ns)*0.05/(l*ns )  ) )  );
    double stdr 	= sqrt(er*(1- (  2*(window*ns)*0.05 /(l*ns ) ) )  );
    BIC_template(segments[i],  BIC_values, densities, densities_r, window, sigma, lambda, foot_print, pi, w, P->threads);   
    double start=-1, rN=0.0 , rF=0.0, rR=0.0, rB=0.0;
    vector<vector<double>> HITS;
    for (int j = 1; j<segments[i]->XN-1; j++){
      if (SCORES){
        double vl 	= BIC_values[j];
        if (std::isnan(double(vl))){
          vl 		= 0;
        }
        FHW_scores<<segments[i]->chrom<<"\t"<<to_string(int(segments[i]->Coordinate(j-1)*ns+segments[i]->start))<<"\t";
        FHW_scores<<to_string(int(segments[i]->Coordinate(j)*ns+segments[i]->start ))<<"\t" <<to_string(vl)<<endl;
      }
      bool HIT = check_hit(BIC_values[j], densities[j], densities_r[j], SC.threshold, ef + CTT*stdf, er + CTT*stdr  );
      if ( HIT ) {
        if (start < 0){
          start = segments[i]->Coordinate(j-1)*ns+segments[i]->start;
        }
        start+=1, rN+=1 , rF+=densities[j], rR+=densities_r[j], rB+=log10( SC.pvalue(BIC_values[j]) + pow(10,-20)) ;	
      } 
      if(not HIT and start > 0 ){
        vector<double> row = {start , segments[i]->Coordinate(j-1)*ns+segments[i]->start, rB/rN , rF/rN, rR/rN  };
        HITS.push_back(row);
        start=-1, rN=0.0 , rF=0.0, rR=0.0, rB=0.0;
      } 		
    }
    for (int j = 0; j < HITS.size();j++){
      segments[i]->bidirectional_bounds.push_back(HITS[j]);	
    }
    segments[i]->bidirectional_bounds 	= merge(segments[i]->bidirectional_bounds, window*0.5);    
  }
  return 1.0;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

/**
 * @brief The refactored template matching algorithm
 * 
 * @param segments  All the chunks on which to run template matching
 * @param P     Joey's ubiquitous parameters hash
 * @param SC    the slice ratio (some preset values)
 * @return double 
 */
double RF_run_global_template_matching(vector<segment*> segments, 
            params * P, slice_ratio SC){
	
  double CTT                    = 5; //filters for low coverage regions, WHY hard coded?!!?

  double ns 			= stod(P->p["-ns"]);
  double window 		= stod(P->p["-pad"])/ns;
  double sigma, lambda, foot_print, pi, w;
  double ct 			= stod(P->p["-bct"]);
  
  sigma 	= stod(P->p["-sigma"])/ns , lambda= ns/stod(P->p["-lambda"]);
  foot_print= stod(P->p["-foot_print"])/ns , pi= stod(P->p["-pi"]), w= stod(P->p["-w"]);
  
  bool SCORES 		= not P->p["-scores"].empty();

  std::cout << "ns: "+ to_string(ns) + " win: " + to_string(window) + " ct: " + to_string(ct)
    << std::endl;
  
  ofstream FHW_scores;
  
  if (SCORES){ FHW_scores.open(P->p["-scores"]); }
  string annotation;
  int prev, prev_start, stop;
  int N;
  int start, center;
 
  std::cout << "segment num: " + to_string(segments.size()) << std::endl;
  for (int i = 0; i < segments.size(); i++){
    std::cout << "i: " + to_string(i) << std::endl;
    double * BIC_values 	= new double[int(segments[i]->XN)];
    double * densities 		= new double[int(segments[i]->XN)];
    double * densities_r 	= new double[int(segments[i]->XN)];

    double l 		=  segments[i]->getXLength(); // maxX-segments[i]->minX;
    double ef 		= segments[i]->fN*( 2*(window*ns)*0.05  /(l*ns ));
    double er 		= segments[i]->rN*( 2*(window*ns)*0.05 /(l*ns ));
    double stdf 	= sqrt(ef*(1- (  2*(window*ns)*0.05/(l*ns )  ) )  );
    double stdr 	= sqrt(er*(1- (  2*(window*ns)*0.05 /(l*ns ) ) )  );
    std::cout << "l: "+ to_string(l) + " ef: " + to_string(ef) + " er: " + to_string(er)
      + " stdf: "+ to_string(stdf) + " stdr: " + to_string(stdr) << std::endl;

    RF_BIC_template(segments[i],  BIC_values, densities, densities_r, window, sigma, lambda, foot_print, pi, w, P->threads);   

  /*
    for (int z = 0; z < int(segments[i]->XN); z++) {
      std::cout << "BIC: " + to_string(BIC_values[z])
         + " den: " + to_string(densities[z])
         + " denR: " + to_string(densities_r[z]) << std::endl;
    }
    */

    double start=-1, rN=0.0 , rF=0.0, rR=0.0, rB=0.0;
    vector<vector<double>> HITS;
    for (int j = 1; j<segments[i]->XN-1; j++){
      if (SCORES){
        double vl 	= BIC_values[j];
        if (std::isnan(double(vl))){
          vl 		= 0;
        }
        FHW_scores<<segments[i]->chrom<<"\t"<<to_string(int(segments[i]->Coordinate(j-1)*ns+segments[i]->start))<<"\t";
        FHW_scores<<to_string(int(segments[i]->Coordinate(j)*ns+segments[i]->start ))<<"\t" <<to_string(vl)<<endl;
      }
      bool HIT = check_hit(BIC_values[j], densities[j], densities_r[j], SC.threshold, ef + CTT*stdf, er + CTT*stdr  );
      if ( HIT ) {
        if (start < 0){
          start = segments[i]->Coordinate(j-1)*ns+segments[i]->start;
        }
        start+=1, rN+=1 , rF+=densities[j], rR+=densities_r[j], rB+=log10( SC.pvalue(BIC_values[j]) + pow(10,-20)) ;	
      } 
      if(not HIT and start > 0 ){
        vector<double> row = {start , segments[i]->Coordinate(j-1)*ns+segments[i]->start, rB/rN , rF/rN, rR/rN  };
        HITS.push_back(row);
        start=-1, rN=0.0 , rF=0.0, rR=0.0, rB=0.0;
      } 		
    }
    for (int j = 0; j < HITS.size();j++){
      segments[i]->bidirectional_bounds.push_back(HITS[j]);	
    }
    segments[i]->bidirectional_bounds 	= merge(segments[i]->bidirectional_bounds, window*0.5);    
  }
  return 1.0;
}

/**
 * @brief Refactored BIC_template function (for understand this damn thing!)
 * 
 */
void RF_BIC_template(segment *data, double *BIC_values, double *densities, double *densities_r, double window,
                     double sigma, double lambda, double foot_print, double pi, double w, int thr)
{
  double vl;
  int NN = int(data->XN);

  int start = 0; int stop = NN;
  int j = start, k = start;
  double N_pos = 0, N_neg = 0;
  double total_density;

  // std::cout << data->write_withData() << std::endl;
  std::cout << data->write_allScalar() << std::endl;

  for (int i = start; i < stop; i++) {
    std::cout << "i: " + to_string(i); 
    // std::cout << "j: " + to_string(j);
    while (j < data->XN and (data->Coordinate(j) - data->Coordinate(i)) < -window) {
      //std::cout << " " + to_string(j);
      N_pos -= data->ForwardCoverage(j);
      N_neg -= data->ReverseCoverage(j);
      j++;
    }
    // std::cout << std::endl;
    // std::cout << "k: " + to_string(k);
    while (k < data->XN and (data->Coordinate(k) - data->Coordinate(i)) < window) {
      // std::cout << " " + to_string(k);
      N_pos += data->ForwardCoverage(k);
      N_neg += data->ReverseCoverage(k);
      k++;
    }
    // std::cout << std::endl;
    std::cout << " j: " + to_string(j);
    std::cout << " k: " + to_string(k) << std::endl;

    int aa = k < data->XN;
    int bb = j < data->XN;
    int cc = N_neg > 0;
    int dd = N_pos > 0;

/*
    std::cout << "aa: " + to_string(aa) << std::endl;
    std::cout << "bb: " + to_string(bb) << std::endl;
    std::cout << "cc: " + to_string(cc) << std::endl;
    std::cout << "dd: " + to_string(dd) << std::endl;
    */

    if (k < data->XN and j < data->XN and k != j and N_neg > 0 and N_pos > 0 and (data->Coordinate(k) - data->Coordinate(j)) > 1.75 * window)
    {
      total_density = (N_pos / (data->Coordinate(k) - data->Coordinate(j))) + (N_neg / (data->Coordinate(k) - data->Coordinate(j)));
      densities[i] = N_pos;
      densities_r[i] = N_neg;

      double mu = (data->Coordinate(k) + data->Coordinate(j)) / 2.;

      BIC_values[i] = RF_BIC3(data->X, j, k, i, N_pos, N_neg,
                           sigma, lambda, foot_print, pi, w);
    }
    else
    {
      BIC_values[i] = 0;
      densities[i] = 0;
      densities_r[i] = 0;
    }
  }
}

double RF_BIC3(double ** X, int j, int k, int i,
	    double N_pos, double N_neg,  double sigma, double lambda, double fp, double pi, double w){

  double N                = N_pos + N_neg;
  double l     = X[0][k] - X[0][j];

  double uni_ll= LOG(pi/  (l))*N_pos + LOG((1-pi)/ (l))*N_neg;


  double pi2      = (N_pos+10000) / (N_neg + N_pos+20000);

  double emg_ll   = 0, p1=0.0,p2=0.0;
  EMG EMG_clf(X[0][i], sigma, lambda, w, pi2  );
  EMG_clf.foot_print      = fp;
  
  for (int i = j; i < k;i++ ){
    p1 = EMG_clf.pdf(X[0][i],1) + (1.0-w)*pi*(1.0/l) , p2 = EMG_clf.pdf(X[0][i],-1) + (1.0-w)*(1.0-pi)*(1.0/l) ;
    if (p1 > 0 and p2 > 0 ){//this should always evalulate!!
      emg_ll+=LOG( p1 )*X[1][i] + LOG( p2 )*X[2][i];
    }else{
    }
  }
  double emg_ratio        = (-2*uni_ll + LOG(N)) / (-2*emg_ll + 20*LOG(N))  ;
  return emg_ratio;
}


