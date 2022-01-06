/**
 * @file BIC.cpp
 * @author Joey Azofeifa 
 * @brief 
 * @version 0.1
 * @date 2016-11-22
 */
#include "model.h"
#include <cmath>
using namespace std;


double gaussian(double x, double mu, double sigma){
   return (1.0/sqrt(2*sigma*M_PI ))*exp(-pow(x-mu,2)/(2*sigma));
}

double BIC3(double ** X, int j, int k, int i,
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

/** The version below was changed in check-in 6197dad63d
 *
 *
double BIC3(double ** X, int j, int k, int i,
            double N_pos, double N_neg,
            double sigma, double lambda, double fp, double pi, double w) {
   int res         = 5;
   double N        = N_pos + N_neg;
   double l        = X[0][k] - X[0][j];
   double pi2      = (N_pos + 10000) / (N_neg + N_pos + 20000);
   double uni_ll   = log(pi2 /   l ) *  N_pos  + log((1 - pi2) /  l ) * N_neg ;
   double MU       = X[0][i];
   double fp_delta = fp / res, best_ll = 0;
   EMG EMG_clf(MU, sigma, lambda, 1.0, pi2 );
   for (int rs = 0 ; rs < res; rs++){
      double emg_ll   = 0, p1 = 0.0, p2 = 0.0;
      EMG_clf.foot_print      = rs*fp_delta;
      for (int iter = j; iter < k; iter++ ) {
         p1 = EMG_clf.pdf(X[0][iter] , 1)  ;
         p2 = EMG_clf.pdf(X[0][iter] , -1) ;

         emg_ll += log( p1 ) *  X[1][iter]  ;
         emg_ll += log( p2 ) *  X[2][iter] ;
      }
      if (emg_ll > best_ll || rs==0){
         best_ll  = emg_ll;
      }
   }

   double arg_bic =  (-2*uni_ll + log(N))/(-2*best_ll + log(N)*4)      ;
   return arg_bic;
}

*/
