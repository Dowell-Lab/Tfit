#include "BIC.h"
#include "model.h"
#include <cmath>
using namespace std;

/** Returns the probability density of a sample point on the gaussian (normal) distribution.
 * @param x The point to sample (as represented in terms of the standard deviation relative to the mean)
 * @param mu mean of the distribution.
 * @param sigma variance of the distribution. 
 * @return 
 */
double gaussian(double x, double mu, double sigma){
   return (1.0/sqrt(2*sigma*M_PI ))*exp(-pow(x-mu,2)/(2*sigma));
}


/** Implements the Bayesian Information Criterion for scoring within the model.
 * As the model optimizes parameters, this function provides a score with which the model can determine
 * the magnitude of improvement.
 * 
 * In effect, the BIC3 function should determine how well a given iteration of the model fits its input data.
 * 
 * The BIC is defined as ln(n)*k-2ln(L), where L is the maximized likelihood function.
 * 
 * This likelihood function is defined as:
 * L=p(x|th, M), where th is the set of parameters that maximize said likelihood function.
 * 
 * In this case, p() is represented by a call to the functions within the EMG class, as EMG
 * encapsulates the distribution fitting code central to tfit's operation.
 * @param X Set of data points on which to compute optimal parameter values.
 * @param j Starting point within the input dataset.
 * @param k Stopping point within the input dataset.
 * @param i Used to index the mu parameter within the input data array.
 * @param N_pos Number of positive reads in the input dataset.
 * @param N_neg Number of negative reads in the input dataset.
 * @param sigma Model parameter sigma
 * @param lambda Model parameter lambda
 * @param fp used to compute a footprint parameter. Presumably specified in base pairs.
 * @param pi Not used.
 * @param w Not used.
 * @return BIC score.
 */
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

