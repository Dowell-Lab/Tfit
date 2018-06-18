#include "FDR.h"
#include "load.h"
#include "model.h"
#include <limits>
#include <iostream>
#include <algorithm>
#include <random>
#include "omp.h"
#include <cmath>
#include <vector>
#include <string>
#include "BIC.h"
#include <math.h> 
using namespace std;

/** Default constructor for the normal class. Does nothing.
 */
normal::normal(){}
/** Full constructor for the normal class. Accepts a parameter X, mean MU, and standard deviation STD.
 * @param X Sample parameter for the distribution.
 * @param MU Mean of the distribution.
 * @param STD Standard deviation of the distribution
 */
normal::normal(double X, double MU, double STD){
  mean=MU, x=X, std=STD, c2 = sqrt(2);
}

/** Returns the cumulative distribution function at the sample point specified.
 * @param x Sample point
 * @return Cumulative distribution at the sample point. (ie. how likely all values less than x combined are).
 */
double normal::cdf(double x){
  double pv = 0.5;
  double Z  =erf( (x-this->mean) / (this->std*c2));
  pv        = 0.5*(1+Z);
  return pv;
}

/** Returns the probability density of the distribution at the specified point.
 * @param x Point at which to sample the distribution.
 * @return Probability density at that point.
 */
double normal::pdf(double x){
  return 1.0 / (sqrt(2*M_PI)*this->std  )*exp(-pow(x-this->mean,2)/(2*pow(this->std,2) ));  
}

/** Default constructor for the exponential class.
 */
exponential::exponential(){};

/** Full constructor for the the exponential class.
 * @param MU Center of the distribution.
 * @param LAMBDA Exponential rate of decay.
 */ 
exponential::exponential(double MU, double LAMBDA){
  this->mu = MU , this->lambda = LAMBDA;
}

/** Returns the probability density of the distribution at the specified point.
 * @param x Point at which to sample the probability density.
 * @return Probability density of the distribution at the specified point.
 */
double exponential::pdf(double x){
  if (x>=mu){
    return this->lambda*exp(-this->lambda*(x-mu) ); 
  }
  return 0.0;
}

/** Default constructor for the pareto class.
 */
pareto::pareto(){};

/** Full constructor for the pareto class.
 * @param MU Center of the distribution
 * @param ALPHA Rate of decay.
 */
pareto::pareto(double MU, double ALPHA){
  this->mu = MU , this->alpha=ALPHA;
};

/** Returns the probability density of the distribution at the specified point.
 * @param x Point at which to sample the probability density.
 * @return Probability density of the distribution at the specified point.
 */
double pareto::pdf(double x){
  if (x>=mu){
    return this->alpha*pow(this->mu,alpha) / pow(x, this->alpha+1);
  }
  return 0.0;
}

/** Returns whether or not the value passed is finite and a number.
 * @param x Value to check.
 * @return Whether or not the value is finite and a number.
 */
bool check_value(float x){
  if (isnan(x) ){
    return false;
  }
  if (not isfinite(x) ){
    return false;
  }
  return true;
}

/** Attempts to train a model on the given input data.
 * @param X Input data on which to train the model.
 * @param mean Starting mean for the model. If the model converges, then this value will be overwritten with the optimized mean.
 * @param sigma Starting standard deviation for the model. If the model converges, then this value will be overwritten with the standard deviation (scaled to base pairs?)
 * @param w Weighting parameter within the model. This is continuously overwritten during training.
 * @param c Exponential decay parameter used within the model. This value is overwritten immediatley.
 * @param EXP Whether or not to use a larger exponential decay parameter during training.
 * @return Whether or not the model converged.
 */
bool EM(vector<vector<double>> X , double & mean, double & sigma, double & w, double & c , bool EXP ){
  mean = 0.3 , sigma = 0.1, w = 0.5;
  if (EXP){
    c=0.5;
  }else{
    c=1.34;
  }
  normal N(1.0, mean,sigma) ;
  exponential E(mean, c);
  pareto P(mean,c);
  
  int t=0, T=300;
  double prevw = -0.1;
  bool converged=false;
  while (t < T and not converged  ){
    double EX=0.0,EX2=0.0, EY=0.0,R1=0.0, R2=0.0,p1=0.0,p2=0.0,r1=0.0,r2=0.0;
    for (int i = 0 ; i < X.size();i++){
      
      double x = X[i][0], y = X[i][1];
      p1= (1-w)*N.pdf(x);
      if (EXP){
	p2=w*E.pdf(x);
      }else{
	p2=w*P.pdf(x);
      }
      if ((p1 + p2) > 0){
	r1= p1 / (p1 + p2) ,r2=p2 / (p1 + p2);
	R1 += r1*y,R2+=r2*y;
	EX +=(r1*x)*y;
	if (EXP){
	  EY +=(r2*x)*y;
	}else{
	  EY += (log(x) - log(N.mean))*r2*y;
	}
	EX2 +=(r1*pow(x-N.mean,2))*y;
      }
    }
    
    N.mean = EX /R1;
    N.std = sqrt(EX2 / R1);
    if (EXP){
      E.lambda = R2/EY ;
    }else{
      P.alpha = R2/EY;
    }
    if( not check_value(N.mean)  or not check_value(N.std)  ){
      converged = false;
      break;
    }
    w = R2 / (R1 + R2);
    if (abs(w-prevw) < pow(10,-6)){
      converged=true;
    } 
    prevw=w, t+=1;
  }
  if(not converged){
    mean = 0.3 , sigma    = 0.05;
  }else{
    
    mean = N.mean , sigma = N.std*15;
  }
  return converged;
} 

/** Default constructor for slice_ratio. Does nothing.
 */
slice_ratio::slice_ratio(){};

/** Full constructor for slice_ratio
 * @param ST Start position (base pairs)
 * @param SP Stop position (base pairs)
 * @param BINS Number of bins to use.
 */
slice_ratio::slice_ratio(double ST, double SP, int BINS){
  this->start = ST, this->stop = SP, this->bins=BINS;
  double step = (stop - start)/double(this->bins);
  for (double c = this->start; c < this->stop ; c+=step){
    vector<double> row  = {c, 0.0}; //step, N, X, X^2
    this->XY.push_back(row);
  }
};

/** Determines the closest slice to the specified sample value.
 * @param x sample value at which to get a slice.
 * @return Index of the closest slice.
 */
int slice_ratio::get_closest(double x){
  int c = 0;
  while (c < this->XY.size() and this->XY[c][0] < x ){
    c+=1;
  }
  if (c > 0){
    return c -1;
  }
  return 0;

}

/** Increments the count of the closest slice to the sample value y.
 * @param y sample value at which to index the slice.
 */
void slice_ratio::insert(double y){
  int c = this->get_closest(y);
  this->XY[c][1]+=1;
}

/** Updates the slice given a specified p-value.
 * @param pval p-value to use.
 */
void slice_ratio::set(double pval){
  this->mean=0.1, this->std=0.1, this->w=0.5, this->c=1.0;
  this->converged = EM(this->XY, this->mean, this->std, this->w, this->c, 1);
  this->norm_all  = normal(0.0, this->mean, this->std);
  double z        = this->norm_all.mean - this->norm_all.std;
  double pv       = this->norm_all.cdf(z);
  while (pv < (1.0 - pval) and z < 100.0){
    z+=0.01;
    pv=this->norm_all.cdf(z);
  }
  threshold=z;
}

/** Performs all of the same operations as set except for updating the mean and convergence value.
 * @param pval p-value to use.
 */
void slice_ratio::set_2(double pval){
  this->norm_all  = normal(0.0, this->mean, this->std);
  double z        = this->norm_all.mean - this->norm_all.std;
  double pv       = this->norm_all.cdf(z);
  while (pv < (1.0 - pval) and z < 100.0){
    z+=0.01;
    pv=this->norm_all.cdf(z);
  }
  threshold=z;
}

/** Computes the p-value based on the values within the slice.
 * @param y The point at which to compute the p-value.
 * @return p-value
 */
double slice_ratio::pvalue(double y){
  double pv = 1.0-this->norm_all.cdf(y);
  return pv;
}

/** Computes a slice_ratio given a set of segments and various other parameters.
 * Use get_slice_pwrapper instead.
 * @depreciated
 * @param segments Set of segments over which to compute the slice_ratio.
 * @param N Number of parameters to estimate.
 * @param CC Threshold (?)
 * @param P Obsolete params object from which to read command line parameters.
 */
slice_ratio get_slice(vector<segment *> segments, int N, double CC, params * P){
  double sigma, lambda, fp, pi, w, window, pval_threshold,ns;
 
  window        = stod(P->p["-pad"]), ns=stod(P->p["-ns"]) ;
  sigma         = stod(P->p["-sigma"])/ns , lambda= ns/stod(P->p["-lambda"]);
  fp            = stod(P->p["-foot_print"])/ns , pi= stod(P->p["-pi"]), w= stod(P->p["-w"]);
  pval_threshold= stod(P->p["-bct"]) ;
  int CN     = segments.size();
  double min_x  = -1 , max_x = -1, n = 0;
  random_device rd;
  mt19937 mt(rd());
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0,1);
  vector<double> XY(N);
  vector<double> CovN(N);
  for (int i = 0 ; i < XY.size(); i++){
    XY[i]=0.0, CovN[i]=0.0;
  }
  #pragma omp parallel for
  for (int n = 0 ; n < N ; n++){
    double U       = distribution(mt);
    double U2      = distribution(mt);
    int NN         = int(U*(CN-1));
    segment * data = segments[NN];
    int c          = U2*int(data->XN);
    int j = c,  k  = c;
    double N_pos = 0 , N_neg =0 ;
    while (j > 0 and (data->X[0][c] - data->X[0][j] )< window){
      N_pos+=data->X[1][j];
      N_neg+=data->X[2][j];
      j--;
    }
    while (k < data->XN and (data->X[0][k] - data->X[0][c] )< window  ){
      N_pos+=data->X[1][k];
      N_neg+=data->X[2][k];
      k++;
    }
    CovN[n] = N_pos + N_neg;
    if (N_pos + N_neg > CC and (data->X[0][k] - data->X[0][j]) > 1.75*window  ){
      
      double val =  BIC3(data->X,  j,  k,  c, N_pos,  N_neg, sigma , lambda, fp , pi, w);
      if (val >0 ){
        XY[n]=val,CovN[n]=N_pos+N_neg;
      }
    }
  }
  string job_name    = P->p["-N"];
  string log_out_dir = P->p["-log_out"];
  ofstream FHW;
  FHW.open(log_out_dir+job_name+"_random_BIC_ratios.csv");
  FHW<<"ratio,N\n";
  for (int i = 0 ; i < XY.size();i++){
    FHW<<to_string(XY[i])+","+to_string(CovN[i]) <<endl;
  }
  FHW.close();
  
  


  //---------------
  for (int n = 0 ; n < XY.size(); n++){
    if (min_x < 0 or XY[n] < min_x ){
      min_x = XY[n];
    }
    if (max_x < 0 or XY[n] > max_x ){
      max_x = XY[n];
    }
  }
  //-------------
  slice_ratio SC(min_x,max_x,400);
  for (int n = 0 ; n < XY.size();n++){

    if (XY[n] > 0.0){
      SC.insert(XY[n] );
    }
  }
  SC.set(pval_threshold);
  return SC;
}

/** Computes a slice_ratio given a set of segments and various other parameters.
 * @param segments Set of segments over which to compute the slice_ratio.
 * @param N Number of parameters to estimate.
 * @param CC Threshold (?)
 * @param pw ParamWrapper object from which to read command line parameters.
 */
slice_ratio get_slice_pwrapper(vector<segment *> segments, int N, double CC, ParamWrapper *pw){
  double sigma, lambda, fp, pi, w, window, pval_threshold,ns;
 
  window        = pw->pad, ns=pw->ns;
  sigma         = pw->sigma/ns , lambda= ns/pw->lambda;
  fp            = pw->footPrint/ns , pi=pw->pi, w=pw->w;
  pval_threshold= pw->llrthresh;
  int CN     = segments.size();
  double min_x  = -1 , max_x = -1, n = 0;
  random_device rd;
  mt19937 mt(rd());
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0,1);
  vector<double> XY(N);
  vector<double> CovN(N);
  for (int i = 0 ; i < XY.size(); i++){
    XY[i]=0.0, CovN[i]=0.0;
  }
  #pragma omp parallel for
  for (int n = 0 ; n < N ; n++){
    double U       = distribution(mt);
    double U2      = distribution(mt);
    int NN         = int(U*(CN-1));
    segment * data = segments[NN];
    int c          = U2*int(data->XN);
    int j = c,  k  = c;
    double N_pos = 0 , N_neg =0 ;
    while (j > 0 and (data->X[0][c] - data->X[0][j] )< window){
      N_pos+=data->X[1][j];
      N_neg+=data->X[2][j];
      j--;
    }
    while (k < data->XN and (data->X[0][k] - data->X[0][c] )< window  ){
      N_pos+=data->X[1][k];
      N_neg+=data->X[2][k];
      k++;
    }
    CovN[n] = N_pos + N_neg;
    if (N_pos + N_neg > CC and (data->X[0][k] - data->X[0][j]) > 1.75*window  ){
      
      double val =  BIC3(data->X,  j,  k,  c, N_pos,  N_neg, sigma , lambda, fp , pi, w);
      if (val >0 ){
        XY[n]=val,CovN[n]=N_pos+N_neg;
      }
    }
  }
  string job_name    = pw->jobName;
  string log_out_dir = pw->logDir;
  ofstream FHW;
  FHW.open(log_out_dir+job_name+"_random_BIC_ratios.csv");
  FHW<<"ratio,N\n";
  for (int i = 0 ; i < XY.size();i++){
    FHW<<to_string(XY[i])+","+to_string(CovN[i]) <<endl;
  }
  FHW.close();
  
  


  //---------------
  for (int n = 0 ; n < XY.size(); n++){
    if (min_x < 0 or XY[n] < min_x ){
      min_x = XY[n];
    }
    if (max_x < 0 or XY[n] > max_x ){
      max_x = XY[n];
    }
  }
  //-------------
  slice_ratio SC(min_x,max_x,400);
  for (int n = 0 ; n < XY.size();n++){

    if (XY[n] > 0.0){
      SC.insert(XY[n] );
    }
  }
  SC.set(pval_threshold);
  return SC;
}
