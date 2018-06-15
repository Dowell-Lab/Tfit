#ifdef USING_ICC
#include <mathimf.h>
#else
#include <math.h>   
#endif
#include <limits>
#include "model_single.h"
#include "load.h"
#include "template_matching.h"
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <random>
#ifdef USING_ICC
#include <mathimf.h>
#else
#include <math.h>   
#endif
#include "model.h"
#include "template_matching.h"
#ifdef USING_ICC
#include <aligned_new>
#endif



/** Default constructor for the NLR class. Use this in conjunction with init().
 */
NLR::NLR(){}

/** Resets internal weights to 0. Not used in the current Tfit codebase.
 */
void NLR::resetSS(){
	EX=0, EX2=0, EY=0, EPI=0, WE=0, WF=0, WR=0,EPIN=0 ;
}

/** Initializes the nonlinear regression classifier.
 * @param data The input segment on which to regress.
 * @param MU Mean parameter for the desired distribution.
 * @param K Scaling factor for the window on which to regress.
 * @param foot_print foot_print parameter for internal EMG model.
 * @param Alpha_1 Used to compute si during set_new_parameters updates.
 * @param Beta_1 Unused
 * @param Alpha_2 Unused
 * @param Beta_2 Unused
 * @param Alpha_3 Unused
 */
void NLR::init(segment * data, double MU , int K, double foot_print, 
	double Alpha_1, double Beta_1, double Alpha_2, double Beta_2, double Alpha_3 ){
	alpha_1=Alpha_1, alpha_2=Alpha_2, alpha_3=Alpha_3, beta_1=Beta_1, beta_2=Beta_2;
	double scale 	= data->SCALE;
	random_device rd;
	mt19937 mt(rd());
	
	uniform_real_distribution<double> dist_lambda_2(1, 500);
	uniform_real_distribution<double> dist_sigma_2(1, 50);
	mu 					= MU;
	si 					= dist_sigma_2(mt)/scale;
	l 					= scale/dist_lambda_2(mt) ;
	pi 					= 0.5;	
	wn = 1.0 / (3*K), wl= 1.0 / (3*K), wr= 1.0 / (3*K);
	bidir 				= EMG(MU, si , l , 1.0 / (3*K), pi);
	bidir.foot_print 	= foot_print;

	forward 			= UNI(MU+(1.0/l ), data->maxX, 1.0 / (3*K), 1, int(data->XN), 0.5);
	reverse 			= UNI(data->minX, MU-(1.0/l ), 1.0 / (3*K), -1, 0,0.5);

}

/** Returns the probability density at the specified point.
 * This function enables each NLR instance to function as a probability distribution.
 * @param x Point at which to get a probability density.
 * @return Probability density sample at x.
 */
double NLR::pdf(double x){
	return bidir.pdf(x,1) + bidir.pdf(x,-1) + forward.pdf(x,1) + forward.pdf(x,-1) + reverse.pdf(x,1) + reverse.pdf(x,-1);
}

/** 
 */
double NLR::addSS(double x, double y, double norm){
	double re, rf,rr, epi, ex, ey, ey2;
	re 	= (bidir.pdf(x,1) + bidir.pdf(x,-1))/norm;
	rf 	= (forward.pdf(x,1) + forward.pdf(x,-1))/norm;
	rr 	= (reverse.pdf(x,1) + reverse.pdf(x,-1))/norm;
	WE+=re*y;
	WF+=rf*y;
	WR+=rr*y;

	epi 	= (bidir.pdf(x,1)  ) / ((bidir.pdf(x,1) )+ (bidir.pdf(x,-1) ) );

	ey 		= max(epi*bidir.EY(x,1)+(1.-epi)*bidir.EY(x,-1),0.);

	ey2 	= epi*bidir.EY2(x,1)+(1.-epi)*bidir.EY2(x,-1);
	ex 		= (x-(ey*epi)  + (ey *(1.-epi)	));
	
	ex 		+= (bidir.foot_print*(1-epi) -bidir.foot_print*epi  );
	
	EPI 	+= epi*y*re;
	EPIN 	+= y*re;
	EY 		+= ey*y*re;
	EX 		+= ex*y*re;
	EX2 	+= max((pow(ex,2) + ey2 - pow(ey,2))*y,0.)*re;
	return  re*y + rf*y + rr*y;
}

/** Updates all parameters within the object given a value N. This function is used internally during regression.
 */
void NLR::set_new_parameters(double N){
	mu 	= EX / WE;
	si 	= sqrt((EX2 - 2*mu*EX + (pow(mu,2)*WE) + 2*beta_1 ) / (WE  + 3 + alpha_1)  );
	l 	= min(WE /EY, 4.0);
	pi 	= (EPI +alpha_3)  / (EPIN + 2*alpha_3);

	wn 	= WE / N;
	wl 	= WR / N;
	wr 	= WF / N;

	bidir.mu 	= mu, bidir.si = si, bidir.l = l, bidir.pi = pi, bidir.w=wn;
	forward.w 	= wr, forward.a = mu + (1.0 / l);
	reverse.w 	= wl, reverse.b = mu - (1.0 / l);

}

/** Prints all fields within the object.
 */
void NLR::print(){
	printf("mu: %f, si: %f,l: %f,pi: %f,w: %f, wf: %f, wr: %f  \n", mu, si, l, pi, wn, wl, wr );
}

/** Full constructor for classifier_single.
 * @param ct Convergence threshold
 * @param mi Maximum iterations
 * @param k K
 * @param TY Type of classification (unused within classifier_single).
 * @param SC Scale (unused within classifier_single).
 * @param Alpha_1 alpha1 model parameter (unused within classifier_single).
 * @param Beta_1 beta1 model parameter (unused within classifier_single).
 * @param Alpha_2 alpha2 model paramter (unused within classifier_single).
 * @param Beta_2 beta2 model parameter (unused within classifier_single).
 * @param Alpha_3 alpha3 model parameter (unused within classifier_single).
 */
classifier_single::classifier_single(double ct, int mi, int k, int TY, double SC,
	double Alpha_1, double Beta_1, double Alpha_2, double Beta_2, double Alpha_3){
	alpha_1=Alpha_1, alpha_2=Alpha_2, alpha_3=Alpha_3, beta_1=Beta_1, beta_2=Beta_2;

	covergence_threshold=ct, max_iterations=mi;
	K = k, type = TY;
	scale 	= SC;
}

/** Default constructor for classifier_single. Does nothing.
 */
classifier_single::classifier_single(){}

/** Fits the model to a single segment of input data.
 * Once fit, the classifier can then be useful in determining higher order model characteristics over an entire
 * dataset rather than a single segment within it.
 * 
 * @param data Segment on which to train.
 * @param foot_print Footprint parameter for the nonlinear regression components used internally.
 * @return log likelihood of the fit. A log likelihood closer to 0 should indicate a better fit.
 */
double classifier_single::fit(segment * data,double foot_print){

	double mu;
	int i;
	double p 	= 0.5;
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> dist_uni(data->minX,data->maxX);
	
	vector<double> mu_seeds;
	for (int i = 0; i < data->parameters.size(); i++){
		mu_seeds.push_back(data->parameters[i][0]);
	}
	components 	= new NLR[K];
	for (int k = 0 ; k < K; k++){
		if (mu_seeds.size()>0){
			i 	= sample_centers(mu_seeds ,  p);
			mu 	= mu_seeds[i];
		}else{
			mu 			= dist_uni(mt);
		}
		components[k].init(data, mu, K, foot_print, 
			alpha_1, beta_1, alpha_2, beta_2, alpha_3);
		components[k].fp 	= foot_print;
	}

	int t 			= 0;
	bool converged  = false;


	double x,y, norm, NNN;
	double prev_ll 	= 0;
	while (t < max_iterations and not converged){
		NNN 	= 0;
		ll 		= 0;
		for (int k = 0; k < K; k++){
			components[k].resetSS();
		}
		//e-step
		for (int i = 0; i < data->XN; i++){
			x=data->X[0][i],y=data->X[1][i];
			norm 	= 0;
			for (int k = 0; k < K;k++){
				norm+=components[k].pdf(x);
			}
			ll+=log(norm)*y;
			for (int k = 0; k < K; k++){
				NNN+=components[k].addSS(x,y, norm);
			}
		}
		//m-step
		for (int k = 0; k < K; k++){
			components[k].set_new_parameters(NNN);
		}
		if (abs(prev_ll - ll) < covergence_threshold){
			converged=true;
		}
		prev_ll=ll;

		t++;
	}


	return ll;
}


