/**
 * @file model.cpp
 * @author Joey Azofeifa
 * @brief Contains \ref UNI (uniform), \ref EMG (full model), \ref NOISE,
 *  \ref component (full set of distros, priors), and \ref classifier (EM convergence
 * info and k number of components) classes.  
 * Note that the EMG class contains both the EMG of initiation and the L of the uniform of 
 * elongation.  Also that the UNI and NOISE are both uniform.  
 * @bug Contains lots of deadends -- as code clearly used to generate from the model and
 * other functions and classes (e.g. NORMAL) that are now missing.
 * @bug HEAVILY overloaded, should be refactored.
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#include "model.h"

#include <math.h>   
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <random>

#include <mpi.h>
#include "omp.h"

#include "load.h"
#include "template_matching.h"

//=============================================
//Helper functions
/**
 * @brief Standard Normal PDF (with a terrible name!)
 * @param x 
 * @return double 
 */
double IN(double x){ //Standard Normal PDF 
	return exp(-pow(x,2)*0.5)/sqrt(2*M_PI);
}
/**
 * @brief Standard Normal CDF (with a terrible name!)
 * 
 * @param x 
 * @return double 
 */
double IC(double x){ //Standard Normal CDF
	return 0.5*(1+erf(x/sqrt(2)));
}
/**
 * @brief Mills Ratio (R)
 * 
 * @param x 
 * @return double 
 */
double R(double x){ //Mills Ratio
	if (x > 4){
		return 1.0 / x;
	}
	double N = IC(x);
	double D = IN(x);
	if (D < pow(10,-15)){ //machine epsilon
		return 1.0 / pow(10,-15);
	}
	return exp(log(1. - N)-log(D));
}
/**
 * @brief a wrapper around isfinite()
 * WHY??
 * 
 * @param x 
 * @return true 
 * @return false 
 */
bool checkNumber(double x){
	if (isfinite(x)){
		return true;
	}
	return false;//wrapper for isfinite()
}
/**
 * @brief wrapper around log(x) that returns nINF if x <=0
 * 
 * @param x 
 * @return double 
 */
double LOG(double x){
	if (x <= 0){
		return nINF;
	}
	return log(x);//will return nINF if x <= 0
}


//=============================================
// Uniform Noise Class
NOISE::NOISE(){} //empty constructor

/**
 * @brief Construct a new NOISE::NOISE object
 * There is also an empty constructor.
 * @param A
 * @param B
 * @param W
 * @param PI
 */
NOISE::NOISE(double A, double B, double W, double PI){
	a=A;
	b=B;
	w=W;
	pi=PI;
}
/**
 * @brief NOISE density function
 * Note that NOISE is just a uniform.
 * 
 * @param x 
 * @param strand 
 * @return double 
 */
double NOISE::pdf(double x, int strand){
	if (strand == 1){
		return (w*pi) / abs(b-a);
	}
	return (w*(1-pi)) / abs(b-a);
}

//=============================================
//Uniform Class
UNI::UNI(){} //empty constructor
/**
 * @brief Construct a new UNI::UNI object
 * 
 * @param start  genomic position of start
 * @param stop   genomic position of stop
 * @param w_i 	 
 * @param strand 
 * @param POS   position 
 * @param Pi   Appears to also be strand info? Appears to overwrite strand? 
 * 
 * @bug Appears to set strand twice. First based on strand then to Pi.
 */
UNI::UNI(double start, double stop, double w_i, int strand, int POS, double Pi){
	a 		= start;
	b 		= stop;
	w 		= w_i;
	st 		= strand;
	pos 	= POS;
	if (st==1){
		pi=1;
	}else{
		pi=0;
	}
	//===================
	//this oversets the constraint that uniform must take either 
	//forward or reverse data points
	pi 		= Pi;
	//===================
	delta_a=0;
	delta_b=0;
	ri_forward=0, ri_reverse=0;
}
/**
 * @brief UNI density function
 * 
 * @param x 
 * @param strand 
 * @return double 
 */
double UNI::pdf(double x, int strand){
	double p;
	if (w==0){
		return 0;
	}

	if ( a<= x and x <=b){

		p= w / abs(b- a);
		p= p*pow(pi, max(0, strand) )*pow(1.-pi, max(0, -strand) );
		return p;
	}
	return 0;
}
/**
 * @brief Print the parameters of the UNI distro.
 * 
 * @return string 
 */
string UNI::print(){
	string text = ("U: " + to_string(a) + "," + to_string(b) 
	+ "," + to_string(w) + "," + to_string(pi));
	return text;
}	

//=============================================
//Exponentially Modified Gaussian class
EMG::EMG(){} //empty constructor
/**
 * @brief Construct a new EMG::EMG object
 * 
 * @param MU 
 * @param SI 
 * @param L 
 * @param W 
 * @param PI 
 */
EMG::EMG(double MU, double SI, double L, double W, double PI ){
	mu 	= MU;
	si 	= SI;
	l  	= L;
	w 	= W;
	pi 	= PI;
	prev_mu = 0 ;
	move_fp = 0;
}
/**
 * @brief Output parameters of this EMG instance.
 * 
 * @return string 
 */
string EMG::print(){
	string text 	= ("N: " + to_string(mu)+","+to_string(si)
		+ "," + to_string(l) + "," + to_string(w) + "," + to_string(pi) + "," + to_string(foot_print) );
	return text;
}
/**
 * @brief EMG density function
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double EMG::pdf(double z, int s ){
        if (w==0){
	  return 0.0;
	}
	if (s==1){
		z-=foot_print;
	}else{
		z+=foot_print;
	}
	double vl 		= (l/2.0)*(s*2*(mu-z) + l*pow(si,2));
	double p;
	if (vl > 100){ //potential for overflow, inaccuracies
		p 			= l*IN((z-mu)/si)*R(l*si - s*((z-mu)/si));
	}else{
		p 			= (l/2)*exp(vl)*erfc((s*(mu-z) + l*pow(si ,2) )/(sqrt(2)*si));
	}
	p     = (l/2)*exp(vl)*erfc((s*(mu-z) + l*pow(si ,2) )/(sqrt(2)*si));
	p     = p*w*pow(pi, max(0, s) )*pow(1-pi, max(0, -s) );
	if (p < pow(10,7) and not isnan(float(p)) ){
	  return p; 
	}
	return 0.0;
}
/**
 * @brief conditional expectation of Y given z_i
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double EMG::EY(double z, int s){
	if (s==1){
		z-=foot_print;
	}else{
		z+=foot_print;
	}
	
	return max(0. , s*(z-mu) - l*pow(si, 2) + (si / R(l*si - s*((z-mu)/si))));
}
/**
 * @brief conditional expectation of Y^2 given z_i
 * 
 * @param z 
 * @param s 
 * @return double 
 */
double EMG::EY2(double z, int s){
	if (s==1){
		z-=foot_print;
	}else{
		z+=foot_print;
	}
	return pow(l,2)*pow(si,4) + pow(si, 2)*(2*l*s*(mu-z)+1 ) + pow(mu-z,2) - ((si*(l*pow(si,2) + s*(mu-z)))/R(l*si - s*((z-mu)/si) )); 
}

//===============================================================================
//functions that help estimate uniform support bounds
/**
 * @brief Get the nearest position. 
 *  This is a helper function for estimating support bounds. 
 * @param data 
 * @param center 
 * @param dist 
 * @return int 
 */
int get_nearest_position(segment * data, double center, double dist){
	int i;

	if (dist < 0 ){
		i=0;
		while (i < (data->XN-1) and (data->X[0][i] -center) < dist){
			i++;
		}
	}else{
		i=data->XN-1;
		while (i >0 and (data->X[0][i] - center) > dist){
			i--;
		}
	}
	return i;
}
/**
 * @brief Get the sum of segment between j and k.
 * 
 * @param data 
 * @param j 
 * @param k 
 * @param st 
 * @return double 
 */
double get_sum(segment * data, int j, int k, int st){
	double S 	= 0;
	for (int i = j; i <k;i++){
		S+=data->X[st][i];
	}
	return S;
}
/**
 * @brief 
 * 
 * @param components 
 * @param data 
 * @param K 
 * @param N 
 */
void update_j_k( component * components,segment * data, int K, double N){
	for (int k = 0 ; k < K;k++){
		if (k > 0){
			components[k].reverse_neighbor 	= &components[k-1]; 
		}
		if (k+1 < K){
			components[k].forward_neighbor 	= &components[k+1];
		}
	}

	double center, w_thresh;

	for (int k = 0; k < K; k++){
		int j 	= k;
		//for the forward
		w_thresh= ( components[k].ALPHA_2 ) / (N + components[k].ALPHA_2*K*3 + K*3 );
		while ( j < K and components[j].forward_neighbor!=NULL and  components[j].forward_neighbor->bidir.w  < pow(10,-2)   ){
			j++;
		}
		double delta 				= (1.0 / components[k].bidir.l);
		components[k].forward.j 	= get_nearest_position(data, components[k].bidir.mu, delta  );
		
		if (j < K and components[j].forward_neighbor!=NULL  ){
			center 						= components[k].bidir.mu + delta;
			components[k].forward.k 	= get_nearest_position(data, center, components[j].forward_neighbor->bidir.mu -center  );
		}else{
			components[k].forward.k 	= get_nearest_position(data, components[k].bidir.mu,data->maxX -components[k].bidir.mu );			
		}
		if (components[k].forward.j >= components[k].forward.k ){
			components[k].forward.j 	= components[k].forward.k;
			components[k].forward.w 	= 0;
		}
		
		j 	= k;
		while (j >=0 and components[j].reverse_neighbor!=NULL and  components[j].reverse_neighbor->bidir.w   < pow(10,-2) ){
			j--;
		}

		if (K>j and j>=0 and components[j].reverse_neighbor!=NULL){
			center 						= components[k].bidir.mu-delta;
			components[k].reverse.j 	= get_nearest_position(data, center, components[j].reverse_neighbor->bidir.mu - center );
		}else{
			components[k].reverse.j 	= get_nearest_position(data, components[k].bidir.mu, data->minX-components[k].bidir.mu );			
		}
		components[k].reverse.k 		= get_nearest_position(data, components[k].bidir.mu, -delta);
		if (components[k].reverse.j >= components[k].reverse.k ){
			components[k].reverse.j 	= components[k].reverse.k;
			components[k].reverse.w 	= 0;
		}		
	}
}
/**
 * @brief 
 * 
 * @param components 
 * @param data 
 * @param K 
 */
void update_l(component * components, segment * data, int K){
	for (int k 	= 0; k < K; k++){
		//forward
		double left_SUM=0, right_SUM=get_sum(data,components[k].forward.j ,components[k].forward.k,1);
		double N 		= left_SUM+right_SUM;
		double null_vl 	= 0.05 / (data->maxX - data->minX  );
		double null_ll 	= N*LOG(1. / (data->maxX - data->minX  ));
		double null_BIC = -2*null_ll + LOG(N);
		double vl 		= 0, w=0;
		double mod_ll, mod_BIC;
		double prev_prev=0, prev=0, current=0;
		double BIC_best = 0;
		int arg_l 	= components[k].forward.k;

		for (int l = components[k].forward.j; l < components[k].forward.k; l++ ){
			left_SUM+=data->X[1][l];
			right_SUM-=data->X[1][l];
			vl 		= 1.0/(data->X[0][l]-data->X[0][components[k].forward.j]);
			w 		= left_SUM/(N);
			mod_ll 	= LOG(w*vl)*left_SUM + LOG(null_vl)*right_SUM ;
			mod_BIC = -2*mod_ll + 5*LOG(N);
			current = null_BIC/mod_BIC;
			if (prev > current and prev > prev_prev and prev > 1.0){
				if (prev > BIC_best){
					BIC_best=prev, arg_l=l;
				}
			}
			prev_prev=prev;
			prev 	= current;			
		}
		components[k].forward.b 	= data->X[0][arg_l];
		//reverse
		arg_l 	= components[k].reverse.j;
		left_SUM=0, right_SUM=get_sum(data,components[k].reverse.j,components[k].reverse.k,2 );
		N 		= left_SUM+right_SUM;
		null_ll 	= N*LOG(1. / (data->maxX - data->minX  ));
		null_BIC = -2*null_ll + LOG(N);
		prev_prev=0, prev=0, current=0,BIC_best = 0;
		for (int l = components[k].reverse.j; l < components[k].reverse.k; l++ ){
			left_SUM+=data->X[2][l];
			right_SUM-=data->X[2][l];
			vl 		= 1.0/(data->X[0][components[k].reverse.k] - data->X[0][l]);
			w 		= right_SUM/(N);
			mod_ll 	= LOG(null_vl)*left_SUM + LOG(w*vl)*right_SUM ;
			mod_BIC = -2*mod_ll + 5*LOG(N);
			current = null_BIC/mod_BIC;
			if (prev > current and prev > prev_prev and prev > 1.0){
				if (prev > BIC_best){	
					BIC_best=prev, arg_l=l;
				}
			}
			prev_prev=prev;
			prev 	= current;
		}
		components[k].reverse.a 	= data->X[0][arg_l];
	}
}

//===============================================================================
//components 
/**
 * @brief Construct a new component::component object
 * Empty constructor.  Sets pointers to NULL.
 * 
 */
component::component(){//empty constructor
	foot_print 			= 0;
	forward_neighbor 	= NULL;
	reverse_neighbor 	= NULL;
} 
/**
 * @brief set the hyperparameters
 * 
 * @param s_0 
 * @param s_1 
 * @param l_0 
 * @param l_1 
 * @param w_0 
 * @param strand_0 
 * @param N 
 * @param K 
 */
void component::set_priors(double s_0, double s_1, 
	double l_0, double l_1, double w_0,double strand_0, double N, int K) { 
	//============================
	//for sigma
	alpha_0 	= 20.46;
	beta_0 		= 10.6;
	//============================
	//for lambda
	alpha_1 	= 20.823;
	beta_1 		= 0.5;
	//==============================
	//for initial length of Uniforms
	alpha_2 	= 1.297;
	beta_2 		= 8260;

	//*****************************************************
	//Priors on parameters for MAP Estimate
	ALPHA_0 = s_0, BETA_0 =s_1; //for sigma
	ALPHA_1 = l_0, BETA_1 =l_1; //for lambda
	ALPHA_2 = w_0; //for weights, Dirichlet
	ALPHA_3 = strand_0; //for strand probs
	//bidir.w 	= (r + ALPHA_2) / (N + ALPHA_2*K*3 + K*3) ;
		
	w_thresh= ( ALPHA_2 ) / (N + ALPHA_2*K*3 + K*3 );
}
/**
 * @brief randomly seed the sigma, pi, lambda ws etc...
 * 
 * @param mu 
 * @param data 
 * @param K 
 * @param scale 
 * @param noise_w 
 * @param termination 
 * @param fp 
 * @param forward_bound 
 * @param reverse_bound 
 */
void component::initialize_bounds(double mu, segment * data , int K, double scale, double noise_w, 
	double termination, double fp, double forward_bound, double reverse_bound){//random seeds...
	foot_print 	= fp;
	EXIT=false;
	if (noise_w>0){
		noise 	= NOISE(data->minX, data->maxX, 
			noise_w, 0.5);
		type 	= 0; 
	}else{

		int complexity=1;
		if (data->strand == "."){
			complexity=3;
		}else{
			complexity=2;
		}

		//====================================
		random_device rd;
		mt19937 mt(rd());
		
		double sigma,lambda, pi_EMG, w_EMG  ;	
		double b_forward,  w_forward;
		double a_reverse,  w_reverse;

		//====================================
		//start sampling
		//for the bidirectional/EMG component
		uniform_real_distribution<double> dist_lambda_2(1, 250);
		uniform_real_distribution<double> dist_sigma_2(1, 250);
		gamma_distribution<double> dist_lengths(1,( (data->maxX-data->minX)/(K)));
		
		sigma 				= dist_sigma_2(mt)/scale;
		lambda 				= scale/dist_lambda_2(mt) ;
		double dist 		=  (1.0/lambda);
		int j 				= get_nearest_position(data, mu, dist);
		int k 				= get_nearest_position(data, mu, forward_bound-mu);
		double pi 			= 0.5;
		if (data->strand == "+"){
			pi 				= 1.0;
		}else if (data->strand=="-") {
			pi 				= 0.;
		}
		forward 			= UNI(mu+(1.0/lambda), data->maxX, 1.0 / (complexity*K), 1, j, pi);
		forward.j=j, forward.k=k;
		

		bidir 				= EMG(mu, sigma, lambda, 1.0 / (complexity*K), 0.5);
		bidir.foot_print 	= fp;
		dist 				=  -(1.0/lambda);
		j 					= get_nearest_position(  data, mu, dist);
		k 					= get_nearest_position(data, mu, reverse_bound-mu);
		reverse 			= UNI(data->minX, mu-(1.0/lambda),  1.0 / (complexity*K) , -1, j,1-pi);
		reverse.j=k, reverse.k=j;
		termination 		= (termination > 1);
		
		type 		= 1;
		
	}
} 
/**
 * @brief Print out parameters of component.
 * 
 */
void component::print(){
	if (type==1){
		string text 	= bidir.print()+ "\n";
		text+=forward.print()+ "\n";
		text+=reverse.print() + "\n";
		cout<<text;
	}else{
		cout<<"NOISE: " << noise.w<<"," <<noise.pi<<endl;
	}
}
/**
 * @brief compute the density at x_i,s_i
 * 
 * @param x 
 * @param st 
 * @return double 
 */
double component::evaluate(double x, int st){
	if (type ==0){ //this is the uniform noise component
		return noise.pdf(x, st);
	}
	if (st==1){
		bidir.ri_forward 	= bidir.pdf(x, st);
		forward.ri_forward 	= forward.pdf(x, st);
		reverse.ri_forward 	= reverse.pdf(x,st);
		return bidir.ri_forward + forward.ri_forward + reverse.ri_forward;
	}
	bidir.ri_reverse 	= bidir.pdf(x, st);
	reverse.ri_reverse 	= reverse.pdf(x, st);
	forward.ri_reverse 	= forward.pdf(x, st);
	return bidir.ri_reverse + reverse.ri_reverse + forward.ri_reverse;
}

/**
 * @brief compute the conditional expectations and add to running total
 *  This is equation #9 in the Azofeifa 2018 paper.
 * 
 * @param x 
 * @param y 
 * @param st 
 * @param normalize 
 */
void component::add_stats(double x, double y, int st, double normalize){
	if (type==0){//noise component
		if (st==1){
			noise.r_forward+=(y*noise.ri_forward/normalize);
			noise.ri_forward=0;
		}else{
			noise.r_reverse+=(y*noise.ri_reverse/normalize);
			noise.ri_reverse=0;
		}

	}else{
		double vl, vl2, vl3;
		if (st==1){
			vl 	= bidir.ri_forward / normalize;
			vl2 = forward.ri_forward/normalize;
			vl3 = reverse.ri_forward/normalize;
			bidir.ri_forward=0, forward.ri_forward=0;
			bidir.r_forward+=(vl*y);
			forward.r_forward+=(vl2*y);
			reverse.r_forward+=(vl3*y);
		
		}else{
			vl 	= bidir.ri_reverse / normalize;
			vl2 = reverse.ri_reverse / normalize;
			vl3 = forward.ri_reverse / normalize;
			bidir.ri_reverse=0, reverse.ri_reverse=0;
			bidir.r_reverse+=(vl*y);

			reverse.r_reverse+=(vl2*y);
			forward.r_reverse+=(vl3*y);
		}
		//now adding all the conditional expectations for the convolution
		if (vl > 0 and y > 0){
			double current_EY 	= bidir.EY(x, st);
			double current_EY2 	= bidir.EY2(x, st);
			double current_EX 	= x-(st*current_EY)-bidir.foot_print*st;
			//	self.C+=max( ((z-self.mu) -E_Y) *r,0)
			// 	self.C+=max((-(z-self.mu) -E_Y)   *r ,0)
			bidir.C+=max((st*(x-bidir.mu) - current_EY  )*vl*y,0.0);
			bidir.ey+=current_EY*vl*y;
			bidir.ex+=current_EX*vl*y;
			bidir.ex2+=(pow(current_EX,2) + current_EY2 - pow(current_EY,2))*vl*y;	
		}
	}
}
/**
 * @brief reset running totals and responsibility terms
 * 
 */
void component::reset(){
	if (type){
		bidir.C=0;
		bidir.ey=0, bidir.ex=0, bidir.ex2=0, bidir.r_reverse=0, bidir.r_forward=0;
		bidir.ri_forward=0, forward.ri_forward=0, forward.ri_reverse=0;
		bidir.ri_reverse=0, reverse.ri_reverse=0, reverse.ri_forward=0;
		forward.r_forward=0, forward.r_reverse=0, reverse.r_reverse=0, reverse.r_forward=0;
		forward.delta_a=0, forward.delta_b=0, reverse.delta_a=0, reverse.delta_b=0;
	}else{
		noise.r_forward=0,noise.r_reverse=0;
		noise.ri_reverse=0,noise.ri_forward=0 ;
		
	}
}
/**
 * @brief used for large responsibility normalization term
 * 
 * @return double 
 */
double component::get_all_repo(){
	if (type==1){
		return bidir.r_forward+bidir.r_reverse+forward.r_forward+reverse.r_reverse;
	}
	return 0.;
}
/**
 * @brief take all the nice sample means and variances
 *  This is equation #10 in Azofeifa 2018.
 * 
 * @param N 
 * @param K 
 */
void component::update_parameters(double N, int K){
	if (type==1){
		//first for the bidirectional
		double r 	= bidir.r_forward + bidir.r_reverse;
		bidir.pi 	= (bidir.r_forward + ALPHA_3) / (r + ALPHA_3*2);
		bidir.w 	= (r + ALPHA_2) / (N + ALPHA_2*K*3 + K*3) ;
		bidir.mu 	= bidir.ex / (r+0.001);

		
		bidir.si 	= pow(abs((1. /(r + 3 + ALPHA_0 ))*(bidir.ex2-2*bidir.mu*bidir.ex + 
			r*pow(bidir.mu,2) + 2*BETA_0  )), 0.5) ;
		if ((r / N) < pow(10,-5) ){ EXIT 	= true; }
		bidir.l 	= min((r+ALPHA_1) / (bidir.ey + BETA_1), 5.);
		bidir.l 	= max(bidir.l, 0.05);
		if (abs(bidir.mu-bidir.prev_mu)< 0.01 ){
			bidir.move_fp 	= true;
		}
		else{
			bidir.prev_mu 	= bidir.mu;
		}
		if (bidir.move_fp){

			bidir.foot_print 	= min( max(bidir.C / (r+0.1),0.0) , 2.5);
		}
		//bidir.foot_print 	= 0.0;
		//now for the forward and reverse strand elongation components
		forward.w 	= (forward.r_forward + ALPHA_2) / (N+ ALPHA_2*K*3 + K*3);
		reverse.w 	= (reverse.r_reverse + ALPHA_2) / (N+ ALPHA_2*K*3 + K*3);
		forward.a 	= bidir.mu  , reverse.b=bidir.mu ;

		//update PIS, this is obviously overwritten if we start the EM seeder with 0/1
		forward.pi 	= (forward.r_forward + 1) / (forward.r_forward + forward.r_reverse+2);
		reverse.pi 	= (reverse.r_forward + 1)/ (reverse.r_forward + reverse.r_reverse+2);
		if (bidir.w==0){
			forward.w 	= 0;
			reverse.w 	= 0;
		}
	}
}

//=========================================================
//sorting functions for the classifier class 
/**
 * @brief Bubble Sort vector of components
 * 
 * @param components 
 * @param K 
 */
void sort_components(component components[], int K){
	bool sorted=true;
	while (sorted){
		sorted=false;
		for (int k = 0; k < K-1; k++){
			if (components[k].bidir.mu > components[k+1].bidir.mu){
				component cp 		= components[k];
				components[k] 		= components[k+1];
				components[k+1] 	= cp;
				sorted 				= true;
			}
		}
	}
}
/**
 * @brief Bubble Sort mu ?
 * 
 * @param X 
 * @return vector<vector<double>> 
 */
vector<vector<double>> sort_mus(vector<vector<double>> X){
	bool changed=true;

	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i][0] > X[i+1][0]){ //sort by starting position
				vector<double> copy 	= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}
/**
 * @brief Bubble sort the data values.
 * 
 * @param X 
 * @param N 
 */
void sort_vector(double X[], int N){
	bool sorted=true;
	while (sorted){
		sorted=false;
		for (int i = 0; i < N-1; i++){
			if (X[i] > X[i+1]){
				double cp = X[i];
				X[i] 		= X[i+1];
				X[i+1] 		= cp;
				sorted 		= true;
			}
		}
	}
}


//=========================================================
//classifier class 
//(most of these constructors are deprecated)
/**
 * @brief Construct a new classifier::classifier object
 * 
 * @param k 
 * @param ct 
 * @param mi 
 * @param nm 
 * @param R_MU 
 * @param alpha_0 
 * @param beta_0 
 * @param alpha_1 
 * @param beta_1 
 * @param alpha_2 
 * @param alpha_3 
 * @param fp 
 */
classifier::classifier(int k, double ct, int mi, double nm,
	double R_MU, double alpha_0, double beta_0,
	double alpha_1, double beta_1, double alpha_2,double alpha_3, double fp){
	foot_print 				= fp;
	K 						= k ;
	seed 					= true;
	convergence_threshold 	= ct;
	max_iterations 			= mi;
	noise_max 				= nm;
	p 						= 0.8;
	last_diff 				= 0;
	r_mu 					= R_MU;

	//=============================
	//hyperparameters
	ALPHA_0=alpha_0, BETA_0=beta_0, ALPHA_1=alpha_1, BETA_1=beta_1;
	ALPHA_2=alpha_2, ALPHA_3=alpha_3;

	move_l = true;
}
/**
 * @brief Construct a new classifier::classifier object
 * 
 * @param k 
 * @param ct 
 * @param mi 
 * @param nm 
 * @param R_MU 
 * @param alpha_0 
 * @param beta_0 
 * @param alpha_1 
 * @param beta_1 
 * @param alpha_2 
 * @param alpha_3 
 * @param MOVE 
 * @param fp 
 */
classifier::classifier(int k, double ct, int mi, double nm,
	double R_MU, double alpha_0, double beta_0,
	double alpha_1, double beta_1, double alpha_2,double alpha_3, bool MOVE, double fp){
	foot_print 				= fp;
	K 						= k ;
	seed 					= true;
	convergence_threshold 	= ct;
	max_iterations 			= mi;
	noise_max 				= nm;
	p 						= 0.8;
	last_diff 				= 0;
	r_mu 					= R_MU;

	//=============================
	//hyperparameters
	ALPHA_0=alpha_0, BETA_0=beta_0, ALPHA_1=alpha_1, BETA_1=beta_1;
	ALPHA_2=alpha_2, ALPHA_3=alpha_3;
	move_l 	= MOVE;
}
/**
 * @brief Construct a new classifier::classifier object
 * 
 * @param ct 
 * @param mi 
 * @param nm 
 * @param R_MU 
 * @param alpha_0 
 * @param beta_0 
 * @param alpha_1 
 * @param beta_1 
 * @param alpha_2 
 * @param alpha_3 
 * @param IP 
 * @param fp 
 */
classifier::classifier(double ct, int mi, double nm,
	double R_MU, double alpha_0, double beta_0,
	double alpha_1, double beta_1, double alpha_2,double alpha_3, vector<vector<double>> IP, double fp){
	
	foot_print 				= fp;
	seed 					= true;
	convergence_threshold 	= ct;
	max_iterations 			= mi;
	noise_max 				= nm;
	p 						= 0.8;
	last_diff 				= 0;
	r_mu 					= R_MU;

	//=============================
	//hyperparameters
	ALPHA_0=alpha_0, BETA_0=beta_0, ALPHA_1=alpha_1, BETA_1=beta_1;
	ALPHA_2=alpha_2, ALPHA_3=alpha_3;
	init_parameters 		= IP;
}

classifier::classifier(){}; 

/**
 * @brief This is the core EM algorithm.
 * 
 * @param data 
 * @param mu_seeds 
 * @param topology 
 * @param elon_move 
 * @return int 
 */
int classifier::fit2(segment * data, vector<double> mu_seeds, int topology,
	 int elon_move ){

	//printf("topology: %d elon_move %d K %d \n", topology, elon_move, K);
	//  printf("K %d \n", K);
	//=========================================================================
	//compute just a uniform model...no need for the EM
	if (K == 0){
		ll 			= 0;
		double 	l 	= (data->maxX-data->minX); //~length
		double pos 	= 0;
		double neg 	= 0;
		// Sum per strand
		for (int i = 0; i < data->XN; i++){
			pos+=data->X[1][i];
			neg+=data->X[2][i];
		}
		double pi 	= pos / (pos + neg); // strand bias?
		// Calculate MLE = -n log (pi/l) where pi is strand (1 -> +; -1 -> -)
		for (int i = 0; i < data->XN; i++){
			if (pi > 0){
				ll+=log(pi / l)*data->X[1][i];
			}
			if (pi < 1){
				ll+=log((1-pi) / l)*data->X[2][i];		
			}
		}
    // Sets the components part of the classifier to a new component, but this is an array?
		components 	= new component[1];
	  // printf("\t l: %9.6f pos: %9.6f neg: %9.6f pi: %9.6f ll: %9.6f \n", l, pos, neg, pi, ll);
		return 1;
	}
	// These are the seeds for the initial location.
	random_device rd;
	mt19937 mt(rd());
	
	int add 	= noise_max>0;
	components 	= new component[K+add];
	//===========================================================================
	//initialize(1) components with user defined hyperparameters
	for (int k = 0; k < K; k++){
		components[k].set_priors(ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3,data->N, K);
	}

	//===========================================================================
	//random seeding, initialize(2), center of pausing components
	int i 	= 0;
	double mu;
	double mus[K];
	for (int k = 0; k < K; k++){
		if (mu_seeds.size()>0  ){
			i 	= sample_centers(mu_seeds ,  p);
			mu 	= mu_seeds[i];
			if (r_mu > 0){
				normal_distribution<double> dist_r_mu(mu, r_mu);
				mu 		= dist_r_mu(mt); // mt is the random number
			}
		}else{
			normal_distribution<double> dist_MU((data->minX+data->maxX)/2., r_mu);
			mu 			= dist_MU(mt); // mt is the random number
		}
		
		mus[k] 	= mu;
		if (mu_seeds.size() > 0  ){
			mu_seeds.erase (mu_seeds.begin()+i);	
		}
	}
	sort_vector(mus, K);
	for (int k = 0; k < K;k++){ //random seeding, initialize(3) other parameters
		components[k].initialize_bounds(mus[k], 
			data, K, data->SCALE , 0., topology,foot_print, data->maxX, data->maxX);
		
	}
	sort_components(components, K);
	for (int k = 0; k < K; k++){
		if (k > 0){
			components[k].reverse_neighbor 	= &components[k-1];
		}else{
			components[k].reverse_neighbor 	= NULL;
		}
		if (k+1 < K){
			components[k].forward_neighbor 	= &components[k+1];
		}else{
			components[k].forward_neighbor 	= NULL;
		}
	}
	// I think this is the noise component ...
	if (add){
		components[K].initialize_bounds(0., data, 0., 0. , noise_max, pi, foot_print, data->minX, data->maxX);
	}

	//===========================================================================
	int t 			= 0; //EM loop ticker
	double prevll 	= nINF; //previous iterations log likelihood
	converged 		= false; //has the EM converged?
	int u 			= 0; //elongation movement ticker
	double norm_forward, norm_reverse,N; //helper variables
	//printf("------------------------------------------\n");
	while (t < max_iterations && not converged){
		//======================================================
		//reset old sufficient statistics
		for (int k=0; k < K+add; k++){
			// components[k].print();
			components[k].reset();
			if (components[k].EXIT){
				converged=false, ll=nINF;
				return 0;
			}
		       
		}
		
		//======================================================
		//E-step, grab all the stats and responsibilities
		ll 	= 0;
		// i -> |D| (Azofeifa 2017 pseudocode) 
		for (int i =0; i < data->XN;i++){
			norm_forward=0;
			norm_reverse=0;
			
			// Equation 7 in Azofeifa 2017: calculate r_i^k
			for (int k=0; k < K+add; k++){ //computing the responsibility terms
				if (data->X[1][i]){//if there is actually data point here...
					norm_forward+=components[k].evaluate(data->X[0][i],1);
				}
				if (data->X[2][i]){//if there is actually data point here...
					norm_reverse+=components[k].evaluate(data->X[0][i],-1);
				}
			}
			if (norm_forward > 0){
				ll+=LOG(norm_forward)*data->X[1][i];
			}
			if (norm_reverse > 0){
				ll+=LOG(norm_reverse)*data->X[2][i];
			}
			
			//now we need to add the sufficient statistics, need to compute expectations
			// Equation 9 in Azofeifa 2017
			for (int k=0; k < K+add; k++){
				if (norm_forward){
					components[k].add_stats(data->X[0][i], data->X[1][i], 1, norm_forward);
				}
				if (norm_reverse){
					components[k].add_stats(data->X[0][i], data->X[2][i], -1, norm_reverse);
				}
			}
		}

		//======================================================
		//M-step, Equation 10 in Azofeifa 2017, Theta_k^(t+1)
		N=0; //get normalizing constant
		for (int k = 0; k < K+add; k++){
			N+=(components[k].get_all_repo());
		}
		
		for (int k = 0; k < K+add; k++){
			components[k].update_parameters(N, K);
		}
		
		if (abs(ll-prevll)<convergence_threshold){
			converged=true;
		}
		if (not isfinite(ll)){
			ll 	= nINF;
			return 0;	
		}
		//======================================================
		//should we try to move the uniform component? e.g. change L
		if (u > 200 ){
			sort_components(components, K);
			//check_mu_positions(components, K);
			if (elon_move){
				update_j_k(components,data, K, N);
				update_l(components,  data, K);
			}
			u 	= 0;
		}

		u++;
		t++;
		prevll=ll;
	}

	return 1;
}
