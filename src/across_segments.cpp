#include <mpi.h>
#include "load.h"
#include "across_segments.h"
#include "model.h"
#include "template_matching.h"
#include "model_single.h"
#include <iostream>
#include <fstream>
#include <map>
#include <time.h>
#include "omp.h"
#include "read_in_parameters.h"
#include "error_stdo_logging.h"
using namespace std;


/** Recursively attempts to return a valid output filename.
 * This function is used so that Tfit doesn't overwrite the results returned from prior runs.
 * @param FILE Base filename
 * @param i initial iterator parameter with which to generate a unique filename. 
 */
string check_file(string FILE, int i){ //don't want to write over an existing file
	string template_file 	= FILE + to_string(i);
	ifstream FH(template_file);
	if (FH){
		FH.close();
		return	check_file(FILE, i+1);
	}
	return template_file;
}

/** Generates a structure containing representations of various command line parameters for use in classification utilizing a depreciated params object.
 * Use make_classifier_struct_free_model_pwrapper instead.
 * @depreciated
 * @param P Params object encapsulating command line arguments.
 * @param data Input reads.
 * @return Map between integers (in this case just indices) and filled out classifier structures.
 */
map<int, vector<classifier> > make_classifier_struct_free_model(params * P, segment * data){

	int min_k 		= stoi(P->p["-minK"]);
	int max_k 		= stoi(P->p["-maxK"]);
	int rounds 		= stoi(P->p["-rounds"]);
	int BDS 		= int(data->centers.size());
	map<int, vector<classifier> > A;
	double noise   	= abs(data->stop - data->start)*stod(P->p["-max_noise"])/(data->fN + data->rN); //expection of binomal(0.05, n) where 0.05 probability mapping noise

	A[0].push_back( classifier(0, stod(P->p["-ct"]), stoi(P->p["-mi"]), noise, 
			stod(P->p["-r_mu"]), stod(P->p["-ALPHA_0"]), stod(P->p["-BETA_0"]), stod(P->p["-ALPHA_1"]), 
			stod(P->p["-BETA_1"]), stod(P->p["-ALPHA_2"]) , stod(P->p["-ALPHA_3"]),0 ));
	double scale 	= stod(P->p["-ns"]);
	if(P->bidir){
		min_k 	= data->counts;
		max_k 	= data->counts;
	}
	for (int k =min_k; k <= max_k; k++){
		for (int r = 0; r < rounds; r++){
			A[k].push_back(classifier(k, stod(P->p["-ct"]), stoi(P->p["-mi"]), stod(P->p["-max_noise"]), 
			stod(P->p["-r_mu"]), stod(P->p["-ALPHA_0"]), stod(P->p["-BETA_0"]), stod(P->p["-ALPHA_1"]), 
			stod(P->p["-BETA_1"]), stod(P->p["-ALPHA_2"]) , stod(P->p["-ALPHA_3"]),0 ));
		}
	
	}
	
	return A;
}

/** Generates a structure containing representations of various command line parameters for use in classification utilizing a ParamWrapper.
 * @param pw ParamWrapper object encapsulating command line arguments.
 * @param data Input reads.
 * @return Map between integers (in this case just indices) and filled out classifier structures.
 */
map<int, vector<classifier> > make_classifier_struct_free_model_pwrapper(ParamWrapper *pw, segment * data){

	int min_k 		= pw->mink;
	int max_k 		= pw->maxk;
	int rounds 		= pw->rounds;
	int BDS 		= int(data->centers.size());
	map<int, vector<classifier> > A;
	double noise   	= abs(data->stop - data->start)*pw->maxNoise/(data->fN + data->rN); //expection of binomal(0.05, n) where 0.05 probability mapping noise

	A[0].push_back( classifier(0, pw->ct, pw->mi, noise, 
			pw->r_mu, pw->alpha0, pw->beta0, pw->alpha1, 
			pw->beta1, pw->alpha2, pw->alpha3, 0));
	double scale 	= pw->ns;
	if(pw->bidir){
		min_k 	= data->counts;
		max_k 	= data->counts;
	}
	for (int k =min_k; k <= max_k; k++){
		for (int r = 0; r < rounds; r++){
			A[k].push_back(classifier(k, pw->ct, pw->mi, noise, pw->r_mu, pw->alpha0, pw->beta0, pw->alpha1, pw->beta1, pw->alpha2, pw->alpha3, 0));
		}
	
	}
	
	return A;
}

/** Full constructor for simple_c_free_mode class.
 * @param FOUND Whether or not the learning process has converged (used as a tag within the object).
 * @param ll Log likelihood of the converged result.
 * @param C Set of components specifying various bidirectional model parameters.
 * @param K Index or ID?
 * @param data Input segment set.
 * @param i Unused.
 * @param forward_N Number of forward strand elements (base pairs; used to set values with respect to call dimensions).
 * @param reverse_N Number of reverse strand elements (see notes for forward_N).
 */
simple_c_free_mode::simple_c_free_mode(bool FOUND, double ll, 
	component C, int K, segment * data, int i, double forward_N, double reverse_N){
	SS[0]=ll, SS[1]=forward_N, SS[2]=reverse_N;
	ID[0]=data->ID, ID[1]=data->start, ID[2]=data->stop, ID[3]=K;
	if (FOUND){
		ID[4]=1; 
	}else{
		ID[4]=0;//didn't converge
	}
	for (int c =0 ; c < 5; c++){
		if (c < data->chrom.size()){
			chrom[c]=data->chrom[c];
		}else{
			chrom[c]='\0';
		}
	}
	chrom[5]='\0';
	if (not FOUND){
		for (int i = 0; i < 12; i++ ){
			ps[i] 	= 0;
		}
	}else{
		ps[0]=C.bidir.mu,ps[1]=C.bidir.si,ps[2]=C.bidir.l, ps[3]=C.bidir.w, ps[4]=C.bidir.pi;
		ps[5]=C.forward.b, ps[6]=C.forward.w, ps[7]=C.forward.pi;
		ps[8]=C.reverse.a, ps[9]=C.reverse.w, ps[10]=C.reverse.pi;
		ps[11]=C.bidir.foot_print;
	}
}

/** Default constructor for simple_c_free_mode.
 */
simple_c_free_mode::simple_c_free_mode(){}

vector<simple_c_free_mode> transform_free_mode(bool FOUND, double ll, component * components, 
	int K, segment * data, int i, double forward_N, double reverse_N) {
	vector<simple_c_free_mode> SC ;
	if (K==0){
		SC.push_back(simple_c_free_mode(false, ll, components[0], K, data, i, forward_N, reverse_N));
	}
	for (int k = 0 ; k < K; k++){
		SC.push_back(simple_c_free_mode(FOUND, ll, components[k], K, data, i,forward_N, reverse_N));
	}
	return SC;

	
}

/** Returns the optimal set of model parameters and results from a given free_mode object.
 * 
 * @param A Set of classifiers to check.
 * @param data Set of segments over which those classifiers (presumably) ran.
 * @param i Unused.
 * 
 * @return Map of ints (used as indices) to simple_c_free_mode objects representing the best fitting models.
 */
map<int, vector<simple_c_free_mode> > get_max_from_free_mode(map<int, vector<classifier> > A, segment * data, int i){
	map<int, vector<simple_c_free_mode> > BEST;
	typedef map<int, vector<classifier> >::iterator it_type_A;
	//get forward and reverse N
	double forward_N =0, reverse_N=0;
	for (int i  = 0 ; i < data->XN;i++){
		forward_N+=data->X[1][i];
		reverse_N+=data->X[2][i];
	}


	for (it_type_A a = A.begin(); a!=A.end(); a++){
		component * best_components;
		double best_ll 	= nINF;
		bool FOUND 		= false;
		int best_k 		= 0;
		for (int r = 0; r < a->second.size(); r++){

			if (A[a->first][r].ll > best_ll){
				best_ll 			= A[a->first][r].ll;
				best_components 	= A[a->first][r].components;
				best_k 				= A[a->first][r].K;
				FOUND 				= true;
			}
		}
		BEST[a->first] 	= transform_free_mode(FOUND, best_ll, best_components, a->first, data, i, forward_N, reverse_N);

	}

	return BEST;
}

/** Runs the model (presumably the core of the behavior seen in the model module) across a set of segments given a depreciated params object.
 * Use run_model_across_free_mode_pwrapper() instead.
 * 
 * @depreciated
 * @param FSI Input segments.
 * @param P Obsolete params object.
 * @param LG Output log file.
 * @return A set of optimal simple_c_free_mode maps as obtained from get_max_from_free_mode().
 */
vector<map<int, vector<simple_c_free_mode> >> run_model_across_free_mode(vector<segment *> FSI, params * P, 
	Log_File * LG){
	vector<map<int, vector<simple_c_free_mode> >> D;
	typedef map<int, vector<classifier> > ::iterator it_type;
	double scale 	= stof(P->p["-ns"]);
	int num_proc 				= omp_get_num_threads();
	int verbose 	= stoi(P->p["-v"]);
	LG->write("(EM) running model across provided intervals...........\n", verbose);
	double N 		= FSI.size();
	double percent 	= 0;
	int elon_move 	= stoi(P->p["-elon"]);
	for (int i = 0 ; i < FSI.size(); i++){
		if ((i / N) > (percent+0.05)){
			LG->write(to_string(int((i / N)*100))+"%,", verbose);
			percent 	= (i / N);
		}

		//first need to populate data->centers
		for (int b = 0 ; b < FSI[i]->bidirectional_bounds.size(); b++){
			double center = FSI[i]->bidirectional_bounds[b][0] +  FSI[i]->bidirectional_bounds[b][1] ;
			center/=2.;
			center-=FSI[i]->start;
			center/=scale;
			FSI[i]->centers.push_back(center);
		}
		segment * data 	= FSI[i];
		map<int, vector<classifier> > A 	= make_classifier_struct_free_model(P, FSI[i]);
		for (it_type k = A.begin(); k!= A.end(); k++){
			int N 	=  k->second.size();
			#pragma omp parallel for num_threads(num_proc)	
			for (int r = 0; r < N; r++ ){
				A[k->first][r].fit2(data, data->centers,0,elon_move);
			}
		}
		D.push_back(get_max_from_free_mode(A, FSI[i], i));
	}
	LG->write("100% done\n", verbose);
	return D;
}

/** Runs the model (presumably the core of the behavior seen in the model module) across a set of segments given a ParamWrapper.
 * 
 * @param FSI Input segments.
 * @param pw Command line arguments used in setting model parameters.
 * @param LG Output log file.
 * @return A set of optimal simple_c_free_mode maps as obtained from get_max_from_free_mode().
 */
vector<map<int, vector<simple_c_free_mode> >> run_model_across_free_mode_pwrapper(vector<segment *> FSI, ParamWrapper *pw, Log_File * LG){
	vector<map<int, vector<simple_c_free_mode> >> D;
	typedef map<int, vector<classifier> > ::iterator it_type;
	double scale 	= pw->ns;
	int num_proc 				= pw->cores;
	int verbose 	= pw->verbose;
	LG->write("(EM) running model across provided intervals...........\n", verbose);
	double N 		= FSI.size();
	double percent 	= 0;
	int elon_move 	= pw->elon;
	for (int i = 0 ; i < FSI.size(); i++){
		if ((i / N) > (percent+0.05)){
			LG->write(to_string(int((i / N)*100))+"%,", verbose);
			percent 	= (i / N);
		}

		//first need to populate data->centers
		for (int b = 0 ; b < FSI[i]->bidirectional_bounds.size(); b++){
			double center = FSI[i]->bidirectional_bounds[b][0] +  FSI[i]->bidirectional_bounds[b][1] ;
			center/=2.;
			center-=FSI[i]->start;
			center/=scale;
			FSI[i]->centers.push_back(center);
		}
		segment * data 	= FSI[i];
		map<int, vector<classifier> > A 	= make_classifier_struct_free_model_pwrapper(pw, FSI[i]);
		for (it_type k = A.begin(); k!= A.end(); k++){
			int N 	=  k->second.size();
			#pragma omp parallel for num_threads(num_proc)	
			for (int r = 0; r < N; r++ ){
				A[k->first][r].fit2(data, data->centers,0,elon_move);
			}
		}
		D.push_back(get_max_from_free_mode(A, FSI[i], i));
	}
	LG->write("100% done\n", verbose);
	return D;
}

/** This computes the average set of model parameters for a given set of segments and an obsolete params object.
 * @depreciated
 * @param segments Input reads.
 * @param P Command line arguments from which to read parameters.
 * @return Set of model parameters.
 */
vector<double> compute_average_model(vector<segment *> segments, params * P){
	//need to compute average model
	double minX 	= 0;
	double maxX 	= 0;
	int test 		= 0;
	double br 		= stod(P->p["-br"]);
	double ns 		= stod(P->p["-ns"]);
	double delta 	= br / ns;
	for (int s = 0 ; s < segments.size(); s++){
		if (segments[s]->N){
			if (segments[s]->maxX > maxX){
				maxX 	= segments[s]->maxX; //this should only evaluate to true once
				test++;
			}
		}
	}
	if (test > 1){
		printf("\nStrange Error in across_segments::compute_average_model\nIgnoring but please consult tFIT contact info\nThank You\n");
	}
	int XN 			= maxX/delta;
	double ** X 	= new double*[3];
	double x 		= 0;
	X[0] 	= new double[XN],X[1] 	= new double[XN],X[2] 	= new double[XN];
	for (int i = 0 ; i < XN;i++){
		X[0][i] 		= x,X[1][i] 		= 0,X[2][i] 		= 0 ;
		x+=delta;
	}
	for (int s = 0 ; s < segments.size(); s++){
		double N 	= 0;
		int j 		= 0;
		if (segments[s]->rN > 1 and segments[s]->fN > 1){
			for (int i = 0 ; i < segments[s]->XN; i++){
				while (j < XN and X[0][j] < segments[s]->X[0][i]){
					j++;
				}
				if (j < XN){
					X[1][j]+=(segments[s]->X[1][i]/segments[s]->fN);
					X[2][j]+=(segments[s]->X[2][i]/segments[s]->rN);
				}
				N+=(segments[s]->X[1][i] + segments[s]->X[2][i] );
			}
		}
	}
	double ll=nINF;
	classifier best_clf;
	for (int r = 0 ; r < 5; r++){
		classifier clf(1, 0.000001, stoi(P->p["-mi"]), 0.3, 
				stod(P->p["-r_mu"]), 10.0, 10.0, 1.0, 
				1.0*segments.size(), 2*segments.size() , stod(P->p["-ALPHA_3"]),0 );
		vector<double> centers 	= {10};
		segment * s 			= new segment("chrX", 0, maxX );
		s->X 					= X;
		s->minX=minX, s->maxX =maxX;
		s->XN 					= XN;
		s->SCALE 				= stod(P->p["-ns"]);
		clf.fit2(s,centers, 0,0);
		if (clf.ll > ll){
			ll 			= clf.ll;
			best_clf 	= clf; 
		}
	}	


	vector<double> parameters(5);


	parameters[0]=	best_clf.components[0].bidir.si*ns;
	parameters[1]=	ns/best_clf.components[0].bidir.l;
	parameters[2]=	best_clf.components[0].bidir.foot_print*ns;
	parameters[3]=	best_clf.components[0].bidir.pi;
	parameters[4]=	best_clf.components[0].bidir.w;


	return parameters;
}	

/** This computes the average set of model parameters for a given set of segments and a ParamWrapper.
 * @param segments Input reads.
 * @param pw Command line arguments from which to read parameters.
 * @return Set of model parameters.
 */
vector<double> compute_average_model_pwrapper(vector<segment *> segments, ParamWrapper *pw){
	//need to compute average model
	double minX 	= 0;
	double maxX 	= 0;
	int test 		= 0;
	double br 		= pw->br;
	double ns 		= pw->ns;
	double delta 	= br / ns;
	for (int s = 0 ; s < segments.size(); s++){
		if (segments[s]->N){
			if (segments[s]->maxX > maxX){
				maxX 	= segments[s]->maxX; //this should only evaluate to true once
				test++;
			}
		}
	}
	if (test > 1){
		printf("\nStrange Error in across_segments::compute_average_model\nIgnoring but please consult tFIT contact info\nThank You\n");
	}
	int XN 			= maxX/delta;
	double ** X 	= new double*[3];
	double x 		= 0;
	X[0] 	= new double[XN],X[1] 	= new double[XN],X[2] 	= new double[XN];
	for (int i = 0 ; i < XN;i++){
		X[0][i] 		= x,X[1][i] 		= 0,X[2][i] 		= 0 ;
		x+=delta;
	}
	for (int s = 0 ; s < segments.size(); s++){
		double N 	= 0;
		int j 		= 0;
		if (segments[s]->rN > 1 and segments[s]->fN > 1){
			for (int i = 0 ; i < segments[s]->XN; i++){
				while (j < XN and X[0][j] < segments[s]->X[0][i]){
					j++;
				}
				if (j < XN){
					X[1][j]+=(segments[s]->X[1][i]/segments[s]->fN);
					X[2][j]+=(segments[s]->X[2][i]/segments[s]->rN);
				}
				N+=(segments[s]->X[1][i] + segments[s]->X[2][i] );
			}
		}
	}
	double ll=nINF;
	classifier best_clf;
	for (int r = 0 ; r < 5; r++){
		classifier clf(1, 0.000001, pw->mi, 0.3, 
				pw->r_mu, 10.0, 10.0, 1.0, 
				1.0*segments.size(), 2*segments.size() , pw->alpha3,0 );
		vector<double> centers 	= {10};
		segment * s 			= new segment("chrX", 0, maxX );
		s->X 					= X;
		s->minX=minX, s->maxX =maxX;
		s->XN 					= XN;
		s->SCALE 				= pw->ns;
		clf.fit2(s,centers, 0,0);
		if (clf.ll > ll){
			ll 			= clf.ll;
			best_clf 	= clf; 
		}
	}	


	vector<double> parameters(5);


	parameters[0]=	best_clf.components[0].bidir.si*ns;
	parameters[1]=	ns/best_clf.components[0].bidir.l;
    printf("Ns value in across_segments: %f\n", ns);
    printf("int ns value: %d\n", (int) ns);
	parameters[2]=	best_clf.components[0].bidir.foot_print*ns;
	parameters[3]=	best_clf.components[0].bidir.pi;
	parameters[4]=	best_clf.components[0].bidir.w;


	return parameters;
}	
