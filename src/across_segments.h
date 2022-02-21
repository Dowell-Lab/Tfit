/**
 * @file across_segments.h
 * @author Joey Azofeifa 
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#ifndef across_segments_H
#define across_segments_H

#include "error_stdo_logging.h"
#include "load.h"
#include "model.h"
#include "read_in_parameters.h"

string check_file(string, int);

map<int, vector<classifier> > make_classifier_struct_free_model(params * P, segment * data);

void run_model_accross_segments(vector<segment*>, 
	params *);
void free_segments(vector<segment*>);

/***
struct simple_c{
	double ll; //loglikelihood of bidirectional(s) 
	
	double noise_ll; //loglikelihood of noise component
	int IDS[4];

	//IDS[0] 	= segment that it belongs to...
	//IDS[1] 	= complexity, number of fitted models
	//IDS[2] 	= number of predicted bidirs
	//IDS[3] 	= bidir prediction, (possible merged segment)
	double ps[13];
	
//	string write_out();

};
***/

/**
 * @brief keeps model fits
 * 
 * @bug this thing is horrible.  it keeps all the parameters in a double array and
 * you just have to *know* which index correponds to which element!
 * 
 */
struct simple_c_free_mode{
	double SS[3];	//log-likelihood, N_forward, N_reverse
	int ID[5] ;  //index of the segment that this belongs,start, stop, converged?
	char chrom[6];

/* 	ps[0]=C.bidir.mu,ps[1]=C.bidir.si,ps[2]=C.bidir.l, ps[3]=C.bidir.w, ps[4]=C.bidir.pi;
		ps[5]=C.forward.b, ps[6]=C.forward.w, ps[7]=C.forward.pi;
		ps[8]=C.reverse.a, ps[9]=C.reverse.w, ps[10]=C.reverse.pi;
		ps[11]=C.bidir.foot_print;
		*/

	double ps[12]; //parameters for the component
	simple_c_free_mode(bool , double, component ,
		int, segment *, int, double, double);
	simple_c_free_mode();

};

vector<map<int, vector<simple_c_free_mode> >> run_model_across_free_mode(vector<segment *> , params *, Log_File * );
vector<double> compute_average_model(vector<segment *> , params * );
	map<int, vector<simple_c_free_mode> > get_max_from_free_mode(map<int, vector<classifier> > A, segment * data, int i);

#endif
