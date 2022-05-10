/**
 * @file load.h
 * @author Joey Azofeifa 
 * @brief Header file for overloaded load.cpp.
 * This is a very overloaded file -- has \ref segment, segment_fits, and \ref node classes
 * Also contains some helper functions.   Should be refactored.
 * 
 * CURRENTLY (May 6, 2022; RDD): Have genomic intervals (Intervals), data (Data), and a
 * class called SetROI (genome level collection) to replace Joey's segment.  Also have
 * ITree to replace the centered interval tree. 
 * 
 * @version 0.1
 */
#ifndef load_H
#define load_H

#include <map>
#include <string>
#include <vector>

#include "read_in_parameters.h"

using namespace std;

class simple_c_free_mode;

class classifier; //forward declare

/**
 * @brief Primary data class which represents a genomic segment of data.
 * Represents both bed files (no data) and bedgraphs (with data).  Also 
 * represents binned, scaled and conditioned data.
 */
class segment{
public:
	string chrom; 	//!< identifier for chromosome 
	int start; //!< genomic coordinates of start
	int stop;  //!< genomic coordinates of stop
	double minX; //!< min coordinate, scaled
	double maxX;   //!< max coordinate, scaled
	string strand; //!< strand information, unspecified = "."

    /** 
	 * @brief vectors of vectors where the inner vector is two doubles
	 * Expected to be [ x y ] with x being a coordinate [0th item] and y [1st item] being a value (reads?)
     */
	vector< vector<double> > forward; //<! corresponds to strand == 1
    /** 
	 * @brief vectors of vectors where the inner vector is two doubles
	 * Expected to be [ x y ] with x being a coordinate [0th item] and y [1st item] being a value (reads?)
     */
	vector< vector<double> > reverse; //<! corresponds to strand == -1

	int ID; //!< Used by MPI code
	int counts; //!< sets min_k in across_segments?

	vector<double> centers; //!< used in the scaling and binning and mu_seeds to fit2 
	vector<vector<double> > bidirectional_bounds; //!< Used by the MPI & slice ratio.

	// Both of below are used in boostrapping code only.  Unclear what they do.
	vector<vector<double>> parameters; 
	map<int, vector<double> > variances;

	/**
	 * @brief This (X) is the smoothed representation of the data.
	 * Vector[0] is coordinate (possibly scaled); [1] is forward (summed for bin)
	 * [2] is reverse (summed for bin).  This is the meat and potatoes data
	 * representation that is used in the EM (fit2).
	 */
	double ** X;  //!< Smoothed data inner is [3] dimensions
	double XN; //!< total number of bins
	double SCALE;  //!< scaling factor

	// I think these are for convenience (calculate once)
	double N;	//!< Total sum of values 

	// I think these are just used for global template matching?
	double fN;	//!< Sum of forward values 
	double rN;	//!< Sum of reverse values 

	vector<vector<double>> fitted_bidirs; //!< mu, sigma, lambda, pi : used in bin()

	// Constructors
	segment(string, int , int);
	segment(string, int , int, int);
	segment(string, int , int, int, string);
	segment();

	// Accessors of commonly used values
	double getXLength ();

	/* FUNCTIONS: */
	string write_interval(); // chr:start-stop,identifier
	string write_allScalar();	// All doubles, strings and ints
	string write_withData();  // includes X
	string write_centers(); 	// Just the centers vector
	string write_bidirectional_bounds(); 	// Just the centers vector

	// bin does the scaling and smoothing of input data (builds X)
	void bin(double, double, bool); // delta, scale, erase
	// add2 appears to add a single data point (coord) to an interval
	void add2(int, double, double); // strand, x, y 

};

/**
 * @brief A node within the interval tree.
 * Interval tree allows for rapid searching based on coordinates.
 */
class node{
public:
	double center;
	int start, stop;
	node * left;  //!< All intervals fully to left of center
	node * right; //!< All intervals fully to the right of center
	vector<segment *> current;	//<! All intervals overlapping center

	void retrieve_nodes(vector<segment * >&);
	void insert_coverage(vector<double>, int);

	// Constructors
	node();	// empty constructor
	node(vector<segment *>);

	/* FUNCTIONS: */
	void searchInterval(int, int, vector<int> &) ;
};

/**
 * @brief The _K_models file is read back into this class. 
 * 
 */
class segment_fits{
public:
	string chrom;
	// At least in model_main, TSS is unused.
	int start, stop, TSS;
	double N;	// N_pos + N_neg
	double N_pos, N_neg;  // counts per strand
	map<int, double> M;	 // <#K, likelihood>
	map<int, string> parameters;  // <#K, parameters_tab_sep_string>
	int model;	// Best model (lowest BIC ratio)
	double BIC_ratio;  // min(null_score / log_likelihood)
	string ID;

	// Constructors
	segment_fits();
	segment_fits(string, int, int, double, double, string);

	/* FUNCTIONS: */
	void identify_best_model(double);
	string write();
};

namespace load{

	vector<segment_fits *> label_tss(string , vector<segment_fits *>   );
	void BIN(vector<segment*>, int, double, bool);

	vector<segment*> load_bedgraphs_total(string, 
		string, string, int , double, string,map<string, int>&,map<int, string>&);

	void write_out_bidirs(map<string , vector<vector<double> > >, string, string, int ,params *, int);
	vector<segment *> load_intervals_of_interest(string,map<int, string>&, params *, bool);

	void collect_all_tmp_files(string , string, int, int );
	vector<segment* > insert_bedgraph_to_segment_joint(map<string, vector<segment *> >  , 
		string , string , string ,int);

	void write_out_models_from_free_mode(map<int, map<int, vector<simple_c_free_mode>  > >,
		params *,int,map<int, string>, int, string &);
	void clear_segments(vector<segment *> );

	vector<segment_fits *> load_K_models_out(string);
	void write_out_bidirectionals_with_penalty(vector<segment_fits*> , params * , int, int );

} // namespace load

#endif
