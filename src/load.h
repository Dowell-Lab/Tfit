/**
 * @file load.h
 * @author Joey Azofeifa, Robin Dowell
 * @brief Header file for overloaded load.cpp.
 * This is a very overloaded file -- has \ref segment, segment_fits, and \ref node classes
 * Also contains some helper functions.   Should be refactored.
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
 * @author Joey Azofeifa, Robin Dowell
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
	double fN;	//!< Sum of forward values 
	double rN;	//!< Sum of reverse values 

	vector<vector<double> > bidirectional_bounds; //!< Used by the MPI
	vector<vector<double>> fitted_bidirs; //!< mu, si, l, pi : used in bin()

	// Constructors
	segment(string, int , int);
	segment(string, int , int, int);
	segment(string, int , int, int, string);
	segment();

	/* FUNCTIONS: */
	// Reporting out (currently unused)
	string write_out();
	// bin does the scaling and smoothing of input data (builds X)
	void bin(double, double, bool); // delta, scale, erase
	// add2 appears to add a single data point (coord) to an interval
	void add2(int, double, double); // strand, x, y 
};

/**
 * @brief A node within the interval tree.
 * Interval tree allows for rapid searching based on coordinates.
 * @author Joey Azofeifa
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
 * @brief Needs documentation.
 * 
 */
class segment_fits{
public:
	string chrom;
	int start, stop, TSS;
	double N;
	double N_pos, N_neg;
	map<int, double> M;
	map<int, string> parameters;
	int model;
	double BIC_ratio;
	string ID;

	// Constructors
	segment_fits();
	segment_fits(string, int, int, double, double, string);

	/* FUNCTIONS: */
	void get_model(double);
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
	void write_out_bidirectionals_ms_pen(vector<segment_fits*> , params * , int, int );

} // namespace load

#endif
