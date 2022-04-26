/**
 * @file load.cpp
 * @author Joey Azofeifa
 * @brief Routines necessary for loading data and segments.
 * This is a very overloaded file -- has \ref segment, segment_fits, and \ref node classes
 * their corresponding code as well as large functions for loading data. 
 * Also contains some helper functions.   Should be refactored.
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#include "load.h"

#include <math.h>   
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "dirent.h"

#include "across_segments.h"
#include "model.h"
#include "model_selection.h"
#include "read_in_parameters.h"
#include "split.h"
#include "template_matching.h"

using namespace std;

/**
 * @brief Constructors: segment class
 * @author Joey Azofeifa 
 *
 * Purpose: create/allocate instances of segment class
 * The segment class contains both strands of data associated
 * with a particular region.  There is an empty constructor option. 
 *
 * @param chr  Chromosome 
 * @param st   Start
 * @param sp   Stop
 * @param Integer identifier (opt)
 * @param STR  Strand (as string) (opt)
 *
 * @bug STR should be a char with only '.' '+' and '-' as valid.
 *
 */
segment::segment(string chr, int st, int sp){
  chrom	= chr;
  start	= st;
  stop	= sp;
  N 		= 0;
  fN 		= 0;
  rN 		= 0;
  minX=st, maxX=sp;
  counts 	= 1;
  XN 		= 0;
  ID 		= 0;
  strand 	= ".";
  
}
// chrom, start, stop, ID
segment::segment(string chr, int st, int sp, int i){
  chrom	= chr;
  start	= st;
  stop	= sp;
  ID 		= i;
  fN 		= 0;
  rN 		= 0;
  N 		= 0;
  minX=st, maxX=sp;
  counts 	= 1;
  XN 		= 0;
  strand 	= ".";
  
}
// chrom, start, stop, ID, strand
segment::segment(string chr, int st, int sp, int i, string STR){
  chrom	= chr;
  start	= st;
  stop	= sp;
  ID 		= i;
  strand 	= STR;

  N 		= 0;
  fN 		= 0;
  rN 		= 0;
  minX=st, maxX=sp;
  counts 	= 1;
  XN 		= 0;
  
}

// empty -- sets defaults
segment::segment(){
  N 		= 0;
  fN 		= 0;
  rN 		= 0;
  counts 	= 1;
  XN 		= 0;
  ID 		= 0;
  strand 	= ".";
  chrom = "";
}

double segment::getXLength () {
	 return maxX-minX;
}

/**
 * @brief Print region identification for segment 
 * @return printable string of format #chr:start-stop,identifier
 */
string segment::write_interval(){
        string  text= ("#" + chrom + ":" + to_string(start) + "-" 
		+ to_string(stop) + "," + to_string(int(N)) );
	return text;
}

/**
 * @brief Print region identificaton and all scalar values in segment
 * 
 * @return printable string with all scalar values in segment
 */
string segment::write_allScalar() {
	std::string identifier = this->write_interval();
  std::string text = (";n: " + to_string(minX) + ",x: " + to_string(maxX) + 
    ":c" + to_string(counts) +  ";XN: " + to_string(XN) + ";SCALE: " + to_string(SCALE)
    + "**" + to_string(N) + ":" + to_string(fN) + "," + to_string(rN));
  return (identifier + " " + strand + " " + text);
}

/**
 * @brief Print allScalar and then contents of X (multiple lines)
 * @return printable string with all scalar and contents of X
 */
string segment::write_withData() {
	std::string allscalar = this->write_allScalar();
  std::string data;
  int numsperline = 15; // this if for convenience and pretty output
  for(int i = 0; i < 3; i++) {  // keeps (index forward reverse)
    data = data + "\ni:" + to_string(i);
    if (XN == 0) {
      data = data + " no data";
    } else {
      for (int j = 0; j < XN; j++) {
        if (j % numsperline == 0) { data += "\n";}
        // trimming decimals to 2 places for pretty printing
        data = data + " " + to_string(X[i][j]).substr(0,std::to_string(X[i][j]).find(".") + 3);
      }
    }
  }
  return (allscalar + data);
}

/**
 * @brief Print out the contents of the centers vector
 * 
 * @return centers vector as string
 */
string segment::write_centers() {
  std::string centers_as_string;
  for (auto & element : centers) {
    centers_as_string = centers_as_string + " " + to_string(element);
  }
  return centers_as_string;
}

/**
 * @brief Given a data point, build segment contents.
 * @author Joey Azofeifa
 * @param strand  here as integer, 1 == forward, -1 == reverse
 * @param x   coordinate
 * @param y   value (read depth assumed)
 * @return (void)
 */
void segment::add2(int strand, double x, double y){
  vector<double> v2(2);
  v2[0] 	= x;
  v2[1] 	= y;
  if (forward.empty() && reverse.empty()){
    minX=x;
    maxX=x;
  }else{
    if (x < minX){
      minX=x;
      start=int(x);
    }
    if (x > maxX){
      maxX=x;
      stop=int(x);
    }
  }
  // [ x y ] where x is the coordinate and y is the value (depth)
  if (strand == 1){
    forward.push_back(v2);
  }else if (strand==-1){
    reverse.push_back(v2);
  }
}


/**
 * @brief For a segment of data, scale and bin (smooth) input data into X vector
 * @author Joey Azofeifa
 * 
 * Note: X is a 3 x XN 2D array of doubles. 
 * X[0] is scaled coordinate
 * X[1] is sum of forward values over bin width (delta)
 * X[2] is sum of reverse values over bin width (delta)
 * 
 * @param delta binning denominator (nts per bin)
 * @param scale scaling constant (sets SCALE)
 * @param erase whether to remove bins with both strands having zero counts
 * @return (void)
 * 
 * @bug This is klunky, long and cryptic.  This function does a LOT of stuff.
 * And reasoning in several steps is unclear.  Needs to be refactored.
 */
void segment::bin(double delta, double scale, bool erase){

  bool debug = false;

  // X is a 3 x XN 2D array of doubles
  // X[0] is scaled coordinate
  // X[1] is sum of forward values over bin width (delta)
  // X[2] is sum of reverse values over bin width (delta)
  X 				= new double*[3];
  SCALE 			= scale;

  int BINS;
  BINS 		= (getXLength())/delta;
  start = minX, stop=maxX;	// Why are we keeping these distinctly?

  for (int j = 0 ; j < 3;j++){
    X[j] 		= new double[BINS];
  }
  N 				= 0;
  fN = 0, rN = 0;
  XN 				= BINS;  // will adjust if erase

  //===================
  //populate bin ranges
  X[0][0] 		= double(minX);
  X[1][0]=0,X[2][0]=0;
	
  for (int i = 1; i < BINS; i++){
    X[0][i] 	= X[0][i-1] + delta;
    X[1][i] 	= 0;
    X[2][i] 	= 0;
  }

  if (debug) {
	  printf("start: %d , stop: %d , bins: %d ,delta: %f, forward: %d, reverse: %d\n", 
			  start, stop, BINS, delta, forward.size(), reverse.size() );
  }

  // ===================
  //BIN forward strand
  int j 	=0;
  for (int i = 0 ; i < forward.size(); i++){
	  // RDD: To me this seems klunky.  Why do it this way?
    while (j < BINS and X[0][j] <=forward[i][0]){
      j++;
    }
    if (j < BINS and forward[i][0]<= X[0][j]){
      X[1][j-1]+=forward[i][1];
      N+=forward[i][1];
      fN+=forward[i][1];
    }
  }

  //BIN reverse strand
  j 	=0;
  for (int i = 0 ; i < reverse.size(); i++){
	  // RDD: To me this seems klunky.  Why do it this way?
    while (j < BINS and X[0][j] <=reverse[i][0]){
      j++;
    }
    if (j < BINS and reverse[i][0]<= X[0][j]){
      X[2][j-1]+=reverse[i][1];
      N+=reverse[i][1];
      rN+=reverse[i][1];
    }
  }

  //===================
  //scale data down for numerical stability
  if (scale){
    for (int i = 0; i < BINS; i ++ ){
      
      X[0][i] 	= (X[0][i]-minX)/scale;
      // X[1][i]/=delta;
      // X[2][i]/=delta;
    }
  }

  //JA: we also want to get rid of those data points that we don't need
  //i.e. the ones where there is no data coverage values on either the 
  //forward or reverse strands
  // RDD: WHY???
  int realN 		= 0;	// number of non-zero bins
  for (int i = 0; i < BINS;i++){
    if (X[1][i]>0 or X[2][i]>0){
      realN++;
    }
  }

  if (erase){  // going to remove the zero bins
    double ** newX 	= new double*[3];
    for (int j=0; j<3;j++){
      newX[j] 	= new double[realN];
    }
    j = 0;
    for (int i = 0; i < BINS; i ++){
      if (X[1][i]>0 or X[2][i]>0){
	newX[0][j] 	= X[0][i];
	newX[1][j] 	= X[1][i];
	newX[2][j] 	= X[2][i];
	j++;
      }
    }
    if (realN!=j){
      printf("WHAT? %d,%d\n", j, realN);
    }
    //clear previous memory
    for (int i = 0; i < 3; i ++){
      delete X[i];
    }
    delete X;
    X 				= newX;
    XN 				= realN;
  }

  // Scaling centers and fitted bidirectionals
  // Q1: Why not make a single scaling -- whereas here we scaled 
  // 	data above and here are scaling these centers/bidir components?
  //
  if (scale){
    if (not centers.empty()){
      for (int i = 0; i < centers.size(); i++){
	centers[i]=(centers[i]-minX)/scale;			
      }
    }
    if (not fitted_bidirs.empty() ){		// need to understand contents of fitted_bidirs
      for (int fb = 0; fb < fitted_bidirs.size(); fb++){
	int center 	= fitted_bidirs[fb][0];
	int std 	= fitted_bidirs[fb][1]*0.5 + (1. /  fitted_bidirs[fb][2]);
		// 1/2* fb[1] + 1/fb[2]], but what are the 1st and 2nd elements?
	int a 		= center - std*3;
	int b 		= center + std*3;
	
	fitted_bidirs[fb][0] = (fitted_bidirs[fb][0] - minX)/scale;
	fitted_bidirs[fb][1] /= scale;
	fitted_bidirs[fb][2] *= scale; 
      }
    }
    
    maxX 			= (getXLength())/scale;
    minX 			=0;
  }
  double S=0;
  for (int i = 0; i < XN; i++){
    S+=X[1][i];
  }
  // Why do we throw away the raw data?
  forward.clear();
  reverse.clear();
}

//================================================================================================
/**
 * @brief Construct a new node::node object
 * @author Joey Azofeifa 
 * 
 * Keep an interval tree for rapid searching for particular genomics coordinates.
 * 
 */
// Empty constructor
node::node(){};

/**
 * @brief Constructor from a set of segments
 * Given a set of n intervals on the number line, we want to construct 
 * a data structure so that we can efficiently retrieve all intervals overlapping
 * another interval or point.
 * Assumes:  segments is a sorted list of intervals (first has smallest start; last largest stop)
 * 
 * @param segments
 */
node::node(vector<segment * > segments ){
  // Initialize the node
	// I believe this center assumes the segments are sorted.
  center 	= (double(segments[0]->start)  + double(segments[segments.size()-1]->stop)) / 2.;
  left=NULL, right=NULL;

  // Will be building sets of segments to left vs right
  vector<segment * > Left;
  vector<segment * > Right;

  for (int i = 0 ; i < segments.size(); i++){
    if (segments[i]->stop < center){ //all intervals completely to the left of the center point
      Left.push_back(segments[i]);
    } else if (segments[i]->start > center){ // all intervals completely to the right of the center point
      Right.push_back(segments[i]);
    } else{	// all intervals overlapping the center point
      current.push_back(segments[i]);
    }
  }
  // Now recursively build the left and right regions:
  if (Left.size() > 0){
    left 	= new node(Left);
	}
  if (Right.size() > 0){
    right 	= new node(Right);
  }
}

/**
 * @brief Add a data point to all of the nodes that include this point on the tree.
 * @author Joey Azofeifa 
 * Appears to assume that you will always be adding data point to
 * the end of the forward/reverse indexed data points (sorted calls?)
 * 
 * @param x a 2D vector [x y] where X is coordinate and y is value
 * @param s strand (as int: 1 is forward; -1 is reverse)
 * @return (void)
 */
void node::insert_coverage(vector<double> x, int s){
  for (int i = 0 ; i < current.size(); i++){
    if (x[0] > current[i]->start and  x[0] < current[i]->stop  ){
      if (s==1){
	      current[i]->forward.push_back(x);
      }else{
	      current[i]->reverse.push_back(x);	
      }
    }
  }	
 
  // Recursively add point to all relevant intervals.
  if (x[0] >= center and right != NULL ){
    right->insert_coverage(x, s);
  }
  if (x[0] <= center and left !=NULL){
    left->insert_coverage(x,  s);
  }
}

/**
 * @brief Effectively counts number of intervals that contain another interval.
 * @author Joey Azofeifa 
 * @param start
 * @param stop
 * @param finds 
 * @return (void)
 */
void node::searchInterval(int start, int stop, vector<int>& finds ){
  for (int i = 0 ; i < current.size(); i++){
    if (stop > current[i]->start and  start < current[i]->stop  ){
      finds.push_back(1);
    }
  }	
  if (start >= center and right != NULL ){
    right->searchInterval(start, stop, finds);
  }
  if (stop <= center and left !=NULL){
    left->searchInterval(start, stop, finds);
  }	
}

/**
 * @brief Collects all segments associated with a node.
 * @author Joey Azofeifa 
 * @param saves
 * @return (void)
 */
void node::retrieve_nodes(vector<segment*> & saves){
  for (int i = 0; i < current.size(); i++){
    saves.push_back(current[i]);
  }
  if (right!= NULL){
    right->retrieve_nodes(saves);
  }
  if (left != NULL){
    left->retrieve_nodes(saves);		
  }
}

//================================================================================================
/**
 * @brief Constructors for segment_fits class
 * @author Joey Azofeifa 
 */
segment_fits::segment_fits(string c, int st, int sp,
			   double n_pos, double n_neg, string id){

  chrom=c, start=st, stop=sp, TSS=0;
  N 	= n_pos + n_neg, N_pos = n_pos, N_neg = n_neg;
  ID 	= id;
  BIC_ratio 	= 0;
}
segment_fits::segment_fits(){
} //empty con

/**
 * @brief Class constructor for loading the K_models formatted file
 * @author Joey Azofeifa 
 * @param c
 * @param st
 * @param sp
 * @param n_pos
 * @param n_neg
 * @param id
 */
void segment_fits::identify_best_model(double penalty){
  typedef map<int, double>::iterator it_type;
  int arg;
  double MIN=INF;
  double score;
  double null_score;

  // We're going to walk through the M map (#K, log_likelihood) 
  for (it_type m 	= M.begin(); m != M.end(); m++){
    if (m->first > 0){
      // So ms_pen is the complexity penalty for additional models!!
      score 		= -2*m->second + log(N)*(penalty*m->first);
    }else{  // For the first entry, store the null_score
    // recall m->second refers to the value (log_likelihood)
      score 		= -2*m->second + log(N) ;		
      null_score 	= -2*m->second + log(N) ;
    }
    if (score < MIN){
      MIN 		= score;
      arg 		= m->first;
      BIC_ratio 	= null_score/score;
    }
  }
  model 	= arg;
}
/**
 * @brief Writes out a set of models (as bed file) that have reasonable
 * parameters.
 * 
 * @return string 
 */
string segment_fits::write(){
  string line 				= "";
  // model is set by identify_best_model to be the index of min model? 
  if (model > 0){
    string forward_N=to_string(N_pos), reverse_N = to_string(N_neg);
    vector<string> params 		= split_by_tab(parameters[model]); // splits parameters
    for (int i = 0 ; i < model; i++){
      vector<string> S 	=  split_by_comma(params[0], ""); 
      double mu 	= stod(split_by_comma(params[0], "")[i]);
      double std 	= stod(split_by_comma(params[1], "")[i]);
      double lam 	= stod(split_by_comma(params[2], "")[i]);
      double pi 	= stod(split_by_comma(params[3], "")[i]);
      double w = stod(split_by_comma(split_by_bar(params[5], "")[i], "")[0]);

      int start = max(mu - (std + lam), 0.0), stop = mu + (std + lam);

      // This is where he decides what are the best hits.  Ironically he's not
      // actually picking K but rather outputting all regions where the parameters are
      // "reasonable" (defined below).  Seems to work OK, but only because most LL end
      // up as INF and therefore the lam ends up INF
      if (std < 5000 and lam < 20000 and w > 0.05 and pi > 0.05 and pi < 0.95) {
        line += chrom + "\t" + to_string(start) + "\t" + to_string(stop) + "\t";
        line += ID + "|";
        line += to_string(BIC_ratio) + "," + to_string(N_pos) + "," + to_string(N_neg) + "\n";
      }
    }
  }
  return line;
}

//================================================================================================
/**
 * @brief Merge segments from the loading_intervals
 * @param segments
 * @param IDS_first
 * @param IDS
 * @param T
 * @return 
 */
vector<segment *> merge_segments(vector<segment *> segments, map<int, string>  IDS_first, map<int, string> & IDS, int & T){
  vector<segment *> new_segments;
	//bubble sort
  bool changed 	= true;
  while (changed){
    changed = false;
    for (int i = 1 ; i < segments.size(); i++){
      if (segments[i-1]->start > segments[i]->start){
	changed 		= true;
	segment * copy 	= segments[i-1];
	segments[i-1] 	= segments[i];
	segments[i] 	= copy;
      }
    }
  }

  int j = 0, N = segments.size(), i =0;
  while (j < N){
    string ID 		= "";
    segment * S 	= new segment(segments[i]->chrom, 
				      segments[j]->start,segments[j]->stop, T, segments[j]->strand );
    while (j < N and segments[j]->start < S->stop and segments[j]->stop > S->start ){
      S->start 	= min(S->start,segments[j]->start);
      S->stop 	= max(S->stop, segments[j]->stop);
      ID+=IDS_first[segments[j]->ID] + ",";
      S->counts+=1;
      j++;
    }
    new_segments.push_back(S);
    IDS[T] 			= ID.substr(0, ID.size()-1);
    T+=1;
  }
  
  return new_segments;
}

/**
 * @brief Helper function: check identifier doesn't have | character?
 * @author Joey Azofeifa 
 * @param INFO
 * @return 
 */
bool check_ID_name(string & INFO){
  bool PASSED 	= true;
  string change 	= "::";
  for (int i = 0; i < INFO.size(); i++){
    if (INFO.substr(i,1)=="|"){
      PASSED 		= false;
      INFO.replace(i, 1, change);
		}
  }
  return PASSED;
}

//================================================================================================
/**
 * @brief Parses a bedgraph into a collection of segments.
 * Also populates information on chromosomes seen within the bedgraph file and 
 * does the data scaling and smoothing.
 * 
 * Assumptions:
 *    Assumes either joint_bedgraph or (forward_strand reverse_strand) are specified (e.g. not empty)
 *    Has a very hard coded 6 character limit on chromsome name (WHY????)
 * 
 * @author Joey Azofeifa 
 * 
 * @param forward_strand Filename of forward strand data 
 * @param reverse_strand Filename of reverse strand data
 * @param joint_bedgraph Filename of joint data (ij)
 * @param BINS how many bases per smoothing -- for bin() function
 * @param scale what is the scaling constant 
 * @param spec_chrom a specified chromosome name, can be "all"
 * 
 * @param chromosomes returns: maps chromsome name to ?? (counter?)
 * @param ID_to_chrom returns: maps index number to chromosome name
 * @return  a vector of segments
 * 
 * @bug This thing has lots of steps and does lots of things, could 
 * stand to be refactored for clarity. 
 */
vector<segment*> load::load_bedgraphs_total(string forward_strand, string reverse_strand, 
		string joint_bedgraph, int BINS, double scale, string spec_chrom, 
		map<string, int>& chromosomes, map<int, string>& ID_to_chrom){

  bool FOUND 	= false;
  if (spec_chrom=="all"){ FOUND = true; }
  bool EXIT = false;   // Indicator for loop management

  vector<string> FILES;	// Keep file names
  int line_number = 0;
  vector<string> lineArray; // Contents of file, split on tab (\t) 
  string line, chrom;
  int start, stop;
  double coverage;

  segment * S =NULL;
  map<string, segment*> G;  // Data associated with a chrom name
  vector<segment*> segments;	// returned variable
 
  if (forward_strand.empty() and reverse_strand.empty()){
    FILES 	= {joint_bedgraph};  // A single (ij) bedgraph with both strand info
  }else if (not forward_strand.empty() and not reverse_strand.empty()){
    FILES 	= {forward_strand, reverse_strand};  // Distinct files per strand
  }
  
  for (int u = 0 ; u < FILES.size(); u++){
	  bool INSERT     = false;	  
	  string prevChrom="";	// What chrom was on the previous line?
	  ifstream FH(FILES[u]) ;
	  if (not FH){ printf("couln't open FILE %s\n", FILES[u].c_str()); }
	  if (EXIT){ break; }

	  // For every line in this file...
	  while (getline(FH, line)){
		  lineArray=string_split(line, '\t');
		  // Have a hard requirement for a four column bedgraph input
		  if (lineArray.size()!=4){
			  EXIT 	= true;
			  printf("\nLine number %d  in file %s was not formatted properly\nPlease see manual\n",line_number, FILES[u].c_str() );
			  break;
		  }
		  line_number++;
	      // Expects: chrom_name \t start \t stop \t coverage \n
		  chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=(stof(lineArray[3]));

		  if (chrom != prevChrom and (chrom==spec_chrom or spec_chrom=="all")  )  {
			  FOUND 		= true;
			  // Why are we restricting chromosome sizes to 6 characters??
			  if (chrom.size()<6){
				  INSERT          = true;
				  FOUND           = true;
			  }
			  if (chrom.size() < 6 and u==0){
				  G[chrom] 	= new segment(chrom, start, stop);
				  INSERT 		= true;
				  FOUND 		= true;
			  } else if(chrom.size() > 6){
				  INSERT 		= false;
			  }
		  }
		  if (FOUND and chrom!= spec_chrom and spec_chrom!= "all"){
			  break;
		  }
		  if (INSERT){
			  if (u==0){ // When u=0 we are either an ij file or positive strand
				  //If an ij file, then coverage sign indicates strand
				  if (coverage > 0) {  
					  for (int xx = start; xx < stop; xx++){
						  G[chrom]->add2(1, double(xx), abs(coverage));
					  }
				  } else {
					  // so zero coverage is always neg strand?!?
					  // What happens if give positive file with zeros!!!?!
					  for (int xx = start; xx < stop; xx++){
						  G[chrom]->add2(-1, xx, abs(coverage));						
					  }
				  }
			  } else {  // If more than one input file, subsequent is neg strand
				  for (int xx = start; xx < stop; xx++){
					  G[chrom]->add2(-1, xx, abs(coverage));	
				  }
			  }
		  }
		  prevChrom=chrom;

	  }
  }
  if (not EXIT) { // EXIT only true if not right format file
	  int c = 1;
	  typedef map<string, segment*>::iterator it_type;

      // For each chromosome in G (each has a single segment?) 
	  for (it_type i = G.begin(); i != G.end(); i++){
		  i->second->bin(BINS, scale, false);	// Scale and smooth data
		  // Building the naming cross referencing: chromosomes, ID_to_chrom
		  if (chromosomes.find(i->second->chrom)==chromosomes.end()){
			  chromosomes[i->second->chrom]=c;
			  ID_to_chrom[c] 	= i->second->chrom;
			  c++;
		  }
	      // Puts this segment into the return collection
		  segments.push_back(i->second);
	  }
  }
  if (not FOUND){
	  segments.clear();
	  printf("couldn't find chromosome %s in bedgraph files\n", spec_chrom.c_str());
  }
  return segments;
}

/**
 * @brief  load a bedfile of intervals
 * @author Joey Azofeifa 
 * @param FILE name of bedfile containing intervals (example: singleregion.bed)
 * @param IDS side effect:  typically given empty and returns this TOO!
 * @param P parameters for this run
 * @param center (initial call in model_run is 0)
 * @return a vector of segments, note these segments are just regions
 * (bedgraph input) -- and have no data associated with them.
 * Also modifies IDS.
 * 
 */
vector<segment*> load::load_intervals_of_interest(string FILE, map<int, string>&  IDS, 
						  params * P, bool center){
  bool debug = false;    // a debugging indicator
  ifstream FH(FILE);
 
  string spec_chrom 	= P->p["-chr"];
  int pad 	        = stoi(P->p["-pad"])+1;

  // if (debug) { printf("\n  spec_chrom: %s, pad: %d\n", spec_chrom.c_str(), pad);}
  
  vector<segment *> G;
  int ct 	= 1;
  map<string, vector<segment * > > GS;
  map<int, string> IDS_first;
  int T 	= 0;
  bool EXIT 		= false;

  if (FH){
    string line, chrom;
    int start, stop;
    int 	i = 0;
    vector<string> lineArray;
    string strand; 
    bool PASSED 	= true;

    // Reading input file line by line
    while(getline(FH, line)){
      // if (debug) { printf("  %s\n", line.c_str()); }
      lineArray=string_split(line, '\t');  // separate on tab
      // Ignore commented lines (start with #) and require at least 3 columns
      if (lineArray[0].substr(0,1)!="#" and lineArray.size()>2){
        if (lineArray.size() > 3){
          // Expects an identifier in column 4
          if (not check_ID_name(lineArray[3]) and PASSED ){
            PASSED 			= false;
            printf("\ninterval id in line: %s, contains a | symbol changing to :: -> %s\n",line.c_str(), 
                lineArray[3].c_str() );
            printf("Will continue to change other occurrences....\n");

          }
          IDS_first[i] 		= lineArray[3]; // identifier for segment
        }else{
          IDS_first[i] 		= "Entry_" + to_string(i+1);	// Provides identifier if not given
        }

        /* This seems like a bug.  In prelim this is  -0.135002,4514,6894
         * but here he's reading it as strand! */
        if (lineArray.size() > 4){
          if (debug) {  printf("  %s\n", lineArray[4].c_str());  }
          strand 		= lineArray[4];
        }else{
          strand 		= ".";
        }
        try{
          if (not center){
            chrom=lineArray[0], start=max(stoi(lineArray[1])-pad, 0), stop=stoi(lineArray[2]) + pad;
          }else{
            int x 	= 	((stoi(lineArray[1]) + stoi(lineArray[2])))/2.;
            start 		= max(x - pad, 0) , stop 	= x + pad;
            chrom=lineArray[0];
          }
        }
        catch(exception& e){
          printf("\n\nIssue with file %s at line %d\nPlease consult manual on file format\n\n",FILE.c_str(), i );
          EXIT=true;
          GS.clear();
          break;
        }
        if (start < stop){
          if (spec_chrom=="all" or spec_chrom==chrom){	    
            segment * S 	= new segment(chrom, start, stop,i,strand);
            GS[S->chrom].push_back(S);
          }
          i++;
        } // Seems like a bug since ignores start > stop input.
      }
    }
  }else{
    printf("couldn't open %s for reading\n", FILE.c_str() );
    EXIT 	= true;
  }
  if (not EXIT){ 
    typedef map<string, vector<segment * > >::iterator it_type;
    IDS 	= IDS_first;
    for (it_type c 	= GS.begin(); c!=GS.end(); c++){
      vector<segment *> m_segs;
      m_segs 	= c->second;
      for (int i = 0 ; i < m_segs.size(); i++){
        G.push_back(m_segs[i]);
      }
    }
  }else{
    G.clear();
  }
  return G;
}

/**
 * @brief 
 * @author Joey Azofeifa 
 * @param A mapping of chrom name to segment array
 * @param forward Filename of forward strand data 
 * @param reverse Filename of reverse strand data
 * @param joint Filename of joint data (ij)
 * @param rank MPI process number
 * @return a vector of segments
 */
vector<segment* > load::insert_bedgraph_to_segment_joint(map<string, vector<segment *> > A , 
    string forward, string reverse, string joint, int rank ){

  bool debug = true;
  map<string, node> NT;
  typedef map<string, vector<segment *> >::iterator it_type_5;

  // Create an interval tree from the existing intervals
  for(it_type_5 c = A.begin(); c != A.end(); c++) {
    NT[c->first] 	= node(c->second);
  }
  int start, stop, N, j;
  double coverage;
  N 	= 0,j 	= 0;
  int strand;
  int o_st, o_sp;
  vector<string> lineArray;
  string chrom, prevchrom, line;
  vector<segment *> segments;
  double center;
  vector<string> FILES;

  if (forward.empty() and reverse.empty()){
    FILES 	= {joint};
  }else if (not forward.empty() and not reverse.empty()) {
    FILES 	= {forward, reverse};
  }
  string FILE;
  for (int i =0; i < FILES.size(); i++){
    FILE=FILES[i];
    ifstream FH(FILE);
    if (FH){
      prevchrom="";
      while (getline(FH, line)){
        lineArray       = string_split(line, '\t');
        if (lineArray.size()==4){
          chrom 		= lineArray[0];
          start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = stod(lineArray[3]);
          if (coverage > 0 and i == 0){
            strand 	= 1;
          }else if (coverage < 0 or i==1){
            strand 	= -1;
          }
          center 	= (stop + start) /2.;
          if (NT.find(chrom)!=NT.end()){
            for (int center_2=start; center_2 < stop; center_2++){
              vector<double> x(2);
              x[0]=double(center_2), x[1] = abs(coverage);
              NT[chrom].insert_coverage(x, strand);
            }

          }
        } else { 
          printf("\n***error in line: %s, not bedgraph formatted\n", line.c_str() );
          segments.clear();
          return segments;
        }
      }
      FH.close();
    }else{
      cout<<"could not open forward bedgraph file: "<<FILE<<endl;
      segments.clear();
      return segments;
    }
  }
  //now we want to get all the intervals and make a vector<segment *> again...
  vector<segment *>NS;
  typedef map<string, node>::iterator it_type_6;
  for (it_type_6 c = NT.begin(); c!=NT.end(); c++){
    c->second.retrieve_nodes(NS);
  }

  return NS;
}


/**
 * @brief   Reads in a _K_models_MLE file
 * @author Joey Azofeifa 
 * @param FILE  (contents produced by write_out_models_from_free_mode)
 * @return  pointer to collection of model fits. 
 */
vector<segment_fits *> load::load_K_models_out(string FILE){
  ifstream FH(FILE);
  string line;
  segment_fits * S = NULL;
  int complexity;
  double ll;
  string chrom;
  int start,stop;
  vector<segment_fits *> segment_fits_all;
  if (FH) {
    while (getline(FH, line)) {
      if (line.substr(0, 1) != "#" and line.substr(0, 1) == ">") {
        // This is the region header line.
        // > ID|chr:start-end|Nforward,Nreverse
        // >ME_36597:0|chr3:57965508-57976184|2572.000000,1885.000000
        if (S != NULL) { segment_fits_all.push_back(S); }
        line = line.substr(1, line.size() - 1); //removes 1st character, ">"

        vector<string> bar_split = split_by_bar(line, "");
        vector<string> comma_split = split_by_comma(bar_split[2], ""); // Nforward Nreverse
        vector<string> colon_split = split_by_colon(bar_split[1], ""); // chr:start-end
        chrom = colon_split[0];
        vector<string> dash_split = split_by_dash(colon_split[1], ""); // start end
        start = stoi(dash_split[0]), stop = stoi(dash_split[1]);
        S = new segment_fits(chrom, start,
                             stop, stod(comma_split[0]), stod(comma_split[1]), bar_split[0]);
      } else if (line.substr(0, 1) == "~" and S != NULL) {
        // This is a model output line.
        // output lines are tab separated fields each comma separated for # models
        line = line.substr(1, line.size() - 1); // removes the ~ first character
        // Example:
        // ~2,-18405.714702	57970016.797935,57974849.280602	69.438995,37.975287	371.734252,408.579080	0.486291,0.586718	121.206766,250.000000	0.617766,0.069621,0.000284|0.310201,0.000341,0.000348	57976184.000000,57976184.000000	57965508.000000,57965508.000000
        // First field is model_complexity,log_likelihood
        vector<string> tab_split = split_by_tab(line);
        vector<string> comma_split = split_by_comma(tab_split[0], "");
        complexity = stoi(comma_split[0]), ll = stod(comma_split[1]);
        S->M[complexity] = ll;
        // subsequent fields were output as: 
        // mus+"\t"+sigmas+"\t"+lambdas+"\t" + pis+"\t" + fps+ "\t" + ws + "\t" + fbs+"\t" +ras;
        string parameters = "";
        // Isn't this just rebuilding the parameters string (e.g. tab_split[1->size]?)
        for (int i = 1; i < tab_split.size(); i++) {
          if (i + 1 < tab_split.size()) {
            parameters += tab_split[i] + "\t";
          } else {
            parameters += tab_split[i];
          } 
        } // for each element tab_split
        S->parameters[complexity] = parameters;
      } // if header or model output line
    } // while each line of the file
    if (S!=NULL){
      segment_fits_all.push_back(S);		
    }
  } else {
    // A file handling error (couldn't open the file)
    printf("couldn't open %s...strange error\n", FILE.c_str() );
  }
  return segment_fits_all;
}

//================================================================================================
vector<vector<double>> bubble_sort_alg(vector<vector<double>> X){
	bool changed 	= true;
	while (changed){
		changed = false;
		for (int i = 1 ; i < X.size(); i++){
			if (X[i-1][0] > X[i][0]){
				changed 				= true;
				vector<double > copy 	= X[i-1];
				X[i-1] 					= X[i];
				X[i] 	= copy;
			}
		}
	}
	return X;
}

/**
 * @brief 
 * @author Joey Azofeifa 
 * @param G
 * @param out_dir
 * @param job_name
 * @param job_ID
 * @param P
 * @param noise
 * @return (void)
 */
void load::write_out_bidirs(map<string , vector<vector<double> > > G, string out_dir, 
			    string job_name,int job_ID, params * P, int noise){
  typedef map<string , vector<vector<double> > >::iterator it_type;
  ofstream FHW;
  FHW.open(out_dir+ job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed");
  FHW<<P->get_header(1);
  int ID 	= 0;
  for (it_type c = G.begin(); c!=G.end(); c++){
    vector<vector<double>> data_intervals 	=  bubble_sort_alg(c->second);
    
    for (int i = 0; i < data_intervals.size(); i++){
      FHW<<c->first<<"\t"<<to_string(int(data_intervals[i][0]))<<"\t"<<to_string(int(data_intervals[i][1]))<<"\tME_"<<to_string(ID)<<"\t";
      FHW<<to_string(data_intervals[i][2] )+"," + to_string(int(data_intervals[i][3] )) + "," + to_string(int(data_intervals[i][4]) )<<endl; 
      ID++;
    }
  }
  FHW.close();
}

/**
 * @brief Writes the _K_models_MLE.tsv file.
 * 
 * @param G 
 * @param P parameters
 * @param job_ID output job name (used as prefix to file name)
 * @param IDS 
 * @param noise 
 * @param file_name  This does NOT appear to be used!
 */
void load::write_out_models_from_free_mode(map<int, map<int, vector<simple_c_free_mode>  > > G, 
	params * P, int job_ID,map<int, string> IDS, int noise, string & file_name){

	//========================================================================================
	//write out each model parameter estimates
	double scale 	= stof(P->p["-ns"]);
	double penalty = stof(P->p["-ms_pen"]) ;
	string out_dir 	= P->p["-o"];
	ofstream FHW;
	file_name 	= out_dir+  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv";
	FHW.open(file_name);
	FHW<<P->get_header(2);
	
	typedef map<int, map<int, vector<simple_c_free_mode>  > >::iterator it_type_1;
	typedef map<int, vector<simple_c_free_mode>  > ::iterator it_type_2;
	typedef vector<simple_c_free_mode>::iterator it_type_3;
	
	typedef map<int, string>::iterator it_type_IDS;
	int IN=0;
	string mus="", sis="", ls="", wEMs="", wPIs="",forward_bs="", forward_ws="",forward_PIs="",reverse_as="", reverse_ws="",reverse_PIs="";
	double w_thresh 	= 0.;
	double ALPHA_2 	= stof(P->p["-ALPHA_2"]);
	string mu, sigma, lambda, pos, neg, ll, pi, w,ra,fb,rw,fw,fp;
	int start;
	string chrom;

	string INFO 	= "";
			
	FHW<<"#ID|chromosome:start-stop|forward strand coverage, reverse strand coverage"<<endl;
	FHW<<"#model complexity,log-likelihood"<<endl;
	FHW<<"#mu_k\tsigma_k\tlambda_k\tpi_k\tfp_k\tw_[p,k],w_[f,k],w_[r,k]\tb_[f,k]\ta_[r,k]"<<endl;
	
	for (it_type_1 s = G.begin(); s!=G.end(); s++){ //iterate over each segment
		FHW<<">" + IDS[s->first]+ "|";
		for (it_type_2 k 	= s->second.begin(); k != s->second.end(); k++){//iterate over each model_complexity
			for (it_type_3 c = k->second.begin(); c!=k->second.end(); c++){
				chrom 		= (*c).chrom;
				INFO 		= chrom + ":" + to_string((*c).ID[1])+"-"+to_string((*c).ID[2]);
				pos 		= to_string((*c).SS[1]);
				neg 		= to_string((*c).SS[2]);
				
			}
		}
		FHW<<INFO<<"|"<<pos+","+neg<<endl;


		for (it_type_2 k 	= s->second.begin(); k != s->second.end(); k++){//iterate over each model_complexity
			string mus 		= "", sigmas="", lambdas="",pis="", ws= ""  ,fbs="",ras="", fws="", rws="", fps="";
			string pos 		= "", neg="";
			string k_header = "~"+ to_string(k->first)+",";
			int NN 			= k->second.size();
			int ii 			= 0;
			for (it_type_3 c = k->second.begin(); c!=k->second.end(); c++){
				chrom 		= (*c).chrom;
				start 		= (*c).ID[1];
				mu 			= to_string((*c).ps[0]*scale + (*c).ID[1] );
				sigma 		= to_string((*c).ps[1]*scale);
				lambda 		= to_string(scale/(*c).ps[2]);
				pi 			= to_string( (*c).ps[4]);
				w  			= to_string( (*c).ps[3]);
				fw 			= to_string( (*c).ps[6]); 
				rw 			= to_string( (*c).ps[9]);
				ra 			= to_string( scale*(*c).ps[8]   + (*c).ID[1] );
				fb 			= to_string( scale*(*c).ps[5]  + (*c).ID[1]);
				fp 			= to_string( scale*(*c).ps[11]  );
				ll 			= to_string((*c).SS[0]);
				if (ii +1 < NN){
					mus+=mu+",";
					sigmas+=sigma+",";
					lambdas+=lambda+",";
					pis+=pi+",";
					ws+=w+ "," + fw+ "," + rw + "|";
					ras+=ra+",";
					fbs+=fb+",";
					fps+=fp+",";
				}else{
					mus+=mu ;
					sigmas+=sigma ;
					lambdas+=lambda ;
					pis+=pi ;
					ws+=w+ "," + fw+ "," + rw ;	
					fbs+=fb;
					ras+=ra;
					fps+=fp ;
				}
				ii++;
			}
			k_header 		+=ll+ "\t";
			FHW<<k_header;
		
			if (k->first>0){
				FHW<<mus+"\t"+sigmas+"\t"+lambdas+"\t" + pis+"\t" + fps+ "\t" + ws + "\t" + fbs+"\t" +ras ;
			}
			FHW<<endl;
		}
	}
	FHW.flush();
}



/**
 * @brief  Writes out the best model as a bed file.
 * @author Joey Azofeifa 
 * @param fits
 * @param P   parameters
 * @param job_ID
 * @param noise
 * @return (void)
 */
void load::write_out_bidirectionals_with_penalty(vector<segment_fits*> fits, params * P, int job_ID, int noise ){
	ofstream FHW;
	FHW.open(P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_bidir_predictions.bed");	
	FHW<<P->get_header(2);
	double penalty 	= stod(P->p["-ms_pen"]);
	for (int i = 0; i < fits.size(); i++){
		fits[i]->identify_best_model(penalty); // applies a complexity penalty and finds min likelihood model
		FHW<<fits[i]->write();
	}
	FHW.flush();
}

//================================================================================================
//
/**
 * @brief 
 * @author Joey Azofeifa 
 * @param segments
 * @param BINS
 * @param scale
 * @param erase
 * @return (void)
 */
void load::BIN(vector<segment*> segments, int BINS, double scale, bool erase){
	for (int i = 0 ; i < segments.size() ; i ++){
		if (segments[i]->forward.size() > 0 or segments[i]->reverse.size() > 0 ){
			segments[i]->bin(BINS, scale, erase);
		}
	}
}

/**
 * @brief 
 * @author Joey Azofeifa 
 * @param dir
 * @param job_name
 * @param nprocs
 * @param job_ID
 * @return (void)
 */
void load::collect_all_tmp_files(string dir, string job_name, int nprocs, int job_ID){
	int c 	= 0;
	time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m_%d_%H_%M", &tstruct);
    string DT 	= buf;
	string OUT 		= dir+ job_name + "-" + to_string(job_ID) +"_" + DT+ ".log";
	ofstream FHW(OUT);
	for (int rank = 0; rank < nprocs; rank++){
		string FILE 	= dir+"tmp_" + job_name+ "-" +to_string(job_ID) + "_" + to_string(rank) + ".log";
		string line;
		ifstream FH(FILE);
		if (FH){
			if (rank!=0){
				FHW<<"=======MPI Call: " + to_string(rank) +"=======\n";
			}
			while (getline(FH, line)){
				if ("#" != line.substr(0,1) or rank==0){
					FHW<<line<<endl;
				}
			}

			FH.close();
			remove( FILE.c_str()) ;
			c++;
		}
	}
}

/**
 * @brief 
 * @author Joey Azofeifa 
 * @param segments
 * @return (void)
 */
void load::clear_segments(vector<segment *> segments){
	for (int i = 0; i < segments.size(); i++){
		if (segments[i]!=NULL){
			delete (segments[i]);
		}
	}
}

/**
 * @brief 
 * @author Joey Azofeifa 
 * @param tss_file
 * @param query_fits
 * @return 
 */
vector<segment_fits *> load::label_tss(string tss_file, vector<segment_fits *> query_fits ){
	vector<segment_fits *> new_fits;
	ifstream FH(tss_file);
	string line;
	map<string, vector<segment *> >G;
	map<string, node> T;
	string chrom,start, stop;
	typedef map<string, vector<segment *>>::iterator it_type;
	if (FH){	
		while (getline(FH, line)){
			vector<string> lineArray 	= split_by_tab(line);
			chrom 	= lineArray[0], start = lineArray[1], stop = lineArray[2];
			segment * S 	= new segment(chrom, stoi(start), stoi(stop));
			G[chrom].push_back(S);
		}
		//make G a node interval tree
		for (it_type c = G.begin(); c!=G.end(); c++){
			T[c->first] 	= node(c->second);
		}
		//label
		for (int i = 0 ; i < query_fits.size(); i++){
			vector<int> FINDS;
			if (T.find(query_fits[i]->chrom) != T.end()){
				T[query_fits[i]->chrom].searchInterval(query_fits[i]->start,query_fits[i]->stop, FINDS);
				if (!FINDS.empty() and query_fits[i]->M[1]!=nINF and query_fits[i]->M[1] > query_fits[i]->M[0]   ){
					query_fits[i]->TSS=1;
					new_fits.push_back(query_fits[i]);
				}
			}
			FINDS.clear();

		}
		
	}else{
		printf("couldn't open tss file %s... weird error\n", tss_file.c_str() );
	}
	return new_fits;	
}


