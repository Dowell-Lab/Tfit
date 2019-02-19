#include <string>
#include <vector>
#include "load.h"
#include "split.h"
#include <iostream>
#include <fstream>
#include <map>
#include "dirent.h"
#include "template_matching.h"
#include "model.h"
#include "read_in_parameters.h"
#include "across_segments.h"
#include "model_selection.h"
#include <cmath>
#include <math.h>
#include <limits>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <random>
#include <exception>

#include <stdio.h>
#include <time.h>
#include <math.h>

using namespace std;
//========================================================================
//The very very important segment class

/** Minimal constructor consisting only of position and chromosome information.
 * Sets all values not corresponding to the given parameters to either 0 or:
 * "." (indicating both positive and negative coverage) for strand
 * 1 for counts.
 * 
 * NOTE: minX and maxX are set to st and sp in addition to start and stop!
 * @param chr Chromosome in the form chrN
 * @param st Starting position in base pairs
 * @param sp Stop position in base pairs.
 */
segment::segment(string chr, int st, int sp) {
   chrom = chr;
   start = st;
   stop  = sp;
   N     = 0;
   fN    = 0;
   rN    = 0;
   minX = st, maxX = sp;
   counts  = 1;
   XN    = 0;
   ID    = 0;
   strand  = ".";
   chrom_ID = 0;

}

/** Slightly less minimal constructor that includes an ID parameter (i).
 * 
 * @param chr Chromosome in the form chrN
 * @param st Starting position in base pairs
 * @param sp Stop position in base pairs.
 * @param i ID for this segment.
 */
segment::segment(string chr, int st, int sp, int i) {
   chrom = chr;
   start = st;
   stop  = sp;
   fN    = 0;
   rN    = 0;
   N     = 0;
   minX = st, maxX = sp;
   counts  = 1;
   XN    = 0;
   ID    = i;
   strand  = ".";
   chrom_ID = 0;

}

/** Full segment constructor that includes an ID parameter in addition to a strand parameter.
 * 
 * @param chr Chromosome in the form chrN
 * @param st Starting position in base pairs
 * @param sp Stop position in base pairs
 * @param i ID for this segment
 * @param STR strand parameter for this segment (in the form +, -, or .)
 */
segment::segment(string chr, int st, int sp, int i, string STR) {
   chrom = chr;
   start = st;
   stop  = sp;
   N     = 0;
   fN    = 0;
   rN    = 0;
   minX = st, maxX = sp;
   counts  = 1;
   XN    = 0;
   ID    = i;
   strand  = STR;
   chrom_ID = 0;

}

/** Default constructor that sets all values to 0 except for the following:
 * counts is set to 1
 * strand is set to "."
 */
segment::segment() {
   N     = 0;
   fN    = 0;
   rN    = 0;
   counts  = 1;
   XN    = 0;
   ID    = 0;
   strand  = ".";
   chrom_ID = 0;
}

/** Returns the segment as a string.
 * @return the segment as a string in the form "#chrom:start-stop,N\n"
 */
string segment::write_out() {
   string text   = ("#" + chrom + ":" + to_string(start) + "-"
                    + to_string(stop) + "," + to_string(int(N)) + "\n");
   return text;
}

/** Pushes the specified vector into the forward or reverse strand vector.
 * This function should also account for the minX and maxX parameters in addition to starting and stopping parameters.
 * @param strand +1 for positive strand, -1 for negative strand. All other values are invalid.
 * @param x The first value of the vector to insert. This is additionally used to alter the min/maxX values.
 * @param y The second value of the vector to insert.
 */
void segment::add2(int strand, double x, double y) {
   vector<double> v2(2);
   v2[0]   = x;
   v2[1]   = y;
   if (forward.empty() && reverse.empty()) {
      minX = x;
      maxX = x;
   } else {
      if (x < minX) {
         minX = x;
         start = int(x);
      }
      if (x > maxX) {
         maxX = x;
         stop = int(x);
      }
   }
   if (strand == 1) {
      forward.push_back(v2);
   } else if (strand == -1) {
      reverse.push_back(v2);
   }
}

/** Sorts the vector passed by the second element of every sub-vector?
 * TODO: look into the arguments of the sort function itself to determine how this is supposed to behave.
 * Original description: "sort vector of vectors by second"
 * @param The vector to be sorted.
 * @return Sorted vector (actually the first parameter)
 */
vector<vector<double>> bubble_sort_by_1(vector<vector<double>> vec) { //sort vector of vectors by second
    sort(vec.begin(), vec.end(),
             [](const  vector<double>& a, const  vector<double>& b) {
     return a[0] < b[0];
   });
   return vec;
}

/** Places data into bins based on various parameters.
 * 
 * @parameter delta
 * @parameter scale
 * @parameter erase
 */
void segment::bin(double delta, double scale, int erase) {
   X         = new double*[3];
   SCALE       = scale;
   int BINS;
   BINS    = (maxX - minX) / delta;
   start = minX, stop = maxX;
   for (int j = 0 ; j < 3; j++) {
      X[j]    = new double[BINS];
   }
   N         = 0;
   fN = 0, rN = 0;
   XN        = BINS;
   //===================
   //populate bin ranges
   X[0][0]     = double(minX);
   X[1][0] = 0, X[2][0] = 0;
   forward     = bubble_sort_by_1(forward);
   reverse     = bubble_sort_by_1(reverse);



   for (int i = 1; i < BINS; i++) {
      X[0][i]   = X[0][i - 1] + delta;
      X[1][i]   = 0;
      X[2][i]   = 0;
   }

   // ===================
   //insert forward strand
   int j   = 0;
   //printf("start: %d , stop: %d , bins: %d ,delta: %f, forward: %d, reverse: %d\n", start, stop, BINS, delta, forward.size(), reverse.size() );
   for (int i = 0 ; i < forward.size(); i++) {
      while (j < BINS and X[0][j] <= forward[i][0]) {
         j++;
      }
      if (j < BINS and forward[i][0] <= X[0][j]) {
         X[1][j - 1] += forward[i][1];
         N += forward[i][1];
         fN += forward[i][1];
      }
   }
   j   = 0;
   //===================
   //insert reverse strand
   for (int i = 0 ; i < reverse.size(); i++) {
      while (j < BINS and X[0][j] <= reverse[i][0]) {
         j++;
      }
      if (j < BINS and reverse[i][0] <= X[0][j]) {
         X[2][j - 1] += reverse[i][1];
         N += reverse[i][1];
         rN += reverse[i][1];
      }
   }
   //===================
   //scale data down for numerical stability
   if (scale) {
      for (int i = 0; i < BINS; i ++ ) {

         X[0][i]   = (X[0][i] - minX) / scale;
         // X[1][i]/=delta;
         // X[2][i]/=delta;
      }
   }
   //we also want to get rid of those data points that we don't need
   //i.e. the ones where there is no data coverage values on either the
   //forward or reverse strands

   int realN     = 0;
   for (int i = 0; i < BINS; i++) {
      if (X[1][i] > 0 or X[2][i] > 0) {
         realN++;
      }
   }
   if (erase) {
      double ** newX  = new double*[3];
      for (int j = 0; j < 3; j++) {
         newX[j]   = new double[realN];
      }
      j = 0;
      for (int i = 0; i < BINS; i ++) {
         if (X[1][i] > 0 or X[2][i] > 0) {
            newX[0][j]  = X[0][i];
            newX[1][j]  = X[1][i];
            newX[2][j]  = X[2][i];
            j++;
         }
      }
      if (realN != j) {
         printf("WHAT? %d,%d\n", j, realN);
      }
      //clear previous memory
      for (int i = 0; i < 3; i ++) {
         delete X[i];
      }
      delete X;
      X         = newX;
      XN        = realN;
   }
   if (scale) {
      if (not centers.empty()) {
         for (int i = 0; i < centers.size(); i++) {
            centers[i] = (centers[i] - minX) / scale;
         }
      }
      if (not fitted_bidirs.empty() ) {
         for (int fb = 0; fb < fitted_bidirs.size(); fb++) {
            int center  = fitted_bidirs[fb][0];
            int std   = fitted_bidirs[fb][1] * 0.5 + (1. /  fitted_bidirs[fb][2]);
            int a     = center - std * 3;
            int b     = center + std * 3;



            fitted_bidirs[fb][0] = (fitted_bidirs[fb][0] - minX) / scale;
            fitted_bidirs[fb][1] /= scale;
            fitted_bidirs[fb][2] *= scale;
         }
      }

      maxX      = (maxX - minX) / scale;
      minX      = 0;
   }
   double S = 0;
   for (int i = 0; i < XN; i++) {
      S += X[1][i];
   }
   forward.clear();
   reverse.clear();
}

//================================================================================================
//interval tree code

/** Default constructor for a node; it does nothing.
 */
node::node() {};

/** Constructor for a node that accepts a vector of segments.
 * Internally, this constructor automatically recursively generates a binary tree
 * based on the segments passed.
 * 
 * @param segments vector of segment pointers.
 */
node::node(vector<segment * > segments ) {
   center  = (double(segments[0]->start)  + double(segments[segments.size() - 1]->stop)) / 2.;
   vector<segment * > Left;
   vector<segment * > Right;
   left = NULL, right = NULL;
   for (int i = 0 ; i < segments.size(); i++) {
      if (segments[i]->stop < center) {
         Left.push_back(segments[i]);
      }
      else if (segments[i]->start > center) {
         Right.push_back(segments[i]);
      }
      else {
         current.push_back(segments[i]);
      }
   }
   if (Left.size() > 0) {
      left  = new node(Left);
   }
   if (Right.size() > 0) {
      right   = new node(Right);
   }
}

/** Inserts coverage data into the node.
 * This process involves pushing all data from a given set of (reads?) into this object's 
 * forward and reverse strand vectors, then automatically passing the data along to its child nodes.
 * This should effectively ensure that all coverage can still be broken down in a binary fashion.
 * 
 * @param x vector of reads
 * @param s Strand in which to insert the vector of reads. +1 for positive strand; -1 for negative strand.
 */
void node::insert_coverage(vector<double> x, int s) {
   for (int i = 0 ; i < current.size(); i++) {
      if (x[0] > current[i]->start and  x[0] < current[i]->stop  ) {
         if (s == 1) {
            current[i]->forward.push_back(x);
         } else {
            current[i]->reverse.push_back(x);
         }
      }
   }

   if (x[0] >= center and right != NULL ) {
      right->insert_coverage(x, s);
   }
   if (x[0] <= center and left != NULL) {
      left->insert_coverage(x,  s);
   }
}

/** Performs a binary search over the tree in order to find any nodes that fit within the specified start and stop values. Based on the output values involved (namely adding unit to the vector of finds), it can be used to find the distance to a given interval.
 * 
 * @param start Start position (base pairs).
 * @param stop Stop position (base pairs).
 * @param finds Vector into which to insert found intervals.
 */
void node::searchInterval(int start, int stop, vector<int>& finds ) {
   for (int i = 0 ; i < current.size(); i++) {
      if (stop > current[i]->start and  start < current[i]->stop  ) {
         finds.push_back(1);
      }
   }
   if (start >= center and right != NULL ) {
      right->searchInterval(start, stop, finds);
   }
   if (stop <= center and left != NULL) {
      left->searchInterval(start, stop, finds);
   }
}

/** This returns the set of all segments in order.
 * @param saves Vector into which to insert segment pointers.
 */
void node::retrieve_nodes(vector<segment*> & saves) {
   for (int i = 0; i < current.size(); i++) {
      saves.push_back(current[i]);
   }
   if (right != NULL) {
      right->retrieve_nodes(saves);
   }
   if (left != NULL) {
      left->retrieve_nodes(saves);
   }
}

//================================================================================================
//a class from loading the K_models formmated file
/** Default constructor for segment_fits. Does nothing.
 */
segment_fits::segment_fits() {
} //empty con
/** Full segment_fits constuctor that sets all parameters except for TSS and BIC_ratio.
 * 
 * @param c Chromosome in the form chrN
 * @param st Starting position (base pairs)
 * @param sp Stop position (base pairs)
 * @param n_pos Number of positive elements in the given segment?
 * @param n_neg Number of negative element sin the given segment?
 * @param id Some form of string identifier for the object.
 */
segment_fits::segment_fits(string c, int st, int sp,
                           double n_pos, double n_neg, string id) {

   chrom = c, start = st, stop = sp, TSS = 0;
   N   = n_pos + n_neg, N_pos = n_pos, N_neg = n_neg;
   ID  = id;
   BIC_ratio   = 0;
}

/** This function iterates over M (a map from integers to doubles) and appears to either score or get log likelihoods for the segment_fits.
 * 
 * @param ms_pen 
 */ 
void segment_fits::get_model(double ms_pen) {
   typedef map<int, double>::iterator it_type;
   int arg;
   double MIN = INF;
   double score;
   double null_score;

   for (it_type m  = M.begin(); m != M.end(); m++) {
      if (m->first > 0) {
         score     = -2 * m->second + log(N) * (ms_pen * m->first);
      } else {
         score     = -2 * m->second + log(N) ;
         null_score  = -2 * m->second + log(N) ;
      }
      if (score < MIN) {
         MIN     = score;
         arg     = m->first;
         BIC_ratio   = null_score / score;
      }
   }
 
 model   = arg;
}

/** Represents the data within the object as a string.
 * 
 * @return The segment_fits as a string.
 */
string segment_fits::write () {
   string line         = "";
   int tag = 0;
   if (model > 0) {
      string forward_N = to_string(N_pos), reverse_N = to_string(N_neg);
      vector<string> params     = split_by_tab(parameters[model], "");
      for (int i = 0 ; i < model; i++) {
         vector<string> S  =  split_by_comma(params[0], "");
         double mu   = stod(split_by_comma(params[0], "")[i]);
         double std  = stod(split_by_comma(params[1], "")[i]);
         double lam  = stod(split_by_comma(params[2], "")[i]);
         double pi   = stod(split_by_comma(params[3], "")[i]);
         double w  = stod(split_by_comma(split_by_bar(params[5], "")[i] , "" )[0] );

// Seems that this might be where we could be losing some bidirectionals...?
// Changed std to 10k from 5k, lam to 30k from 20k, w to 0.01 from 0.05 -- believe this increases the number of calls, but for the best...?
          
         int start   = max(mu - (std + lam), 0.0), stop = mu + (std + lam);
         if (std  < 5000 and lam < 20000 and w > 0.05 and pi > 0.05 and pi < 0.95 ) {
            line += chrom + "\t" + to_string(start) + "\t" + to_string(stop);
            line += "\t" + ID + "." + to_string(++tag) + "|" + to_string(BIC_ratio);
            line += "\t" + to_string(lround(mu)) + "\t" + to_string(lround(std)) + "\t" + to_string(w) + "\t" + to_string(lround(lam)) + "\t" + to_string(lround(N_pos)) + "\t" + to_string(lround(N_neg)) + "\n";
         } 
      }
   }
   return line;
}

//================================================================================================
/**merge segments from loading_intervals
 * @param segments vector of segments to merge
 * @param IDS_first map from ints to strings of IDs where the ID value is presumably first.
 * @param IDS map from ints to strings in which to place revised ID mappings.
 * @param T int representing the index of the last element in IDS.
 * @return vector of merged segment pointers.
 */
vector<segment *> merge_segments(vector<segment *> segments, map<int, string>  IDS_first, map<int, string> & IDS, int & T) {
   vector<segment *> new_segments;
   //bubble sort
   int changed   = 1;
   while (changed) {
      changed = 0;
      for (int i = 1 ; i < segments.size(); i++) {
         if (segments[i - 1]->start > segments[i]->start) {
            changed     = 1;
            segment * copy  = segments[i - 1];
            segments[i - 1]   = segments[i];
            segments[i]   = copy;
         }
      }
   }

   int j = 0, N = segments.size(), i = 0;
   while (j < N) {
      string ID     = "";
      segment * S   = new segment(segments[i]->chrom,
                                  segments[j]->start, segments[j]->stop, T, segments[j]->strand );
      while (j < N and segments[j]->start < S->stop and segments[j]->stop > S->start ) {
         S->start  = min(S->start, segments[j]->start);
         S->stop   = max(S->stop, segments[j]->stop);
         ID += IDS_first[segments[j]->ID] + ",";
         S->counts += 1;
         j++;
      }
      new_segments.push_back(S);
      IDS[T]      = ID.substr(0, ID.size() - 1);
      T += 1;
   }

   return new_segments;
}
/** Checks to determine whether or not a given input string matches various rules. This is likely used
 * when reading segments from disk.
 * @param INFO Input string that gets transformed into a different representation?
 * @return whether or not the check passed.
 */
int check_ID_name(string & INFO) {
   int PASSED  = 1;
   string change   = "::";
   for (int i = 0; i < INFO.size(); i++) {
      if (INFO.substr(i, 1) == "|") {
         PASSED    = 0;
         INFO.replace(i, 1, change);
      }
   }
   return PASSED;
}

//================================================================================================
//LOADING from file functions...need to clean this up...

/** Reads all elements from a bedgraph and returns a set of segments.
 * 
 * @param forward_strand Forward strand bed3 file. This can be left blank if "joint_bedgraph" is specified.
 * @param reverse_strand Reverse strand bed3 file. This can be left blank if "joint_bedgraph" is specified.
 * @param joint_bedgraph Bed6 file containing both forward and reverse strand reads. This can take the place of both "forward_strand" and "reverse_strand."
 * @param BINS Number of bins to use in segmenting the input bedgaph file. This should affect how various data is mapped and searched for internally.
 * @param scale Scaling parameter used in generating individual segments. 
 * @param spec_chrom User-specified chromosome. This enables Tfit to efficiently run on a subset of a whole dataset.
 * @param chromosomes Map in which to place chromosome to index mappings.
 * @param ID_to_chrom Map in which to place index to chromosome mappings.
 * @return Set of segments read from the specified input file(s).
 */
vector<segment*> load::load_bedgraphs_total(string forward_strand,
      string reverse_strand, string joint_bedgraph, int BINS, double scale, string spec_chrom, map<string, int>& chromosomes
      , map<int, string>& ID_to_chrom) {
   int FOUND   = 0;
   if (spec_chrom == "all") {
      FOUND   = 1;
   }
   map<string, segment*>   G;
   vector<segment*> segments;
   vector<string> FILES;
   if (forward_strand.empty() and reverse_strand.empty()) {
      FILES   = {joint_bedgraph};
   } else if (not forward_strand.empty() and not reverse_strand.empty()) {
      FILES   = {forward_strand, reverse_strand};
   }

   string line, chrom;
   int start, stop;
   double coverage;
   vector<string> lineArray;
   segment * S = NULL;
   int EXIT    = 0;
   int line_number = 0;
   for (int u = 0 ; u < FILES.size(); u++) {
      int INSERT     = 0;
      string prevChrom = "";
      ifstream FH(FILES[u]) ;
      if (not FH ) {
         printf("couln't open FILE %s\n", FILES[u].c_str());
      }
      if (EXIT) {
         break;
      }
      while (getline(FH, line)) {
         lineArray=splitter(line, "\t");
         //lineArray = string_split(line, '\t');
         if (lineArray.size() != 4) {
            EXIT  = 1;
            printf("\nLine number %d  in file %s was not formatted properly\nPlease see manual\n", line_number, FILES[u].c_str() );
            break;
         }
         line_number++;
         chrom = lineArray[0], start = stoi(lineArray[1]), stop = stoi(lineArray[2]), coverage = (stof(lineArray[3]));
         if (chrom != prevChrom and (chrom == spec_chrom or spec_chrom == "all")  )  {
            FOUND     = 1;
            if (chrom.size() < 6) {
               INSERT          = 1;
               FOUND           = 1;
            }
            if (chrom.size() < 6 and u == 0) {
               G[chrom]  = new segment(chrom, start, stop );
               INSERT    = 1;
               FOUND     = 1;
            } else if (chrom.size() > 6) {
               INSERT    = 0;
            }
         }
         if (FOUND and chrom != spec_chrom and spec_chrom != "all") {
            break;
         }
         if (INSERT) {
            if (u == 0) {
               if (coverage > 0) {
                  for (int xx = start; xx < stop; xx++) {
                     G[chrom]->add2(1, double(xx), abs(coverage));
                  }
               } else {
                  for (int xx = start; xx < stop; xx++) {
                     G[chrom]->add2(-1, double(xx), abs(coverage));
                  }
               }
            } else {
               for (int xx = start; xx < stop; xx++) {
                  G[chrom]->add2(-1, double(xx), abs(coverage));
               }
            }
         }
         prevChrom = chrom;

      }
   }
   if (not EXIT) {
      int c = 1;
      typedef map<string, segment*>::iterator it_type;
      for (it_type i = G.begin(); i != G.end(); i++) {
         i->second->bin(BINS, scale, 0);
         if (chromosomes.find(i->second->chrom) == chromosomes.end()) {
            chromosomes[i->second->chrom] = c;
            ID_to_chrom[c]  = i->second->chrom;
            c++;
         }
         segments.push_back(i->second);
      }
   }
   if (not FOUND) {
      segments.clear();
      printf("couldn't find chromosome %s in bedgraph files\n", spec_chrom.c_str());
   }
   return segments;
}

/** Loads a set of segments from a given input training file. This file should indicate groundtruth
 * regions of transcription from which Tfit can generate a useful model. This function uses a depreciated params object. Please use load_intervals_of_interest_pwrapper instead.
 * @depreciated
 * @param FILE Input training file path.
 * @param IDS Map of integers to strings in which to place ID to (chromosome?) mappings.
 * @param P Obsolete params object from which to read user parameters such as desired chromosome, etc.
 * @param center Whether or not the input reads file indicates the center of the given set of distributions instead of entire regions of transcription.
 * @return Set of segments representing training intervals.
 */
vector<segment*> load::load_intervals_of_interest(string FILE, map<int, string>&  IDS,
      params * P, int center) {
   ifstream FH(FILE);

   string spec_chrom   = P->p["-chr"];
   int pad           = stoi(P->p["-pad"]) + 1;

   vector<segment *> G;
   int ct  = 1;
   map<string, vector<segment * > > GS;
   map<int, string> IDS_first;
   int T   = 0;
   int EXIT    = 0;
   if (FH) {
      string line, chrom;
      int start, stop;
      int   i = 0;
      vector<string> lineArray;
      string strand;
      int PASSED  = 1;

      while (getline(FH, line)) {
         //lineArray=splitter(line, "\t");
         lineArray = string_split(line, '\t');
         if (lineArray[0].substr(0, 1) != "#" and lineArray.size() > 2) {
            if (lineArray.size() > 3) {
               if (not check_ID_name(lineArray[3]) and PASSED ) {
                  PASSED      = 0;
                  printf("\ninterval id in line: %s, contains a | symbol changing to :: -> %s\n", line.c_str(), lineArray[3].c_str() );
                  printf("Will continue to change other occurrences....\n");

               }
               IDS_first[i]    = lineArray[3];
            } else {
               IDS_first[i]    = "Entry_" + to_string(i + 1);
            }
            if (lineArray.size() > 4) {
               strand    = lineArray[4];
            } else {
               strand    = ".";
            }
            try {
               if (not center) {
                  chrom = lineArray[0], start = max(stoi(lineArray[1]) - pad, 0), stop = stoi(lineArray[2]) + pad;
               } else {
                  int x   =   ((stoi(lineArray[1]) + stoi(lineArray[2]))) / 2.;
                  start     = max(x - pad, 0) , stop  = x + pad;
                  chrom = lineArray[0];
               }
            }
            catch (exception& e) {
               printf("\n\nIssue with file %s at line %d\nPlease consult manual on file format\n\n", FILE.c_str(), i );
               EXIT = 1;
               GS.clear();
               break;
            }
            if (start < stop) {
               if (spec_chrom == "all" or spec_chrom == chrom) {
                  segment * S   = new segment(chrom, start, stop, i, strand);
                  GS[S->chrom].push_back(S);
               }
               i++;
            }
         }
      }
   } else {
      printf("couldn't open %s for reading\n", FILE.c_str() );
      EXIT  = 1;
   }
   if (not EXIT) {
      typedef map<string, vector<segment * > >::iterator it_type;
      IDS   = IDS_first;
      for (it_type c  = GS.begin(); c != GS.end(); c++) {
         vector<segment *> m_segs;
         m_segs  = c->second;
         for (int i = 0 ; i < m_segs.size(); i++) {
            G.push_back(m_segs[i]);
         }
      }
   } else {
      G.clear();
   }
   return G;
}

/** Loads a set of segments from a given input training file. This file should indicate groundtruth
 * regions of transcription from which Tfit can generate a useful model. 
 * @depreciated
 * @param FILE Input training file path.
 * @param IDS Map of integers to strings in which to place ID to (chromosome?) mappings.
 * @param pw ParamWrapper object from which to read user parameters such as desired chromosome, etc.
 * @param center Whether or not the input reads file indicates the center of the given set of distributions instead of entire regions of transcription.
 * @return Set of segments representing training intervals.
 */
vector<segment*> load::load_intervals_of_interest_pwrapper(string FILE, map<int, string>&  IDS,
      ParamWrapper *pw, int center) {
   ifstream FH(FILE);

   string spec_chrom   = pw->chromosome;
   int pad           = pw->pad + 1;

   vector<segment *> G;
   int ct  = 1;
   map<string, vector<segment * > > GS;
   map<int, string> IDS_first;
   int T   = 0;
   int EXIT    = 0;
   if (FH) {
      string line, chrom;
      int start, stop;
      int   i = 0;
      vector<string> lineArray;
      string strand;
      int PASSED  = 1;

      while (getline(FH, line)) {
         //lineArray=splitter(line, "\t");
         lineArray = string_split(line, '\t');
         if (lineArray[0].substr(0, 1) != "#" and lineArray.size() > 2) {
            if (lineArray.size() > 3) {
               if (not check_ID_name(lineArray[3]) and PASSED ) {
                  PASSED      = 0;
                  printf("\ninterval id in line: %s, contains a | symbol changing to :: -> %s\n", line.c_str(), lineArray[3].c_str() );
                  printf("Will continue to change other occurrences....\n");

               }
               IDS_first[i]    = lineArray[3];
            } else {
               IDS_first[i]    = "Entry_" + to_string(i + 1);
            }
            if (lineArray.size() > 4) {
               strand    = lineArray[4];
            } else {
               strand    = ".";
            }
            try {
               if (not center) {
                  chrom = lineArray[0], start = max(stoi(lineArray[1]) - pad, 0), stop = stoi(lineArray[2]) + pad;
               } else {
                  int x   =   ((stoi(lineArray[1]) + stoi(lineArray[2]))) / 2.;
                  start     = max(x - pad, 0) , stop  = x + pad;
                  chrom = lineArray[0];
               }
            }
            catch (exception& e) {
               printf("\n\nIssue with file %s at line %d\nPlease consult manual on file format\n\n", FILE.c_str(), i );
               EXIT = 1;
               GS.clear();
               break;
            }
            if (start < stop) {
               if (spec_chrom == "all" or spec_chrom == chrom) {
                  segment * S   = new segment(chrom, start, stop, i, strand);
                  GS[S->chrom].push_back(S);
               }
               i++;
            }
         }
      }
   } else {
      printf("couldn't open %s for reading\n", FILE.c_str() );
      EXIT  = 1;
   }
   if (not EXIT) {
      typedef map<string, vector<segment * > >::iterator it_type;
      IDS   = IDS_first;
      for (it_type c  = GS.begin(); c != GS.end(); c++) {
         vector<segment *> m_segs;
         m_segs  = c->second;
         for (int i = 0 ; i < m_segs.size(); i++) {
            G.push_back(m_segs[i]);
         }
      }
   } else {
      G.clear();
   }
   return G;
}

/** This inserts a set of reads into a set of preexisting segments per-chromosome.
 * @param A Chromosome-indexed map of segment vectors in which to insert the given data.
 * @param forward Forward strand bed3 file (may be left blank if "joint" is specified.)
 * @param reverse Reverse strand bed3 file (may be left blank if "joint" is specified.)
 * @param joint Combined bed6 file with both forward and reverse strand reads. (May be used in place of "forward" and "reverse")
 * @param rank Rank parameter obtained from the MPI runtime.
 * @return Vector of processed segments.
 */
vector<segment* > load::insert_bedgraph_to_segment_joint(map<string, vector<segment *> > A ,
      string forward, string reverse, string joint, int rank ) {
   map<string, node> NT;
   typedef map<string, vector<segment *> >::iterator it_type_5;
   for (it_type_5 c = A.begin(); c != A.end(); c++) {
      NT[c->first]  = node(c->second);
   }
   int start, stop, N, j;
   double coverage;
   N   = 0, j   = 0;
   int strand;
   int o_st, o_sp;
   vector<string> lineArray;
   string chrom, prevchrom, line;
   vector<segment *> segments;
   double center;
   vector<string> FILES;
   if (forward.empty() and reverse.empty()) {
      FILES   = {joint};
   } else if (not forward.empty() and not reverse.empty()) {
      FILES   = {forward, reverse};
   }
   string FILE;
   for (int i = 0; i < FILES.size(); i++) {
      FILE = FILES[i];
      ifstream FH(FILE);
      if (FH) {
         prevchrom = "";
         while (getline(FH, line)) {
            //lineArray   = splitter2(line, "\t");
            lineArray       = string_split(line, '\t');
            if (lineArray.size() == 4) {
               chrom     = lineArray[0];
               start = stoi(lineArray[1]), stop = stoi(lineArray[2]), coverage = stod(lineArray[3]);
               if (coverage > 0 and i == 0) {
                  strand  = 1;
               } else if (coverage < 0 or i == 1) {
                  strand  = -1;
               }
               center  = (stop + start) / 2.;
               if (NT.find(chrom) != NT.end()) {
                  for (int center_2 = start; center_2 < stop; center_2++) {
                     vector<double> x(2);
                     x[0] = double(center_2), x[1] = abs(coverage);
                     NT[chrom].insert_coverage(x, strand);
                  }

               }
            }
            else {
               printf("\n***error in line: %s, not bedgraph formatted\n", line.c_str() );
               segments.clear();
               return segments;
            }
         }
         FH.close();

      } else {
         cout << "could not open forward bedgraph file: " << FILE << endl;
         segments.clear();
         return segments;
      }
   }
   //now we want to get all the intervals and make a vector<segment *> again...
   vector<segment *>NS;
   typedef map<string, node>::iterator it_type_6;
   for (it_type_6 c = NT.begin(); c != NT.end(); c++) {
      c->second.retrieve_nodes(NS);
   }

   return NS;
}

/** Returns a set of segment_fits from a given input model bidir file.
 * This should be used in conjunction with the bidir module to provide useful model inputs.
 * @param FILE Input file from which to read fits.
 * @return A set of segment_fits representing individual segment models.
 */
vector<segment_fits *> load::load_K_models_out(string FILE) {
   ifstream FH(FILE);
   string line;
   segment_fits * S = NULL;
   int complexity;
   double ll;
   string chrom;
   int start, stop;
   vector<segment_fits *> segment_fits_all;
   if (FH) {
      while (getline(FH, line)) {
         if (line.substr(0, 1) != "#" and line.substr(0, 1) == ">") {
            if (S != NULL) {
               segment_fits_all.push_back(S);
            }
            line              = line.substr(1, line.size() - 1);

            vector<string> bar_split    = split_by_bar(line, "");
            vector<string> comma_split    = split_by_comma(bar_split[2], "");
            vector<string> colon_split    = split_by_colon(bar_split[1], "");
            vector<string> dash_split     = split_by_dash(colon_split[1], "");
            chrom = colon_split[0], start   = stoi(dash_split[0]) , stop = stoi(dash_split[1]) ;
            S     = new segment_fits(chrom, start,
                                     stop, stod(comma_split[0]), stod(comma_split[1]), bar_split[0] );

         } else if (line.substr(0, 1) == "~" and S != NULL) {
            line  = line.substr(1, line.size() - 1);
            vector<string> tab_split    = split_by_tab(line, "");
            vector<string> comma_split    = split_by_comma(tab_split[0], "");
            complexity  = stoi(comma_split[0]), ll  = stod(comma_split[1]);
            S->M[complexity] = ll;
            string parameters   = "";
            for (int i = 1; i < tab_split.size(); i++) {
               if (i + 1 < tab_split.size()) {
                  parameters += tab_split[i] + "\t";
               } else {
                  parameters += tab_split[i];
               }
            }
            S->parameters[complexity]   = parameters;

         }
      }
      if (S != NULL) {
         segment_fits_all.push_back(S);
      }
   } else {
      printf("couldn't open %s...strange error\n", FILE.c_str() );
   }
   return segment_fits_all;
}


//================================================================================================
//WRITE out to file functions

/** Sorts a set of sets of values.
 * 
 * @param X Values to sort
 * @return Sorted set
 */
vector<vector<double>> bubble_sort_alg(vector<vector<double>> X) {
   int changed   = 1;
   while (changed) {
      changed = 0;
      for (int i = 1 ; i < X.size(); i++) {
         if (X[i - 1][0] > X[i][0]) {
            changed         = 1;
            vector<double > copy  = X[i - 1];
            X[i - 1]          = X[i];
            X[i]  = copy;
         }
      }
   }
   return X;
}

/** Dumps a set of bidirectional predictions and fits to the specified output file. This function uses the obsolete params class and is thus depreciated. Please use write_out_bidirs_pwrapper instead.
 * 
 * @depreciated
 * @param G Chromosome-indexed map of (what appear to be?) model parameters.
 * @param out_dir Directory in which to write all of the bidirectional output files.
 * @param job_name Job name prefix with which to name the output files.
 * @param job_ID Job ID number likely obtained from the MPI runtime.
 * @param P Command line arguments encapsulated in an obsolete params object.
 */
void load::write_out_bidirs(map<string , vector<vector<double> > > G, string out_dir,
                            string job_name, int job_ID, params * P, int noise) {
   typedef map<string , vector<vector<double> > >::iterator it_type;
   ofstream FHW;
   FHW.open(out_dir);
   FHW << P->get_header(1);
   int ID  = 0;
   for (it_type c = G.begin(); c != G.end(); c++) {
      vector<vector<double>> data_intervals   =  bubble_sort_alg(c->second);

      for (int i = 0; i < data_intervals.size(); i++) {
         FHW << c->first << "\t" << to_string(int(data_intervals[i][0])) << "\t" << to_string(int(data_intervals[i][1])) << "\tBIDIR_" << to_string(ID) << "\t";
         FHW << to_string(data_intervals[i][2] ) + "," + to_string(int(data_intervals[i][3] )) + "," + to_string(int(data_intervals[i][4]) ) << endl;
         ID++;
      }
   }
   FHW.close();
}

/** Dumps a set of bidirectional predictions and fits to the specified output file.
 * 
 * @param G Chromosome-indexed map of (what appear to be?) model parameters.
 * @param out_dir Directory in which to write all of the bidirectional output files.
 * @param job_name Job name prefix with which to name the output files.
 * @param job_ID Job ID number likely obtained from the MPI runtime.
 * @param pw Command line arguments encapsulated in a ParamWrapper.
 */
void load::write_out_bidirs_pwrapper(map<string , vector<vector<double> > > G, string out_dir,
                            string job_name, int job_ID, ParamWrapper *pw, int noise) {
   typedef map<string , vector<vector<double> > >::iterator it_type;
   ofstream FHW;
   FHW.open(out_dir);
   FHW << pw->getHeader(1);
   int ID  = 0;
   for (it_type c = G.begin(); c != G.end(); c++) {
      vector<vector<double>> data_intervals   =  bubble_sort_alg(c->second);

      for (int i = 0; i < data_intervals.size(); i++) {
         FHW << c->first << "\t" << to_string(int(data_intervals[i][0])) << "\t" << to_string(int(data_intervals[i][1])) << "\tPRELIM_" << to_string(ID) << endl;
         ID++;
      }
   }
   FHW.close();
}

/** Dumps a set of fit model parameters to a file. This function uses a depreciated params object. Please use write_out_models_from_free_mode_pwrapper instead.
 * 
 * @depreciated
 * @param G Map between (segment ids or indices?) and simple_c_free_mode models.
 * @param P Obsolete params object containing user command line options.
 * @param job_ID Job identifier possibly obtained from the MPI runtime.
 * @param IDS Map of ids to chromosome names.
 * @param noise Unused.
 * @param file_name Output file path. This name is modified by the function so that it is possible to refer back to the output file.
 */
void load::write_out_models_from_free_mode(map<int, map<int, vector<simple_c_free_mode>  > > G,
      params * P, int job_ID, map<int, string> IDS, int noise, string & file_name) {

   //========================================================================================
   //write out each model parameter estimates
   double scale  = stof(P->p["-ns"]);
   double penality = stof(P->p["-ms_pen"]) ;
   string dir  = P->p["-l"];
   ofstream FHW;
   file_name   = dir +  P->p["-N"] + "-" + to_string(job_ID) +  "_K_models_MLE.tsv";
   FHW.open(file_name);
   FHW << P->get_header(2);

   typedef map<int, map<int, vector<simple_c_free_mode>  > >::iterator it_type_1;
   typedef map<int, vector<simple_c_free_mode>  > ::iterator it_type_2;
   typedef vector<simple_c_free_mode>::iterator it_type_3;

   typedef map<int, string>::iterator it_type_IDS;
   int IN = 0;
   string mus = "", sis = "", ls = "", wEMs = "", wPIs = "", forward_bs = "", forward_ws = "", forward_PIs = "", reverse_as = "", reverse_ws = "", reverse_PIs = "";
   double w_thresh   = 0.;
   double ALPHA_2  = stof(P->p["-ALPHA_2"]);
   string mu, sigma, lambda, pos, neg, ll, pi, w, ra, fb, rw, fw, fp;
   int start;
   string chrom;

   string INFO   = "";

   FHW << "#ID|chromosome:start-stop|forward strand coverage, reverse strand coverage" << endl;
   FHW << "#model complexity,log-likelihood" << endl;
   FHW << "#mu_k\tsigma_k\tlambda_k\tpi_k\tfp_k\tw_[p,k],w_[f,k],w_[r,k]\tb_[f,k]\ta_[r,k]" << endl;

   for (it_type_1 s = G.begin(); s != G.end(); s++) { //iterate over each segment
      FHW << ">" + IDS[s->first] + "|";
      for (it_type_2 k  = s->second.begin(); k != s->second.end(); k++) { //iterate over each model_complexity
         for (it_type_3 c = k->second.begin(); c != k->second.end(); c++) {
            chrom     = (*c).chrom;
            INFO    = chrom + ":" + to_string((*c).ID[1]) + "-" + to_string((*c).ID[2]);
            pos     = to_string((*c).SS[1]);
            neg     = to_string((*c).SS[2]);

         }
      }
      FHW << INFO << "|" << pos + "," + neg << endl;


      for (it_type_2 k  = s->second.begin(); k != s->second.end(); k++) { //iterate over each model_complexity
         string mus    = "", sigmas = "", lambdas = "", pis = "", ws = ""  , fbs = "", ras = "", fws = "", rws = "", fps = "";
         string pos    = "", neg = "";
         string k_header = "~" + to_string(k->first) + ",";
         int NN      = k->second.size();
         int ii      = 0;
         for (it_type_3 c = k->second.begin(); c != k->second.end(); c++) {
            chrom     = (*c).chrom;
            start     = (*c).ID[1];
            mu      = to_string((*c).ps[0] * scale + (*c).ID[1] );
            sigma     = to_string((*c).ps[1] * scale);
            lambda    = to_string(scale / (*c).ps[2]);
            pi      = to_string( (*c).ps[4]);
            w       = to_string( (*c).ps[3]);
            fw      = to_string( (*c).ps[6]);
            rw      = to_string( (*c).ps[9]);
            ra      = to_string( scale * (*c).ps[8]   + (*c).ID[1] );
            fb      = to_string( scale * (*c).ps[5]  + (*c).ID[1]);
            fp      = to_string( scale * (*c).ps[11]  );
            ll      = to_string((*c).SS[0]);
            if (ii + 1 < NN) {
               mus += mu + ",";
               sigmas += sigma + ",";
               lambdas += lambda + ",";
               pis += pi + ",";
               ws += w + "," + fw + "," + rw + "|";
               ras += ra + ",";
               fbs += fb + ",";
               fps += fp + ",";
            } else {
               mus += mu ;
               sigmas += sigma ;
               lambdas += lambda ;
               pis += pi ;
               ws += w + "," + fw + "," + rw ;
               fbs += fb;
               ras += ra;
               fps += fp ;
            }
            ii++;
         }
         k_header    += ll + "\t";
         FHW << k_header;

         if (k->first > 0) {
            FHW << mus + "\t" + sigmas + "\t" + lambdas + "\t" + pis + "\t" + fps + "\t" + ws + "\t" + fbs + "\t" + ras ;
         }
         FHW << endl;
      }
   }
   FHW.flush();
}

/** Dumps a set of fit model parameters to a file.
 * 
 * @depreciated
 * @param G Map between (segment ids or indices?) and simple_c_free_mode models.
 * @param pw ParamWrapper object containing user command line options.
 * @param job_ID Job identifier possibly obtained from the MPI runtime.
 * @param IDS Map of ids to chromosome names.
 * @param noise Unused.
 * @param file_name Output file path. This name is modified by the function so that it is possible to refer back to the output file.
 */
void load::write_out_models_from_free_mode_pwrapper(map<int, map<int, vector<simple_c_free_mode>  > > G,
      ParamWrapper *pw, int job_ID, map<int, string> IDS, int noise, string & file_name) {

   //========================================================================================
   //write out each model parameter estimates
   double scale  = pw->ns;
   double penality = pw->penalty;
   string dir  = pw->logDir;
   ofstream FHW;
   file_name   = dir + pw->jobName + "-" + to_string(job_ID) +  "_MLE.tsv";
   FHW.open(file_name);
   FHW << pw->getHeader(2);

   typedef map<int, map<int, vector<simple_c_free_mode>  > >::iterator it_type_1;
   typedef map<int, vector<simple_c_free_mode>  > ::iterator it_type_2;
   typedef vector<simple_c_free_mode>::iterator it_type_3;

   typedef map<int, string>::iterator it_type_IDS;
   int IN = 0;
   string mus = "", sis = "", ls = "", wEMs = "", wPIs = "", forward_bs = "", forward_ws = "", forward_PIs = "", reverse_as = "", reverse_ws = "", reverse_PIs = "";
   double w_thresh   = 0.;
   double ALPHA_2  = pw->alpha2;
   string mu, sigma, lambda, pos, neg, ll, pi, w, ra, fb, rw, fw, fp;
   int start;
   string chrom;

   string INFO   = "";

   FHW << "#ID|chromosome:start-stop|forward strand coverage, reverse strand coverage" << endl;
   FHW << "#model complexity,log-likelihood" << endl;
   FHW << "#mu_k\tsigma_k\tlambda_k\tpi_k\tfp_k\tw_[p,k],w_[f,k],w_[r,k]\tb_[f,k]\ta_[r,k]" << endl;

   for (it_type_1 s = G.begin(); s != G.end(); s++) { //iterate over each segment
      FHW << ">" + IDS[s->first] + "|";
      for (it_type_2 k  = s->second.begin(); k != s->second.end(); k++) { //iterate over each model_complexity
         for (it_type_3 c = k->second.begin(); c != k->second.end(); c++) {
            chrom     = (*c).chrom;
            INFO    = chrom + ":" + to_string((*c).ID[1]) + "-" + to_string((*c).ID[2]);
            pos     = to_string((*c).SS[1]);
            neg     = to_string((*c).SS[2]);

         }
      }
      FHW << INFO << "|" << pos + "," + neg << endl;


      for (it_type_2 k  = s->second.begin(); k != s->second.end(); k++) { //iterate over each model_complexity
         string mus    = "", sigmas = "", lambdas = "", pis = "", ws = ""  , fbs = "", ras = "", fws = "", rws = "", fps = "";
         string pos    = "", neg = "";
         string k_header = "~" + to_string(k->first) + ",";
         int NN      = k->second.size();
         int ii      = 0;
         for (it_type_3 c = k->second.begin(); c != k->second.end(); c++) {
            chrom     = (*c).chrom;
            start     = (*c).ID[1];
            mu      = to_string((*c).ps[0] * scale + (*c).ID[1] );
            sigma     = to_string((*c).ps[1] * scale);
            lambda    = to_string(scale / (*c).ps[2]);
            pi      = to_string( (*c).ps[4]);
            w       = to_string( (*c).ps[3]);
            fw      = to_string( (*c).ps[6]);
            rw      = to_string( (*c).ps[9]);
            ra      = to_string( scale * (*c).ps[8]   + (*c).ID[1] );
            fb      = to_string( scale * (*c).ps[5]  + (*c).ID[1]);
            fp      = to_string( scale * (*c).ps[11]  );
            ll      = to_string((*c).SS[0]);
            if (ii + 1 < NN) {
               mus += mu + ",";
               sigmas += sigma + ",";
               lambdas += lambda + ",";
               pis += pi + ",";
               ws += w + "," + fw + "," + rw + "|";
               ras += ra + ",";
               fbs += fb + ",";
               fps += fp + ",";
            } else {
               mus += mu ;
               sigmas += sigma ;
               lambdas += lambda ;
               pis += pi ;
               ws += w + "," + fw + "," + rw ;
               fbs += fb;
               ras += ra;
               fps += fp ;
            }
            ii++;
         }
         k_header    += ll + "\t";
         FHW << k_header;

         if (k->first > 0) {
            FHW << mus + "\t" + sigmas + "\t" + lambdas + "\t" + pis + "\t" + fps + "\t" + ws + "\t" + fbs + "\t" + ras ;
         }
         FHW << endl;
      }
   }
   FHW.flush();
}

/** Writes a set of bidirectional model fits utilizing a penalty value passed from the command line (the -ms_pen parameter at the time of this writing). This function uses a depreciated params object. Please use write_out_bidirectionals_ms_pen_pwrapper instead.
 * 
 * @depreciated
 * @param fits Set of models to write.
 * @param P Obsolete params object from which to read command line parameters.
 * @param job_ID Job identifier most likely obtained from the MPI runtime.
 * @param noise Unused.
 */
void load::write_out_bidirectionals_ms_pen(vector<segment_fits*> fits, params * P, int job_ID, int noise ) {
   ofstream FHW;
   FHW.open(P->p["-o"]);
   FHW << P->get_header(2);
   double penality   = stod(P->p["-ms_pen"]);
   for (int i = 0; i < fits.size(); i++) {
      fits[i]->get_model(penality);
      FHW << fits[i]->write();
   }
   FHW.flush();
}

/** Writes a set of bidirectional model fits utilizing a penalty value passed from the command line (the -ms_pen parameter at the time of this writing).
 * 
 * @param fits Set of models to write.
 * @param pw ParamWrapper object from which to read command line parameters.
 * @param job_ID Job identifier most likely obtained from the MPI runtime.
 * @param noise Unused.
 */
void load::write_out_bidirectionals_ms_pen_pwrapper(vector<segment_fits*> fits, ParamWrapper * pw, int job_ID, int noise ) {
   ofstream FHW;
   FHW.open(pw->outputDir);
   FHW << pw->getHeader(2);
   FHW << "#chrom\tstart\tend\tbidir\tmu\tsigma\tomega\tlamdba\tcov_pos\tcov_neg" << endl;
   double penality   = pw->penalty;
   for (int i = 0; i < fits.size(); i++) {
      fits[i]->get_model(penality);
      FHW << fits[i]->write();
   }
   FHW.flush();
}

//================================================================================================
//misc.

/** Bins a set of segments.
 * 
 * @param segments Set of segments to bin.
 * @param BINS Number of bins to use.
 * @param scale Scaling factor to use in the binning process.
 * @param erase Whether or not to delete existing data after binning (?)
 */
void load::BIN(vector<segment*> segments, int BINS, double scale, int erase) {
   for (int i = 0 ; i < segments.size() ; i ++) {
      if (segments[i]->forward.size() > 0 or segments[i]->reverse.size() > 0 ) {
         segments[i]->bin(BINS, scale, erase);
      }
   }
}

/** Deletes all temporary files used in the fitting and model output process.
 * 
 * @param dir Directory in which temporary files were created.
 * @param job_name User-specified name for the job.
 * @param nprocs Number of workers obtained from the MPI runtime.
 * @param job_ID Job identifier likely obtained from the MPI runtime.
 */
void load::collect_all_tmp_files(string dir, string job_name, int nprocs, int job_ID) {
   int c   = 0;
   time_t     now = time(0);
   struct tm  tstruct;
   char       buf[80];
   tstruct = *localtime(&now);
   // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
   // for more information about date/time format
   strftime(buf, sizeof(buf), "%m%d_%H%M", &tstruct);
   string DT     = buf;
   string OUT    = dir + job_name + "-" + to_string(job_ID) + "_" + DT + ".log";
   ofstream FHW(OUT);
   for (int rank = 0; rank < nprocs; rank++) {
      string FILE   = dir + "tmp_" + job_name + "-" + to_string(job_ID) + "_" + to_string(rank) + ".log";
      string line;
      ifstream FH(FILE);
      if (FH) {
         if (rank != 0) {
            FHW << "=======MPI Call: " + to_string(rank) + "=======\n";
         }
         while (getline(FH, line)) {
            if ("#" != line.substr(0, 1) or rank == 0) {
               FHW << line << endl;
            }
         }

         FH.close();
         remove( FILE.c_str()) ;
         c++;
      }
   }
}

/** Deletes all segments in the specified vector.
 * 
 * @param segments Set of segments to delete.
 */
void load::clear_segments(vector<segment *> segments) {
   for (int i = 0; i < segments.size(); i++) {
      if (segments[i] != NULL) {
         delete (segments[i]);
      }
   }
}

/** Loads a TSS input file.
 * @param tss_file Input TSS training file.
 * @param query_fits Set of segment_fits to match against the given TSS inputs.
 * @return New set of segment_fits with TSS values inserted.
 */
vector<segment_fits *> load::label_tss(string promoterTSS, vector<segment_fits *> query_fits ) {
   vector<segment_fits *> new_fits;
   ifstream FH(promoterTSS);
   string line;
   map<string, vector<segment *> >G;
   map<string, node> T;
   string chrom, start, stop;
   typedef map<string, vector<segment *>>::iterator it_type;
   if (FH) {
      while (getline(FH, line)) {
         vector<string> lineArray  = split_by_tab(line, "");
         chrom   = lineArray[0], start = lineArray[1], stop = lineArray[2];
         segment * S   = new segment(chrom, stoi(start), stoi(stop));
         G[chrom].push_back(S);
      }
      //make G a node interval tree
      for (it_type c = G.begin(); c != G.end(); c++) {
         T[c->first]   = node(c->second);
      }
      //label
      for (int i = 0 ; i < query_fits.size(); i++) {
         vector<int> FINDS;
         if (T.find(query_fits[i]->chrom) != T.end()) {
            T[query_fits[i]->chrom].searchInterval(query_fits[i]->start, query_fits[i]->stop, FINDS);
            if (!FINDS.empty() and query_fits[i]->M[1] != nINF and query_fits[i]->M[1] > query_fits[i]->M[0]   ) {
               query_fits[i]->TSS = 1;
               new_fits.push_back(query_fits[i]);
            }
         }
         FINDS.clear();

      }

   } else {
      printf("Couldn't open bidirectional file %s. Check that file exists and is in standard BED format (chr, start, end).\n", promoterTSS.c_str() );
   }
   return new_fits;
}
















