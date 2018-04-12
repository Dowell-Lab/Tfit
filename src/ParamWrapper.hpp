#ifndef PARAMWRAPPER_HPP
#define PARAMWRAPPER_HPP
#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <map>
#include <vector>

#include <sstream>

//using namespace std;

/* This class is similar to the ParamWrapper used by FStitch. It sgoal is to make it easier to
 * change, add, and replace command line options in Tfit. 
 */
 
class ParamWrapper
{
private:
public:
    bool exit;
    bool verbose;
    
    //Use a set of booleans like read_in_parameters to represent modules for ease of porting:
    bool select;
    bool bidir;
    bool model;
    
    //Set of parameters:
    std::string module;
    std::string forwardStrand;
    std::string reverseStrand;
    std::string mergedStrand;
    
    std::string jobName;
    std::string outputDir;
    std::string logDir;
    
    std::string regionsOfInterest;
    
    //Set of bidir parameters with specific default settings:
    bool mle;
    int br; //This is undocumented in read_in_parameters;
    bool elon;
    double pad; //This is likely padding for ... something.
    std::string promoterTSS;
    std::string chromosome; //Default=all.
    int llrthresh; //LLR threshold specified by -bct default=1
    int ns;
    int fdr;
    
    double lambda; //default=200
    double sigma; //default=10
    double pi; //default=0.5
    double w; //default=0.5
    
    int footPrint; //default=86
    
    //Set of model parameters with specific default settings:
    int mink; //default=1
    int maxk; //default=1
    int rounds; //default=5
    double ct; //default=0.0001
    int mi; //default=2000
    int r_mu; //Default=0; undocumented.
    
    //Additional bleeding edge parameters for the model module:
    bool experimentalValsSpecified;
    double alpha0; //default=1
    double beta0; //default=1
    double alpha1; //prior 1 for lambda. No recommendation or default given.
    double beta1; //prior 2 for lambda. No recommendation or default given.
    double alpha2; //symmetric prior on mixing weights (higher values=more thorough attempt to ... "components of equal mixing weights". default=100
    double alpha3; //symmetric prior on the strand bias (higher value=more thorough attempt to find bidirectional events with equal strand bias. default=100
    
    //Additional parameters that came up when modifying existing code to support ParamWrappers:
    std::string scores;
    int penalty;
    double maxNoise;
    
    //Methods:
    ParamWrapper();
    ParamWrapper(int argc, char **argv);
    void printUsage();
    
    void display(int nodes, int cores);
    std::string getHeader(int id);
};
#endif 
