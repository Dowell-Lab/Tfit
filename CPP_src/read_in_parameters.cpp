#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include "read_in_parameters.h"
#include <stdio.h>
#include <ctype.h>
#include <cctype>
#include <stdio.h>
#include <time.h>
#include "split.h"
#ifdef USING_ICC
#include <aligned_new>
#endif
using namespace std;

params::params(){
	p["-N"] 		= "EMG";
	p["-v"] 		= "1";
	
	p["-h"] 		= "";
	p["--help"] 	= "";

	p["-config"] 	= "";
	p["-i"] 		= "";
	p["-j"] 		= "";
	p["-k"] 		= "";
	p["-tss"] 		= "";
	p["-o"] 		= "";
	p["-q"] 		= "";
	p["-log_out"] 	= "";

	p["-merge"] 	= "0";
	p["-pad"] 		= "0";
	p["-br"] 		= "300";
	p["-ns"] 		= "100";
	p["-minK"] 		= "1";
	p["-maxK"] 		= "3";
	p["-rounds"] 	= "10";
	p["-ct"] 		= "0.0001";
	p["-bct"] 		= "0.95";
	p["-ms_pen"] 	= "1";
	p["-MLE"] 		= "0";
	p["-select"] 	= "0";
	p["-max_noise"] = "0.05";
	p["-chr"] 		= "all";
	p["-elon"] 		= "0";
	p["-mi"] 		= "2000";
	p["-r_mu"] 		= "0";
	
	p["-ALPHA_0"] 	= "1";
	p["-BETA_0"] 	= "1";
	p["-ALPHA_1"] 	= "1";
	p["-BETA_1"] 	= "1";
	p["-ALPHA_2"] 	= "1";
	p["-ALPHA_3"] 	= "1";
		
	
	N 				= 0;
	module 			= "";
	EXIT 			= 0;
	bidir 			= 0;
	model 			= 0;
	select 			= 0;

}
bool is_decimal(const std::string& s){
   if(s.empty() || std::isspace(s[0]) || std::isalpha(s[0])) return false ;
   char * p ;
   strtod(s.c_str(), &p) ;
   return (*p == 0) ;
}
bool is_integer(const std::string& s){
	std::string::const_iterator it = s.begin();
	while (it != s.end() && std::isdigit(*it)) ++it;
	return !s.empty() && it == s.end();
}

bool is_path(string FILE){
	ifstream FH(FILE);
	if (FH){
		return true;
	}else{
		return false;
	}

}

vector<string> params::validate_parameters(){
	vector<string> errors;
	for (int i = 0; i < 12; i++){
		if (i <8 and not  is_integer(p[isIntGroup[i]])  ){
			string line = "User provided input for (" + string(isIntGroup[i]) + ") "  ;
			line+= + "'"+string(p[isIntGroup[i]])+ "'"+ " is not integer valued";
			errors.push_back(line);
		}
		if (not  is_decimal(p[isDecGroup[i]])){
			string line = "User provided input for (" + string(isDecGroup[i]) + ") ";
			line+= + "'"+ string(p[isDecGroup[i]] ) + "'" + " is not decimal valued";
			errors.push_back(line);
		}
		if (i < 8){
			string path_FILE 	= p[isPathGroup[i]];
			if (isPathGroup[i] == "-o" and !path_FILE.empty() and path_FILE.substr(path_FILE.size()-1, 1) != "/"  ){
				p[isPathGroup[i]] 	= path_FILE+ "/";
			}
		}
		if (i < 8 and not is_path(p[isPathGroup[i]])){

			if (!(bidir and isPathGroup[i] == "-k" ) and (!((bidir or model ) and isPathGroup[i]=="-tss")) and (!((bidir or model ) and isPathGroup[i]=="-q"))    ){
				string line = "User provided input for (" + string(isPathGroup[i]) + ") ";
				line+= + "'"+ string(p[isPathGroup[i]]) + "'"  + " path does not exist";
				errors.push_back(line);
			}
		}
	}
	return errors;
}

const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m/%d/%Y %X", &tstruct);

    return buf;
}
void params::help(){
	string header 	= "";
	header+="--------------------------------------------------------------------------------------\n";
	header+="                        transcriptional inference (tINF)                 \n";
	printf("%s\n",header.c_str() );
	printf("                   ....description of application modules....             \n");
	printf("                                  (critical)                            \n\n");
	
	printf("bidir     : must be provided immediatly following the application call \"EMGU\"\n");
	printf("              to scan for bidirectional signal across provided dataset, uses\n");
	printf("              a poisson noise background model to assess enrichment and \n");
	printf("              quick moment estimation to asses divergent transcription shape \n");
	printf("model     : must be provided immediatly following the application call \"EMGU\"\n");
	printf("              to perform maximum likelihood or a-posteriori parameter inference\n");
	printf("              recommended for accruacy; especially for point estimate on\n");
	printf("              RNA polymerase II loading position needed for TF ID-ing \n");
	printf("select    : must be provided immediatly following the application call \"EMGU\"\n");
	printf("              model selection is performed via penalized bayesian information\n");
	printf("              criteria. To set the penality, we consider an ROC curve over signal\n");
	printf("              with no moment estimate prediction (TN) and bidirectionals found near TSS (TP) \n");
	
	printf("\n\n");
	header="";
	header+="                ....description of non-default parameters....          \n";
	header+="                                 (critical)                             \n";
	printf("%s\n",header.c_str() );
	printf("-N        : Name of Job; default is \"EMG\"\n");
	printf("-i        : /path/to/forward/strand/bedgraph/file\n");
	printf("              this file should be bedgraph formatted\n");
	printf("              chromosome[tab]start[tab]stop[tab]coverage[newline]\n");
	printf("-j        : /path/to/reverse/strand/bedgraph/file\n");
	printf("              this file should be bedgraph formatted\n");
	printf("              chromosome[tab]start[tab]stop[tab]coverage[newline]\n");
	printf("-k        : /path/to/interval/file\n");
	printf("              this file should be formatted as a bed file\n");
	printf("              chromosome[tab]start[tab]stop[newline]\n");
	printf("              for each provided interval we will fit a mixture\n");
	printf("              from -minK components to -maxK components\n");
	printf("-o        : /path/to/output/directory/\n");
	printf("              the \"bidir\" module will create a file in this\n");
	printf("              directory {-N}_prelim_bidir_hits.bed\n");
	printf("              the \"model\" module will create two files in this\n");
	printf("              directory {-N}_K_models_MLE.tsv and\n");
	printf("              {-N}_bidirectional_hits_intervals.bed\n");
	printf("-log_out  : /path/to/directory/where/tmp/log_files/will/be/stored\n");
	printf("              as the application runs, this directory will house tmp_{-N}.log\n");
	printf("              log files for each MPI process (and they will be removed at the end)\n");
	printf("              these will be updated and will give application progress\n");
	printf("\n");
	printf("                    ....description of default parameters....          \n");	
	printf("                                (non-critical)                         \n");
	printf("\n");


	printf("-chr      : (chromosome ID; i.e. chr1) specific chromosome to run on (default is \"all\")\n");
	printf("-merge    : (boolean integer) will merge overlaping intervals and run model on joint\n");
	printf("               recommended for bidirectional de novo not recommended for evaluating \n");
	printf("               gene intervals, default = 0\n");
	
	printf("-elon     : (boolean integer) adjust support of elongation component, (default=0)\n");
	printf("              useful only when fitting to FStitch[1] or groHMM[2] output intervals\n");
	printf("-pad      : (positive integer) each provided interval will be extended\n");
	printf("              in both the five-prime and three-prime direction (default=1000)\n");
	printf("-MLE      : (boolean integer) specific to the bidir module, will perform parameter\n");
	printf("              inference via EM (highly recommended for accuracy)\n");
	printf("-select   : (boolean integer) specific to the bidir and model module will perform\n");
	printf("              ROC analysis and run the selection module\n");
	printf("-ms_pen   : (positive floating) penality term in BIC criteria for model selection\n");
	printf("              (default = 1)\n");
	
	printf("\n");
	printf("                    ....description of default parameters....          \n");	
	printf("                        (non-critical, advanced useage)                 \n");
	printf("\n");
	printf("-minK     : (positive integer) minimum number of model components to try\n");	
	printf("              (default=1)\n");
	printf("-maxK     : (positive integer) maximum number of model components to try\n");
	printf("              (default=5)\n");
	printf("-rounds   : (positive integer) number of random initializations to the\n");
	printf("              EM algorithm (default=10)\n" );	                  
	printf("-mi       : (positive integer) maximum number of iterates to the\n");
	printf("              EM algorithm (default=2000)\n" );	                  
	printf("-ct       : (positive decimal) EM log-likelihood convergence threshold\n");
	printf("              (default=0.0001)\n" );	                  
	printf("-ALPHA_0  : hyperparameter (1) for the Normal Inverse Wishart prior for loading variance (sigma)\n" );	                  
	printf("              (default=1; weak)\n" );	                  
	printf("-BETA_0   : hyperparameter (2) for the Normal Inverse Wishart fprior for loading variance (sigma)\n" );	                  
	printf("              (default=1; weak)\n" );	                  
	printf("-ALPHA_1  : hyperparameter (1) for the Gamma prior for initiating length (lambda)\n" );	                  
	printf("              (default=1; weak)\n" );	                  
	printf("-BETA_1   : hyperparameter (2) for the Gamma prior for initiating length (lambda)\n" );	                  
	printf("              (default=1; weak)\n" );	                  
	printf("-ALPHA_2  : hyperparameter (symmetric) for the Dirichlet prior for component mixing weights\n" );	                  
	printf("              (default=1; weak)\n" );	                 
	printf("-ALPHA_3  : hyperparameter (symmetric) for the Beta prior for strand bias\n" );	                  
	printf("              (default=1; weak)\n" );	                 
	 
	
	printf("\n\n");
	printf("-config   :  all parameters may be specified in a config file with     \n");	
	printf("              flag = value [newline] syntax, parameters specified on the command line\n" );	                 
	printf("              following the -config flag will override any parameters given in the config file\n" );	                 
	
	
	printf("\nQuestions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu\n" );
		
	printf("--------------------------------------------------------------------------------------------\n");

	
	

	
	



	
	


	
	



}
void params::display(int nodes, int cores){
	bool MLE 	= stoi(p["-MLE"]);
	bool SELECT = stoi(p["-select"]);
	string header 	= "";
	header+="----------------------------------------------------------------\n";
	header+="             transcriptional inference (tINF)                   \n";
	if (bidir){
	header+="               ...bidirectional scanner...     \n";
	}if (bidir and MLE){
	header+="             ....coupled to mixture model....     \n";
	}if (bidir and select){
	header+="       ....coupled to BIC penality optimization....     \n";
	}
	if (model){
	header+="               ...running mixture model...                      \n";
	}
	if (select){
	header+="            ....BIC penality optimization....                      \n";		
	}
	printf("%s\n",header.c_str() );
	printf("-N         : %s\n", p["-N"].c_str()  );
	printf("-i         : %s\n", p["-i"].c_str()  );
	printf("-j         : %s\n", p["-j"].c_str()  );
	printf("-MLE       : %s\n",  p["-MLE"].c_str());
	printf("-select    : %s\n",  p["-select"].c_str());
	if (model or MLE){
	printf("-k         : %s\n", p["-k"].c_str()  );
	}

	if (select){
		printf("-q         : %s\n", p["-q"].c_str());
		printf("-tss       : %s\n", p["-tss"].c_str()  );
	}
	printf("-o         : %s\n", p["-o"].c_str()  );
	printf("-log_out   : %s\n", p["-log_out"].c_str()  );
	printf("-chr       : %s\n", p["-chr"].c_str());	
	printf("-br        : %s\n", p["-br"].c_str());	
	printf("-rounds    : %s\n", p["-rounds"].c_str()  );
	printf("-elon      : %s\n", p["-elon"].c_str()  );
	printf("-pad       : %s\n", p["-pad"].c_str()  );
	if (!model){
	printf("-bct       : %s\n", p["-bct"].c_str());
	}
	printf("-minK      : %s\n", p["-minK"].c_str());
	printf("-maxK      : %s\n", p["-maxK"].c_str());
	printf("-threads   : %d\n",  cores);
	printf("-MPI_np    : %d\n",  nodes);
	printf("Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu\n" );
	printf("----------------------------------------------------------------\n" );
	
}


string params::get_header(int ID){
	string header = "";
	header+="#----------------------------------------------------\n";
	header+="#Date Time    : "+currentDateTime()+"\n";
	header+="#-N           : "+p["-N"]+"\n";
	header+="#-i           : "+p["-i"]+"\n";
	header+="#-j           : "+p["-j"]+"\n";
	header+="#-k           : "+p["-k"]+"\n";
	header+="#-o           : "+p["-o"]+"\n";
	header+="#-ns          : "+p["-ns"]+"\n";
	header+="#-br          : "+p["-br"]+"\n";
	header+="#-bct         : "+p["-bct"]+"\n";
	header+="#-elon        : "+p["-elon"]+"\n";
	header+="#-pad         : "+p["-pad"]+"\n";
	header+="#-minK        : "+p["-minK"]+"\n";
	header+="#-maxK        : "+p["-maxK"]+"\n";
	header+="#-ms_pen      : "+p["-ms_pen"]+"\n";
	header+="#-mi          : "+p["-mi"]+"\n";
	header+="#-ct          : "+p["-ct"]+"\n";
	header+="#-rounds      : "+p["-rounds"]+"\n";
	header+="#-ALPHA_0     : "+p["-ALPHA_0"]+"\n";	
	header+="#-BETA_0      : "+p["-BETA_0"]+"\n";	
	header+="#-BETA_1      : "+p["-BETA_1"]+"\n";	
	header+="#-ALPHA_2     : "+p["-ALPHA_2"]+"\n";	
	header+="#-ALPHA_3     : "+p["-ALPHA_3"]+"\n";	
	header+="#----------------------------------------------------\n";
	
	return header;
}



void split_config_line(string line, string& param, string& value){
	bool   S=false;
	for (int i = 0 ; i < line.size(); i++){
		if (line.substr( i, 1) == "="  ){
			S 	= true;
		}
		else if ( not S and   !isspace(line[i])   ){
			param+=line[i];
		}else if(S and line.substr( i, 1) == "#" ){
			break;
		}else if (S and !isspace(line[i])   ){
			value+=line[i];
		}
	}
}

void fill_in_config_file(string FILE, params * P, int rank){
	ifstream FH(FILE);
	if (FH){
		string line;
		while (getline(FH, line)){
			if (line.substr(0,1)!= "#" and line.substr(0,1)!= "~" and !line.empty()){
				string param="", value="";
				split_config_line(line, param, value);
				if (!param.empty() and param.substr(0,1)!= "-" and P->p.find(param)==P->p.end()){
					if (rank == 0){
						printf("option in config file %s is not valid\n", param.c_str() );
					}
					P->EXIT=1;
				}else{
					P->p[param]=value;
				}
			
			}
		}
	}else{
		if (rank==0){
			printf("could not open config file %s\n",FILE.c_str() );
		}
		P->EXIT 	= 1;
	}
}

void fill_in_options(char* argv[],params * P, int rank){
	bool bidir 		= P->bidir;
	bool model 		= P->model;
	string F 		= "";
	char * COM 		= "-";
		
	while (*argv){
		if ((*argv)[0] == COM[0]){
			F 	= string(*argv); 
			if ( P->p.find(F) ==P->p.end() ){
				if (rank == 0){
					printf("Unknown user option: %s\n", F.c_str() );
				}
				P->EXIT 	= 1;
				break;	
			}
			if (F.substr(0,2) == "-h" or F.substr(0,7)=="--help" ){
				P->EXIT=true;
				P->help();
				break;
			}
		
		}
		else if (not F.empty()) {
			if ((model or bidir or select) && P->p.find(F) !=P->p.end()){
				P->p[F]=string(*argv);
				if (F=="-config"){
					fill_in_config_file(P->p[F], P, rank );
				}		
				F="";
			}
		}
		argv++;
	}
}




int read_in_parameters( char* argv[], params * P, int rank){	
	string userModParameter = "";
	argv = ++argv;
	if (not *argv){
		if (rank==0){
			printf("No Model found, please specify either bidir, model or select\n");
		}
		P->EXIT = 1;
		return 1;
	}

	if (not P->EXIT){
		string F 	= *argv;
		if ((F.size()==2 and F.substr(0,2) == "-h") or (F.size()==6 and F.substr(0,6)=="--help") ){

			P->EXIT=true;
			P->help();
			return 1;
		}
		if (F.size()==5 and F.substr(0,5) == "bidir") {
			P->bidir 	= 1;
	 	}
		else if (F.size()==5 and  F.substr(0,5) == "model") {
			P->model 	= 1;
		}
		else if(F.size() == 6 and F.substr(0,6)=="select"){
			P->select 	= 1;
		}
		else{
			if (rank == 0){
				printf("couldn't understand user provided option: %s\n",F.c_str() );
			}
			P->EXIT = 1;
		}
		if (not P->EXIT){
			argv = ++argv;
			fill_in_options(argv, P, rank);
		}
	}
	if (not P->EXIT){
		vector<string> errors = P->validate_parameters();
		if (!errors.empty()){
			P->EXIT=1;
		}
		if (rank == 0){
			for (int i = 0 ; i < errors.size(); i++){
				printf("%s\n",errors[i].c_str() );
			}
		}
	
	}
	return 1;
}



















