#include "ParamWrapper.hpp" 


void ParamWrapper::printUsage()
{
    printf("Usage: TFit modulename [arguments]\n");
    printf("Where modulename is one of the following:\n");
    printf("\tbidir - This module searches the genome for areas resembling bidirectional transcription\n");
    printf("\t\tby comparing a fixed template mixture model to a noise model by a log-likelihood ratio score.\n");
    printf("\tmodel - This module attempts to generate an optimal set of parameters per region instead \n");
    printf("\t\tof using a fixed set of parameters for the entire genome\n");
    
    printf("\n\nWhere [arguments] is one or more of the following for the bidir module:\n");
    printf("\t-i\tForward bedgraph file\n");
    printf("\t-j\tReverse bedgraph file.\n");
    printf("\t-ij\tBoth forward and reverse bedgraph file. This parameter may be used in place of -i and -j if reads are in one bedgraph file.\n");
    printf("\t-N\tJob name.\n");
    printf("\t-o\tOutput directory. If it does not exist, it will be created.\n");
    printf("\t-log_out\tLog file output directory. If it does not exist, it will be created.\n");
    printf("\nAdditional (optional) parameters for the bidir module:\n");
    printf("\t-tss\tTranscription model path. Models are provided for hg19 and mm10 in the annotations directory of this project.\n");
    printf("\t-chr\tRun bidir only on the specified chromosome by name. The default is, \"all\"\n");
    printf("\t-bct\tLLR threshold. Default=1\n");
    printf("\nIf -tss is not specified, the values normally specified by the model file can be specified explicitly using the following:\n");
    printf("\t-lambda\tThis is the entry length parameter for the EMG density function. default=200 bp\n");
    printf("\t-sigma\tThis is the variance parameter for the EMG density function. default=10bp\n");
    printf("\t-pi\tThis is the strand bias parameter for the EMG density function. default=0.5\n");
    printf("\t-w\tThis is the pausing probability parameter for the EMG denisty function. default=0.5\n");
    
    printf("\n\nWhere [arguments] is one or more of the following for the model module:\n");
    printf("\t-i\tForward bedgraph file\n");
    printf("\t-j\tReverse bedgraph file\n");
    printf("\t-ij\tBoth forward and reverse bedgraph file. This parameter may be used in place of -i and -j if reads are in one bedgraph file.\n");
    printf("\t-k\tBedgraph file containing a set of regions of interest.\n");
    printf("\t-N\tJob name.\n");
    printf("\t-o\tOutput directory. If it does not exist, it will be created.\n");
    printf("\t-log_out\tLog file output directory. If it does not exist, it will be created.\n");
    printf("\nAdditional (optional) parameters for the model module:\n");
    printf("\t-mink\tMinimum number of finite mixtures to consider. default=1\n");
    printf("\t-maxk\tMaximum number of finite mixtures to consider. default=1\n");
    printf("\t-rounds\tNumber of random seeds to use in the model. default=5\n");
    printf("\t-ct\tConvergence threshold after which processing stops. default=0.0001\n");
    printf("\t-mi\tMaximum number of model iterations after which processing stops. default=2000\n");
    printf("The model module currently has experimental support for parameter inference. Ie. it can attempt to estimate lambda,\n");
    printf("sigma, pi, and w via conjugate priors. These values are specified using the following parameters:\n");
    printf("\t-ALPHA_0\tPrior parameter 1 for sigma. Recommended value=1\n");
    printf("\t-BETA_0\tPrior parameter 2 for sigma. Recommended value=1\n");
    printf("\t-ALPHA_1\tPrior parameter 1 for lambda.\n");
    printf("\t-BETA_1\tPrior parameter 2 for lambda.\n");
    printf("\t-ALPHA_2\tSymmetric prior on mixing weights. Higher values=stronger attempt to find components of equal mixing weights.\n");
    printf("\t\t\tRecommended value=100\n");
    printf("\t-ALPHA_3\tSymmetric prior on the strand bias. Higher values=stronger attempt to find bidirectional events with equal strand bias.\n");
    printf("\t\t\tRecommended value=100\n");
    
    printf("-elon     : (boolean integer) adjust support of elongation component, (default=0)\n");
	printf("              useful only when fitting to FStitch[1] or groHMM[2] output intervals\n");
}

ParamWrapper::ParamWrapper()
{
    this->module="";
    //This constructor sets blank or default values.
    this->forwardStrand="";
    this->reverseStrand="";
    this->mergedStrand="";
    this->jobName="";
    this->outputDir="";
    this->logDir="";
    this->regionsOfInterest="";
    this->promoterTSS="";
    this->chromosome="";
    this->llrthresh=1;
    this->lambda=200;
    this->sigma=10;
    this->pi=0.5;
    this->w=0.5;
    this->mink=1;
    this->maxk=1;
    this->rounds=5;
    this->ct=0.0001;
    this->mi=2000;
    this->ns=100;
    this->experimentalValsSpecified=false;
    this->alpha0=1;
    this->beta0=1;
    //There are no default values for this parameter.
    this->alpha1=0;
    this->beta1=0;
    this->alpha2=100;
    this->alpha3=100;
    this->br=25;
    this->fdr=0;
}

ParamWrapper::ParamWrapper(int argc, char **argv)
{
    int i;
    map<string, string> paramMap;
    char *prevCmd;
    map<string, string>::iterator it;
    
    this->exit=false;
    this->verbose=false;
    this->module="";
    //This constructor sets blank or default values.
    this->forwardStrand="";
    this->reverseStrand="";
    this->mergedStrand="";
    this->jobName="";
    this->outputDir="";
    this->logDir="";
    this->regionsOfInterest="";
    this->promoterTSS="";
    this->chromosome="";
    this->llrthresh=1;
    this->lambda=200;
    this->sigma=10;
    this->pi=0.5;
    this->w=0.5;
    this->mink=1;
    this->maxk=1;
    this->rounds=5;
    this->ct=0.0001;
    this->mi=2000;
    this->experimentalValsSpecified=false;
    this->alpha0=1;
    this->beta0=1;
    this->ns=100;
    //There are no default values for this parameter.
    this->alpha1=1;
    this->beta1=1;
    this->alpha2=100;
    this->alpha3=100;
    this->br=25;
    this->pad=2000;
    this->footPrint=86;
    this->fdr=0; //This is an undocumented parameter.
    
    if(argc==1)
    {
        this->printUsage();
        this->exit=true;
        return;
    }
    
    else if(argc==2)
    {
        printf("Error: no arguments specified for module %s\n", argv[1]);
        this->printUsage();
        this->exit=true;
        return;
    }
    
    this->module=string(argv[1]);
    
    if(this->module!="bidir" && this->module!="select" && this->module!="model")
    {
        printf("Invalid module specification: %s\n", argv[1]);
        this->exit=true;
        return;
    }
    
    this->bidir=this->module=="bidir";
    this->model=this->module=="model";
    this->select=this->module=="select";
    
    for(i=2;i<argc;i++)
    {
        if(argv[i][0]=='-' && strlen(argv[i])!=1)
        {
            prevCmd=argv[i];
            //Check to ensure that the argument is supposed to be able to work without additional parameters:
            if(!strcmp(prevCmd, "-v"))
            {
                this->verbose=true;
            }
            
            else if((!strcmp(prevCmd, "-h"))||(!strcmp(prevCmd, "--help")))
            {
                this->printUsage();
                this->exit=true;
                return;
            }
        }
        
        else
        {
            paramMap.insert(pair<string, string>(prevCmd, argv[i]));
            prevCmd=NULL;
        }
    }
    
    for(it=paramMap.begin();it!=paramMap.end();it++)
    {
        if(it->first=="i")
        {
            this->forwardStrand=it->second;
        }
        
        else if(it->first=="-j")
        {
            this->reverseStrand=it->second;
        }
        
        else if(it->first=="-ij")
        {
            this->mergedStrand=it->second;
        }
        
        else if(it->first=="-N")
        {
            this->jobName=it->second;
        }
        
        else if(it->first=="-ns")
        {
            this->ns=atoi(it->second.c_str());
        }
        
        else if(it->first=="-o")
        {
            this->outputDir=it->second;
        }
        
        else if(it->first=="-fdr")
        {
            this->fdr=atoi(it->second.c_str());
        }
        
        else if(it->first=="-log_out")
        {
            this->logDir=it->second;
        }
        
        else if(it->first=="-tss")
        {
            this->promoterTSS=it->second;
        }
        
        else if(it->first=="-chr")
        {
            this->chromosome=it->second;
        }
        
        else if(it->first=="-bct")
        {
            this->llrthresh=atoi(it->second.c_str());
        }
        
        else if(it->first=="-lambda")
        {
            this->lambda=atof(it->second.c_str());
        }
        
        else if(it->first=="-sigma")
        {
            this->sigma=atof(it->second.c_str());
        }
        
        else if(it->first=="-pi")
        {
            this->pi=atof(it->second.c_str());
        }
        
        else if(it->first=="-w")
        {
            this->w=atof(it->second.c_str());
        }
        
        else if(it->first=="-mink")
        {
            this->mink=atoi(it->second.c_str());
        }
        
        else if(it->first=="-pad")
        {
            this->pad=atof(it->second.c_str());
        }
        
        else if(it->first=="-maxk")
        {
            this->maxk=atoi(it->second.c_str());
        }
        
        else if(it->first=="-rounds")
        {
            this->rounds=atoi(it->second.c_str());
        }
        
        else if(it->first=="-ct")
        {
            this->ct=atof(it->second.c_str());
        }
        
        else if(it->first=="-mi")
        {
            this->mi=atoi(it->second.c_str());
        }
        
        else if(it->first=="-ALPHA_0")
        {
            this->experimentalValsSpecified=true;
            this->alpha0=atof(it->second.c_str());
        }
        
        else if(it->first=="-BETA_0")
        {
            this->experimentalValsSpecified=true;
            this->beta0=atof(it->second.c_str());
        }
        
        else if(it->first=="-foot_print")
        {
            this->footPrint=atoi(it->second.c_str());
        }
        
        else if(it->first=="-ALPHA_1")
        {
            this->experimentalValsSpecified=true;
            this->alpha1=atof(it->second.c_str());
        }
        
        else if(it->first=="-BETA_1")
        {
            this->experimentalValsSpecified=true;
            this->beta1=atof(it->second.c_str());
        }
        
        else if(it->first=="-ALPHA_2")
        {
            this->experimentalValsSpecified=true;
            this->alpha2=atof(it->second.c_str());
        }
        
        else if(it->first=="-ALPHA_3")
        {
            this->experimentalValsSpecified=true;
            this->alpha3=atof(it->second.c_str());
        }
        
        else if(it->first=="-MLE")
        {
            this->mle=atoi(it->second.c_str());
        }
        
        else if(it->first=="-elon")
        {
            this->elon=atoi(it->second.c_str());
        }
    }
    
    if((i=="" || j=="") && ij=="")
    {
        printf("Either i or j empty. Please specify i or j.\n");
        this->exit=true;
    }
    
    //Todo: add more sanity checks.
}

//Functions implemented to reach feature parity with read_in_parameters:
void ParamWrapper::display(int nodes, int cores){
	//We have no analogs, since these parameters aren't in the repo documentation:
    //bool MLE 	= stoi(p["-MLE"]);
	//bool SELECT = stoi(p["-select"]);
	string header 	= "";
	header+="----------------------------------------------------------------\n";
	header+="             transcriptional inference (tINF)                   \n";
	if (bidir){
	header+="               ...bidirectional scanner...     \n";
	}//if (this->module=="bidir" and MLE){
	//header+="             ....coupled to mixture model....     \n";
	//}
	if (bidir and select){
	header+="       ....coupled to BIC penalty optimization....     \n";
	}
	if (model){
	header+="               ...running mixture model...                      \n";
	}
	if (select){
	header+="            ....BIC penalty optimization....                      \n";		
	}
	printf("%s\n",header.c_str() );
	printf("-N         : %s\n", this->jobName.c_str());
	if (not p["-ij"].empty()){
		printf("-ij        : %s\n", this->mergedStrand.c_str());
	}else{
		printf("-i         : %s\n", this->forwardStrand.c_str());
		printf("-j         : %s\n", this->reverseStrand.c_str());
	}
	//if (model or MLE){
	//printf("-k         : %s\n", p["-k"].c_str()  );
	//}
	if (this->promoterTSS!=""){
		printf("-tss       : %s\n", promoterTSS.c_str()  );
	}
	printf("-o         : %s\n", this->outputDir.c_str()  );
	printf("-log_out   : %s\n", this->logDir.c_str()  );
	printf("-MLE       : %d\n", this->mle);
	printf("-chr       : %s\n", this->chromosome.c_str());	
	printf("-br        : %d\n", this->br);	
	printf("-rounds    : %d\n", this->rounds);
	if (this->elon){
		printf("-elon      : %d\n", this->elon);
	}
	printf("-pad       : %d\n", this->pad);
	if (!model){
	printf("-bct       : %s\n", this->llrthresh.c_str());
	}
	if (model){
		printf("-minK      : %d\n", this->mink);
		printf("-maxK      : %d\n", this->maxk);
	}
	printf("-threads   : %d\n",  cores);
	printf("-MPI_np    : %d\n",  nodes);
	printf("\nQuestions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu\n" );
    printf("\nRevisions made by michael[dot]gohde[at]colorado[dot]edu\n");
	printf("----------------------------------------------------------------\n" );
	
}

string ParamWrapper::get_header(int ID){
	string header = "";
	header+="#----------------------------------------------------\n";
	header+="#Date Time    : "+currentDateTime()+"\n";
	header+="#-N           : "+this->jobName+"\n";
	if (p["-ij"].empty()){
		header+="#-i           : "+this->forwardStrand+"\n";
		header+="#-j           : "+this->reverseStrand+"\n";
	}else{
		header+="#-ij          : "+this->mergedStrand+"\n";
	}
	if (ID!=1){
		header+="#-k           : "+this->regionsOfInterest+"\n";
	}
	header+="#-o           : "+this->outputDir+"\n";
	header+="#-br          : "+this->br+"\n";
	if (ID==1){
		header+="#-bct         : "+this->llrthresh+"\n";
		header+="#-pad         : "+this->pad+"\n";
	}
	if (ID!=1){
		header+="#-elon        : "+this->elon+"\n";
		header+="#-minK        : "+this->mink+"\n";
		header+="#-maxK        : "+this->maxk+"\n";
		header+="#-mi          : "+this->mi+"\n";
		header+="#-ct          : "+this->ct+"\n";
		header+="#-rounds      : "+this->rounds+"\n";
	}
	if (ID!=1){
		header+="#-ALPHA_0     : "+this->alpha0+"\n";	
		header+="#-BETA_0      : "+this->beta0+"\n";	
		header+="#-BETA_1      : "+this->beta1+"\n";	
		header+="#-ALPHA_2     : "+this->alpha2+"\n";	
		header+="#-ALPHA_3     : "+this->alpha3+"\n";	
	}
	if (ID==1){
		header+="#-sigma       : "+ this->sigma +"\n";
		header+="#-lambda      : "+ this->lambda+"\n";
		header+="#-foot_print  : "+ this->footPrint+"\n";
		header+="#-pi          : "+ this->pi+"\n";
		header+="#-w           : "+ this->w+"\n";
	}
	header+="#----------------------------------------------------\n";
	return header;
}
