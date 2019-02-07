/** \brief ParamWrapper parses and preprocesses command line arguments.
 * 
 * The ParamWrapper class was implemented to enable a greater degree of abstraction within tfit
 * with respect to how command line arguments are passed and parsed. Under the previous model,
 * each function that made use of values gathered from the command line would search through a dictionary
 * of argument-value pairs. This meant that any changes to the command line specification of Tfit would
 * require functions throughout the codebase to be revised.
 */
#include "ParamWrapper.hpp" 

/** Prints a usage statement including all currently supported command line arguments.
 */
void ParamWrapper::printUsage()
{
    printf("\n==================================================================================================================\n");
    printf("Transcription Fit (Tfit)\n\n");
    printf("A tool for modeling RNAPII behavior and annotating bidirectionals in nascent transcription data.\n\n");
    printf("For a full description, questions, and error reporting, see:\n");
    printf("\tDowell Lab @ https://github.com/Dowell-Lab/Tfit\n\n");
    printf("==================================================================================================================\n");
    printf("Usage: TFit [modulename] [arguments]\n\n");
    printf("Where modulename is one of the following:\n");
    printf("\tprelim - This module searches the genome for regions resembling bidirectional transcription (active RNAPII)\n");
    printf("\t\tby comparing a fixed template mixture model to a noise model by a log-likelihood ratio score.\n");
    //printf("\tbidir_old - This module implements the functionality seen in versions of Tfit used in publications\n");
    printf("\tmodel - This module models RNAPII behavior given the segment regions obtained from the bidir module.\n");
    printf("==================================================================================================================\n\n");
    printf("Where [arguments] is one or more of the following for the `bidir` module:\n\n");
    printf("Required arguments:\n");
    printf("\t-f   \t--forward   \t<FILE.pos.bedGraph> \t Forward (positive) bedGraph file (only required if -bg unused).\n");
    printf("\t-r   \t--reverse   \t<FILE.neg.bedGraph> \t Reverse (negative) bedGraph file (only required if -bg unused).\n");
    printf("\t-bg  \t--bedgraph  \t<FILE.cat.bedGraph> \t Concatenated forward and reverse bedgraph file.\n");
    printf("\t-N   \t--jobName   \t<MYJOB>             \t Job name for log output (no extensions).\n"); 
    printf("\t-o   \t--output    \t<PRELIMHITS.bed>    \t Output file name (.bed) with directory extension.\n");
    printf("\t-l   \t--logOut    \t</project/logs/>    \t Log file output directory.\n");
    printf("\nAdditional (optional) arguments for the bidir module:\n");
    //This only has an effect (now option -bd) if the -mle is run. We will not give this option yet
    //printf("\t-tss\tTranscription model path. Models are provided for hg19 and mm10 in the annotations directory of this project.\n");
    printf("\t-s   \t--segment   \t<SEGFILE.bed>       \t BED file that specifies sample regions of interest (e.g. FStitch output)\n");
    printf("\t-chr \t--chromosome\t<chrX>              \t Run bidir only on the specified chromosome. Default = all\n");
    printf("\t-n   \t--threads   \t<integer>           \t Number of threads to run the job on. Default=1\n");
    //printf("\t-bct\tLLR threshold. Default=1\n");
    //printf("\nIf -tss is not specified, the values normally specified by the model file can be specified explicitly using the following:\n");
    //printf("\t-lambda\tThis is the entry length parameter for the EMG density function. default=200 bp\n");
    //printf("\t-sigma\tThis is the variance parameter for the EMG density function. default=10bp\n");
    //printf("\t-pi\tThis is the strand bias parameter for the EMG density function. default=0.5\n");
    //printf("\t-w\tThis is the pausing probability parameter for the EMG denisty function. default=0.5\n");
    //printf("\t-scores\tSome form of score output file. This parameter is presently undocumented.\n");
    //printf("\t-r_mu\tSome classification parameter. Default=0. This parameter is presently undocumented.\n");
    //printf("\t-ms_pen\tPenalty term for model selection. Default=1.\n");
    //printf("\t-max_noise\tMaximum noise threshold. Default=0.05. This parameter is presently undocumented.\n");
    //printf("\t-fdr\tGenerate a likelihood score distribution on the input data. This parameter has yet to be fully documented and tested.\n");
    printf("==================================================================================================================\n\n");
    printf("Where [arguments] is one or more of the following for the `model` module:\n\n");
    printf("Required arguments:\n");
    printf("\t-f   \t--forward   \t<FILE.pos.bedGraph> \t Forward (positive) bedGraph file (only required if -bg unused).\n");
    printf("\t-r   \t--reverse   \t<FILE.neg.bedGraph> \t Reverse (negative) bedGraph file (only required if -bg unused).\n");
    printf("\t-bg  \t--bedgraph  \t<FILE.cat.bedGraph> \t Concatenated forward and reverse bedgraph file.\n");
    printf("\t-N   \t--jobName   \t<MYJOB>             \t Job name for log and model output files (no extensions).\n"); 
    printf("\t-o   \t--output    \t<FINALHITS.bed>     \t Output file name (.bed extension).\n");
    printf("\t-l   \t--logOut    \t</projects/logs/>   \t Log file output directory.\n");
    printf("\nAdditional (optional) parameters for the model module:\n");
    printf("\t-s   \t--segment   \t<PRELIMHITS.bed>    \t BED file that specifies sample regions of interest (e.g. bidir output)\n");
    printf("\t-bd  \t--bidirs    \t<FILE.bed>          \t Regions of bidirectionals for fine-tuning parameters.\n");
    printf("\t-chr \t--chromosome\t<chrX>              \t Run model only on the specified chromosome. Default = all\n");
    printf("\t-n   \t--threads   \t<integer>           \t Number of threads to run the job on. Default=1\n");
    printf("==================================================================================================================\n\n");
    printf("\t-h   \t--help      \t Usage information\n");
    printf("\t--version        \t Print version.\n");
 
/* None of these seem to work appropriately and are also very confusing for the typical user. Instead, I have made it such that similar 
 * to the bidir module, we can input a "training" to fine-tune these parameters (to an extent... some are still hard-coded and I'm not sure
 * where). That said, it seems to work "ok" and we can still make the config file an option for fine-tuning maybe.
 */
    printf("\t-mink\tMinimum number of finite mixtures to consider. default=1\n");
    printf("\t-maxk\tMaximum number of finite mixtures to consider. default=3\n");
    //printf("\t-rounds\tNumber of random seeds to use in the model. default=5\n");
    printf("\t-ct\tConvergence threshold after which processing stops. default=0.0001\n");
    //printf("\t-mi\tMaximum number of model iterations after which processing stops. default=2000\n");
// This is not currently true. It looks like it tries an then defaults back to 2000
    //printf("The model module currently has experimental support for parameter inference. Ie. it can attempt to estimate lambda,\n");
    //printf("sigma, pi, and w via conjugate priors. These values are specified using the following parameters:\n");
    //printf("\t-ALPHA_0\tPrior parameter 1 for sigma. Recommended value=1\n");
    //printf("\t-BETA_0\tPrior parameter 2 for sigma. Recommended value=1\n");
    //printf("\t-ALPHA_1\tPrior parameter 1 for lambda.\n");
    //printf("\t-BETA_1\tPrior parameter 2 for lambda.\n");
    //printf("\t-ALPHA_2\tSymmetric prior on mixing weights. Higher values=stronger attempt to find components of equal mixing weights.\n");
    //printf("\t\t\tRecommended value=100\n");
    //printf("\t-ALPHA_3\tSymmetric prior on the strand bias. Higher values=stronger attempt to find bidirectional events with equal strand bias.\n");
    //printf("\t\t\tRecommended value=100\n");
    
    //printf("-elon     : (boolean integer) adjust support of elongation component, (default=0)\n");
	//printf("              useful only when fitting to FStitch[1] or groHMM[2] output intervals\n");
}

/** Sets all values to their defaults as per the old read_in_parameters codebase.
 */
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
    this->maxk=3;
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
    this->debug=false;
    this->cores=1;
    this->version=false;
}

/** Parses the arguments passed to tfit and attempts to store them internally.
 * @param argc The number of command line arguments (first parameter to main())
 * @param argv The command line arguments array (second parameter to main())
 */
ParamWrapper::ParamWrapper(int argc, char **argv)
{
    int i;
    std::map<std::string, std::string> paramMap;
    char *prevCmd;
    std::map<std::string, std::string>::iterator it;
    
    //NOTE: Parameters were changed to the defaults seen in read_in_parameters.cpp, rather than the documented defaults.
    this->debug=false;
    this->exit=false;
    this->verbose=true;//false;
    this->module="";
    //This constructor sets blank or default values.
    this->forwardStrand="";
    this->reverseStrand="";
    this->mergedStrand="";
    this->jobName="EMG"; //"";
    this->outputDir="";
    this->logDir="";
    this->regionsOfInterest="";
    this->promoterTSS="";
    this->chromosome="all";//"";
    this->llrthresh=0.95;//1;
    //These parameters were modified given the behavior of stock Tfit:
    this->lambda=2000;
    this->sigma=123;
    this->pi=0.5;
    this->w=0.9;
    this->mink=1;
    this->maxk=3;//1;
    this->rounds=10;//5;
    this->ct=0.0001;
    this->mi=2000;
    this->experimentalValsSpecified=false;
    this->alpha0=1;
    this->beta0=1;
    this->ns=100;
    //There are no default values for this parameter.
    this->alpha0=1;
    this->beta0=1;
    this->alpha1=1;
    this->beta1=1;
    this->alpha2=1;//100;
    this->alpha3=1;//100;
    this->br=25;
    this->pad=500; // 2000 was default -- seemed to be too large for most regions provided and results in a lot of type I and II error
    this->footPrint=100;
    this->fdr=0; //This parameter appears to change how this software computes prior distributions in the bidir module. 
                 //When set to 0, Tfit will use a shortcut and/or precomputed model. 
    this->scores=""; //This is another undocumented parameter.
    this->r_mu=0; //yet another undocumented parameter.
    this->penalty=1; //There's documentation, but only in read_in_parameters.
    this->maxNoise=0.05; //This seems to only be used in across_segments.
    this->mle=0; //This parameter runs the model module after bidir, IIRC.
    this->elon=0;
    this->cores=1;
    this->version=false;
    
    
    if(argc==1)
    {
        printf("Error: no arguments specified for module %s\n", argv[1]);
        this->printUsage();
        this->exit=true;
        return;
    }
    
    else if(argc==2)
    {
        this->printUsage();
        this->exit=true;
        return;
    }
    
    this->module=std::string(argv[1]);
    
    if(this->module!="prelim" && this->module!="select" && this->module!="model" && this->module!="bidir_old")
    {
        printf("Invalid module specification: %s\n", argv[1]);
        this->exit=true;
        return;
    }
    
    this->prelim=this->module=="prelim";
    this->bidirOld=this->module=="bidir_old";
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
            
            else if(!strcmp(prevCmd, "-debug"))
            {
                this->debug=true;
            }
            
            else if((!strcmp(prevCmd, "-v"))||(!strcmp(prevCmd, "--version")))
            {
                this->version=true;
                return;
            }

        }
        
        else
        {
            paramMap.insert(std::pair<std::string, std::string>(prevCmd, argv[i]));
            prevCmd=NULL;
        }
    }
    
    for(it=paramMap.begin();it!=paramMap.end();it++)
    {
        if(it->first=="-f" || it->first=="--forward" || it->first=="-pos")
        {
            this->forwardStrand=it->second;
        }
        
        else if(it->first=="-r" || it->first=="--reverse" || it->first=="-neg")
        {
            this->reverseStrand=it->second;
        }
        
        else if(it->first=="-bg" || it->first=="--bedgraph" || it->first=="--combined")
        {
            this->mergedStrand=it->second;
        }
        
        else if(it->first=="-N" || it->first=="--jobname")
        {
            this->jobName=it->second;
        }
        
        else if(it->first=="-ns")
        {
            this->ns=atoi(it->second.c_str());
        }
        
        else if(it->first=="-o" || it->first=="--output")
        {
// For sake of user friendliness, we're just going to run the bidir and model moduels separately, so we don't need this anymore
            //if(it->second.back()!='/')
            //{
                //if(this->verbose)
                //{
                //    printf("Appending slash to end of output path.\n");
                //}
                
                this->outputDir=it->second;
            //}
            //
            //else
            //{
            //    this->outputDir=it->second;
            //}
        }
        
        else if(it->first=="-fdr" || it->first=="-FDR")
        {
            this->fdr=atoi(it->second.c_str());
        }
        
        else if(it->first=="--logOut" || it->first=="-l")
        {
            if(it->second.back()!='/')
            {
                if(this->verbose)
                {
                    printf("Appending slash to end of log output path.\n");
                }
                
                this->logDir=it->second+"\n";
            }
            
            else
            {
                this->logDir=it->second;
            }
        }
        
        else if(it->first=="-ms_pen")
        {
            this->penalty=atof(it->second.c_str());
        }
        
        else if(it->first=="-bd" || it->first=="--bidirs")
        {
            this->promoterTSS=it->second;
        }
        
        else if(it->first=="-chr" || it->first=="--chromosome")
        {
            this->chromosome=it->second;
        }
        
        else if(it->first=="-bct")
        {
            this->llrthresh=atoi(it->second.c_str());
        }
        
        else if(it->first=="-max_noise")
        {
            this->maxNoise=atof(it->second.c_str());
        }
        
        else if(it->first=="-lambda")
        {
            this->lambda=atof(it->second.c_str());
        }
        
        else if(it->first=="-scores")
        {
            this->scores=it->second;
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
            printf("Reading footprint parameter.\n");
            this->footPrint=atof(it->second.c_str());
            //Added for debugging purposes.
            printf("this->footprint is %f\n", atof(it->second.c_str()));
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
        
        else if(it->first=="-mle")
        {
            this->mle=atoi(it->second.c_str());
        }
        
        else if(it->first=="-elon")
        {
            this->elon=atoi(it->second.c_str());
        }
        
        else if(it->first=="-s" || it->first=="--segment")
        {
            printf("Setting regions of interest to %s\n", it->second.c_str());
            this->regionsOfInterest=it->second;
        }
        
        //yet another undocumented parameter!
        else if(it->first=="-r_mu")
        {
            this->r_mu=atoi(it->second.c_str());
        }
        
        else if(it->first=="--threads" || it->first=="-n")
        {
            this->cores=atoi(it->second.c_str());
        }
    }
    
    if((this->forwardStrand=="" || this->reverseStrand=="") && this->mergedStrand=="")
    {
        printf("Either arguments bgr or bgf empty or specified files do not exist. Please specify both forward AND revese strand if the bedGraph has not been concatenated.\n");
        this->exit=true;
    }
    
    /* This has now been taken care of earlier:
    if(this->outputDir!="")
    {
        if(this->outputDir[this->outputDir.length()-2]!='/')
        {
            this->outputDir=this->outputDir+"/";
        }
    }*/
    
    //Todo: add more sanity checks.
}

//Functions implemented to reach feature parity with read_in_parameters:

/** Displays a summary of relevant command line arguments.
 * 
 * @param nodes Number of nodes on which Tfit is to run. This is obtained from the MPI runtime.
 * @param cores Number of cores on which Tfit is to run. This is obtained from the MPI runtime.
 */
void ParamWrapper::display(int nodes, int cores){
    //We have no analogs, since these parameters aren't in the repo documentation:
    //bool MLE 	= stoi(p["-MLE"]);
    //bool SELECT = stoi(p["-select"]);
    std::string header 	= "";
    header+="----------------------------------------------------------------\n";
    header+="             Transcription Fit (Tfit)                   \n";
    if (prelim){
    }
    //if (this->module=="bidir" and MLE){
    //header+="             ....coupled to mixture model....     \n";
    //}
    if (prelim and select){
    header+="       ....coupled to BIC penalty optimization....     \n";
    }
    if (model){
    header+="               ...running mixture model...                      \n";
    }
    if (select){
    header+="            ....BIC penalty optimization....                      \n";		
    }
    printf("%s\n",header.c_str() );
    printf("Job Name          : %s\n", this->jobName.c_str());
    if (this->mergedStrand!=""){
        printf("bedGraph          : %s\n", this->mergedStrand.c_str());
    }else{
        printf("Forward bedGraph : %s\n", this->forwardStrand.c_str());
        printf("Reverse bedGraph : %s\n", this->reverseStrand.c_str());
    }
    //if (model or MLE){
    //printf("-k         : %s\n", p["-k"].c_str()  );
    //}
    if (this->promoterTSS!=""){
        printf("Bidirectionals    : %s\n", this->promoterTSS.c_str());
    }
    printf("Output Dir        : %s\n", this->outputDir.c_str()  );
    printf("Log Out Dir       : %s\n", this->logDir.c_str()  );
    printf("-MLE       : %d\n", this->mle);
    printf("-chr       : %s\n", this->chromosome.c_str());	
    printf("-br        : %d\n", this->br);	
    printf("-rounds    : %d\n", this->rounds);
    if (this->elon){
        printf("-elon      : %d\n", this->elon);
    }
    printf("-pad       : %lf\n", this->pad);
    if (!model){
    printf("-bct       : %d\n", this->llrthresh);
    }
    if (model){
        printf("-minK      : %d\n", this->mink);
        printf("-maxK      : %d\n", this->maxk);
    }
    printf("--threads  : %d\n",  cores);
    printf("--MPI_np    : %d\n",  nodes);
    printf("\nQuestions/Bugs? https://github.com/Dowell-Lab/Tfit" );
    printf("\nRevisions made by michael[dot]gohde[at]colorado[dot]edu\n");
    printf("----------------------------------------------------------------\n" );
    
}

/** Returns the date and time as a C++ string.
 * 
 * @return Date and time as a string.
 */
const std::string cdt() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m/%d/%Y %X", &tstruct);

    return buf;
}

/** Returns a summary of command line arguments.
 * 
 * @param ID likely specifies either execution mode (such as the running module) or is obtained from MPI.
 * @return Pretty-printed string representing command line arguments.
 */
std::string ParamWrapper::getHeader(int ID){
    std::ostringstream os;
    os<<"#----------------------------------------------------\n";
    os<<"#Date Time      : "<<cdt()<<"\n";
    os<<"#-N --jobName   : "<<this->jobName<<"\n";
    if (this->mergedStrand==""){
        os<<"#-f --forward   : "<<this->forwardStrand<<"\n";
        os<<"#-r --reverse   : "<<this->reverseStrand<<"\n";
    }else{
        os<<"#-bg --bedgraph : "<<this->mergedStrand<<"\n";
    }
    if (ID!=1){
        os<<"#-s --segment   : "<<this->regionsOfInterest<<"\n";
    }
    os<<"#-o --output      : "<<this->outputDir<<"\n";
    os<<"#-l --logOut      : "<<this->logDir<<"\n";    
    os<<"#-br            : "<<this->br<<"\n";
    if (ID==1){
        os<<"#-bct           : "<<this->llrthresh<<"\n";
        os<<"#-pad           : "<<this->pad<<"\n";
    }
    if (ID!=1){
        os<<"#-elon          : "<<this->elon<<"\n";
        os<<"#-minK          : "<<this->mink<<"\n";
        os<<"#-maxK          : "<<this->maxk<<"\n";
        os<<"#-mi            : "<<this->mi<<"\n";
        os<<"#-ct            : "<<this->ct<<"\n";
        os<<"#-rounds        : "<<this->rounds<<"\n";
    }
    if (ID!=1){
        os<<"#-ALPHA_0     : "<<this->alpha0<<"\n";
        os<<"#-BETA_0      : "<<this->beta0<<"\n";
        os<<"#-BETA_1      : "<<this->beta1<<"\n";
        os<<"#-ALPHA_2     : "<<this->alpha2<<"\n";
        os<<"#-ALPHA_3     : "<<this->alpha3<<"\n";
    }
    if (ID==1){
        os<<"#-sigma       : "<< this->sigma <<"\n";
        os<<"#-lambda      : "<< this->lambda<<"\n";
        os<<"#-foot_print  : "<< this->footPrint<<"\n";
        os<<"#-pi          : "<< this->pi<<"\n";
        os<<"#-w           : "<< this->w<<"\n";
    }
    os<<"#----------------------------------------------------\n";
    return os.str();
}
