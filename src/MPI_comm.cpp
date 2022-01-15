/**
 * @file MPI_comm.cpp
 * @author Joey Azofeifa
 * @brief 
 * @version 0.1
 * @date 2016-05-20
 * 
 */
#include "MPI_comm.h"

#include <dirent.h>
#include <stddef.h>
#include <sys/types.h>

#include <vector>

#include <mpi.h>

#include "across_segments.h"
#include "load.h"
#include "read_in_parameters.h"
#include "split.h"

using namespace std;

vector<segment *> MPI_comm::slice_segments(vector<segment *> segments, int rank, int nprocs){
  
  int count 	= segments.size()/ nprocs;
  if (count==0){
    count 	= 1;
  }
  int start 	= rank * count;
  int stop 	= min(start + count, int(segments.size()));
  
  if (rank==(nprocs-1)){
    stop 	= int(segments.size());
  }
  if (start >= stop){
    start 	= stop;
  }
  vector<segment *> new_segs(segments.begin() + start, segments.begin()+stop);
  
  return new_segs;	
}

struct bounds{
public:
  double lower_upper[5];
};
int MPI_comm::gather_all_bidir_predicitions(vector<segment *> all, 
					    vector<segment *> segments , 
					    int rank, int nprocs, string out_file_dir, string job_name, int job_ID, params * P, int noise){
  
  map<string , vector<vector<double> > > G;
  map<string , vector<vector<double> > > A;
  //insert data from root
  int N 	= all.size();
  int S;
  int count 	= N/ nprocs;
  
  if (count==0){
    count 	= 1;
  }
  bounds B;
  MPI_Datatype mystruct;
  
  int blocklens[1]={5};
  MPI_Datatype old_types[1] = {MPI_DOUBLE}; 
  MPI_Aint displacements[1];
  displacements[0] 	= offsetof(bounds, lower_upper);
  MPI_Type_create_struct( 1, blocklens, displacements, old_types, &mystruct );
  
  MPI_Type_commit( &mystruct );
  int NN = 0;
  
  if (rank == 0){
    
    for (int i = 0 ; i < segments.size(); i++){
      for (int j = 0; j < segments[i]->bidirectional_bounds.size();j++){
	NN+=1;
       
	G[segments[i]->chrom].push_back(segments[i]->bidirectional_bounds[j]);
      }
    }
    
    
    for (int j =1; j < nprocs; j++){
      
      int start 	= j * count;
      int stop 	= min(start + count, int(N ));	
      
      if (j==nprocs-1){
	stop 	= N;
      }
      if (start >= stop){
	start 	= stop;
      }
      
      
      
      for (int i = 0; i < (stop-start) ; i++){
	MPI_Recv(&S, 1, MPI_INT, j, i, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	NN+=S;	
       
	G[all[ start+i ]->chrom] 	= vector<vector<double> >(S);
	for (int u = 0; u < S; u++ ){
	  MPI_Recv(&B, 1, mystruct, j, u, MPI_COMM_WORLD,MPI_STATUS_IGNORE);				  
	  for (int l = 0 ; l < 5;l++){
	    G[all[ start+i ]->chrom][u].push_back(B.lower_upper[l]);
	  }
	}
      }
      
      
    }
    
    
  }else  {
   
    int IS = segments.size();
    for (int i = 0;  i < segments.size();i++ ){
      int S 	= segments[i]->bidirectional_bounds.size();
      
      MPI_Ssend(&S, 1, MPI_INT, 0,i,MPI_COMM_WORLD );
      
      for (int u=0; u < segments[i]->bidirectional_bounds.size(); u++){				
	bounds B;
	
	for (int l = 0 ; l < 5; l++){
	  B.lower_upper[l] = segments[i]->bidirectional_bounds[u][l];
	}
	
	
	MPI_Send(&B, 1, mystruct, 0, u, MPI_COMM_WORLD);
	
      }
      
    }
    
    
  }
  cout<<"-------------------------"<<endl;
  if (rank==0 and not out_file_dir.empty()){
    load::write_out_bidirs(G, out_file_dir, job_name, job_ID, P, noise);
  }
  
  return NN;
}



struct simple_seg_struct{
	char chrom[6];
	char strand[2];
	int st_sp[4]; //first->start, second->stop
};


/* Function: MPI_comm::send_out_single_fit_assignments
 *
 * Purpose: 
 *
 * Args:
 *   FSI  a list of intervals/regions
 *   rank   MPI process identifier
 *   nprocs MPI total processes available
 *
 * Assumptions:
 *
 * Returns: a mapping of X to segments
 */
map<string, vector<segment *> > MPI_comm::send_out_single_fit_assignments(vector<segment *> FSI, int rank, int nprocs ){
  bool debug = false;
	map<string, vector<segment *> > GG;
	int N 		= FSI.size();
	int count 	= N / nprocs;
	if (count == 0){ count 	= 1; }
	simple_seg_struct sss;
	MPI_Datatype mystruct;
	
	int blocklens[3]={6,2,4};
	MPI_Datatype old_types[3] = {MPI_CHAR,MPI_CHAR, MPI_INT}; 
	MPI_Aint displacements[3];
	displacements[0] 	= offsetof(simple_seg_struct, chrom);
	displacements[1] 	= offsetof(simple_seg_struct, strand);
	displacements[2] 	= offsetof(simple_seg_struct, st_sp);
	
	MPI_Type_create_struct( 3, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );

	vector<simple_seg_struct> runs;
	int start, stop;
	int S;
	if (rank == 0){
		//first send out the number you are going to send
		for (int j =0; j < nprocs; j++){
			start 	= j*count;
			stop 	= (j+1)*count;
			stop 		= min(stop, N);
			if (j == nprocs-1){
				stop 	= N;
			}
			if (start > stop){
				stop 	= start; //in this case there are more nodes than segments...so
			}
			S 		= stop-start;
			if (j > 0){
				MPI_Send(&S, 1, MPI_INT, j,1, MPI_COMM_WORLD);
			}
			//now send out structs
			int u 	= 0;
			for (int i =start; i< stop; i++){
				simple_seg_struct SSS;
				for (int c = 0; c < 6; c++){
					if (c <FSI[i]->chrom.size() ){
						SSS.chrom[c] 	= FSI[i]->chrom[c];
					}else{
						SSS.chrom[c] 	= '\0';
					}
				}
				SSS.chrom[5]	= '\0';
				SSS.strand[0] 	= FSI[i]->strand[0];
				SSS.strand[1] 	= '\0';
				SSS.st_sp[0] 	= FSI[i]->start;
				SSS.st_sp[1] 	= FSI[i]->stop;
				SSS.st_sp[2] 	= FSI[i]->ID;
				SSS.st_sp[3] 	= FSI[i]->counts;
				if (j >0){
					MPI_Send(&SSS, 2, mystruct, j, u, MPI_COMM_WORLD  );
				}else{
					runs.push_back(SSS);
				}
				u++;
			}
		}

	}else{
		MPI_Recv(&S, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		for (int u = 0; u < S; u++){
			MPI_Recv(&sss, 2, mystruct, 0,u,MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			runs.push_back(sss);
		}	
	}
	//now convert to GG type
	for (int i = 0; i < runs.size(); i++){
		segment * ns 	= new segment(runs[i].chrom, runs[i].st_sp[0], runs[i].st_sp[1], runs[i].st_sp[2], runs[i].strand);
		ns->counts 		= runs[i].st_sp[3];
		GG[ns->chrom].push_back(ns) ;
	}

	typedef map<string, vector<segment *> >::iterator it_type;

	return GG;
}


vector<double> MPI_comm::send_out_parameters(vector<double> parameters, int rank, int nprocs){
	vector<double> new_parameters;
	double * P 	= new double[5];
	if (rank==0){
		for (int i = 0 ; i < 5;i++){
			P[i] 	= parameters[i];
		}
		for (int j = 1 ; j < nprocs;j++){
			MPI_Ssend(&P[0], 5, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
		}

	}else{
		MPI_Recv(&P[0], 5, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
	}
	for (int k = 0 ; k < 5; k++){
		new_parameters.push_back(P[k]);
	}

	return new_parameters;

}



int MPI_comm::get_job_ID(string path, string job_ID, int rank, int nprocs){
	//string OUT = dir+"EMGU-" + to_string(job_ID) +"_" + DT+ ".log";
	//tmp_EMGU-0_2.log
	int current 			= 1;
	string current_file 	= "";
	vector<string> line_array;
	int jN 	= job_ID.size();
	bool FOUND = false;
	if (rank==0){
		DIR * dirFile = opendir( path.c_str() );
		if ( dirFile ){
			struct dirent* hFile;
			errno = 0;
			while (( hFile = readdir( dirFile )) != NULL ){
				current_file 	= hFile->d_name;
				if (current_file.substr(0,jN) == job_ID){
					if (current_file.substr(jN, 1) =="-"){
						int t 		= 1;
						while (current_file.substr(jN+t,1) != "_" and t < current_file.size() ){
							t++;
						}
						if (t < current_file.size()){
							current 	= max(current, stoi( current_file.substr(jN+1,t) ) );	
						}else{
							printf("WHOOOO, MPI_comm getting job_ID\n");

						}
						FOUND 		= true;
					}
				}
				
				if (current_file.substr(0, jN + 4 ) =="tmp_" +job_ID){

					if (current_file.substr(jN+4,1) =="-"){
						int t 		= 1;
						while (current_file.substr(jN+4+t,1)   != "_" and t < current_file.size() ){
							t++;
						}
						if (t < current_file.size()){
							current 	= max(current, stoi( current_file.substr(jN+5,t) ) );
						}else{
							printf("WHOOOO, MPI_comm getting job_ID\n");
						}
						FOUND 		= true;
					}
					
				}

			}	
			closedir( dirFile );
		}else{
			cout<<"couldn't open temp log directory"<<endl;
		}
		if (FOUND){
			current++;
		}
		for (int j =1; j < nprocs; j++ ){
			MPI_Send(&current, 1, MPI_INT, j,1, MPI_COMM_WORLD);		
		}


	}else{
		MPI_Recv(&current, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
	}
	return current;
}

	
map<int, map<int, vector<simple_c_free_mode>  > > MPI_comm::gather_all_simple_c_free_mode(vector<map<int, vector<simple_c_free_mode> >> FITS, 
	int rank, int nprocs){
	

	simple_c_free_mode sc_fm;
	MPI_Datatype mystruct;
	
	int blocklens[4]={3,5, 6, 12};
	MPI_Datatype old_types[4] = {MPI_DOUBLE, MPI_INT, MPI_CHAR, MPI_DOUBLE}; 
	MPI_Aint displacements[4];
	displacements[0] 	= offsetof(simple_c_free_mode, SS);
	displacements[1] 	= offsetof(simple_c_free_mode, ID);
	displacements[2] 	= offsetof(simple_c_free_mode, chrom);
	displacements[3] 	= offsetof(simple_c_free_mode, ps);
	
	
	MPI_Type_create_struct( 4, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );


	vector<simple_c_free_mode> recieved;

	typedef map<int, vector<simple_c_free_mode> >::iterator model_it;
	if (rank==0){
		//first want to receive the number of segments each child processes ran
		for (int j = 0; j < nprocs; j++){
			int S =0;
			if (j > 0){
				MPI_Recv(&S, 1, MPI_INT, j, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for (int s = 0; s < S; s++){
					MPI_Recv(&sc_fm, 1, mystruct,j,s+1,  MPI_COMM_WORLD,MPI_STATUS_IGNORE );
					recieved.push_back(sc_fm);
				}
				
			}else{
				for (int s = 0; s < FITS.size(); s++){
					for (model_it k = FITS[s].begin(); k!=FITS[s].end(); k++){
						for (int sc= 0; sc < k->second.size();sc++ ){
							recieved.push_back( k->second[sc] );
						}
					}
				}					
			}
		}

	}else{

		int S 	= 0;
		for (int s = 0; s < FITS.size(); s++){
			for (model_it k = FITS[s].begin(); k!=FITS[s].end(); k++){
				for (int sc= 0; sc < k->second.size();sc++ ){
					S++;
				}
			}
		}
		MPI_Ssend(&S, 1, MPI_INT, 0,0, MPI_COMM_WORLD);
		S=0;
		for (int s = 0; s < FITS.size(); s++){
			for (model_it k = FITS[s].begin(); k!=FITS[s].end(); k++){
				for (int sc= 0; sc < k->second.size();sc++ ){
					simple_c_free_mode sc_FM 	= k->second[sc];
					MPI_Ssend(&sc_FM, 1, mystruct, 0, S+1, MPI_COMM_WORLD);	
					S++;
				}
			}
		}
		
	}
	//printf("Rank: %d,%d\n",rank, recieved.size());
	map<int, map<int, vector<simple_c_free_mode>  > > G;
	typedef vector<simple_c_free_mode>::iterator it_type_fm;

	for (it_type_fm sc = recieved.begin(); sc!=recieved.end(); sc++){
		G[(*sc).ID[0]][(*sc).ID[3]].push_back(*sc);
	}

	return G;

}


void MPI_comm::wait_on_root(int rank, int nprocs){
	int S 	= 0;
	if (rank==0){
		for (int j = 1; j < nprocs; j++){
			MPI_Ssend(&S, 1, MPI_INT, j,0, MPI_COMM_WORLD);

		}
	}else{
		MPI_Recv(&S, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
	}
}


map<string, vector<segment *> > MPI_comm::convert_segment_vector(vector<segment *> FSI){
	map<string, vector<segment *>> GG;
	for (int s = 0 ; s < FSI.size(); s++){
		GG[FSI[s]->chrom].push_back(FSI[s]);
	}
	return GG;
}










