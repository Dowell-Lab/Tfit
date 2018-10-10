##Author: 		Joey Azofeifa
#model_across is essentially the file handler for model.
#model.py is the heavy lifiting of "math"
#model_across takes model and parallelizes it
#
#
import multiprocessing as mp
import model
import os
import load
def checkFileExists(FILE, i):
	if os.path.exists(FILE + str(i)):
		return checkFileExists(FILE, i+1)
	return FILE + str(i)

def wrapper_fit_function(X, k=2, 
max_iterations=200, convergence_thresh=0.0001,
move_uniform=0):
	clf 	= model.EMGU(max_ct=convergence_thresh, max_it=max_iterations, K=k, bayes=False, noise=True, 
			noise_max=0.1, moveUniformSupport=0, cores=1, seed=True)
	clf.fit(X)
	return clf.ll , clf.rvs, clf.converged, clf.resets,clf

	
##Important stuff right here##
#https://www1.udel.edu/biology/rosewc/kaap686/notes/EMG%20analysis.pdf
#link might have some explination of EMG or be totally wrong

def run(D, bic, rounds, max_k, 
	standardize, convergence_thresh,
			max_iterations, move_uniform, 
			write_out_dir, specific_chromosome):
	#if BIC is 0 don't perform model selection and output each model 
	#from 1  to max_k individually
	FILE 	= write_out_dir+"EMG_model_fits_" + specific_chromosome + "_"
	FILE 	= checkFileExists(FILE, 0)
	FHW 	= open(FILE, "w")
	for d in D:
		d.X[:,0]-=min(d.X[:,0])
		d.X[:,0]/=standardize #devide everything by 10 
		if move_uniform == 0: #lets parrallelize the rest of this
			models 	= list()
			for k in range(max_k):
				output 	= mp.Queue()#multiprocessing function
				#defining a function within a for loop
				def wrapper_fit_function_pp(X, k, output,
					max_iterations=max_iterations, convergence_thresh=convergence_thresh,
					move_uniform=move_uniform):
					#(line 270) 
					clf 	= model.EMGU(max_ct=convergence_thresh, max_it=max_iterations, K=k, bayes=False, noise=True, 
							noise_max=0.1, moveUniformSupport=0, cores=4)
					clf.fit(X)
					output.put((clf.ll , clf.rvs, clf.converged, clf.resets, clf))
				#runs EMGU in mulitprocessing here 
				processes 		= [mp.Process(target=wrapper_fit_function_pp, args=(d.X, k,output) )for i in range(rounds)]
				for p in processes:
				    p.start()
				for p in processes:
					p.join()
				keepers 	= [output.get() for p in range(rounds)]
				models.append(min(keepers))
		else:

			models 	= [max([ wrapper_fit_function(d.X, k,
				max_iterations=max_iterations,
				convergence_thresh=convergence_thresh,
				move_uniform=move_uniform ) for r in range(rounds)]) for k in range(max_k+1)]
			

		FHW.write("#"+d.print_info())
		for ll,model, converged, resets,clf in models:
			FHW.write("~"+str(len(model)) + "," + str(ll) + "," + str(converged) + "," + str(resets)+ "\n")
			model_txt = "\n".join([ m.__str__() for m in model])
			FHW.write(model_txt+"\n")
		FHW.flush()

if __name__ == "__main__":
	
	X	=  load.grab_specific_region("chr1",836835, 843549,
		pos_file="/Users/jackdempsey/Desktop/Tfit_All/Tfit/examples/test_j_neg.BedGraph",
		neg_file="/Users/jackdempsey/Desktop/Tfit_All/Tfit/examples/test_j_neg.BedGraph",
		SHOW 	=False, bins=300)
		 	
	#stuff = wrapper_fit_function(X, k=2, 
	#	max_iterations=200, convergence_thresh=0.0001,
	#	move_uniform=0)
	
	#D 			= load.formatted_file(formatted_file, bins,specific_chromosome)
	#test = run(D , bic, rounds = 5 , max_k =2 , standardize = 10, convergence_thresh  = 0.0001, max_iterations = 5 , move_uniform = 0 , write_out_dir= "/Users/jackdempsey/Desktop/Tfit_All/Tfit/" , specific_chromosome= "CHR!")
	
	
	
	
	
