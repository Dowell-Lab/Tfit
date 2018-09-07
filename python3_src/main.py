#================================================================================================
#Author: 		Joey Azofeifa
#Creation Date: 05/20/2015
#Last Update:   06/30/2015
#Input Parameters
#(-i) annotation file: <strand>\t<chrom>\t<start>\t<stop>\t<info>\n
#(-j) forward strand genome coverage file: <chrom>\t<start>\t<stop>\t<coverage>\n
#(-k) reverse strand genome coverage file: <chrom>\t<start>\t<stop>\t<coverage>\n
#(-c) specific chromosome; default set to run against all
#(-d) density threshold; default set to 0.1, i.e. one coverage value per 10 basepairs
#(-b) max number of model components; default set to 1
#(-r) binning resolution; default set to 10
#(-s) scaling resolution; default set to 200
#(-p) pad; default set to 0
#
#Edited and translated to python3 Summer/Fall 2018 by Jack Dempsey
#
#If you wish to run the hardcoded function, you may just execute:$python3 main.py
#
#If you wish to have more specific then you may with the added arguements:
#:$python3 main.py runModel -i /direction/to/input/bedgraph -wo /write/out/directory -sc hg19
#
#
#================================================================================================

import sys, load, parse_argv, model_across
def hardcode():
	SI 	= True
	if SI: #single isoform genes that overlap a single FStitch call
		#===================
		#union of forward and reverse strand FStitch calls
		#BELOW is hardcoding
		#====================
		#FILES: 
		root 			= "/Users/jackdempsey/Desktop/Tfit_All/Tfit/"
		refseqFILE 		= root+"annotation_files/hg19_TSS.bed"
		FS_forward		= root+"examples/SRR1105737.pos.sorted.BedGraph"#need output from fstitch .bed files
		FS_reverse  	= root+"examples/SRR1105737.minus.sorted.BedGraph"#need output from fstitch
		forward_file_bg = root+"examples/SRR1105737.pos.sorted.BedGraph"
		reverse_file_bg = root+"examples/SRR1105737.minus.sorted.BedGraph"
		write_out 		= root+"examples/single_isoform_FStitch.tsv"
		#====================
		#====================
		#INPUT PARAMETERS
		bins 	= 300


		RF 		= load.gene_annotations(refseqFILE)
		FS 		= load.FStitch_annotations(FS_forward, FS_reverse, merge=True)
		#print("RF: ", RF)# "FS: ", FS)
		filtered= load.filter_single_overlaps(FS, RF)
		load.bedGraphFile(forward_file_bg, reverse_file_bg,
			filtered, SHOW=False, test=False,
			write_out=write_out)



def run(argv):

	aw 		= parse_argv.run(argv)
	if aw.exit:
		aw.printErrors()
		print("exiting...")
		return False
	if aw.mod == "formatData":
		if aw.mod2=="FStitchSingleIsoform":
			#====================================
			refseqFILE 		= aw.G["-ref"][0]
			FS_forward 		= aw.G["-ffs"][0]
			FS_reverse 		= aw.G["-rfs"][0]
			forward_file_bg = aw.G["-fbg"][0]
			reverse_file_bg = aw.G["-rbg"][0]
			write_out_dir 	= aw.G["-wo"][0]
			pad 			= float(aw.G["-pad"][0])
			#====================================
			RF 		= load.gene_annotations(refseqFILE)
			FS 		= load.FStitch_annotations(FS_forward, FS_reverse, merge=True)
			filtered= load.filter_single_overlaps(FS, RF)
			load.bedGraphFile(forward_file_bg, reverse_file_bg,
				filtered, SHOW=False, test=False,
				write_out=write_out_dir)

		elif aw.mod2=="RefSeqOnly":
			print("RefSeq only option not supported yet")

		elif aw.mod2=="FStitchMerged":
			FS_forward 		= aw.G["-ffs"][0]
			FS_reverse 		= aw.G["-rfs"][0]
			forward_file_bg = aw.G["-fbg"][0]
			reverse_file_bg = aw.G["-rbg"][0]
			write_out_dir 	= aw.G["-wo"][0]
			pad 			= float(aw.G["-pad"][0])
			FS 				= load.FStitch_annotations(FS_forward, FS_reverse,
				merge=True, pad=pad)
			load.bedGraphFile(forward_file_bg, reverse_file_bg,
				FS, SHOW=False, test=False,
				write_out=write_out_dir)



	else: #run MODEL!
		#====================================
		formatted_file 		= aw.G["-i"][0]
		write_out_dir 		= aw.G["-wo"][0]
		max_k 				= int(aw.G["-k"][0])
		rounds 				= int(aw.G["-it"][0])
		bins 				= int(aw.G["-b"][0])
		specific_chromosome = aw.G["-sc"][0]
		bic 				= float(aw.G["-bic"][0])
		standardize 		= float(aw.G["-st"][0])
		convergence_thresh 	= float(aw.G["-mc"][0])
		max_iterations 		= int(aw.G["-mt"][0])
		move_uniform 		= int(aw.G["-mu"][0])
		#====================================
		D 					= load.formatted_file(formatted_file, bins,specific_chromosome)
		model_across.run(D, bic, rounds,
			max_k, standardize, convergence_thresh,
			max_iterations,move_uniform, write_out_dir, specific_chromosome)







if __name__ == "__main__":

	#example input: runModel -i /Users/jackdempsey/Desktop/Tfit_All/Tfit/examples/SRR1105737.pos.sorted.BedGraph -wo /Users/jackdempsey/Desktop/Tfit_All/Tfit/examples/tests -sc hg19
	if(len(sys.argv)>1):
		run(sys.argv)
	else:
		print("No Arguments Were given, Running Hardcoded version \n")
		hardcode()

	pass
