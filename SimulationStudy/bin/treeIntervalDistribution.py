## A short script to extract the distribution of tree widths from a treeSeq file

import pyslim
import argparse
import pandas as pd

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--trees", "-t", 
			required = True,
			dest = "trees",
			type = str, 
			help = "The treeSeq file")
	parser.add_argument("--output", "-o", 
			required = True,

			dest = "output",
			type = str, 
			help = "The name of the output file (extension not added - do that yourself!)")

	args = parser.parse_args()
	
	ts = pyslim.load( args.trees )
	count = 0
	
	output = open( args.output, "w")

	output.write( "tree_no,start,end,dist\n")
	
	for i in ts.trees():
		count += 1

		outputLine = [ count, i.interval[0], i.interval[1], i.interval[1] - i.interval[0] + 1 ]	

		
		output.write(",".join( map( str, outputLine ) )+ "\n")

	output.close()
		
if __name__ == "__main__":
	main()