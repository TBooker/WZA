## A short script to extract sites from the neutral chromosome from BayPass files
import argparse


def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--csv", 
			required = True,
			dest = "csv",
			type = str, 
			help = "The name of CSV file")
	parser.add_argument("--bayPass", 
			type = str,
			dest = "bayPass", 
			help = "The name of the BayPass file?",
			default = "none")			
	parser.add_argument("--output", 
			type = str,
			dest = "output", 
			help = "The name of the output file?",
			default = "none")			
	args = parser.parse_args()

	lineCount_1 = 0
	for i in open(args.csv):
		if i.startswith("position"): continue
		lineCount_1 += 1
		pos = float(i.strip().split(',')[0] )

		if pos >= 8000000:
			neutralChromLine = lineCount_1
			break

	output_file = open(args.output, "w")

	lineCount_2 = 0

	for i in open(args.bayPass):
		lineCount_2 += 1
		if lineCount_2 >= neutralChromLine:
			output_file.write(i)
	output_file.close()

main()
