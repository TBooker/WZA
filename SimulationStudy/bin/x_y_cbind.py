import pandas as pd
import argparse

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--x",
			required = True,
			dest = "file_x",
			type = str,
			help = "The x file")

	parser.add_argument("--y",
			required = True,
			dest = "file_y",
			type = str,
			help = "The y file - i.e. the one that will be pasted to the x-file")

	parser.add_argument("--output",
			required = True,
			dest = "output",
			type = str,
			help = "What name do you want to give to the output file? [DON'T GIVE FILE EXTENSION, THE SCRIPT MAKES TWO OUTPUT FILES")

	parser.add_argument("--retain",
			required = False,
			dest = "retain",
			nargs="+",
			type = str,
			help = "Give the columns from y that you want to add to x")


	args = parser.parse_args()

	gea_data = pd.read_csv( args.file_x, engine = "python")

	bayPass_data = pd.read_csv( args.file_y ,sep = r"\s+", engine = "python")

	if args.retain == None:
		cols = ["BF(dB)"]

	for i in cols:
		gea_data[i] = bayPass_data[i]

	gea_data.to_csv(args.output,
					index = False)
main()
