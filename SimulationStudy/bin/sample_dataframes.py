import pandas as pd
import sys, glob

all_data = []
for i in glob.glob(sys.argv[1]+"/*.csv"):
	rep = i.split("/")[-1].split("_")[0]
	temp = pd.read_csv(i) 
	temp = temp[temp.maf > 0.05]
	temp["rep"] = rep
	all_data.append(temp)

#pd.concat(all_data).sample(frac = 0.1, replace = False).to_csv(sys.argv[2],index = False)
pd.concat(all_data).to_csv(sys.argv[2],index = False)
