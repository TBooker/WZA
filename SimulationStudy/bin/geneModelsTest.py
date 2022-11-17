import numpy as np
## Model "Genes" by examining GEA within 10,000bp windows

geneStarts = np.array( [49500+ (i*100000) for i in range(10)] )

geneEnds = geneStarts + 999

for s, e in zip(geneStarts, geneEnds):
	print(s, e)
print("\n")

geneStarts = np.array( [ i for i in range(0,1000000,10000) ][1::2] ) - 500

geneEnds = geneStarts + 999

for s, e in zip(geneStarts, geneEnds):
	print(s, e)
