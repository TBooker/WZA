import pandas as pd
import numpy as np

def genomeMakerWithNeutral(genome_size = 1000000, gene_number = 10):


	geneStarts = np.array( [ i for i in range(0,1000000,5000) ][1::4] )

	geneEnds = geneStarts + 9999

	return( pd.DataFrame([geneStarts, geneEnds]).transpose() )

genome = genomeMakerWithNeutral()

windows = [[49500+ (i*100000),50499 + (i*100000)] for i in range(10)]

print(genome)

print(windows)
