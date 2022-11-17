# A script to implement the spline fitting business

import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

roller = 50
wza  = pd.read_csv(sys.argv[1])

print(wza.shape)
# remove null Z values - they won't help us

wza_t =  wza[ ~wza.Z.isnull() ].reset_index()


wza_s = wza_t.sort_values('SNPs')

#print(list(wza_s.Z))
rolled_Z_vars = wza_s.Z.rolling(window = roller).var()
masking_array = ~rolled_Z_vars.isnull()
rolled_Z_sd = np.sqrt(rolled_Z_vars)[masking_array]
rolled_Z_means = wza_s.Z.rolling(window = roller).mean()[masking_array]
rolled_mean_SNP_number = wza_s.SNPs.rolling(window = roller).mean()[masking_array]


# Generating weights for polynomial function with degree =2 - standard deviation
sd_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_sd, deg = 2)
# Generating model with the given weights - the model as a function
sd_polynomial_model = np.poly1d(sd_weights)

# Generating weights for polynomial function with degree =2 - mean
mean_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_means, deg = 2)
# Generating model with the given weights - the model as a function
mean_polynomial_model = np.poly1d(mean_weights)


xp = np.linspace(1,wza["SNPs"].max(),100)



sd_predictions = sd_polynomial_model( wza["SNPs"] )
mean_predictions = mean_polynomial_model( wza["SNPs"] )

print(wza["Z"][987])
print( mean_predictions[987])
print( sd_predictions[987])

wza_p_values = []
# I think I have to loop as I don't think stats.cdf allows vectorised input
for i in range(wza.shape[0]):
    p_val = 1- norm.cdf(wza["Z"][i], loc = mean_predictions[i], scale = sd_predictions[i])
    print(i,wza["SNPs"][i], wza["Z"][i],  mean_predictions[i], sd_predictions[i], p_val)
    wza_p_values.append(p_val)


wza["Z_pVal"] = wza_p_values

print( wza.sort_values('Z_pVal'))

plt.scatter(rolled_mean_SNP_number, rolled_Z_sd )
plt.plot(xp, sd_polynomial_model(xp), "-", c = "r")
plt.show()

#print( norm.cdf(wza["Z"][987], loc = mean_predictions[987], scale = sd_predictions[987]) )

wza.to_csv("TEMP",index = False)
sys.exit()
