import pandas as pd
import numpy as np

# Simulate a recombination rate map as a Markov Chain
def run_mc_sim( P, num_iters = 1000 ):

    num_states = P.shape[0]

    states = np.zeros( num_iters, dtype = np.int8 )

 # initialize variable for first state
    states[0] = np.random.randint(0,num_states-1)

    for t in range(1,num_iters):

    # probability vector to simulate next state X_{t+1}
        p  = P[int(states[t-1]), ]
        rng = np.random.default_rng()
        state_sample = np.min(np.where(  rng.multinomial(1, p, size=1)[0] ==1 ) )

    ## draw from multinomial and determine state
        states[t] =  int(state_sample)
    return states
transitionMatrix = np.array([[0.3, 0.7, 0.0],
                [0.15, 0.7, 0.15],
                [0.00, 0.7, 0.3]])


num_chains  = 1
num_iterations = 1000


wza_0 = pd.read_csv("/home/booker/work/GEA/simulations/SimulationStudy/s0.0136/recombination_map/BC_Map_output_i100.wza.csv")
wza_1 = pd.read_csv("/home/booker/work/GEA/simulations/SimulationStudy/s0.0136/recombination_map/BC_Map_output_i1000.wza.csv")
wza_2 = pd.read_csv("/home/booker/work/GEA/simulations/SimulationStudy/s0.0136/BC_Map.wza.csv")
print(list(wza_0))
ranger = [i for i in range(1000)]

def gimmeScores(wza_0, wza_1, wza_2, stat):
    scores = np.array([ wza_0[wza_0.rep == i][stat],
                            wza_1[wza_1.rep == i][stat],
                            wza_2[wza_2.rep == i][stat] ] )
    return(scores)


output_dfs = []
for i in range(1,30):
    recMap = run_mc_sim( transitionMatrix )

    if len( set( [wza_0[wza_0.rep == i].shape[0],
                    wza_1[wza_1.rep == i].shape[0] ,
                    wza_2[wza_2.rep == i].shape[0]]) ) != 1: continue

    WZA_scores = gimmeScores(wza_0, wza_1, wza_2, "empR_Z")[recMap,ranger]
    WZA_kendall = gimmeScores(wza_0, wza_1, wza_2, "kendall_Z")[recMap,ranger]
    WZA_RDA_scores = gimmeScores(wza_0, wza_1, wza_2, "RDA_Z")[recMap,ranger]
    RDA_scores = gimmeScores(wza_0, wza_1, wza_2, "RDA_pVal")[recMap,ranger]
    TC_scores = gimmeScores(wza_0, wza_1, wza_2, "top_candidate_p_TC2")[recMap,ranger]
    RDA_TC_scores = gimmeScores(wza_0, wza_1, wza_2, "RDA_top_candidate_p_TC2")[recMap,ranger]

    temp_df = pd.DataFrame( { "empR_Z":WZA_scores,
                    "kendall_Z":WZA_kendall,
                    "RDA_Z":WZA_RDA_scores,
                    "top_candidate_p_TC2":TC_scores,
                    "RDA_top_candidate_p_TC2":RDA_TC_scores,
                    "LA":wza_0[wza_0.rep == i]["LA"],
                    "kendall_rank":wza_1[wza_1.rep == i]["kendall_rank"],
                    "RDA_rank":wza_1[wza_1.rep == i]["RDA_rank"],
                    "gene":wza_1[wza_1.rep == i]["gene"],
                    "rep":i} )

    output_dfs.append( temp_df )


pd.concat( output_dfs ).to_csv("recombinationRateDatasets.csv")
