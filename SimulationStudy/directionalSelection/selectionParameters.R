#########
######### Directional Selection
#########

## 40% Fitness difference
s = 0.0055
nLoci = 12
(1+s*-7)^nLoci /(1+s*7)^nLoci 


## 40% Fitness difference
s = 0.003
nLoci = 12
(1+s*-7)^nLoci /(1+s*7)^nLoci 

## 25% Fitness difference
s = 0.0017
nLoci = 12
(1+s*-7)^nLoci /(1+s*7)^nLoci 

## 10% Fitness difference
s = 0.00065
nLoci = 12
(1+s*-7)^nLoci /(1+s*7)^nLoci 

## 5% Fitness difference
s = 0.0003
nLoci = 12
(1+s*-7)^nLoci /(1+s*7)^nLoci 

## 1% Fitness difference
s = 0.000065
nLoci = 12
(1+s*-7)^nLoci /(1+s*7)^nLoci 


#########
######### Stabilising Selection
#########

gaussianFitness <- function(phen, opt, Vs){
  exp( (-1*( phen - opt)^2)/(2*Vs)) 
}

## 40% Fitness difference
Vs = 192
gaussianFitness(-7,7,Vs)/ gaussianFitness(7,7,Vs)

## 10% Fitness difference
Vs = 900
gaussianFitness(-7,7,Vs)/ gaussianFitness(7,7,Vs)

