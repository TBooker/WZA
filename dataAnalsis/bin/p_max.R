## An R implementation of Mike's Pmax

pMaxTest <- function( p, setSize){
  numTests = length(p)
  sortedP = sort(p, decreasing = T)
  pThreshold = sortedP[1 + numTests - setSize] 
  sum( dbinom( setSize:numTests, numTests, pThreshold))
}

# Point estimate of p from a zero-truncated binomial distribution (Rider 1955)
truncBinomPointEst <- function(x_raw, f_raw, k){
  n = max(x_raw)
  x = x_raw[k+1:n]
  f = f_raw[k+1:n]
  T0 = sum(f)  
  T1 = sum(x*f)  
  T2 = sum(x*x*f)  
## This is just a point estimate, but bootstrapping could be used to obtain a CI 
  return ( (T2 - k*T1)/( (n-1)*T1 - (k - 1)*n*T0))
}

# Do a sequenctial p-value correction on the p-values for each set of tests
estimate_binomial_p <- function(M){
  FiveWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x, 5)), "holm" )
  FourWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x, 4)), "holm" )
  ThreeWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x,3)), "holm" )
  TwoWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x,2)), "holm" )
  OneWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x,1)), "holm" )


  FiveWayHits <- sum(FiveWay<0.05)
  FourWayHits <- sum(FourWay<0.05)
  ThreeWayHits <- sum(ThreeWay<0.05)
  TwoWayHits <- sum(TwoWay<0.05)
  OneWayHits <- sum(OneWay<0.05)

  genes = c(0,OneWayHits,TwoWayHits, ThreeWayHits, FourWayHits, FiveWayHits)

  x_raw = c(0,1,2,3,4,5)

  return ( truncBinomPointEst(x_raw,genes,1) )
}

par(mfrow = c(4,1))

### Here's a set of silly simulations of the neutral case:

## Every gene has an effect pulled from the normal distribution
## any convergence will be by chance and should get 
## corrected out - but remember that the iterative pValue thing is not figure out

neutral = c()
for (i in 1:50){
  
  pop1 <- 1-rank(c(rnorm(10000,0,1)))/10000
  pop2 <- 1-rank(c(rnorm(10000,0,1)))/10000
  pop3 <- 1-rank(c(rnorm(10000,0,1)))/10000
  pop4 <- 1-rank(c(rnorm(10000,0,1)))/10000
  pop5 <- 1-rank(c(rnorm(10000,0,1)))/10000
  
  Matrix <- as.matrix( cbind(pop1,pop2,pop3, pop4, pop5) )
  neutral[i] = estimate_binomial_p( Matrix )
  
}

hist(neutral, xlim = c(0,1), main = "Distribution of estimated p from neutral simulations")

######## MILD CONVERGENCE
## Now there are now 1000 genes in each species genome that 
## has a Z score that is pulled from a distribution that has a mean of 2

mild = c()
for (i in 1:50){
  
  pop1 <- 1-rank(c(rnorm(1000,2,1),rnorm(9000,0,1)))/10000
  pop2 <- 1-rank(c(rnorm(1000,2,1),rnorm(9000,0,1)))/10000
  pop3 <- 1-rank(c(rnorm(1000,2,1),rnorm(9000,0,1)))/10000
  pop4 <- 1-rank(c(rnorm(1000,2,1),rnorm(9000,0,1)))/10000
  pop5 <- 1-rank(c(rnorm(1000,2,1),rnorm(9000,0,1)))/10000
  
  Matrix <- as.matrix( cbind(pop1,pop2,pop3, pop4, pop5) )
  mild[i] = estimate_binomial_p( Matrix )
  
}

hist(mild, xlim = c(0,1), main = "Distribution of estimated p from mild convergence simulations")


######## MEDIUM CONVERGENCE
## As before, but there are 3,000 genes pulled

medium = c()
for (i in 1:50){
  
  pop1 <- 1-rank(c(rnorm(3000,2,1),rnorm(7000,0,1)))/10000
  pop2 <- 1-rank(c(rnorm(3000,2,1),rnorm(7000,0,1)))/10000
  pop3 <- 1-rank(c(rnorm(3000,2,1),rnorm(7000,0,1)))/10000
  pop4 <- 1-rank(c(rnorm(3000,2,1),rnorm(7000,0,1)))/10000
  pop5 <- 1-rank(c(rnorm(3000,2,1),rnorm(7000,0,1)))/10000
  
  Matrix <- as.matrix( cbind(pop1,pop2,pop3, pop4, pop5) )
  medium[i] = estimate_binomial_p( Matrix )
  
}

hist(medium, xlim = c(0,1), main = "Distribution of estimated p from medium convergence simulations")





######## NO CONVERGENCE
## 
inferno = c()
for (i in 1:1){
  
  pop1 <- sample(1-rank(c(rnorm(100,10,1),rnorm(9900,0,1)))/10000)
  pop2 <- sample(1-rank(c(rnorm(100,10,1),rnorm(9900,0,1)))/10000)
  pop3 <- sample(1-rank(c(rnorm(100,10,1),rnorm(9900,0,1)))/10000)
  pop4 <- sample(1-rank(c(rnorm(100,10,1),rnorm(9900,0,1)))/10000)
  pop5 <- sample(1-rank(c(rnorm(100,10,1),rnorm(9900,0,1)))/10000)
  
  Matrix <- as.matrix( cbind(pop1,pop2,pop3, pop4, pop5) )
  inferno[i] = estimate_binomial_p( Matrix )
  
}

hist(inferno, xlim = c(0,1), main = "Distribution of estimated p from inferno convergence simulations")

