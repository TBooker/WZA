rm(list= ls())


first_set <- paste("gene",0:11,sep = "")
second_set <- paste("gene",49-(0:11),sep = "")

output_data <- list()
count = 0

for (i in 1:100){ 
  file <- paste("/media/booker/HOWDY/GEA/G.1.1/BC_map/",i,"_0.3_0.0001_100.csv", sep = '')
  if (file.exists(file)){
    count = count +1
    rep <- read.csv(file)
    rep$stat <- -log10(rep$pop_k_tau_p_value)
    
    start = rep[rep$gene%in% first_set,]
    labels_1 <- rep("1",nrow(start))
    end = rep[rep$gene%in% second_set,]
    labels_2 <- rep("2",nrow(end))
    output_data[[count]] <- list(start$stat,labels_1,end$stat,labels_2)
  }
  else{
    next
  }
}

## Do within chromosome test first
chromResultsMat <- matrix(nrow = length(output_data), ncol = 3)

for (i in 1:length(output_data)){
  chromResultsMat[i,1] = i
  slice <- output_data[[i]]
  dataSet <- c(slice[[1]], slice[[3]])
  labels <- c(slice[[2]], slice[[4]])
  fit <- lm(dataSet ~ labels)
  chromResultsMat[i,2] <- anova(fit)["Residuals","Mean Sq"] ## Within group variance
  chromResultsMat[i,3] <- anova(fit)["labels","Mean Sq"] ## Between group variance
}


chromResultsDF <- data.frame( chromResultsMat )
names(chromResultsDF) <- c("id","within","between")
chromResultsDF$ratio <- chromResultsDF$between/chromResultsDF$within
head(chromResultsDF)



## Do between chromosome test second
betweenResultsMat <- matrix(nrow = length(output_data), ncol = 3)
for (i in 1:length(output_data)){
  betweenResultsMat[i,1] = i
  slice_1 <- output_data[[sample(1:length(output_data), 1)]]
  slice_2 <- output_data[[sample(1:length(output_data), 1)]]
  dataSet <- c(slice_1[[1]], slice_2[[3]])
  labels <- c(slice_1[[2]], slice_2[[4]])
  fit <- lm(dataSet ~ labels)
  betweenResultsMat[i,2] <- anova(fit)["Residuals","Mean Sq"] ## Within group variance
  betweenResultsMat[i,3] <- anova(fit)["labels","Mean Sq"] ## Between group variance
}


betweenResultsDF <- data.frame( betweenResultsMat )
names(betweenResultsDF) <- c("id","within","between")
betweenResultsDF$ratio <- betweenResultsDF$between/betweenResultsDF$within

head(betweenResultsDF)
library(ggplot2)

ggplot(data = betweenResultsDF, aes(x = ratio))+
  geom_histogram(aes(fill = "between simulation"), alpha = 0.5)+
  geom_histogram(data = chromResultsDF, aes(x = ratio,fill = "within simulation"), alpha = 0.5)
  