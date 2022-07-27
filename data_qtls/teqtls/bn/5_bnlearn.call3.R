#!/usr/bin/env Rscript 

library(bnlearn)
library(GenABEL)
#read input file
args = commandArgs(trailingOnly=TRUE)

input_dir=args[1]
output=args[2]

getBest <- function(posteriors){
  max_model <- max(posteriors)
  if(posteriors[1] == max_model){return("m1")}
  if(posteriors[2] == max_model){return("m2")}
  if(posteriors[3] == max_model){return("m3")}
}

getBNs <- function(input){

  d = as.data.frame(data.table::fread(input, stringsAsFactor=FALSE, head=TRUE))
  
  d <- d[,-1]
  
  col <- names(d)
  triplet <- paste(col,collapse=";")
  
  colnames(d)=c("V","G","T")
  #d$V[which(d$V == "./.")] <- NA 
  d <- d[!is.na(d$V),]
  d$V <- as.numeric(as.character(d$V))
  d$G <- as.numeric(as.character(d$G))
  d$"T" <- as.numeric(as.character(d$"T")) 
  d$V = (d$V - mean(d$V))/sd(d$V)
  
  #In statistics, @@standardization@@ is the process of putting different variables on the same scale. This process allows you to compare scores between different types of variables. Typically, to standardize variables, you calculate the mean and standard deviation for a variable. Then, for each observed value of the variable, you subtract the mean and divide by the standard deviation.
  
  # V = variant
  # T = Transposable ELEMENT
  # G = Gene
  
  #list of 4 models to test
  m=rep("", 3)
  m[1]="[V][T|V][G|T]"	#V -> T -> G # Scenario where the variants affects the TE that affects the G
  m[2]="[V][G|V][T|G]"	#V -> G -> T # Scenario where the variant affects the Gene and the Gene affects the TE
  m[3]="[V][T|V][G|V]"	#V -> G, V-> T # Scenario where the variants affects independently the TE and Gene.
  
  #calcuate P(D|G) where G=1,2,3
  loglik=rep(0, 3)
  for (i in 1:3) {
    net = model2network(m[i])
    loglik[i]=score(net, d, type="bge")
  }
  
  #assume constant prior over the 4 network configurations, i.e. P(G) = 0.25, and compute posterior
  prior=rep(1/3, 3)
  posteriors = exp(loglik - max(loglik)) * prior # scores are in log scale. substraction means division!!!
  posteriors = posteriors / sum(posteriors) # Division by the sum so that all of the probabilities add to 1!!
  
  #write output
  #cat ("M1=", signif(posteriors[1],3), " M2=", signif(posteriors[2],3), "M3=", signif(posteriors[3],3), "\n")
  best <- getBest(posteriors)
  
  df <- data.frame("variant" = col[1],
                   "gene" = col[2],
                   "te" = col[3],
                   "triplet" = triplet,
                   "L1" = loglik[1],
                   "L2" = loglik[2],
                   "L3" = loglik[3],
                   "M1" = signif(posteriors[1],3),
                   "M2" = signif(posteriors[2],3),
                   "M3" = signif(posteriors[3],3),
                   "max" = max(posteriors),
                   "best" = getBest(posteriors)
                   )
  
  return (df)
}


#input_dir <- "/scratch/lykoskou/GTExTE/data_qtls/teqtls/triplets/Adipose_Subcutaneous"

FILES <- Sys.glob(paste0(input_dir,"/*_triplet_content.txt"))
#print(FILES)
dataFrame <- data.frame()
for (rrr in 1:length(FILES))
{
  if (rrr %% 100 == 0) { cat ("  * Processed",rrr,"files\n")}
  
  input <- FILES[rrr]
  
  dataFrame <- rbind(dataFrame, getBNs(input))  
  
}

data.table::fwrite(dataFrame, output, quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")


