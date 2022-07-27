#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(GenABEL))


bed_file <- args[1]
cov_file <- args[2]
out_file <- args[3]

# Verbose
cat("* Residualizing expression data.
* Bed file: [",bed_file, "]\n* Covariates file: [", cov_file,"]\n* out file: [",out_file,"]\n")

#Reading files
cat("\t* Reading bed file.")
bed <- as.data.frame(data.table::fread(paste("zcat ",bed_file,sep=""),header=T,stringsAsFactors=F,sep="\t"))
cat("\t* Reading covariates")
cov <- as.data.frame(data.table::fread(paste("cat ",cov_file,sep=""),header=T,stringsAsFactors=F,sep="\t"))

# Verbose
cat("Number of phenotypes to residualize: ",nrow(bed),"\n")
cat("Numver of covariates to be used: ",nrow(cov),"\n")

# Processing covariates, removing first column containing pc info, ordering based on sample names in bed file and transposing.
rownames(cov) <- cov$SampleID 
cov <- cov[,-1]
cov <- cov[names(bed[,-c(1:6)])]
cov.t <- as.matrix(t(cov))
cov.t <- scale(cov.t, center=T,scale=F)

# Main loop for residualization
results <- foreach(rrr=1:nrow(bed), .combine = rbind, .multicombine = T) %do% {
    cat(rrr,"-",nrow(bed),"\r")

    x <- as.data.frame(bed[rrr,])
    Y <- as.numeric(x[,-c(1:6)])

    lm_res <- lm(formula = Y ~ cov.t)
    
    # Getting the residuals and adding the mean.
    exp <- residuals(lm_res) + summary(lm_res)$coefficients[1]
    df_exp <- data.frame("exp" = exp)
    rownames(df_exp) <- names(bed[,-c(1:6)])

    return(data.frame("#chr" = x$"#chr",
                      "start" = x$start, 
                      "end" = x$end,
                      "gene" = x$id,
                      "info" = x$info,
                      "strand" = x$strand,
                      t(df_exp)))
}
# If you have sample names starting with an integer, fucking R puts an X in front. The line below removes the Xs.
colnames(results) <- gsub("X","",colnames(results))

colnames(results) <- gsub(".chr","#chr",colnames(results))

#Write results in file.
cat("Writing results in file.\n")
write.table(results,paste0(out_file,sep=""),col.names=T,row.names=F,sep="\t",quote=F)

cat("Done :)\n")
