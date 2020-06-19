library(boot)
#DataFolder = "../../../../../../../Volumes/disk1/home/acevedo.ashley/Closed/Ebola/SetB/RoNi/v2_output_combined/"
#DataFolder = "../../../../Desktop/Pradeep/Q20ThresholdFiles/InVitroTranscription/Uncapped/"

#Load Q20 threshold count data
Data = read.table(paste(DataFolder,"Q20threshold.txt",sep=""), header=FALSE, col.names=c("POSITIONS", "REFERENCE", "A", "C", "G", "T"), colClasses=c("integer","character","integer","integer","integer","integer"))
attach(Data)

COVERAGE = A + C + G + T

#Inititate vectors to hold data in a more easily parseable format
Positions = c()
Variants = c()
VariantCounts = c()
Coverage = c()

#Vector/list to help identify possible variants and extract their counts
VARIANTS = c("A","C","G","T")
Counts = list(A, C, G, T)

#Convert count data into a more easily parseable fromat (with variant types and their counts explicitly specified)
for (i in 1:length(POSITIONS)){
	for (j in 1:length(VARIANTS)){
		if (VARIANTS[j] != REFERENCE[i]){
			Positions = c(Positions, POSITIONS[i])
			Variants = c(Variants, paste(REFERENCE[i],VARIANTS[j],sep=""))
			VariantCounts = c(VariantCounts, Counts[[j]][i])
			Coverage = c(Coverage, COVERAGE[i])
		}
	}
}

Frequencies = VariantCounts/Coverage

#Likelihood function for the probability of the observed data (muliplication of probabilities of observing the data at all nonsense positions) given a proposed rate
LikelihoodFunction = function(mu, counts, trials, LB, UB){
	if (mu > UB | mu < LB){
		return(0)
	} else {
		Probabilities = dbinom(counts,size=trials,prob=mu)
		LogLikelihood = sum(log(Probabilities))
		return(-LogLikelihood)
	}
}

OUTFILENAME = "MutationRates_Q20MLE.txt"

#Open and write header for the output file
cat(paste("Type","MLE","SE",sep="\t"), file=paste(DataFolder,OUTFILENAME,sep=""), sep="\n")

VariantTypes = c("AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG")

#Perform maximum likelihood estimation for each type of variant
for (TYPE in VariantTypes){

	InitialCounts = VariantCounts[(Variants==TYPE)&(Coverage>100000)]
	InitialCoverage = Coverage[(Variants==TYPE)&(Coverage>100000)]
	

	#Estimate the mutation rate from the data to use as a starting point for optimization of the likelihood function
	muEstimate = sum(InitialCounts)/sum(as.numeric(InitialCoverage))

	#Identify outliers using a binomial test - variants where the frequencies are significantly different from the mean are removed from further analysis
	repeat{
		RefinedCounts = c()
		RefinedCoverage = c()
		for (i in 1:length(InitialCounts)){
			if (binom.test(InitialCounts[i], InitialCoverage[i], muEstimate, alternative=c("greater"))$p.value > .05){
				RefinedCounts = c(RefinedCounts, InitialCounts[i])
				RefinedCoverage = c(RefinedCoverage, InitialCoverage[i])
			}
		}
		muEstimate = sum(RefinedCounts)/sum(as.numeric(RefinedCoverage))
		if (length(RefinedCounts) != length(InitialCounts)){
			InitialCounts = RefinedCounts
			InitialCoverage = RefinedCoverage
		} else {break}
	
	}
	
	#Bootstrapping function to provide maximum likelihood estimate and a bootstrapped standard error
	MLEBootstrapFunction = function(d, i, start){
		MaximumLikelihoodEstimate = optim(start, LikelihoodFunction, counts=d[i,]$counts , trials=d[i,]$trials , method="Brent", upper=start*10, lower=start/10, UB=start*10, LB=start/10)
		return(MaximumLikelihoodEstimate$par)
	}

	#Combine count data by column into a dataframe to allow sampling of variant and total counts as pairs
	SamplingMatrix = data.frame(counts=RefinedCounts, trials=RefinedCoverage)
	
	#Non-parametric bootstrap of count data in maximum likelihood function
	MLEBootstrap = boot(SamplingMatrix, statistic=MLEBootstrapFunction, R=1000, start = muEstimate)
	
	#t0 equals the population estimate
	#t is the vector of bootstrap estimates
	cat(paste(TYPE, MLEBootstrap$t0, sd(MLEBootstrap$t), sep="\t"), file=paste(DataFolder,OUTFILENAME,sep=""), sep="\n", append=TRUE)

}	

detach(Data)