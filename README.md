# Cavefish2017
Scripts and data from "The Importance of Selection in the Evolution of Blindness in Cavefish"

## Figures

This directory contains `R` scripts used to analyze the model in the paper and create descriptive figures.

## Simulations

This directory contains `R` and `Python` scripts for finite-population simulations. The scripts contain comments describing their function and usage.

## Data

This directory contains summarized results of the finite-population simulations.

```R
library(reshape2)
a <- read.csv("drift_sims_log_connected_u6.csv")
qq = acast(a,m~s,value.var="q1",mean)
write.table(qq,"drift-2500.txt")
qq = acast(a,m~s,value.var="q2",mean)
write.table(qq,"drift-5000.txt")
qq = acast(a,m~s,value.var="q3",mean)
write.table(qq,"drift-7500.txt")
qq = acast(a,m~s,value.var="q",mean)
write.table(qq,"drift-10000.txt")
```

### KAlleles

This directory contains summarized results of the k-allele simulations after 10,000 generates.
The following statistics are reported:

 * q: average cave-adaptive allele frequency
 * a: average fitness in the population
 * f: average phenotype in the population
 * l: average genetic load in the population

Each file contains a table in which the each cell corresponds to a specific migration (row) x selection (column) parameter set. Row labels and column labels are in quotes.
Each cell is the average of 100 simulations.
The file names contain both the population size (n=1000), the number of loci (k = 1,2,4,6,12),
and the statistic reported.

The following code was used to generated the summary files from the simulation output.

```R
library(data.table)
library(reshape2)

# Read the simulation results, which are really, really big
d = fread("/tmp/k-allele-sim.csv")

# For each k, extract the results from generation 10000
# For each statistic, construct a migration x selection table.
# Save results
for(kk in c(1,2,4,6,12)) {
	a = subset(d, k == kk & g == 10000)
	for(v in c("q","a","f","l")) {
		qq = acast(a,m~s,value.var=v,mean)
		fn = sprintf("drift-n1000-k%02d-%s.txt",kk,v)
		write.table(qq,fn)
	}
}
```
