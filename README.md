# Cavefish2017
Scripts and data from "The Importance of Selection in the Evolution of Blindness in Cavefish"

## Cavefish.nb

Mathematica notebook verifying the analysis of our deterministic model.

## Figures

This directory contains `R` scripts used to analyze the model in the paper and create descriptive figures.

## Simulations

This directory contains `R` and `Python` scripts for finite-population simulations. The scripts contain comments describing their function and usage.

## Data

This directory contains summarized results of the finite-population simulations. 

### KAlleles

This directory contains summarized results of the k-allele simulations after 10,000 generations.
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

### n1000

This directory contains summarized results of the single-allele drift simulations with n = 1000. Allele frequencies were recorded after 2500, 5000, 7500, and 10000 generations (when simulations were run for a total of 10000 generations), or after 1,250,000, 2,500,000, 3,500,000, and 5,000,000 generations (when simulations were run for a total of 5,000,000 generations). Each file is named with the number of generations. Files labeled with mk are results from simulations where the connection between the surface and cave populations was intermittant, with the probability of switching between connected and disconnected is 0.1 and is a markov process. The statistics are reported as described above.

The following code was used to generate the summary files for fully connect population from the simulation output. Other outputs were similarly generated.

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

### n100

This directory contains summarized results of the single-allele drift simulations as described above, but with n = 100.
