#!/bin/sh

# This script ties k-allele-sim.R and k-allele-sim-list.R together

RSCRIPT="Rscript --vanilla"

(echo ""; sleep 1; $RSCRIPT k-allele-sim-list.R) |
	parallel --col-sep " " $RSCRIPT k-allele-sim.R |
	gzip > k-allele-sim.csv.gz
