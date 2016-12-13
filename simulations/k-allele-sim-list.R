#!/usr/bin/Rscript --vanilla

# This script generates a list of simulations to run

# Rscript --vanilla k-allele-sim-list.R | gzip > k-allele-jobs.txt.gz
# zcat k-allele-jobs.txt.gz | parallel --col-sep " " Rscript --vanilla k-allele-sim.R

r = 100

G = 10000
s = 10**seq(-6,2,0.1)
m = 10**seq(-8,0,0.1)
h = c(0.5)
# N = c(1000,100)
# k = c(1,2,4,6,12,18)
# u = c(1e-6,1e-5)
# Q = c(0.01,0.1)
N = c(1000)
k = c(1,2,4,6,12,18)
u = c(1e-6)
Q = c(0.01)

params = list(N=N,k=k,h=h,u=u,Q=Q,s=s,m=m)

sim = sprintf("G=%s",G)
for(p in names(params)) {
	ss = sprintf("%s=%s",p,params[[p]])
	if(p == "Q") {
		ss = paste0(ss,sprintf(" q0=%s",params[[p]]))
	}
	sim = paste(rep(sim,each=length(ss)),ss,sep=" ")
}

for(x in sim) {
	x = rep(x,each=r)
	x = paste0(x,collapse="\n")
	cat(x,"\n")
}
