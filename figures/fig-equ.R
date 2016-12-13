#!/usr/bin/Rscript --vanilla

# This figure describes how the locations of equilibria change as
# s increases from near zero to near 1. It contains circles that
# mark the true equilibria and lines that denote approximations.

# setup figure
library(RColorBrewer)
library(extrafont)

scriptname = sub("^--file=",'',grep("^--file=",commandArgs(),value=TRUE)[1])
pdfname = sub("\\.[^.]*$", ".pdf", scriptname)

cz <- brewer.pal(8, "Set1")

# calculate the roots of delta-q
h = function(m,s,u,qn=sqrt(u)) {
	A = -s
	B = s*(qn*m*(1-u)-m*(1-u)+1)
	C = -m*(1-u)-u
	D = qn*m*(1-u)+u
	r = polyroot(c(D,C,B,A))
	isR = sapply(Im(r),function(x) isTRUE(all.equal.numeric(x,0)))
	ifelse(isR,Re(r),NA)
}

# default parameters
s = 0.1
m = 0.01
u = 1e-6
qn = sqrt(u/0.01)

s = 10^seq(-2,0,0.0001)
smin = 0.01
smax = 1

# start figure
pdf(pdfname, width=3.46, height=3.46,family="Helvetica")
par(mai=c(0.55, 0.55, 1/72, 1/72), mgp=c(1.7,0.6,0),
	las=1,cex.axis=1.2^-1)

# figure out plotting range in q-space
qmin = h(m,smin,u,qn)[1]
qmax = h(m,smax,u,qn)[3]
q = seq(qmin,qmax,(qmax-qmin)/10000)

# calculate a line representing the true equilibria in the graph
sh = (m*(q-qn)-u*(1-q))/(q*q*(1-q-m*(1-qn)))
shl = log10(sh)
# estimate the location of circles so that they are evenly
# spaced along the line
qd = q[-1]-q[-length(q)]
sd = (shl[-1]-shl[-length(shl)])/(log10(smax)-log10(smin))
ld = sqrt(qd*qd+sd*sd)
sumd = sum(ld)
cumd = cumsum(ld)
seqd = seq(0.0,sumd,sumd/60)
w = sapply(seqd, function(x) which(cumd >= x)[1])
hd = seqd-(cumd[w]-ld[w])
xx = q[w]
yy = shl[w]
aa = q[w+1]-xx
bb = shl[w+1]-yy
a = hd/sqrt(1+bb*bb/(aa*aa))
b = a*bb/aa*(log10(smax)-log10(smin))

# Calculate approximations for q-hat given s
s = sh
qha = ((m+u)-sqrt((m+u)^2-4*(1-m)*s*(u+m*qn)))/(2*(1-m)*s)
qhbu = 0.5*(1-m*(1-qn)-sqrt((1-m*(1-qn))^2-4*m/s))
qhbl = ((m+u)+sqrt((m+u)^2-4*(1-m)*s*(u+m*qn)))/(2*(1-m)*s)
qhc = 0.5*(1-m*(1-qn)+sqrt((1-m*(1-qn))^2-4*m/s))

# Cleanup bad approximations
llim = 4*m/(1-m*(1-qn))^2
ulim = (m+u)^2/(4*(1-m)*(m*qn+u))

qha[s > ulim] = NA
qhbu[s < llim | s > ulim] = NA
qhbl[s < llim | s > ulim] = NA
qhc[s < llim ] = NA

# Setup plot
y = cbind(qha,qhbl,qhbu,qhc)
matplot(y,log10(s),col=cz,type="n",xlim=c(0,1),pch=19,
	xlab="q",ylab="s",yaxt="n")

# Plot the estimate range of three equilibria
abline(log10(llim),0,lty=2,col="gray25")
abline(log10(ulim),0,lty=2,col="gray25")

# Plot approximations
matlines(y,log10(s),col=cz[c(1,2,3,4)],lty=1,
	xlim=c(0,1),pch=19,lwd=1.5)

# Plot truth
matpoints(xx+a,yy+b,pch=21,cex=1.2^-2,col="#000000CF")

# Plot axes
p = c(1,0.3,0.1,0.03,0.01)
axis(2,at=log10(p),p)

legend("bottomright",c(
	expression(hat(q)),
	expression(hat(q)["a,1"]),
	expression(hat(q)["b,1"]),
	expression(hat(q)["b,2"]),
	expression(hat(q)["c,2"])
	),
  cex=1.2^-2, pt.cex=1,
  col=c("black",cz[1:4]), lty=c(NA,1,1,1,1),
  lwd=c(1.5,2,2,2,2), pch=c(21,NA,NA,NA,NA)
)

# finalize figure
invisible(dev.off())
embed_fonts(pdfname, options="-DPDFSETTINGS=/prepress")
