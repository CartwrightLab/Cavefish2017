#!/usr/bin/Rscript --vanilla

# script setup
library(RColorBrewer)
library(extrafont)
scriptname = sub("^--file=",'',grep("^--file=",commandArgs(),value=TRUE)[1])
pdfname = sub("\\.[^.]*$", ".pdf", scriptname)

cz <- brewer.pal(8, "Set1")[-6]
#cz = "gray25"

# setup parameters
s = 0.1
m = 0.01
u = 1e-6
Q = 0.01
q = seq(0,1,0.001)
h = 0.0

# calculate the change in allele frequency (delta-q)
dq = function(qq,s) {
	q0 = qq
	qq = ((1+s)*qq*qq + qq*(1-qq)) / ((1+s)*qq*qq + 2*qq*(1-qq)+(1-qq)*(1-qq))
	qq = qq*(1-m)+m*Q
	qq = qq+(1-qq)*u
	qq-q0
}

# calculate the roots of delta-q
dqroots = function(s) {
  A = -(1-2*h)*s
  B = s*(Q*m*(1-u)-m*(1-u)+1+h*(-3+m*(1-2*Q)*(1-u)-u))
  C = -m*(1-u)*(1+h*(1-2*Q)*s)-u+h*s*(1+u)
  D = Q*m*(1-u)+u
  r = polyroot(c(D,C,B,A))
  isR = sapply(Im(r),function(x) isTRUE(all.equal.numeric(x,0)))
  rr = ifelse(isR & 0 <= Re(r) & Re(r) <= 1,Re(r),NA)
  rr[which(!is.na(rr))]
}

# for several values of s, calculate dq for many values of q
# also estimate the roots.
mm = c()
rr = list()
for(s in c(0,0.05,0.10,0.15,0.20,0.25,0.30)) {
	d = dq(q,s)
	mm = cbind(mm,d)
	r = dqroots(s)
	rr[[as.character(s)]] = r
}

# begin figure generation.
pdf(pdfname, width=7, height=3.46,family="Helvetica")
par(mai=c(0.6, 0.6, 0.1, 0.1), mgp=c(1.7,0.6,0), mfcol=c(1,2),
	las=1,cex.axis=1.2^-1)

# setup figure
yl = c(-0.010,0.015)
matplot(q,mm,type="n",lty=1,col=cz,lwd=1.5,xlab="q",
	ylab=expression(Delta*q),ylim=yl)

# plot the y=0 axis
abline(h=0,lty=2)
# plot delta-q
matlines(q,mm,type="l",lty=1,col=cz,lwd=1.5,
	xlab="q",ylim=yl)
# plot location of equilibria
for(i in seq_along(rr)) {
	col = cz[i]
	for(r in rr[[i]]) {
		y = yl[1] + (yl[2]-yl[1])*c(-0.02,0.06)
		lines(c(r[1],r[1]),y,col=col,lwd=2,lty=1)
	}
}

# Same as above, but in a smaller area
yl = c(-0.00025*2/3,0.00025)
matplot(q,mm,type="n",lty=1,col=cz,lwd=1.5,xlab="q",
	ylab=expression(Delta*q), yaxt="n",
	ylim=yl,xlim=c(0,0.05))
axis(2,at=c(-0.0001,0.000,0.0001,0.0002),labels=c("-0.0001","0.0000","0.0001","0.0002"))
abline(h=0,lty=2)
matlines(q,mm,type="l",lty=1,col=cz,lwd=1.5,
	xlab="q")
for(i in seq_along(rr)) {
	col = cz[i]
	for(r in rr[[i]]) {
		y = yl[1] + (yl[2]-yl[1])*c(-0.02,0.06)
		lines(c(r[1],r[1]),y,col=col,lwd=2,lty=1)
	}
}

# finalize figure
invisible(dev.off())
embed_fonts(pdfname, options="-DPDFSETTINGS=/prepress")
