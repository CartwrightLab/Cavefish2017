#!/usr/bin/Rscript --vanilla
#setup figure
scriptname = sub("^--file=",'',grep("^--file=",commandArgs(),value=TRUE)[1])
pdfname = sub("\\.[^.]*$", ".pdf", scriptname)

library(extrafont)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

bcz = brewer.pal(10,"RdYlBu")
keylevels = seq(0,1,length.out=length(bcz)+1)
keytext = keylevels

linecol = "#000000FF"

# find the final population frequency
# additive model with k alleles
f = function(m,s,k,u,Q,q0) {
  B = -s + m*s - 2*k*m*s + s*u - 2*k*s*u - m*s*u + 2*k*m*s*u
  C = -2*k*m + s - m*s + 2*k*m*Q*s - 2*k*u + 2*k*m*u - s*u + 2*k*s*u + m*s*u - 2*k*m*Q*s*u
  D = 2*k*m*Q + 2*k*u - 2*k*m*Q*u
  r = polyroot(c(D,C,B))
  isR = sapply(Im(r),function(x) isTRUE(all.equal.numeric(x,0)))
  rr = ifelse(isR & 0 <= Re(r) & Re(r) <= 1,Re(r),NA)
  if(sum(is.na(rr)) > 0)
    return(rr[which(!is.na(rr))])
  if(q0 < rr[2])
    return(rr[1])
  return(rr[3])
}

# returns the value of m that makes q an equilibrium value
mequ = function(s,k,u,Q,q=1/2) {
  ((-1 + q)*(2*k*u + q*s*(1 + (-1 + 2*k)*u)))/((-((-1 + q)*q*s) + 2*k*(q - Q)*(1 + q*s))*(-1 + u))
}

# find m* = sup { m : q_infinity >= qh }
mstar = function(ess,k,u,Q,qh=1/2) {
  ret = rep(0,length(ess))
  for(j in seq_along(ess)) {
    s = ess[j]
    m = mequ(s,k,u,Q,qh)
    ret[j] = m
  }
  ret
}

sstar = function(m,N) {
  (1+4*m*N)/(2*N*(1-m))
}

# Figure Layout
width = 7.08
keyheight = (7.08-2*(34/72))/length(bcz)+50/72
height = width+keyheight
h1 = 0
h2 = keyheight/height
h3 = h2+(1-h2)/2
h4 = 1
w1 = 0
w2 = 0.5
w3 = 1
scrmat = rbind(
  c(w1,w2,h3,h4),
  c(w2,w3,h3,h4),
  c(w1,w2,h2,h3),
  c(w2,w3,h2,h3),
  c(w1,w3,h1,h2),
  c(w1,w3,h2,h4)
)

# pretty axes helper function
adjustaxs = function(x) {
  xe = sapply(log10(x), function(u) as.expression(bquote(10^.(u))))
  ifelse(abs(log10(x)) > 3, xe, x)  
}

# start figure
pdf(pdfname, width=width, height=height,family="Helvetica")
par(mgp=c(1.5,0.6,0),las=1,cex.axis=1.2^-2)
split.screen(scrmat)
par(cex = 1)

u = 1e-6
Q = 0.01

tmv = c("Inf1",1,"Inf12",12)

tml = list(
  c(expression(k == 1),expression(N == infinity),expression(t == infinity)),
  c(expression(k == 1),expression(N == 1000),expression(t == 10000)),
  c(expression(k == 12),expression(N == infinity),expression(t == infinity)),
  c(expression(k == 12),expression(N == 1000),expression(t == 10000))
)

# key
screen(5)
par(mai = c(40,34,10,34)/72)
plot.new()
plot.window(ylim=c(0,1), xlim=range(keylevels), xaxs="i", yaxs="i")
rect(keylevels[-length(keylevels)], 0, keylevels[-1L], 1, col = bcz)
x = keytext
box()
axis(1,keylevels,x)
text(keylevels[1]+0.0*(keylevels[2]-keylevels[1]),0.5,"Sightedness",col="white",cex=1.2^1,font=2,pos=4)
text(keylevels[11]+0.0*(keylevels[2]-keylevels[1]),0.5,"Blindness",col="white",cex=1.2^1,font=2,pos=2)
title(xlab="q")

# axes labels
screen(6)
par(mai = c(43,42,43,34)/72)
par(mgp=c(2,0.6,0))
title(xlab="Migration (m)", ylab="Selection (s)")

# plot panels
for(tm in seq_along(tmv) ) {
  screen(tm)
  # each panel has different margins
  if(tm == 1) {
    par(mai=c(2, 38, 38, 2)/72 )
  } else if(tm == 2) {
    par(mai=c(2, 2, 38, 38)/72 )
  } else if(tm == 3) {
    par(mai=c(38, 38, 2, 2)/72 )
  } else {
    par(mai=c(38, 2, 2, 38)/72 )
  }

  m = rev(10^-seq(0,8,0.1))
  s = rev(10^-seq(-2,6,0.1))
  ss = (6+log10(s))/8
  mm = 1+log10(m)/8

  k = tmv[tm]

  # Either calculate analytic results, or read a simulation results file  
  if(grepl("^Inf",k)) {
    k = as.numeric(sub("^Inf","",k))
    o = outer(m,s,Vectorize(function(x,y) f(x,y,k,u,Q,Q)))
  } else {
    k = as.numeric(k)
    fn = sprintf("../data/kalleles/drift-n1000-k%02d-q.txt",k)
    o = as.matrix(read.table(fn,check.names=F))
  }

  # setup data for contour plot
  fz = o
  xsel = seq_along(m)
  ysel = seq_along(s)
  fz = fz[xsel,ysel]

  fx = seq(0, 1, length.out = nrow(fz))
  fy = seq(0, 1, length.out = ncol(fz))
  
  mlim = log10(range(m[xsel]))
  slim = log10(range(s[ysel]))


  # estimate m* (aka s*) for k=1 and k=k
  mstar_1 = (log10(mstar(s,1,u,Q,0.5))-mlim[1])/(mlim[2]-mlim[1])
  mstar_k = (log10(mstar(s,k,u,Q,0.5))-mlim[1])/(mlim[2]-mlim[1])

  x = rev(10^-seq(0,8,1))
  xx = (log10(x)-mlim[1])/(mlim[2]-mlim[1])
  y = rev(10^-seq(-2,6,1))
  yy = (log10(y)-slim[1])/(slim[2]-slim[1])

  # begin plot
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
  .filled.contour(fx,fy,fz,levels=keylevels,col=bcz)

  # draw m* (aka s*) lines
  lines(mstar_1,ss,lwd=2,col="white",lty=3)
  lines(mstar_k,ss,lwd=3,col="white",lty=1)
  
  # labels
  hght <- 1.75*strheight("|M")
  wdth1 <- strwidth("N")-strwidth("k")
  wdth3 <- strwidth("N")-strwidth("t")
  txt = tml[[tm]]

  ww = c(wdth1,0,wdth3)

  text(0.65+ww,0+hght*(1.0*length(txt)-(seq_along(txt)-1)),
    txt,pos=4,col="white",cex=1.2,font=2)
  text(0,1,letters[tm],adj=c(-0.5,1.5),col="white",cex=1.2^2,font=2)

  # pretty axes
  box()
  x = 10^seq(-8,0,1)
  xx = (log10(x)-mlim[1])/(mlim[2]-mlim[1])
  x = adjustaxs(x)
  y = 10^seq(-6,2,1)
  yy = (log10(y)-slim[1])/(slim[2]-slim[1])
  y = adjustaxs(y)

  if(tm == 2 || tm == 4) {
    x[1] = ""
  }
  if(tm == 1 || tm == 2) {
    axis(3,xx,x)
    y[1] = ""
  } else {
    axis(1,xx,x)
  }
  if(tm == 1 || tm == 3) {
    axis(2,yy,y)
  } else {
    axis(4,yy,y)
  }
  #title(ylab="s",xlab="m")
}

# finalize figure
invisible(dev.off())
embed_fonts(pdfname, options="-DPDFSETTINGS=/prepress")
