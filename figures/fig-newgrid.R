#!/usr/bin/Rscript --vanilla

# setup figure
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
f = function(m,s,u,Q,h,q0) {
  A = -(1-2*h)*s
  B = s*(Q*m*(1-u)-m*(1-u)+1+h*(-3+m*(1-2*Q)*(1-u)-u))
  C = -m*(1-u)*(1+h*(1-2*Q)*s)-u+h*s*(1+u)
  D = Q*m*(1-u)+u
  r = polyroot(c(D,C,B,A))
  isR = sapply(Im(r),function(x) isTRUE(all.equal.numeric(x,0)))
  rr = ifelse(isR & 0 <= Re(r) & Re(r) <= 1,Re(r),NA)
  if(sum(is.na(rr)) > 0)
    return(rr[which(!is.na(rr))])
  if(q0 < rr[2])
    return(rr[1])
  return(rr[3])
}

# returns the value of m that makes q an equilibrium value
mequ = function(s,h,u,Q,q=1/2) {
  m1 = ((-1 + q)*(u + q*s*(q + h*(1 - 2*q + u))))/
   ((q - Q + q*(q - q*Q + h*(-1 + q)*(-1 + 2*Q))*s)*(-1 + u))
   m1
}

# find the valid roots of delta q
rroots = function(m,s,u,Q,h) {
  A = -(1-2*h)*s
  B = s*(Q*m*(1-u)-m*(1-u)+1+h*(-3+m*(1-2*Q)*(1-u)-u))
  C = -m*(1-u)*(1+h*(1-2*Q)*s)-u+h*s*(1+u)
  D = Q*m*(1-u)+u
  
  V = c(D,C,B,A)

  r = polyroot(V)
  isR = sapply(Im(r),function(x) isTRUE(all.equal.numeric(x,0)))
  rr = ifelse(isR & 0 <= Re(r) & Re(r) <= 1,Re(r),NA)
  rr[which(!is.na(rr))]
}

# find the valid values of m that make the determinant 0
droots = function(s,u,Q,h) {
  pow = function(x,y) x**y

  E = s*pow(s*(1 + h*(-1 + u)) + u,2)*(-4*u + 8*h*u + pow(h,2)*s*pow(1 + u,2))

  D = -2*s*(-1 + u)*(2*pow(h,4)*(-1 + 2*Q)*pow(s,3)*u*(-1 + pow(u,2)) - 2*(Q*s*(s + 3*s*u + (9 - 5*u)*u) + u*(3*u + s*(-5 - 3*s + 5*u))) + 
      h*(12*pow(u,2) + pow(s,2)*(-1 - u*(22 + u) + Q*(9 + u*(14 + u))) + s*u*(-39 + u*(26 + u) - Q*(-63 + u*(38 + u)))) + 
      pow(h,3)*pow(s,2)*(3 - s*(1 + u)*(-4 - u*(1 + u) + Q*(5 + (-2 + u)*u)) + u*(18 - u*(5 + 4*u) + 4*Q*(-9 + u*(5 + 2*u)))) + 
      pow(h,2)*s*(3*(13 - 3*u)*u - pow(s,2)*(1 + u)*(2 + u) + 4*s*u*(3 + u*(2 + u)) + 
         Q*(pow(s,2)*(1 + u)*(3 + u) + 2*u*(-27 + u*(18 + u)) - s*(9 + u*(-16 + u*(11 + 4*u))))))
  C = s*pow(-1 + u,2)*(12*(-1 + 2*h)*u + s*(1 + 9*(2 - 3*Q)*Q - 22*u + 2*Q*(2 + 9*Q)*u + pow(-1 + Q,2)*pow(u,2) + 
         2*h*(3 + Q*(-9 + u)*(5 + u) + u*(28 + u) - 2*pow(Q,2)*(-27 + u*(18 + u))) + 
         pow(h,2)*(-15 - 2*u*(9 + u) + 4*Q*(27 + pow(u,2)) + 4*pow(Q,2)*(-27 + u*(18 + u)))) + 
      pow(h,2)*pow(s,3)*(6 + u*(6 + u) + 2*pow(h,2)*pow(1 - 2*Q,2)*(-1 + 3*pow(u,2)) - 2*Q*(9 + u*(8 + u)) + pow(Q,2)*(13 + u*(10 + u)) - 
         2*h*(-1 + 2*Q)*(-6 - u*(7 + 3*u) + Q*(5 + u*(4 + 3*u)))) + 
      2*pow(s,2)*(-6*(-1 + Q)*(Q + (-1 + Q)*u) + h*(3 + 5*u*(4 + u) + pow(Q,2)*(3 + u)*(7 + 5*u) - 2*Q*(14 + u*(21 + 5*u))) + 
         pow(h,3)*(-1 + 2*Q)*(9 + (5 - 6*u)*u + 2*Q*(-9 + u*(5 + 6*u))) - pow(h,2)*(2 + u*(11 + 10*u) + 2*pow(Q,2)*u*(15 + 11*u) - 2*Q*(8 + u*(19 + 16*u)))))
  B = -2*s*pow(-1 + u,3)*(-2 + 2*pow(h,4)*pow(-1 + 2*Q,3)*pow(s,3)*u - 2*pow(-1 + Q,2)*pow(s,2)*(3*Q + (-1 + Q)*u) + 
      (-1 + Q)*s*(1 + 9*Q + (-1 + Q)*u) - pow(h,3)*pow(s - 2*Q*s,2)*(-5 + 4*u - 8*Q*u + s*(-4 + 3*Q + 3*(-1 + Q)*u)) + 
      pow(h,2)*(-1 + 2*Q)*s*(-9 + pow(Q,2)*s*(3*(-4 + s) + (-20 + s)*u) + s*(-4 - 8*u + s*(2 + u)) + Q*(2*(9 + u) + s*(15 - 5*s + 28*u - 2*s*u))) + 
      h*(4 + (-1 + Q)*pow(s,2)*(3 + 6*u + 2*Q*(-13 + 12*Q - 11*u + 8*Q*u)) - s*(2 + u + Q*(-34 - 5*u + 4*Q*(9 + u)))))
  A = pow(1 - Q + h*(-1 + 2*Q),2)*pow(s,2)*(1 + s*(2*h*pow(1 - 2*Q,2) - 4*(-1 + Q)*Q + pow(h,2)*pow(1 - 2*Q,2)*s))*pow(-1 + u,4)

  r = polyroot(c(E,D,C,B,A))
  isR = sapply(Im(r),function(x) isTRUE(all.equal.numeric(x,0)))
  rr = ifelse(isR & 0 <= Re(r) & Re(r) <= 1,Re(r),NA)
  rr[which(!is.na(rr))]
}

# find m* = sup { m : q_infinity >= qh }
mline_low = function(h,ess,u,Q,qh=1/2) {
  ret = rep(0,length(ess))
  for(j in seq_along(ess)) {
    s = ess[j]
    m = mequ(s,h,u,Q,qh)
    dm = rroots(m=m,s=s,h=h,u=u,Q=Q)
    if(length(dm)>1) {
      rm = droots(s=s,h=h,u=u,Q=Q)
      m = max(rm[rm < m])
    }
    ret[j] = m
  }
  ret
}

# find m* = inf { m : q_infinity <= qh && there is a sing root }
mline_high = function(h,ess,u,Q,qh=1/2) {
  ret = rep(0,length(ess))
  for(j in seq_along(ess)) {
    s = ess[j]
    m = mequ(s,h,u,Q,qh)
    dm = rroots(m=m,s=s,h=h,u=u,Q=Q)
    if(length(dm)>1) {
      rm = droots(s=s,h=h,u=u,Q=Q)
      m = min(rm[rm > m])
    }
    ret[j] = m
  }
  ret
}

# find s* = inf { s : q invades deterministically }
sstar = function(m,N,u) {
  -2+(1 + 4*N)/(2*(1 - m)*N*(1 - u))
}

# wright 1931 provides an equation for the mean value of the interaction of m + u + v
sstar_wright = function(m,N,u,Q) {
  log(m*(1-Q)/(u+m*Q))/(4*N)
}

# Setup the figure regions
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

# helper function to make pretty axes
adjustaxs = function(x) {
  xe = sapply(log10(x), function(u) as.expression(bquote(10^.(u))))
  ifelse(abs(log10(x)) > 3, xe, x)  
}

# setup figure
pdf(pdfname, width=width, height=height,family="Helvetica")
par(mgp=c(1.5,0.6,0),las=1,cex.axis=1.2^-2)
split.screen(scrmat)
par(cex = 1)

h = 0
u = 1e-6
Q = sqrt(u/0.01)

tmv = c("Inf","10000","5000000","10000mk")

tml = list(
  c("Constant",expression(N == infinity), expression(t == infinity)),
  c("Constant",expression(N == 1000), expression(t == 10000)),
  c("Constant",expression(N == 1000), expression(t == 5 %*% 10^6)),
  c("Episodic",expression(N == 1000), expression(t == 10000))
)

# Plot the key
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

# Plot axis labels.
screen(6)
par(mai = c(43,42,43,34)/72)
par(mgp=c(2,0.6,0))
title(xlab="Migration (m)", ylab="Selection (s)")

# Plot the Individual panels
for(tm in seq_along(tmv) ) {
  screen(tm)
  # Margins vary for each panel
  if(tm == 1) {
    par(mai=c(2, 38, 38, 2)/72 )
  } else if(tm == 2) {
    par(mai=c(2, 2, 38, 38)/72 )
  } else if(tm == 3) {
    par(mai=c(38, 38, 2, 2)/72 )
  } else {
    par(mai=c(38, 2, 2, 38)/72 )
  }

  # parameter ranges
  m = rev(10^-seq(0,8,0.05))
  s = rev(10^-seq(-2,6,0.05))
  ss = (6+log10(s))/8
  mm = 1+log10(m)/8

  # Either calculate analytic results, or read a simulation results file
  if(tmv[tm] == "Inf") {
    o = outer(m,s,Vectorize(function(x,y) f(x,y,u,Q,h,Q)))
  } else {
    fn = sprintf("../data/n1000/drift-%s.txt",tmv[tm])
    o = as.matrix(read.table(fn,check.names=F))
  }

  # Prepare contour map
  fz = o
  xsel = seq_along(m)
  ysel = seq_along(s)
  fz = fz[xsel,ysel]

  fx = seq(0, 1, length.out = nrow(fz))
  fy = seq(0, 1, length.out = ncol(fz))
  
  mlim = log10(range(m[xsel]))
  slim = log10(range(s[ysel]))

  # calculate m* lines.
  mstar_low = (log10(mline_low(h,s,u,Q,0.5))-mlim[1])/(mlim[2]-mlim[1])
  mstar_high = (log10(mline_high(h,s,u,Q,0.5))-mlim[1])/(mlim[2]-mlim[1])

  x = rev(10^-seq(0,8,1))
  xx = (log10(x)-mlim[1])/(mlim[2]-mlim[1])
  y = rev(10^-seq(-2,6,1))
  yy = (log10(y)-slim[1])/(slim[2]-slim[1])

  # begin plot
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
  .filled.contour(fx,fy,fz,levels=keylevels,col=bcz)

  # draw m* (aka s*) lines
  lines(mstar_high,ss,lwd=2,col="white",lty=2)
  lines(mstar_low,ss,lwd=2,col="white")

  # Draw lines such that they all share a dot at the intersection
  if(tm == 3) {
    # 4 N m == 1
    gridx = (log10(1/(4*1000))-mlim[1])/(mlim[2]-mlim[1])
    # 2 N s == 1
    gridy = (log10(1/(2*1000))-slim[1])/(slim[2]-slim[1])

    lines(c(gridx,gridx),c(gridy,0),lwd=2,lty=3,col="white")
    lines(c(gridx,gridx),c(gridy,1),lwd=2,lty=3,col="white")    
    lines(c(gridx,0),c(gridy,gridy),lwd=2,lty=3,col="white")
    lines(c(gridx,1),c(gridy,gridy),lwd=2,lty=3,col="white")
    hh = strheight("|M")
    text(1-0.5*hh,gridy+0.6*hh,expression(2*N*s == 1),pos=2,col="white",offset=0)
    text(gridx-0.6*hh,1-0.6*hh,expression(4*N*m == 1),pos=2,col="white",srt=90,offset=0)
  }

  if(tm %in% c(2,4)) {
    # m Q == u
    gridx = (log10(u/Q)-mlim[1])/(mlim[2]-mlim[1])
    # 2 N s == 1
    # use gridy to match the dots from the other panel
    gridy = (log10(1/(2*1000))-slim[1])/(slim[2]-slim[1])    
    lines(c(gridx,gridx),c(gridy,0),lwd=2,lty=3,col="white")
    lines(c(gridx,gridx),c(gridy,1),lwd=2,lty=3,col="white")
    hh = strheight("|M")
    text(gridx-0.6*hh,1-0.6*hh,expression(Q*m == u),pos=2,col="white",srt=90,offset=0)    
  }

  # Label graph
  hght <- 1.75*strheight("|M")
  wdth <- strwidth("N")-strwidth("t")

  txt = tml[[tm]]
  if(length(txt)==2) {
    ww = c(0,wdth)
  } else {
    ww = c(0,0,wdth) 
  }
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
