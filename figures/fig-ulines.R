#!/usr/bin/Rscript --vanilla
# see fig-hlines.R for a description
scriptname = sub("^--file=",'',grep("^--file=",commandArgs(),value=TRUE)[1])
pdfname = sub("\\.[^.]*$", ".pdf", scriptname)

library(extrafont)
library(RColorBrewer)

bcz = brewer.pal(5,"Set1")

pow = function(x,y) x**y

# find the valid, values of m that make the determinant 0
droots = function(s,u,Q,h) {
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

# find the roots
r = function(m,s,u,Q,h) {
  A = -(1-2*h)*s
  B = s*(Q*m*(1-u)-m*(1-u)+1+h*(-3+m*(1-2*Q)*(1-u)-u))
  C = -m*(1-u)*(1+h*(1-2*Q)*s)-u+h*s*(1+u)
  D = Q*m*(1-u)+u
  
  polyroot(c(D,C,B,A))
}

rroots = function(m,s,u,Q,h) {
  A = -(1-2*h)*s
  B = s*(Q*m*(1-u)-m*(1-u)+1+h*(-3+m*(1-2*Q)*(1-u)-u))
  C = -m*(1-u)*(1+h*(1-2*Q)*s)-u+h*s*(1+u)
  D = Q*m*(1-u)+u
  
  V = c(D,C,B,A)
  #V = V[1:max(which(V != 0))]

  r = polyroot(V)
  isR = sapply(Im(r),function(x) isTRUE(all.equal.numeric(x,0)))
  rr = ifelse(isR & 0 <= Re(r) & Re(r) <= 1,Re(r),NA)
  rr[which(!is.na(rr))]
}

d = function(m,s,u,Q,h) {
	A = -(1-2*h)*s
	B = s*(Q*m*(1-u)-m*(1-u) + 1+h*(-3+m*(1-2*Q)*(1-u)-u))
  C = -m*(1-u)*(1+h*(1-2*Q)*s)-u+h*s*(1+u)
	D = Q*m*(1-u)+u

	det = 18*A*B*C*D - 4*B*B*B*D + B*B*C*C - 4*A*C*C*C - 27*A*A*D*D
	det
}

deltaq = function(m,s,u,Q,h,q) {
  A = -(1-2*h)*s
  B = s*(Q*m*(1-u)-m*(1-u) + 1+h*(-3+m*(1-2*Q)*(1-u)-u))
  C = -m*(1-u)*(1+h*(1-2*Q)*s)-u+h*s*(1+u)
  D = Q*m*(1-u)+u

  del = A*q*q*q+B*q*q+C*q+D
  del
}

# returns the value of m that makes q an equilibrium value
mequ = function(s,h,u,Q,q=1/2) {
	m1 = ((-1 + q)*(u + q*s*(q + h*(1 - 2*q + u))))/
   ((q - Q + q*(q - q*Q + h*(-1 + q)*(-1 + 2*Q))*s)*(-1 + u))
   m1
}

# returns the value of m that makes p and equilibrium phenotype value
mequp = function(s,h,u,Q,p=1/2) {
  if(-1+2*h == 0) {
    q = p
  } else {
    q = (-h + sqrt(p - 2*p*h + h*h))/(1 - 2*h)
  } 
  mequ(s,h,u,Q,q)
}

ulines1 = function(h,ess,yue,Q,qh=1/2) {
  ret = matrix(0,length(ess),length(yue))
  for(i in seq_along(yue)) {
    u = yue[i]
    for(j in seq_along(ess)) {
      s = ess[j]
      m = mequ(s,h,u,Q,qh)
      dm = rroots(m=m,s=s,h=h,u=u,Q=Q)
      if(length(dm)>1) {
        rm = droots(s=s,h=h,u=u,Q=Q)
        m = max(rm[rm < m])
      }
      ret[j,i] = m
    }
  }
  ret
}

ulines2 = function(h,ess,yue,Q,ph=1/2) {
  ret = matrix(0,length(ess),length(yue))
  for(i in seq_along(yue)) {
    u = yue[i]
    for(j in seq_along(ess)) {
      s = ess[j]
      m = mequp(s,h,u,Q,ph)
      dm = rroots(m=m,s=s,h=h,u=u,Q=Q)
      if(length(dm)>1) {
        rm = droots(s=s,h=h,u=u,Q=Q)
        m = max(rm[rm < m])
      }
      ret[j,i] = m
    }
  }
  ret
}

adjustaxs = function(x) {
  xe = sapply(log10(x), function(u) as.expression(bquote(10^.(u))))
  ifelse(abs(log10(x)) > 3, xe, x)  
}

em = 10^seq(-8,0,0.1)
ess = 10^seq(-7,4,0.1)

width = 7.08
height = width*(2/3)

pdf(pdfname, width=width, height=height,family="Helvetica")
par(mgp=c(1.8,0.6,0),mai=c(0.55,0.55,0.1,0.1), las=1,cex.axis=1.2^-2)
par(mfrow=c(2,2))
par(cex = 1)

yue = c(1e-7,1e-6,1e-5,1e-4,1e-3)
Q = 0.01
aitch = c(1/2,0,1)

for(page in seq_along(aitch)) {
	h = aitch[page]

	####################################################################################
	# Threshold for 'b' to be most common allele

	hl = ulines1(h,ess,rev(yue),Q,1/2)
	matplot(log10(hl),log10(ess/hl),type="l",col=rev(bcz),lty=1,lwd=2,
	  xlab="m",ylab="s/m",xaxt="n",yaxt="n",xlim=c(-6,0),ylim=c(0,2)
	)
	text(-6,2,"a",cex=1.2^1,font=2,adj=c(0.5,1))

	x = 10^seq(-8,0,1)
	axis(1,log10(x),adjustaxs(x))
	y = c(1,3,10,30,100)
	axis(2,log10(y),adjustaxs(y))

	leg = sapply(log10(yue), function(u) {as.expression(bquote(10^.(u)))})

	legend(-5.25,2.2,legend=leg,ncol=2,fill=bcz,bty="n")
	text(-5.25,1.75,"u",cex=1.2)

	text(0,0,expression(q[infinity] >= 0.5),cex=1.2,adj=c(1,-0.2))

	####################################################################################
	# Threshold for blindness to be the most common phenotype

	hl = ulines2(h,ess,rev(yue),Q,1/2)
	matplot(log10(hl),log10(ess/hl),type="l",col=rev(bcz),lty=1,lwd=2,
	  xlab="m",ylab="s/m",xaxt="n",yaxt="n",xlim=c(-6,0),ylim=c(0,2)
	)
	text(-6,2,"b",cex=1.2^1,font=2,adj=c(0.5,1))

	x = 10^seq(-8,0,1)
	axis(1,log10(x),adjustaxs(x))
	y = c(1,3,10,30,100)
	axis(2,log10(y),adjustaxs(y))

	text(0,0,expression(a[infinity] >= 0.5),cex=1.2,adj=c(1,-0.2))

	####################################################################################
	# Threshold for 'b' to become "fixed"

	hl = ulines1(h,ess,rev(yue),Q,0.99)
	matplot(log10(hl),log10(ess/hl),type="l",col=rev(bcz),lty=1,lwd=2,
	  xlab="m",ylab="s/m",,xaxt="n",yaxt="n",xlim=c(-6,0),ylim=c(1.5,5)
	)
	text(-6,5,"c",cex=1.2^1,font=2,adj=c(0.5,1))

	x = 10^seq(-8,0,1)
	axis(1,log10(x),adjustaxs(x))
	y = c(1,10,100,1000,10000,100000)
	axis(2,log10(y),adjustaxs(y))

	text(0,1.5,expression(q[infinity] >= 0.99),cex=1.2,adj=c(1,-0.2))

	####################################################################################
	# Threshold for blindness to become "fixed"

	hl = ulines2(h,ess,rev(yue),Q,0.99)
	matplot(log10(hl),log10(ess/hl),type="l",col=rev(bcz),lty=1,lwd=2,
	  xlab="m",ylab="s/m",,xaxt="n",yaxt="n",xlim=c(-6,0),ylim=c(1.5,5)
	)
	text(-6,5,"d",cex=1.2^1,font=2,adj=c(0.5,1))

	x = 10^seq(-8,0,1)
	axis(1,log10(x),adjustaxs(x))
	y = c(1,10,100,1000,10000,100000)
	axis(2,log10(y),adjustaxs(y))

	text(0,1.5,expression(a[infinity] >= 0.99),cex=1.2,adj=c(1,-0.2))
}
invisible(dev.off())
embed_fonts(pdfname, options="-DPDFSETTINGS=/prepress")
