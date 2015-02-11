#1D continuous space model for examining stochastic evolutionary processes during range expansion
# Author: Ben Phillips


where<-Sys.info()["sysname"]
if (where=="Darwin"){
  setwd("~/evo-dispersal/evostoch/")
  #system(paste("R CMD SHLIB Density1D.c"))
  system(paste("R CMD SHLIB PointMetrics1D.c"))
}
if (where=="Linux"){
  setwd("/scratch/jc227089/evo-dispersal/evostoch/")
  #system(paste("R CMD SHLIB Density1D.c"))
  system(paste("R CMD SHLIB -L/usr/lib64/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4/ PointMetrics1D.c"))
}
dyn.load("PointMetrics1D.so")




# Initialises a matrix with n individuals
# individuals have a location (spX, spY) and are 
# initialised with phenotype for D and H
init.inds<-function(n, spX, Hmean, Dmean, h2H, h2D, VPH, VPD){
	X<-runif(n, -spX, spX) #position
	D<-rnorm(n, Dmean, sqrt(h2D*VPD)) # individual sd
	H<-rnorm(n, Hmean, sqrt(h2H*VPH)) #Habitat match
	DP<-D+rnorm(n, 0, sqrt((1-h2D)*VPD))
	HP<-H+rnorm(n, 0, sqrt((1-h2H)*VPH))
	cbind(X, D, H, DP, HP)		
}


# hassel commins population growth where b=1 (contest competition)
# generates expected fecundity
#hass.comm<-function(lambda, K, N){
#	lambda/(1+(lambda-1)*N/K)
#}

# Beverton Holt population growth (same as Hass Comins above)
bev.holt<-function(N, R, a){
  #a<-(R0-1)/K
  R/(1+a*N)
}

# gets density, mean and (if evovar==TRUE) variance in breeding values
metrics<-function(popmatrix, evovar=TRUE){
	if (is.null(popmatrix)) return(NULL)
	out<-.Call("metrics", R_X=popmatrix[,"X"], R_H=popmatrix[,"H"],
	R_D=popmatrix[,"D"], R_n=nrow(popmatrix), R_ev=as.numeric(evovar))
	colnames(out)<-c("Nx", "MeanD", "MeanH", "SDD", "SDH")
	if (evovar==FALSE) out<-out[,1:3]
	out
}

# gets density, mean and (if evovar==TRUE) variance in breeding values
# differs from above in that calculates through space and then approximates back to the individual
ca.metrics<-function(popmatrix, evovar=TRUE){
	if (is.null(popmatrix)) return(NULL)
  bins<-floor(min(popmatrix[,"X"])):ceiling(max(popmatrix[,"X"]))
  cn<-c("Nx", "MeanD", "MeanH", "SDD", "SDH")
  out.sp<-.Call("sum_metrics", R_X=popmatrix[,"X"], R_H=popmatrix[,"H"],
             R_D=popmatrix[,"D"], R_n=nrow(popmatrix), R_bins=bins, R_nbins=length(bins),
             R_ev=as.numeric(evovar), R_bw=1)
  if (evovar==FALSE) {out.sp<-out.sp[,1:3]; cn<-cn[1:3]}
  out.ind<-c()
  for (ii in 1:ncol(out.sp)){
  	out.ind<-cbind(out.ind, approx(x=bins, y=out.sp[,ii], xout=popmatrix[,"X"])$y)
  }
  colnames(out.ind)<-cn
  matrix(out.ind, ncol=length(cn), dimnames=list(NULL, cn))
}

# gets density, mean and (if evovar==TRUE) variance in breeding values
# differs from above in that bw is adjustable and metrics are scored to spatial bins rather than individuals
sum.metrics<-function(popmatrix, evovar=TRUE, bw=1){
	if (is.null(popmatrix)) return(NULL)
  bins<-floor(min(popmatrix[,"X"])):ceiling(max(popmatrix[,"X"]))
  out<-.Call("sum_metrics", R_X=popmatrix[,"X"], R_H=popmatrix[,"H"],
             R_D=popmatrix[,"D"], R_n=nrow(popmatrix), R_bins=bins, R_nbins=length(bins),
             R_ev=as.numeric(evovar), R_bw=bw)
  colnames(out)<-c("Nx", "MeanD", "MeanH", "SDD", "SDH")
  if (evovar==FALSE) out<-out[,1:3]
  out<-cbind(X=bins, out)
  out
}

#Finds a mate for each individual
mate<-function(popmatrix, dens){
	if (is.null(popmatrix)) return(NULL)
	out<-.Call("mate", R_X=popmatrix[,"X"], R_dens=dens, R_n=nrow(popmatrix),
		R_bw=1)
	out
}


# Function defining the strength of stabilising selection
h.surv<-function(n=1, H, Opt=0, k=2){
  n*exp(-k*(H-Opt)^2)
}

VE<-function(h, VP){
	(1-h)*VP
}

Allee<-function(dens, theta=0.5){
	dens/(dens+theta)
}


# reproduces individuals based on density, kills parents
# kills offspring based on habitat
# kills dispersing offspring based on dispersal cost
# disperses individuals based on inherited P
repro.disp<-function(popmatrix, a, R0, h2H, h2D, VPH, VPD, evovar=TRUE, Allee.fx=FALSE){
	if (is.null(popmatrix)) return(NULL)
	if (nrow(popmatrix)==1) return(NULL) #cases where only one individual left
	SDVED<-sqrt(VE(h2D, VPD)) # define environmental variances
	SDVEH<-sqrt(VE(h2H, VPH))
	# calculate traits through space
	mets<-metrics(popmatrix, evovar)
	#print('mets OK')
	# reproduction
	hab<-h.surv(1, popmatrix[,"HP"])
	if (Allee.fx) {hab<-hab*Allee(dens=mets[,"Nx"])}
	offs<-rpois(nrow(popmatrix), hab*bev.holt(mets[,"Nx"], R0, a)) # surviving offspring
	if (sum(offs)<=1) return(NULL) #catch extinction
	m<-mate(popmatrix, mets[,"Nx"]) #find a mate
	m[m==-99]<-NA # catch lonely ones
	# calculate new G and P
	if (evovar==FALSE){
		SDD<-sqrt(h2D*VPD)
		SDH<-sqrt(h2H*VPH)
	} else {
		SDD<-mets[,"SDD"]
		SDH<-mets[,"SDH"]
	}
	#get midparent values
	popmatrix[,"D"]<-0.5*(popmatrix[,"D"]+popmatrix[m,"D"])
	popmatrix[,"H"]<-0.5*(popmatrix[,"H"]+popmatrix[m,"H"])
	
	inds<-1:length(offs) #vector indexes for offspring
	inds<-rep(inds, times=offs)
	popmatrix<-popmatrix[inds,] #replace parents with offspring. Offspring inherit location and P from parents
	if (nrow(popmatrix)==1) return(NULL) #cases where only one individual left
	SDD<-SDD[inds]
	SDH<-SDH[inds]
	

	
	#Add segregation and environmental variances
	popmatrix[,"D"]<-rnorm(nrow(popmatrix), popmatrix[,"D"], sqrt(0.5*SDD^2))
	popmatrix[,"DP"]<-popmatrix[,"D"]+rnorm(nrow(popmatrix), 0, SDVED)
	popmatrix[,"H"]<-rnorm(nrow(popmatrix), popmatrix[,"H"], sqrt(0.5*SDH^2))
	popmatrix[,"HP"]<-popmatrix[,"H"]+rnorm(nrow(popmatrix), 0, SDVEH)
	
	# Disperse offspring
	disp<-rnorm(nrow(popmatrix), mean=popmatrix[,"X"], sd=exp(popmatrix[,"D"])) #disperse offspring
  	#surv<-rbinom(n=length(disp), size=1, prob=d.surv(abs(disp-popmatrix[,"X"])))
  	popmatrix[,"X"]<-disp
	popmatrix
}

# calculates each individual's fitness
w<-function(popmatrix, a, R0){
	mets<-ca.metrics(popmatrix, evovar=FALSE)
	hab<-h.surv(1, popmatrix[,"HP"])
	hab*bev.holt(mets[,"Nx"], R0, a)
}


# plots population density and mean trait values through space
plotter.mean<-function(popmatrix, a, R0, filename=NULL, H.init=0, Hsd.init=0.06, D.init=log(4), Dsd.init=0.06, ...){
	popmatrix<-popmatrix[order(popmatrix[,"X"]),]
	summ<-sum.metrics(popmatrix, TRUE, ...)
	summ[,'SDH']<-summ[,'SDH']^2
	summ[,'SDD']<-summ[,'SDD']^2
  	if (!is.null(filename)) pdf(file=filename)
	par(mfrow=c(3, 2), mar=c(2,5,1,1), cex.lab=1.4, oma=c(5, 0, 0, 0))
	
	tmp.yl<-c(min(summ[,"MeanH"], H.init), max(summ[,"MeanH"], H.init))
  	plot(summ[,"X"], summ[,"MeanH"], xlab="", ylab=quote(Mean~italic(z[w])~value), col="darkorange", type="l", ylim=tmp.yl)
  	arrows(x0=0, y0=tmp.yl[1], y1=tmp.yl[2], col='grey50', length=0)
  	arrows(x0=min(summ[,"X"]), y0=H.init, x1=max(summ[,"X"]), col='grey50', length=0)  	
  	
  	tmp.yl<-c(min(summ[,"SDH"], Hsd.init), max(summ[,"SDH"], Hsd.init))
  	plot(summ[,"X"], summ[,"SDH"], xlab="", ylab=quote(Genetic~variance~of~italic(z[w])), col="darkorange", type="l", ylim=tmp.yl)
	arrows(x0=0, y0=tmp.yl[1], y1=tmp.yl[2], col='grey50', length=0)
  	arrows(x0=min(summ[,"X"]), y0=Hsd.init, x1=max(summ[,"X"]), col='grey50', length=0)  	  	
	
	tmp.yl<-c(min(summ[,"MeanD"], D.init), max(summ[,"MeanD"], D.init))
	plot(summ[,"X"], summ[,"MeanD"], xlab="", ylab=quote(Mean~italic(z[d])~value), col="darkgreen", type="l", ylim=tmp.yl)
  	arrows(x0=0, y0=tmp.yl[1], y1=tmp.yl[2], col='grey50', length=0)
  	arrows(x0=min(summ[,"X"]), y0=D.init, x1=max(summ[,"X"]), col='grey50', length=0)  	
	
	tmp.yl<-c(min(summ[,"SDD"], Dsd.init), max(summ[,"SDD"], Dsd.init))
	plot(summ[,"X"], summ[,"SDD"], xlab="", ylab=quote(Genetic~variance~of~italic(z[d])), col="darkgreen", type="l", ylim=tmp.yl)
	arrows(x0=0, y0=tmp.yl[1], y1=tmp.yl[2], col='grey50', length=0)
  	arrows(x0=min(summ[,"X"]), y0=Dsd.init, x1=max(summ[,"X"]), col='grey50', length=0)
	
	plot(summ[,"X"], summ[,"Nx"], xlab="", ylab="Population density", type="l")
	arrows(x0=0, y0=min(summ[,"Nx"]), y1=max(summ[,"Nx"]), col='grey50', length=0)
	
	wvec<-w(popmatrix, a, R0)
	plot(popmatrix[,"X"], wvec, xlab="", ylab="Fitness", type="l")
	arrows(x0=0, y0=min(wvec), y1=max(wvec), col='grey50', length=0)
	
	mtext(text="Location", outer=TRUE, side=1, line=3)
	if (!is.null(filename)) dev.off()
}



	
# iterates pop over ngens and collects data
mother<-function(n, spX, a, R0, ngens, Hmean, Dmean, h2H, h2D, VPH, VPD, evovar, Allee.fx=FALSE, Plot){
  pop<-init.inds(n=n, spX=spX, Hmean=Hmean, Dmean=Dmean,
	h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD)
  for(ii in 1:ngens){
  	if (Plot==TRUE) plotter.mean(pop, a, R0)
    pop<-repro.disp(popmatrix=pop, a=a, R0=R0, 
    	h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD, evovar=evovar, Allee.fx=Allee.fx)
    if (is.null(pop)) break
  }
  # Get summary stats
  xlim<-range(pop[,"X"])
  mets<-metrics(pop, evovar)
  sset<-pop[,"X"] %in% xlim
  out<-cbind(xlim, mets[sset,])
  out<-list(parameters=list(ngens=ii, a=a, R0=R0, n=n, spX=spX,
            ngens=ngens, Hmean=Hmean, Dmean=Dmean, h2H=h2H, h2D=h2D, VPH=VPH,
            VPD=VPD, evovar=evovar), out=out, pop=pop)
  out
}

# iterates pop over ngens and collects data
# takes a popmatrix as an argument, so same (or different) pop can be used each time
mother.fed<-function(popmatrix, n, spX, a, R0, ngens, Hmean, Dmean, h2H, h2D, VPH, VPD, evovar, Allee.fx=FALSE, Plot){
  pop<-popmatrix
  for(ii in 1:ngens){
  	if (Plot==TRUE) plotter.mean(pop, a, R0)
    pop<-repro.disp(popmatrix=pop, a=a, R0=R0, 
    	h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD, evovar=evovar, Allee.fx=Allee.fx)
    if (is.null(pop)) return(NULL)
  }
  # Get summary stats
  xlim<-range(pop[,"X"])
  mets<-metrics(pop, evovar)
  sset<-pop[,"X"] %in% xlim
  out<-cbind(pop[sset,"X"], mets[sset,])
  out<-list(parameters=list(ngens=ii, a=a, R0=R0, n=n, spX=spX,
            ngens=ngens, Hmean=Hmean, Dmean=Dmean, h2H=h2H, h2D=h2D, VPH=VPH,
            VPD=VPD, evovar=evovar), out=out, pop=pop)
  out
}

#from the web, for plotting
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

############################## Functions for 2D Version ###########################

# Initialises a matrix with n individuals
# individuals have a location (spX, spY) and are 
# initialised with phenotype for D and H
init.inds2D<-function(n, spX, Hmean, Dmean, h2H, h2D, VPH, VPD, spY=spX){
	X<-runif(n, -spX, spX) #position
	Y<-runif(n, -spY, spY) #position
	D<-rnorm(n, Dmean, sqrt(h2D*VPD)) # individual sd
	H<-rnorm(n, Hmean, sqrt(h2H*VPH)) #Habitat match
	DP<-D+rnorm(n, 0, sqrt((1-h2D)*VPD))
	HP<-H+rnorm(n, 0, sqrt((1-h2H)*VPH))
	cbind(X, Y, D, H, DP, HP)		
}

# gets density, mean and (if evovar==TRUE) variance in breeding values, 2D case
metrics2D<-function(popmatrix, evovar=TRUE){
	if (is.null(popmatrix)) return(NULL)
	out<-.Call("metrics2D", R_X=popmatrix[,"X"], R_Y=popmatrix[,"Y"], R_H=popmatrix[,"H"],
	R_D=popmatrix[,"D"], R_n=nrow(popmatrix), R_ev=as.numeric(evovar))
	colnames(out)<-c("Nx", "MeanD", "MeanH", "SDD", "SDH")
	if (evovar==FALSE) out<-out[,1:3]
	out
}

#Finds a mate for each individual
mate2D<-function(popmatrix, dens){
	if (is.null(popmatrix)) return(NULL)
	out<-.Call("mate2D", R_X=popmatrix[,"X"], R_Y=popmatrix[,"Y"], R_dens=dens, R_n=nrow(popmatrix),
		R_bw=1)
	out[out==-99]<-NA # catch lonely ones
	out
}

# reproduces individuals based on density, kills parents
# kills offspring based on habitat
# kills dispersing offspring based on dispersal cost
# disperses individuals based on inherited P
repro.disp2D<-function(popmatrix, a, R0, h2H, h2D, VPH, VPD, evovar=TRUE, Allee.fx=FALSE){
	if (is.null(popmatrix)) return(NULL)
	if (nrow(popmatrix)==1) return(NULL) #cases where only one individual left
	SDVED<-sqrt(VE(h2D, VPD)) # define environmental variances
	SDVEH<-sqrt(VE(h2H, VPH))
	# calculate traits through space
	mets<-metrics2D(popmatrix, evovar)
	#print('mets OK')
	# reproduction
	hab<-h.surv(1, popmatrix[,"HP"])
	if (Allee.fx) {hab<-hab*Allee(dens=mets[,"Nx"])}
	offs<-rpois(nrow(popmatrix), hab*bev.holt(mets[,"Nx"], R0, a)) # surviving offspring
	if (sum(offs)<=1) return(NULL) #catch extinction
	m<-mate2D(popmatrix, mets[,"Nx"]) #find a mate
	# calculate new G and P
	if (evovar==FALSE){
		SDD<-sqrt(h2D*VPD)
		SDH<-sqrt(h2H*VPH)
	} else {
		SDD<-mets[,"SDD"]
		SDH<-mets[,"SDH"]
	}
	#get midparent values
	popmatrix[,"D"]<-0.5*(popmatrix[,"D"]+popmatrix[m,"D"])
	popmatrix[,"H"]<-0.5*(popmatrix[,"H"]+popmatrix[m,"H"])
	
	inds<-1:length(offs) #vector indexes for offspring
	inds<-rep(inds, times=offs)
	popmatrix<-popmatrix[inds,] #replace parents with offspring. Offspring inherit location and P from parents
	if (nrow(popmatrix)==1) return(NULL) #cases where only one individual left
	SDD<-SDD[inds]
	SDH<-SDH[inds]
	

	
	#Add segregation and environmental variances
	popmatrix[,"D"]<-rnorm(nrow(popmatrix), popmatrix[,"D"], sqrt(0.5*SDD^2))
	popmatrix[,"DP"]<-popmatrix[,"D"]+rnorm(nrow(popmatrix), 0, SDVED)
	popmatrix[,"H"]<-rnorm(nrow(popmatrix), popmatrix[,"H"], sqrt(0.5*SDH^2))
	popmatrix[,"HP"]<-popmatrix[,"H"]+rnorm(nrow(popmatrix), 0, SDVEH)
	
	# Disperse offspring
	disp<-rnorm(nrow(popmatrix), mean=0, sd=exp(popmatrix[,"D"])) #disperse offspring
	dirn<-runif(nrow(popmatrix), 0, 2*pi)
	popmatrix[,"X"]<-popmatrix[,"X"]+disp*cos(dirn)
	popmatrix[,"Y"]<-popmatrix[,"Y"]+disp*sin(dirn)
	popmatrix
}

# iterates pop over ngens and collects data
mother2D<-function(n, spX, a, R0, ngens, Hmean, Dmean, h2H, h2D, VPH, VPD, evovar, Allee.fx=FALSE, Plot){
  pop<-init.inds2D(n=n, spX=spX, Hmean=Hmean, Dmean=Dmean,
	h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD)
  for(ii in 1:ngens){
  	if (Plot==TRUE) plotter.mean(pop, a, R0)
  	cat("Calculating generation ", ii, " of ", ngens, "\n")
    pop<-repro.disp2D(popmatrix=pop, a=a, R0=R0, 
    	h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD, evovar=evovar, Allee.fx=Allee.fx)
    if (is.null(pop)) break
  }
  # Get summary stats
  out<-list(parameters=list(ngens=ii, a=a, R0=R0, n=n, spX=spX,
            ngens=ngens, Hmean=Hmean, Dmean=Dmean, h2H=h2H, h2D=h2D, VPH=VPH,
            VPD=VPD, evovar=evovar), pop=pop)
  out
}

# iterates pop over ngens and collects data
# takes a popmatrix as an argument, so same (or different) pop can be used each time
mother.fed2D<-function(popmatrix, n, spX, a, R0, ngens, Hmean, Dmean, h2H, h2D, VPH, VPD, evovar, Allee.fx=FALSE, Plot){
  pop<-popmatrix
  for(ii in 1:ngens){
  	if (Plot==TRUE) plotter.mean(pop, a, R0)
  	cat("Calculating generation ", ii, " of ", ngens, "\n")
    pop<-repro.disp2D(popmatrix=pop, a=a, R0=R0, 
    	h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD, evovar=evovar, Allee.fx=Allee.fx)
    if (is.null(pop)) return(NULL)
  }
  # Get summary stats
  xlim<-range(pop[(pop[,"Y"]<1 & pop[,"Y"]>-1),"X"])
  ylim<-range(pop[(pop[,"X"]<1 & pop[,"X"]>-1),"Y"])
  #mets<-metrics2D(pop, evovar)
  #sset.x<-pop[,"X"] %in% xlim
  #sset.y<-pop[,"Y"] %in% xlim
  #out<-cbind(pop[sset,"X"], mets[sset,])
  out<-cbind(xlim, ylim)
  
  out<-list(parameters=list(ngens=ii, a=a, R0=R0, n=n, spX=spX,
            ngens=ngens, Hmean=Hmean, Dmean=Dmean, h2H=h2H, h2D=h2D, VPH=VPH,
            VPD=VPD, evovar=evovar), out=out, pop=pop)
  out
}

# Plots population as XY plot with heat colours denoting dispersal breeding values
pop.plot2D<-function(popmatrix, trait, circ=TRUE, ...){
	clrs<-switch(trait,
		D=rev(heat.colors(50)),
		H=topo.colors(50)
	)
	cpal<-switch(trait,
		D=clrs[cut(popmatrix[,trait], 50)],
		H=clrs[cut(popmatrix[,trait], 50)]
	)
	radius<-mean(c(abs(quantile(popmatrix[,"X"], probs=c(0.005, 0.995))), 
						abs(quantile(popmatrix[,"Y"], probs=c(0.005, 0.995)))))
	#xylim<-range(c(popmatrix[,"X"], popmatrix[,"Y"]))
	plot(popmatrix[,"X"], popmatrix[,"Y"],
		col=cpal,
		bty="l",
		#xlim=xylim,
		#ylim=xylim, 
		...)
	if (circ==TRUE){
		symbols(0, 0, circles=radius, 
			inches=FALSE, 
			add=TRUE)
	}
}
