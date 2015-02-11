#To run locally


source("evostochFunctions.R")

library(SDMTools)


R0<-4
Nstar<-10

spX<-0
R0<-R0
a<-(R0-1)/Nstar
n<-20
ngens<-30
Hmean<-0
Dmean<-log(4)
h2H<-0.3
h2D<-h2H
VPH<-0.2
VPD<-0.2
nreps<-10


pop<-init.inds2D(n, spX, Hmean, Dmean, h2H, h2D, VPH, VPD)

samp.freq<-5
genset<-ngens/samp.freq
out<-vector(mode="list", length=genset)

for (ii in 1:genset){
	out[[ii]]<-mother.fed2D(pop, n=n, spX=spX, a=a, R0=R0, ngens=samp.freq, Hmean=Hmean, 
		Dmean=Dmean, h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD, evovar=TRUE, Plot=FALSE)
	cat("Sampling...\n")
	pop<-out[[ii]]$pop
}
save(out, file="2DDemoRun.RData")


#### Plotting ####

plotset<-c(1, 3, 6) # which list elements to plot
pdf(file="progression.pdf", width=12, height=length(plotset)*6)
	par(mfrow=c(length(plotset), 2))
	xylim<-range(c(out[[genset]]$pop[,"X"], out[[genset]]$pop[,"Y"]))
	leg.pnts<-cbind(x=rep(c(xylim[1], xylim[1]*0.9), 2),
					y=rep(c(xylim[2]*0.95, xylim[2]*0.7), each=2))
	ct<-FALSE
	
	for (ii in plotset){
		if (ii==genset) ct<-TRUE
		pop.plot2D(out[[ii]]$pop, trait="D",
			circ=ct,
			xlim=xylim,
			ylim=xylim,
			xlab="X",
			ylab="Y",
			main=paste("Generation", ii*samp.freq))
		if (ii==genset){
			rect(xylim[1], -1, xylim[2], 1, col="grey70", border=NA)
			rect(-1, xylim[1], 1, xylim[2], col="grey70", border=NA)
		}
		
		legend.gradient(leg.pnts, 
			rev(heat.colors(50)), 
			limits=round(range(out[[ii]]$pop[,"D"]),2),
			title=quote(b[d]))
		
		pop.plot2D(out[[ii]]$pop, trait="H",
			circ=ct,
			xlim=xylim,
			ylim=xylim,
			xlab="X",
			ylab="Y",
			main=paste("Generation", ii*samp.freq))
		if (ii==genset){
			#rect(xylim[1], -1, xylim[2], 1, col="grey70", border=NA)
			#rect(-1, xylim[1], 1, xylim[2], col="grey70", border=NA)
		}
			
		legend.gradient(leg.pnts, 
			topo.colors(50), 
			limits=round(range(out[[ii]]$pop[,"H"]),2),
			title=quote(b[w]))
	}
dev.off()

temp<-out[[genset]]$pop
temp<-subset(temp, (temp[,"Y"]<1 & temp[,"Y"]>-1))
temp<-temp[,colnames(temp)!="Y"]
plotter.mean(temp, a=out[[genset]]$parameters$a, R0=out[[genset]]$parameters$R0,
	filename="samplePop2D.pdf", H.init=0, Hsd.init=0.06, D.init=log(4), Dsd.init=0.06)


