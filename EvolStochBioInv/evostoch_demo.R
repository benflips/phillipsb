# Demonstration of a run of the evolutionary stochasticity model
# Author: Ben Phillips

# This runs a number of replicates with and without evolution and saves them to a file

# note that plotting functions are also found in the functions source file

# to run in the wd in which the functions file is located.

source("evostochFunctions.R")



# parameter values
spX<-0
R0<-20
a<-(R0-1)/50 #Nstar=50
n<-20
ngens<-30
Hmean<-0
Dmean<-log(4)
h2H<-0.3
h2D<-h2H
VPH<-0.2
VPD<-0.2
nreps<-1


out<-vector(mode="list", length=2*nreps)
pop<-init.inds(n, spX, Hmean, Dmean, h2H, h2D, VPH, VPD)
for (ii in 1:nreps){
	print(ii)
	repl<-mother.fed(pop, n=n, spX=spX, a=a, R0=R0, ngens=ngens, Hmean=Hmean, 
		Dmean=Dmean, h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD, evovar=TRUE, Plot=FALSE)
	out[[ii]]<-repl
}
#Turn evolution off
h2H<-0
h2D<-h2H
pop[,"H"]<-Hmean
pop[,"D"]<-Dmean

for (ii in (nreps+1):(2*nreps)){
	print(ii)
	repl<-mother.fed(pop, n=n, spX=spX, a=a, R0=R0, ngens=ngens, Hmean=Hmean, 
		Dmean=Dmean, h2H=h2H, h2D=h2D, VPH=VPH, VPD=VPD, evovar=TRUE, Plot=FALSE)
	out[[ii]]<-repl
}

  	
save(out, file="evostoch_varR0.RData")
