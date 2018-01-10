
source("gamete_trajectories_hs.R")
source("survival_probs.R")

popsize = 100
chrlen = 1000000
rate = .01
tstart = 20
seldmi = 0
selanc = 0
domdmi = 0
starting.freqs = rev(c(0,.5,.5,0))
source.freqs = c(1,0,0,0)
recrate = 1e-8
dmisep = .5
npts = 20
maxlen = 1
pop = 1
Ls = 1


#------------Survival Probabilities------------------------
compute.psivec = function(popsize,rate,v,dmisep,recrate,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs){
    recomb = dmisep*recrate
    WF.det.list = simulate.trajecs.hs(popsize,rate,recomb,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs)
    psivec = survival.probs(WF.det.list,v,w,tstart)
    write.table(psivec,"psivec.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
}

#-------------probability of ancestry until kth step--------
compute.bk = function() {
    system(paste("python get_tracts.py",rate,tstart,npts,maxlen,pop,Ls))


}



nDist = as.numeric(paste(unlist(read.table("nDist.txt"))))
