



#-----simulation parameters-------
popsize = 2000
chrlen = 1000000
rate = .01
tstart = 50
dmipos = c(.3,.7)
seldmi = .2
selanc = 0
domdmi = 0
#---------------------------------
recomb = (.7-.3)*chrlen*1e-8
starting.freqs = rev(c(0,.5,.5,0))
source.freqs = c(1,0,0,0)
reps = 20
#AB Ab aB ab

#------------Frequency Trajectories Hybrid Swarm--------------
source("gamete_trajectories_hs.R")
make.fitnessmatrix(seldmi,selanc,domdmi)
WF.det.list = simulate.trajecs.hs(popsize,rate,recomb,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs)
WF.sto.list = WF.trajec.reps.hs(popsize,rate,recomb,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs,reps)
plot.WF.versus.deterministic.hs(WF.det.list,WF.sto.list)
