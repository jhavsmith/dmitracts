#-----simulation parameters-------
popsize = 100
chrlen = 1000000
rate = .1
tstart = 20
dmipos = c(.3,.7)
seldmi = 0
selanc = 0
domdmi = 0
#---------------------------------
recomb = (.7-.3)*chrlen*1e-8
starting.freqs = rev(c(0,.5,.5,0))
source.freqs = c(1,0,0,0)
reps = 20
#AB Ab aB ab

#------------Survival Probabilities--------------
source("gamete_trajectories_hs.R")
source("survival_probs.R")

#pdf("surv_probs_LD.pdf")
seldmi = 0
dmipos = c(.3,.7)
recomb = (dmipos[2]-dmipos[1])*chrlen*1e-8
WF.det.list = simulate.trajecs.hs(popsize,rate,recomb,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs)
WF.sto.list = WF.trajec.reps.hs(popsize,rate,recomb,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs,reps)
survival.time.matrix.WF = plot.survival.prob.versus.gens(WF.det.list,WF.sto.list,dmipos,chrlen,recomb)

seldmi = .3
dmipos = c(.3,.7)
recomb = (dmipos[2]-dmipos[1])*chrlen*1e-8
WF.det.list = simulate.trajecs.hs(popsize,rate,recomb,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs)
WF.sto.list = WF.trajec.reps.hs(popsize,rate,recomb,tstart,seldmi,selanc,seldom,source.freqs,starting.freqs,reps)
survival.time.matrix.WF = plot.survival.prob.versus.gens.again(WF.det.list,WF.sto.list,dmipos,chrlen,recomb)
#dev.off()
