#---------------Gamete frequency trajectory simulations-----------------

popsize = 100000      # number of individuals per deme
m = 0.01         # migration rate between species
msub = 0.01      # migration rate among demes within species
t = 1000        # number of generations. Also initial state for t.
t.prime = 110   # next state for t
reps = 20       # number of WF simulations
w = 100000      # basepairs from locus A to the junction
v = 100000      # basepairs from the junction to locus B
r = 1e-8*(w+v)  # recombination rate between A and B
beta = .1       # fitness matrix parameters
beta1 = .1
s = .1
alpha = 0
alpha1 = 0
starting.freqs = rev(c(0,.5,.5,0)) # starting freqs for the hybrid swarm.
source.freqs = c(1,0,0,0)          # infinite source population gamete frequencies.
#AB Ab aB ab

#------------Frequency Trajectories Hybrid Swarm--------------
source("gamete_trajectories_hs.R")
make.fitnessmatrix(alpha,alpha1,beta,beta1,s)
WF.det.list = simulate.trajecs.hs(popsize,m,r,t,alpha,alpha1,beta,beta1,s,source.freqs,starting.freqs)
WF.sto.list = WF.trajec.reps.hs(popsize,m,r,t,alpha,alpha1,beta,beta1,s,source.freqs,starting.freqs,reps)
plot.WF.versus.deterministic.hs(WF.det.list,WF.sto.list)
