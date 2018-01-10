
source("tracts_js.R")

rate = .01
tstart = 5
bins = 20
initial.probs = 'stationary'

mig = make.continuous.mig(rate,tstart)
tract.dists = run.tract.length.functions.with.edges(mig, bins, initial.probs)
sum.tract.dist = tract.dists[[1]]+tract.dists[[2]]+tract.dists[[3]]
rel.tract.dist = sum.tract.dist/sum(sum.tract.dist)
rel.tract.dist
