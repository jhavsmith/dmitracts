
library(zipfR)

make.pulse.mig = function(rate,tstart) {
    mig = matrix(nrow=tstart+2,ncol=2)
    mig[1:length(mig)] = 0
    mig[nrow(mig),] = c(1-rate,rate)
    return(mig)
}

make.continuous.mig = function(rate,tstart) {
    mig = matrix(nrow=tstart+2,ncol=2)
    mig[,1] = c(0,0,rep(0,tstart))
    mig[,2] = c(0,0,rep(rate,tstart))
    mig[nrow(mig),] = c(1,0)
    return(mig)
}

make.Q.matrix = function(mig, states) {
    state.labs = apply(states,1,paste,collapse='.')
    index.grid = make.index.grid(state.labs)
    t.matrix = matrix(nrow=length(state.labs),ncol=length(state.labs))
    index.cells = expand.grid(c(1:length(state.labs)),c(1:length(state.labs)))
    dimnames(t.matrix) = list(state.labs,state.labs)
    tau.list = make.tau.dict(mig,states)
    tau.vec = tau.list$tau.vec
    tau.keys = tau.list$tau.keys
    pb = txtProgressBar(1,nrow(index.grid),1,style=3)
    for (i in 1:nrow(index.grid)) {
        tau.pop = which(tau.keys[,2]==index.grid[i,4])
        tau.t =  which(tau.keys[,1]==index.grid[i,3])
        tau.sum = tau.vec[intersect(tau.pop,tau.t)]
        transition.rate = sum(tau.sum[1:(min(index.grid[i,c(1,3)])-2)])
        t.matrix[index.cells[i,1],index.cells[i,2]] = transition.rate
        setTxtProgressBar(pb, i)
    }
    diag(t.matrix) = 0
    diag(t.matrix) = -rowSums(t.matrix)
    close(pb)
    return(t.matrix)
}

make.index.grid = function(states) {
    state.grid = expand.grid(unlist(states),unlist(states))
    index.grid = matrix(nrow=nrow(state.grid),ncol=4)
    row = paste(unlist(strsplit(paste(unlist(state.grid[,1])),split='.',fixed=TRUE)))
    col = paste(unlist(strsplit(paste(unlist(state.grid[,2])),split='.',fixed=TRUE)))
    index.grid[,1] = as.numeric(row[seq(1,length(row),2)])
    index.grid[,2] = as.numeric(row[seq(2,length(row),2)])
    index.grid[,3] = as.numeric(col[seq(1,length(col),2)])
    index.grid[,4] = as.numeric(col[seq(2,length(col),2)])
    return(index.grid)
}



make.tau.dict = function(mig, states) {
    totmig = mig[,1]+mig[,2]
    tau.vec = c(1:(sum(states[,1])-nrow(states)))
    tau.keys = matrix(nrow=length(tau.vec),ncol=3)
    index = 1
    t.vec = states[,1]
    pop.vec = states[,2]
    for (i in 1:length(t.vec)) {
        temp.vec = c(1:(t.vec[i]-1))
        t.index = 1
        for (tau in temp.vec) {    #MAKE THIS AN EMPTY ARRAY = 1
            if ((t.vec[i]-1) < (tau+1)) {
                prod = 1
            } else {
                prod = prod((1-totmig)[(tau+1):(t.vec[i]-1)])
            }
            temp.vec[t.index] = totmig[t.vec[i]]*prod
            t.index = t.index+1
            index = index+1
        }
        tau.keys[c((index-length(temp.vec)):(index-1)),1] = states[i,1]
        tau.keys[c((index-length(temp.vec)):(index-1)),2] = states[i,2]
        tau.keys[c((index-length(temp.vec)):(index-1)),3] = c(1:length(temp.vec))
        tau.vec[(index-length(temp.vec)):(index-1)] = temp.vec
    }
    tau.list = list("tau.vec"=tau.vec,"tau.keys"=tau.keys)
    return(tau.list)
}

make.P.matrix = function(Q.matrix) {
    q = max(-diag(Q.matrix))
    P.matrix = Q.matrix/q
    diag(P.matrix) = 1-(-diag(Q.matrix)/q)
    return(P.matrix)
}

get.stationary.vec = function(P.matrix) {
    eigen_vect = eigen(t(P.matrix))$vectors[,1]
    stationary.vec = eigen_vect/sum(eigen_vect)
    return(stationary.vec)
}

get.bn.vec = function(P.matrix, stationary.vec, states, maxrate) {
    recover()
    init.probs = stationary.vec
    init.probs[which(states[,2]==1)] = 0
    new = init.probs %*% P.matrix
    newrest = new[which(states[,2]==1)]
    newrest = newrest/sum(newrest)
    shortmat = t(P.matrix[which(states[,2]==1),which(states[,2]==1)])
    escapes = 1 - sum(shortmat)
    maxmorgans = 100
    nit = as.integer(maxmorgans*maxrate)
    bn.vec = rep(0,nit)
    for (i in 1:nit) {
        escapes
        prob.diffs = probs.b - probs.a
        bn.vec[i] = prob.diffs[length(prob.diffs)]
        probs.a = probs.b
    }
    if (length(which(bn.vec==0))==0) {
        print("bn.vec needs to be longer")
    }
    if ((length(which(bn.vec==0))>0)) {
        bn.vec = bn.vec[1:max(which(bn.vec<1e-16))]
    }
    return(bn.vec)
}

get.Z = function(bn.vec, big.t, L) {
    Lambda = length(bn.vec)
    k.vec = c(1:Lambda)
    Z = L + sum((bn.vec*k.vec)/big.t)
    return(Z)
}

get.phi.x = function(Q.zero, bn.vec, x) {
    Lambda = length(bn.vec)
    erlang.k = dgamma(x,rate=1/Q.zero,shape=c(1:Lambda))
    phi.vec = bn.vec*erlang.k
    phi.x = sum(phi.vec)
    return(phi.x)
}

get.phi.x.i = function(bn.vec, x, big.t, L, Z) {
    Lambda = length(bn.vec)
    erlang.k = dgamma(x,rate=1/big.t,shape=c(1:Lambda))
    phi.vec = bn.vec*erlang.k
    phi.x = ((L-x)/Z) * sum(phi.vec)
    return(phi.x)
}

get.phi.x.e = function(bn.vec, x, big.t, Z) {
    Lambda = length(bn.vec)
    k.vec = c(1:Lambda)
    igamma = exp(Igamma(k.vec,(big.t*x),lower=FALSE,log=TRUE) - lfactorial(k.vec-1))
    phi.x = (2/Z)*sum(bn.vec * igamma)
    return(phi.x)
}

get.phi.x.f = function(bn.vec, x, big.t, L, Z) {
    Lambda = length(bn.vec)
    k.vec = c(1:Lambda)
    igamma.a = exp(log(L) + Igamma(k.vec,(big.t*L),lower=FALSE,log=TRUE) - lfactorial(k.vec-1))
    igamma.b = exp(Igamma((k.vec+1),(big.t*L),lower=FALSE,log=TRUE) - lfactorial(k.vec-1))
    igamma.sum = igamma.a + igamma.b
    if ((L-x)==0) {
        dirac = 1
    } else {
        dirac = 0
    }
    phi.x = (dirac/Z) * sum(bn.vec * igamma.sum)
    return(phi.x)
}

get.tract.length.distr = function(Q.zero, bn.vec, x.vec) {
    tract.length.distr = c(1:length(x.vec))
    for (i in 1:length(x.vec)) {
        tract.length.distr[i] = get.phi.x(Q.zero, bn.vec, x.vec[i])
    }
    return(tract.length.distr)
}

get.tract.length.distr.with.edges = function(Q.zero, bn.vec, bins, mig) {
    big.t = nrow(mig)
    L = bins
    Z = get.Z(bn.vec, big.t, L)
    x.vec = seq(1,L,1)
    tract.length.distr.i = c(1:L)
    tract.length.distr.e = c(1:L)
    tract.length.distr.f = c(1:L)
    tract.length.distr.g = c(1:L)
    for (i in 1:length(x.vec)) {
        tract.length.distr.i[i] = get.phi.x.i(bn.vec, x.vec[i], big.t, L, Z)
        tract.length.distr.e[i] = get.phi.x.e(bn.vec, x.vec[i], big.t, Z)
        tract.length.distr.f[i] = get.phi.x.f(bn.vec, x.vec[i], big.t, L, Z)
        tract.length.distr.g[i] = get.phi.x(Q.zero, bn.vec, x.vec[i])
    }
    tract.length.distr.list = list("tract.length.distr.i"=tract.length.distr.i,
                                   "tract.length.distr.e"=tract.length.distr.e,
                                   "tract.length.distr.f"=tract.length.distr.f,
                                   "tract.length.distr.g"=tract.length.distr.g)
    return(tract.length.distr.list)
}

run.tract.length.functions = function(mig, bins, initial.probs) {
    Q.matrix = make.Q.matrix(mig)
    maxrate = max(-diag(Q.matrix))
    Q.zero = 1/maxrate
    P.matrix = make.P.matrix(Q.matrix)
    if (initial.probs[1]=="stationary") {
        stationary.vec = get.stationary.vec(P.matrix)
    }
    if (initial.probs[1]!="stationary") {
        stationary.vec = initial.probs
    }
    bn.vec = get.bn.vec(P.matrix,stationary.vec)
    tract.length.distr = get.tract.length.distr(Q.zero, bn.vec, bins)
    return(tract.length.distr)
}

run.tract.length.functions.with.edges = function(mig, bins, initial.probs) {
    states = which(mig!=0,arr.ind = T)
    states = states[order(states[,1]),]
    Q.matrix = make.Q.matrix(mig,states)
    maxrate = max(-diag(Q.matrix))
    Q.zero = 1/maxrate
    P.matrix = make.P.matrix(Q.matrix)
    if (initial.probs[1]=="stationary") {
        stationary.vec = get.stationary.vec(P.matrix)
    }
    if (initial.probs[1]!="stationary") {
        stationary.vec = initial.probs
    }
    bn.vec = get.bn.vec(P.matrix, stationary.vec, states, maxrate)
    tract.length.distr.list = get.tract.length.distr.with.edges(Q.zero, bn.vec, bins, mig)
    return(tract.length.distr.list)
}

run.tract.length.functions.with.selection = function(mig, bins, initial.probs) {
    Q.matrix = make.Q.matrix(mig)
    Q.zero = 1/max(-diag(Q.matrix))
    P.matrix = make.P.matrix(Q.matrix)
    if (initial.probs[1]=="stationary") {
        stationary.vec = get.stationary.vec(P.matrix)
    }
    if (initial.probs[1]!="stationary") {
        stationary.vec = initial.probs
    }
    bn.vec = get.bn.vec(P.matrix,stationary.vec)
    tract.length.distr = get.tract.length.distr.with.selection(Q.zero, bn.vec, bins)
    return(tract.length.distr)
}
