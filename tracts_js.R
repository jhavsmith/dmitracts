      
library(zipfR)
        
make.pulse.mig = function(rate,tstart) {
    mig = matrix(nrow=tstart+2,ncol=2)
    mig[1:length(mig)] = 0 
    mig[nrow(mig),] = c(1-rate,rate)
    return(mig)
}   

make.continuous.mig = function(rate,tstart) {
    mig = matrix(nrow=tstart+2,ncol=2)
    mig[,1] = c(0,0,rep(1-rate,tstart))
    mig[,2] = c(0,0,rep(rate,tstart))
    mig[nrow(mig),] = c(1,0)
    return(mig)
}            
                       
make.Q.matrix = function(mig) {
    m.vec.a = mig[,2] 
    m.vec.b = mig[,1]  
    t.matrix = matrix(nrow=nrow(mig)*2,ncol=nrow(mig)*2)
    dim.names = expand.grid(c(1:nrow(mig)),c(1:2))
    index.cells = expand.grid(c(1:(nrow(mig)*2)),c(1:(nrow(mig)*2)))
    index.grid = expand.grid(c(1:nrow(mig)),c(1:2),c(1:nrow(mig)),c(1:2))
    dim.names = apply(dim.names,1,paste,collapse=".")
    dimnames(t.matrix) = list(dim.names,dim.names)
    names(index.grid) = c("t.0","p.0","t.1","p.1")     
    pb = txtProgressBar(1,nrow(index.grid),1,style=3)                                                         
    for (i in 1:nrow(index.grid)) {   
        t.matrix[index.cells[i,1],index.cells[i,2]] = transition.rates(m.vec.a, m.vec.b, index.grid[i,1], index.grid[i,3], index.grid[i,2], index.grid[i,4])       
        setTxtProgressBar(pb, i)                                                             
    }                                      
    diag(t.matrix) = 0                     
    diag(t.matrix) = -rowSums(t.matrix)   
    close(pb)                               
    return(t.matrix)                        
}     
                                                                                                                                                                                                                                                                                                                                                               
transition.rates = function(m.vec.a, m.vec.b, t, t.prime, p, p.prime) {
    if (p==1) {
        if (p != p.prime) {
            m.vec = m.vec.a
        }
        if (p == p.prime) {
            m.vec = 1-m.vec.a
        } 
        no.m.vec = 1-m.vec               
    }   
    if (p==2) {
        if (p != p.prime) {
            m.vec = m.vec.b  
        }
        if (p == p.prime) {
            m.vec = 1-m.vec.b
        }    
        no.m.vec = 1-m.vec 
    }             
    tau.vec = c(1:(min(t, t.prime)-1))
    for (tau in 1:(min(t, t.prime)-1)) {  
         tau.vec[tau] = prod(no.m.vec[(tau+1):(t.prime-1)])
    }                        
    transition.rate = m.vec[t.prime]*sum(tau.vec) 
#    selection.transition.rate = 
#    transition.rate = neutral.transition.rate * selection.transition.rate
    return(transition.rate) 
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
                                                 
get.bn.vec = function(P.matrix, stationary.vec) {
    bn.vec = c(1:100000)
    init.probs = stationary.vec[1:(length(stationary.vec)/2)]
    sum.probs = sum(stationary.vec[((length(stationary.vec)/2)+1):length(stationary.vec)])
    probs.a = c(init.probs,sum.probs)   
    P.matrix.a = matrix(nrow=((nrow(P.matrix)/2)+1),ncol=((nrow(P.matrix)/2)+1))
    tr.range = c(1:(nrow(P.matrix)/2))
    P.matrix.a[tr.range,tr.range] = P.matrix[tr.range,tr.range] 
    P.matrix.a[c(1:(ncol(P.matrix.a)-1)),ncol(P.matrix.a)] = rowSums(P.matrix[c(1:(nrow(P.matrix)/2)),((ncol(P.matrix)/2)+1):ncol(P.matrix)])
    P.matrix.a[((nrow(P.matrix)/2)+1),] = c(rep(0,(ncol(P.matrix.a)-1)),1)
    for (i in 1:length(bn.vec)) {
        probs.b = probs.a %*% P.matrix.a
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
                                                                           
get.tract.length.distr.with.edges = function(Q.zero,bn.vec, x.vec, t.vec) {
    big.t = length(t.vec)
    L = length(x.vec)
    Z = get.Z(bn.vec, big.t, L)
    tract.length.distr.i = c(1:length(x.vec))
    tract.length.distr.e = c(1:length(x.vec))
    tract.length.distr.f = c(1:length(x.vec))
    tract.length.distr.g = c(1:length(x.vec))
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
    Q.zero = 1/max(-diag(Q.matrix))         
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
    tract.length.distr.list = get.tract.length.distr.with.edges(Q.zero, bn.vec, bins)
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




