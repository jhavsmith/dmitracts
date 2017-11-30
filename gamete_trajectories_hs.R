
make.population.hs = function(starting.freqs) {
    gamete.x = starting.freqs
    AB.x = gamete.x[1]; Ab.x = gamete.x[2]
    aB.x = gamete.x[3]; ab.x = gamete.x[4]
    zygote.x = matrix(nrow=3,ncol=3)
    zygote.x[1,1] = AB.x*AB.x                 #AABB
    zygote.x[1,2] = 2*AB.x*Ab.x               #AABb
    zygote.x[1,3] = Ab.x*Ab.x                 #AAbb
    zygote.x[2,1] = 2*AB.x*aB.x               #AaBB
    zygote.x[2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
    zygote.x[2,3] = 2*Ab.x*ab.x               #Aabb
    zygote.x[3,1] = aB.x*aB.x                 #aaBB
    zygote.x[3,2] = 2*aB.x*ab.x               #aaBb
    zygote.x[3,3] = ab.x*ab.x                 #aabb
    pop.list = list("zygote.x"=zygote.x,"gamete.x"=gamete.x)
    return(pop.list)
}

make.fitnessmatrix = function(seldmi,selanc,seldom) {
    if (domdmi==0) {
        h0=1
        h1=1
    }
    if (domdmi==1) {
        h0=1/4
        h1=1/2
    }
    if (domdmi==2) {
        h0=0
        h1=1/2
    }
    if (domdmi==3) {
        h0=0
        h1=0
    }
    zygote.w = matrix(nrow=3,ncol=3)
    zygote.w[1,1] = 1 - seldmi      #AABB
    zygote.w[1,2] = 1 - seldmi*h1   #AABb
    zygote.w[1,3] = 1                #AAbb
    zygote.w[2,1] = 1 - seldmi*h1   #AaBB
    zygote.w[2,2] = 1 - seldmi*h0   #AaBb
    zygote.w[2,3] = 1 - (selanc/2)  #Aabb
    zygote.w[3,1] = 1               #aaBB
    zygote.w[3,2] = 1 - (selanc/2)  #aaBb
    zygote.w[3,3] = 1 - selanc      #aabb
    return(zygote.w)
}





recombination.and.selection.hs = function(pop.list,zygote.w,r) {
    #getting the zygote fitnesses
    AABB.w = zygote.w[1,1]; AABb.w = zygote.w[1,2]; AAbb.w = zygote.w[1,3]
    AaBB.w = zygote.w[2,1]; AaBb.w = zygote.w[2,2]; Aabb.w = zygote.w[2,3]
    aaBB.w = zygote.w[3,1]; aaBb.w = zygote.w[3,2]; aabb.w = zygote.w[3,3]
    zygote.x = pop.list$zygote.x
    gamete.x = pop.list$gamete.x
    AB.x = gamete.x[1]; Ab.x = gamete.x[2]
    aB.x = gamete.x[3]; ab.x = gamete.x[4]
    #marginalizing the gamete fitnesses
    AB.w = (AABB.w*AB.x)+(AABb.w*Ab.x)+(AaBB.w*aB.x)+(AaBb.w*ab.x)
    Ab.w = (AABb.w*AB.x)+(AAbb.w*Ab.x)+(AaBb.w*aB.x)+(Aabb.w*ab.x)
    aB.w = (AaBB.w*AB.x)+(AaBb.w*Ab.x)+(aaBB.w*aB.x)+(aaBb.w*ab.x)
    ab.w = (AaBb.w*AB.x)+(Aabb.w*Ab.x)+(aaBb.w*aB.x)+(aabb.w*ab.x)
    gamete.w = c(AB.w, Ab.w, aB.w, ab.w)
    #mean fitness
    mean.w = sum(zygote.w*zygote.x)
    #LD
    D = AB.x*ab.x - Ab.x*aB.x
    #gamete frequencies after recombination and selection
    selection.a = (gamete.w[c(1,4)]*gamete.x[c(1,4)] - r*zygote.w[2,2]*D)/mean.w
    selection.b = (gamete.w[c(2,3)]*gamete.x[c(2,3)] + r*zygote.w[2,2]*D)/mean.w
    pop.list$gamete.x = c(selection.a[1],selection.b[1],selection.b[2],selection.a[2])
    return(pop.list)
}

migration.hs = function(pop.list,source.freqs,m) {
    gametes.x = pop.list$gamete.x
    gametes = (1-(2*m))*gametes.x + m*source.freqs + m*rev(source.freqs)
    pop.list$gamete.x = gametes
    return(pop.list)
}

make.zygotes.hs = function(pop.list,N) {
    zygote.x = pop.list$zygote.x
    gamete.x = pop.list$gamete.x
    AB.x = gamete.x[1]
    Ab.x = gamete.x[2]
    aB.x = gamete.x[3]
    ab.x = gamete.x[4]
    zygote.x[1,1] = AB.x*AB.x                 #AABB
    zygote.x[1,2] = 2*AB.x*Ab.x               #AABb
    zygote.x[1,3] = Ab.x*Ab.x                 #AAbb
    zygote.x[2,1] = 2*AB.x*aB.x               #AaBB
    zygote.x[2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
    zygote.x[2,3] = 2*Ab.x*ab.x               #Aabb
    zygote.x[3,1] = aB.x*aB.x                 #aaBB
    zygote.x[3,2] = 2*aB.x*ab.x               #aaBb
    zygote.x[3,3] = ab.x*ab.x                 #aabb
    pop.list$zygote.x = zygote.x
    return(pop.list)
}

make.zygotes.WF.hs = function(pop.list,N) {
    zygote.x = pop.list$zygote.x
    gamete.x = pop.list$gamete.x
    gamete.vec = rmultinom(1,(2*N),prob=gamete.x)/(2*N)
    AB.x = gamete.vec[1]
    Ab.x = gamete.vec[2]
    aB.x = gamete.vec[3]
    ab.x = gamete.vec[4]
    zygote.x[1,1] = AB.x*AB.x                 #AABB
    zygote.x[1,2] = 2*AB.x*Ab.x               #AABb
    zygote.x[1,3] = Ab.x*Ab.x                 #AAbb
    zygote.x[2,1] = 2*AB.x*aB.x               #AaBB
    zygote.x[2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
    zygote.x[2,3] = 2*Ab.x*ab.x               #Aabb
    zygote.x[3,1] = aB.x*aB.x                 #aaBB
    zygote.x[3,2] = 2*aB.x*ab.x               #aaBb
    zygote.x[3,3] = ab.x*ab.x                 #aabb
    pop.list$zygote.x = zygote.x
    pop.list$gamete.x = gamete.vec
    return(pop.list)
}

simulate.trajecs.hs = function(N,m,r,t,seldmi,selanc,seldom,source.freqs,starting.freqs) {
    pop.list = make.population.hs(starting.freqs)
    zygote.w = make.fitnessmatrix(seldmi,selanc,seldom)
    gamete.trajecs = matrix(nrow=(t+1),ncol=4)     # gamete freqs are ordered AB Ab aB ab
    gamete.trajecs[1,] = pop.list$gamete.x
    allele.trajecs = matrix(nrow=(t+1),ncol=2)     # allele freqs are for A and b
    allele.trajecs[1,] = c(sum(pop.list$gamete.x[c(1,2)]),sum(pop.list$gamete.x[c(2,4)]))
    for (j in 2:(t+1)) {
        pop.list = recombination.and.selection.hs(pop.list,zygote.w,r)
        pop.list = migration.hs(pop.list,source.freqs,m)
        pop.list = make.zygotes.hs(pop.list,N)
        gamete.trajecs[j,]  = pop.list$gamete.x
        allele.trajecs[j,1:2] = c(sum(pop.list$gamete.x[c(1,2)]),sum(pop.list$gamete.x[c(2,4)]))
    }
    trajec.list = list('AB.matrix'=gamete.trajecs[,1],'Ab.matrix'=gamete.trajecs[,2],
                       'aB.matrix'=gamete.trajecs[,3],'ab.matrix'=gamete.trajecs[,4],
                       'A.matrix'=allele.trajecs[,1],'b.matrix'=allele.trajecs[,2],
                       'zygote.w'=zygote.w)
    return(trajec.list)
}

make.zygotes.hs = function(pop.list,N) {
    zygote.x = pop.list$zygote.x
    gamete.x = pop.list$gamete.x
    AB.x = gamete.x[1]
    Ab.x = gamete.x[2]
    aB.x = gamete.x[3]
    ab.x = gamete.x[4]
    zygote.x[1,1] = AB.x*AB.x                 #AABB
    zygote.x[1,2] = 2*AB.x*Ab.x               #AABb
    zygote.x[1,3] = Ab.x*Ab.x                 #AAbb
    zygote.x[2,1] = 2*AB.x*aB.x               #AaBB
    zygote.x[2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
    zygote.x[2,3] = 2*Ab.x*ab.x               #Aabb
    zygote.x[3,1] = aB.x*aB.x                 #aaBB
    zygote.x[3,2] = 2*aB.x*ab.x               #aaBb
    zygote.x[3,3] = ab.x*ab.x                 #aabb
    pop.list$zygote.x = zygote.x
    pop.list$gamete.x = gamete.x
    return(pop.list)
}

simulate.trajecs.WF.hs = function(N,m,r,t,seldmi,selanc,seldom,source.freqs,starting.freqs,reps) {
    pop.list = make.population.hs(starting.freqs)
    zygote.w = make.fitnessmatrix(seldmi,selanc,seldom)
    gamete.trajecs = matrix(nrow=(t+1),ncol=4)     # gamete freqs are ordered AB Ab aB ab
    gamete.trajecs[1,] = pop.list$gamete.x
    allele.trajecs = matrix(nrow=(t+1),ncol=2)     # allele freqs are for A and b
    allele.trajecs[1,] = c(sum(pop.list$gamete.x[c(1,2)]),sum(pop.list$gamete.x[c(2,4)]))
    for (j in 2:(t+1)) {
        pop.list = recombination.and.selection.hs(pop.list,zygote.w,r)
        pop.list = migration.hs(pop.list,source.freqs,m)
        pop.list = make.zygotes.WF.hs(pop.list,N)
        gamete.trajecs[j,]  = pop.list$gamete.x
        allele.trajecs[j,1:2] = c(sum(pop.list$gamete.x[c(1,2)]),sum(pop.list$gamete.x[c(2,4)]))
    }
    trajec.list = list('AB.matrix'=gamete.trajecs[,1],'Ab.matrix'=gamete.trajecs[,2],
                       'aB.matrix'=gamete.trajecs[,3],'ab.matrix'=gamete.trajecs[,4],
                       'A.matrix'=allele.trajecs[,1],'b.matrix'=allele.trajecs[,2],
                       'zygote.w'=zygote.w)
    return(trajec.list)
}

WF.trajec.reps.hs = function(N,m,r,t,seldmi,selanc,seldom,source.freqs,starting.freqs,reps) {
    AB.matrix = matrix(nrow=t+1,ncol=reps)
    Ab.matrix = matrix(nrow=t+1,ncol=reps)
    aB.matrix = matrix(nrow=t+1,ncol=reps)
    ab.matrix = matrix(nrow=t+1,ncol=reps)
    A.matrix  = matrix(nrow=t+1,ncol=reps)
    b.matrix  = matrix(nrow=t+1,ncol=reps)
    pb = txtProgressBar(1,reps,1,style=3)
    for (b in 1:reps) {
        trajec.list.WF = simulate.trajecs.WF.hs(N,m,r,t,seldmi,selanc,seldom,source.freqs,starting.freqs,reps)
        AB.matrix[,b] = trajec.list.WF[[1]]
        Ab.matrix[,b] = trajec.list.WF[[2]]
        aB.matrix[,b] = trajec.list.WF[[3]]
        ab.matrix[,b] = trajec.list.WF[[4]]
        A.matrix[,b]  = trajec.list.WF[[5]]
        b.matrix[,b]  = trajec.list.WF[[6]]
        setTxtProgressBar(pb, b)
    }
    close(pb)
    zygote.w = trajec.list.WF[[7]]
    WF.rep.list = list('AB.matrix'=AB.matrix,'Ab.matrix'=Ab.matrix,
                       'aB.matrix'=aB.matrix,'ab.matrix'=ab.matrix,
                       'A.matrix' =A.matrix, 'b.matrix' =b.matrix,
                       'zygote.w'=zygote.w)
    return(WF.rep.list)
}

find.equilibrium.freqs.hs = function(pops,N,m,msub,r,seldmi,selanc,seldom,source.freqs) {
    pop.list = make.population.hs(starting.freqs)
    zygote.w = make.fitnessmatrix(seldmi,selanc,seldom)
    gamete.trajecs = matrix(nrow=(t+1),ncol=4)     # gamete freqs are ordered AB Ab aB ab
    gamete.trajecs[1,] = pop.list$gamete.x
    allele.trajecs = matrix(nrow=(t+1),ncol=2)     # allele freqs are for A and b
    allele.trajecs[1,] = c(sum(pop.list$gamete.x[c(1,2)]),sum(pop.list$gamete.x[c(2,4)]))
    for (j in 2:(t+1)) {
        pop.list = recombination.and.selection.hs(pop.list,zygote.w,r)
        pop.list = migration.hs(pop.list,source.freqs,m)
        pop.list = make.zygotes.WF.hs(pop.list,N)
        gamete.trajecs[j,]  = pop.list$gamete.x
        allele.trajecs[j,1:2] = c(sum(pop.list$gamete.x[c(1,2)]),sum(pop.list$gamete.x[c(2,4)]))
        frequency.diffs = abs(allele.trajecs[j,]-allele.trajecs[(j-1),])
        if (length(which(frequency.diffs<1e-5))==4) {
            break
        }
    }
    eq.freq.list = list('AB.eq.a'=gamete.trajecs[j,1],'Ab.eq.a'=gamete.trajecs[j,2],
                        'aB.eq.a'=gamete.trajecs[j,3],'ab.eq.a'=gamete.trajecs[j,4],
                        'A.eq.a'=allele.trajecs[j,1],'b.eq.a'=allele.trajecs[j,2],
                        'zygote.w'=zygote.w)
    return(eq.freq.list)
}

plot.WF.versus.deterministic.hs = function(WF.det.list,WF.sto.list) {
    library('RColorBrewer')
    colors = brewer.pal(4,'Paired')
    par(mar=c(5,5,4,2))
    plot(WF.sto.list[[1]][,1],col=colors[1],typ='l',main='Gamete Frequencies',
         lwd=2,ylim=c(0,1),xlab='Generation',ylab='Frequency',cex.lab=1.8,cex.main=1.8,cex.axis=1.8)
    lines(WF.sto.list[[2]][,1],col=colors[2],lwd=1)
    lines(WF.sto.list[[3]][,1],col=colors[3],lwd=1)
    lines(WF.sto.list[[4]][,1],col=colors[4],lwd=1)
    for (i in 2:reps) {
        lines(WF.sto.list[[1]][,i],col=colors[1],lwd=1)
        lines(WF.sto.list[[2]][,i],col=colors[2],lwd=1)
        lines(WF.sto.list[[3]][,i],col=colors[3],lwd=1)
        lines(WF.sto.list[[4]][,i],col=colors[4],lwd=1)
    }
    lines(WF.det.list[[1]],col=colors[1],lwd=3,lty=2)
    lines(WF.det.list[[2]],col=colors[2],lwd=3,lty=1)
    lines(WF.det.list[[3]],col=colors[3],lwd=3,lty=2)
    lines(WF.det.list[[4]],col=colors[4],lwd=3,lty=2)
    legend("topright",legend=c("AB","Ab","aB","ab"),col=c(colors),lwd=5,cex=1.8)
}
