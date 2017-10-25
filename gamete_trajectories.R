
make.populations = function(source.freqs,pops) {
    gamete.x.a = source.freqs
    AB.x = gamete.x.a[1]; Ab.x = gamete.x.a[2]
    aB.x = gamete.x.a[3]; ab.x = gamete.x.a[4]
    zygote.x.a = matrix(nrow=3,ncol=3)
    zygote.x.a[1,1] = AB.x*AB.x                 #AABB
    zygote.x.a[1,2] = 2*AB.x*Ab.x               #AABb
    zygote.x.a[1,3] = Ab.x*Ab.x                 #AAbb
    zygote.x.a[2,1] = 2*AB.x*aB.x               #AaBB
    zygote.x.a[2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
    zygote.x.a[2,3] = 2*Ab.x*ab.x               #Aabb
    zygote.x.a[3,1] = aB.x*aB.x                 #aaBB
    zygote.x.a[3,2] = 2*aB.x*ab.x               #aaBb
    zygote.x.a[3,3] = ab.x*ab.x                 #aabb
    zygote.x.b = apply(apply(zygote.x.a,2,rev),1,rev)
    gamete.x.b = rev(source.freqs)
    if (pops > 1) {
        zygote.x.a = rbind(zygote.x.a,matrix(nrow=3*(pops-1),ncol=3))
        zygote.x.b = rbind(zygote.x.b,matrix(nrow=3*(pops-1),ncol=3))
        gamete.x.a = rbind(gamete.x.a,matrix(nrow=(pops-1),ncol=4))
        gamete.x.b = rbind(gamete.x.b,matrix(nrow=(pops-1),ncol=4))
        zygote.index = seq(4,pops*3,3)
        zygote.index = cbind(zygote.index,seq(6,pops*3,3))
        for (i in 1:(pops-1)) {
            zygote.x.a[(zygote.index[i,1]:zygote.index[i,2]),] = zygote.x.a[(1:3),]
            zygote.x.b[(zygote.index[i,1]:zygote.index[i,2]),] = zygote.x.b[(1:3),]
            gamete.x.a[(i+1),] = gamete.x.a[1,]
            gamete.x.b[(i+1),] = gamete.x.b[1,]
        }
    }
    pop.list = list("zygote.x.a"=zygote.x.a,"zygote.x.b"=zygote.x.b,
                    "gamete.x.a"=gamete.x.a,"gamete.x.b"=gamete.x.b)
    return(pop.list)
}

recombination.and.selection = function(pop.list,zygote.w,pops,r) {
    #getting the zygote fitnesses
    AABB.w = zygote.w[1,1]; AABb.w = zygote.w[1,2]; AAbb.w = zygote.w[1,3]
    AaBB.w = zygote.w[2,1]; AaBb.w = zygote.w[2,2]; Aabb.w = zygote.w[2,3]
    aaBB.w = zygote.w[3,1]; aaBb.w = zygote.w[3,2]; aabb.w = zygote.w[3,3]
    #index over both species
    for (i in 1:2) {
        zygote.x = pop.list[[i]]
        gamete.x = pop.list[[(i+2)]]
        zygote.index = seq(1,pops*3,3)
        zygote.index = cbind(zygote.index,seq(3,pops*3,3))
        if (pops==1) {
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
            if (i==1) {
                gamete.x = c(selection.a[1],selection.b[1],selection.b[2],selection.a[2])
            }
            if (i==2) {
                gamete.x = c(selection.a[1],selection.b[1],selection.b[2],selection.a[2])
            }
        }
        if (pops>1) {
            #index over populations within each species
            for (j in 1:pops) {
                AB.x = gamete.x[j,1]; Ab.x = gamete.x[j,2]
                aB.x = gamete.x[j,3]; ab.x = gamete.x[j,4]
                zygote.pop.x = zygote.x[zygote.index[j,1]:zygote.index[j,2],]
                #marginalizing the gamete fitnesses
                AB.w = (AABB.w*AB.x)+(AABb.w*Ab.x)+(AaBB.w*aB.x)+(AaBb.w*ab.x)
                Ab.w = (AABb.w*AB.x)+(AAbb.w*Ab.x)+(AaBb.w*aB.x)+(Aabb.w*ab.x)
                aB.w = (AaBB.w*AB.x)+(AaBb.w*Ab.x)+(aaBB.w*aB.x)+(aaBb.w*ab.x)
                ab.w = (AaBb.w*AB.x)+(Aabb.w*Ab.x)+(aaBb.w*aB.x)+(aabb.w*ab.x)
                gamete.w = c(AB.w, Ab.w, aB.w, ab.w)
                #mean fitness
                mean.w = sum(zygote.w*zygote.pop.x)
                #LD
                D = AB.x*ab.x - Ab.x*aB.x
                #gamete frequencies after recombination and selection
                selection.a = ((c(AB.w,ab.w)*c(AB.x,ab.x)) - r*zygote.w[2,2]*D)/mean.w
                selection.b = ((c(Ab.w,aB.w)*c(Ab.x,aB.x)) + r*zygote.w[2,2]*D)/mean.w
                gamete.x[j,] = c(selection.a[1],selection.b[1],selection.b[2],selection.a[2])
            }
        }
        pop.list[[i+2]] = gamete.x
    }
    return(pop.list)
}

migration = function(pop.list,source.freqs,m,msub,pops) {
    gametes.a = pop.list[[3]]
    gametes.b = pop.list[[4]]
    if (pops==1) {
        sp.a = (1-m-msub)*gametes.a + m*gametes.b + msub*source.freqs
        sp.b = (1-m-msub)*gametes.b + m*gametes.a + msub*rev(source.freqs)
        gametes.a = sp.a
        gametes.b = sp.b
    }
    if (pops==2) {
        hybrid.pop.a = (1-m-msub)*gametes.a[pops,] + m*gametes.b[pops,] + msub*gametes.a[(pops-1),]
        hybrid.pop.b = (1-m-msub)*gametes.b[pops,] + m*gametes.a[pops,] + msub*gametes.b[(pops-1),]
        periph.pop.a = (1-(2*msub))*gametes.a[1,] + msub*(source.freqs + gametes.a[pops,])
        periph.pop.b = (1-(2*msub))*gametes.b[1,] + msub*(rev(source.freqs) + gametes.b[pops,])
        gametes.a = rbind(periph.pop.a,hybrid.pop.a)
        gametes.b = rbind(periph.pop.b,hybrid.pop.b)
    }
    if (pops>2) {
        hybrid.pop.a = (1-m-msub)*gametes.a[pops,] + m*gametes.b[pops,] + msub*gametes.a[(pops-1),]
        hybrid.pop.b = (1-m-msub)*gametes.b[pops,] + m*gametes.a[pops,] + msub*gametes.b[(pops-1),]
        periph.pop.a = (1-(2*msub))*gametes.a[1,] + msub*(source.freqs + gametes.a[pops,])
        periph.pop.b = (1-(2*msub))*gametes.b[1,] + msub*(rev(source.freqs) + gametes.b[pops,])
        sub.pops.a = gametes.a
        sub.pops.b = gametes.b
        for (i in 2:(pops-1)) {
            sub.pops.a[i,] = (1-(2*msub))*gametes.a[i,] + msub*(gametes.a[i-1,] + gametes.a[i+1,])
            sub.pops.b[i,] = (1-(2*msub))*gametes.b[i,] + msub*(gametes.b[i-1,] + gametes.b[i+1,])
        }
        gametes.a = rbind(periph.pop.a,sub.pops.a[2:(pops-1),])
        gametes.a = rbind(gametes.a,hybrid.pop.a)
        gametes.b = rbind(periph.pop.b,sub.pops.b[2:(pops-1),])
        gametes.b = rbind(gametes.b,hybrid.pop.b)
    }
    pop.list[[3]] = gametes.a
    pop.list[[4]] = gametes.b
    return(pop.list)
}

make.zygotes = function(pop.list,pops,N) {
    for (i in 1:2) {
        zygote.x = pop.list[[i]]
        gamete.x = pop.list[[(i+2)]]
        zygote.index = seq(1,pops*3,3)
        zygote.index = cbind(zygote.index,seq(3,pops*3,3))
        if (pops==1) {
            AB.x = gamete.x[1]; Ab.x = gamete.x[2]
            aB.x = gamete.x[3]; ab.x = gamete.x[4]
            zygote.x[1,1] = AB.x*AB.x                 #AABB
            zygote.x[1,2] = 2*AB.x*Ab.x               #AABb
            zygote.x[1,3] = Ab.x*Ab.x                 #AAbb
            zygote.x[2,1] = 2*AB.x*aB.x               #AaBB
            zygote.x[2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
            zygote.x[2,3] = 2*Ab.x*ab.x               #Aabb
            zygote.x[3,1] = aB.x*aB.x                 #aaBB
            zygote.x[3,2] = 2*aB.x*ab.x               #aaBb
            zygote.x[3,3] = ab.x*ab.x                 #aabb
        }
        if (pops>1) {
            for (j in 1:pops) {
                AB.x = gamete.x[j,1]; Ab.x = gamete.x[j,2]
                aB.x = gamete.x[j,3]; ab.x = gamete.x[j,4]
                index.a = zygote.index[j,1]
                index.b = zygote.index[j,2]
                zygote.x[c(index.a:index.b),][1,1] = AB.x*AB.x                 #AABB
                zygote.x[c(index.a:index.b),][1,2] = 2*AB.x*Ab.x               #AABb
                zygote.x[c(index.a:index.b),][1,3] = Ab.x*Ab.x                 #AAbb
                zygote.x[c(index.a:index.b),][2,1] = 2*AB.x*aB.x               #AaBB
                zygote.x[c(index.a:index.b),][2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
                zygote.x[c(index.a:index.b),][2,3] = 2*Ab.x*ab.x               #Aabb
                zygote.x[c(index.a:index.b),][3,1] = aB.x*aB.x                 #aaBB
                zygote.x[c(index.a:index.b),][3,2] = 2*aB.x*ab.x               #aaBb
                zygote.x[c(index.a:index.b),][3,3] = ab.x*ab.x                 #aabb
            }
        }
        pop.list[[i]] = zygote.x
    }
    return(pop.list)
}

make.zygotes.WF = function(pop.list,pops,N) {
    for (i in 1:2) {
        zygote.x = pop.list[[i]]
        gamete.x = pop.list[[(i+2)]]
        zygote.index = seq(1,pops*3,3)
        zygote.index = cbind(zygote.index,seq(3,pops*3,3))
        if (pops==1) {
            gamete.vec = rmultinom(1,(2*N),prob=gamete.x)/(2*N)
            AB.x = gamete.vec[1]; Ab.x = gamete.vec[2]
            aB.x = gamete.vec[3]; ab.x = gamete.vec[4]
            zygote.x[1,1] = AB.x*AB.x                 #AABB
            zygote.x[1,2] = 2*AB.x*Ab.x               #AABb
            zygote.x[1,3] = Ab.x*Ab.x                 #AAbb
            zygote.x[2,1] = 2*AB.x*aB.x               #AaBB
            zygote.x[2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
            zygote.x[2,3] = 2*Ab.x*ab.x               #Aabb
            zygote.x[3,1] = aB.x*aB.x                 #aaBB
            zygote.x[3,2] = 2*aB.x*ab.x               #aaBb
            zygote.x[3,3] = ab.x*ab.x                 #aabb
            gamete.x = gamete.vec
        }
        if (pops>1) {
            for (j in 1:pops) {
                gamete.vec = rmultinom(1,(2*N),prob=gamete.x[j,])/(2*N)
                AB.x = gamete.vec[1]; Ab.x = gamete.vec[2]
                aB.x = gamete.vec[3]; ab.x = gamete.vec[4]
                index.a = zygote.index[j,1]
                index.b = zygote.index[j,2]
                zygote.x[c(index.a:index.b),][1,1] = AB.x*AB.x                 #AABB
                zygote.x[c(index.a:index.b),][1,2] = 2*AB.x*Ab.x               #AABb
                zygote.x[c(index.a:index.b),][1,3] = Ab.x*Ab.x                 #AAbb
                zygote.x[c(index.a:index.b),][2,1] = 2*AB.x*aB.x               #AaBB
                zygote.x[c(index.a:index.b),][2,2] = 2*AB.x*ab.x + 2*Ab.x*aB.x #AaBb
                zygote.x[c(index.a:index.b),][2,3] = 2*Ab.x*ab.x               #Aabb
                zygote.x[c(index.a:index.b),][3,1] = aB.x*aB.x                 #aaBB
                zygote.x[c(index.a:index.b),][3,2] = 2*aB.x*ab.x               #aaBb
                zygote.x[c(index.a:index.b),][3,3] = ab.x*ab.x                 #aabb
                zygote.x = zygote.x
                gamete.x[j,] = gamete.vec
            }
        }
        pop.list[[i]] = zygote.x
        pop.list[[i+2]] = gamete.vec
    }
    return(pop.list)
}

simulate.trajecs = function(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs) {
    pop.list = make.populations(source.freqs,pops)
    zygote.w = make.fitnessmatrix(alpha,alpha1,beta,beta1,s)
    gamete.trajecs = matrix(nrow=(t+1),ncol=pops*8)
    gamete.trajecs[1,1:4] = pop.list[[3]]
    gamete.trajecs[1,5:8] = pop.list[[4]]
    allele.trajecs = matrix(nrow=(t+1),ncol=pops*4)
    allele.trajecs[1,1:2] = c(sum(pop.list[[3]][c(1,2)]),sum(pop.list[[3]][c(2,4)]))
    allele.trajecs[1,3:4] = c(sum(pop.list[[4]][c(1,2)]),sum(pop.list[[4]][c(2,4)]))
    for (j in 2:(t+1)) {
        pop.list = recombination.and.selection(pop.list,zygote.w,pops,r)
        pop.list = migration(pop.list,source.freqs,m,msub,pops)
        pop.list = make.zygotes(pop.list,pops,N)
        gamete.trajecs[j,(1:(pops*8/2))]  = pop.list[[3]]
        gamete.trajecs[j,5:8] = pop.list[[4]]
        allele.trajecs[j,1:2] = c(sum(pop.list[[3]][c(1,2)]),sum(pop.list[[3]][c(2,4)]))
        allele.trajecs[j,3:4] = c(sum(pop.list[[4]][c(1,2)]),sum(pop.list[[4]][c(2,4)]))
    }
    trajec.list = list('AB.matrix.a'=gamete.trajecs[,1],'Ab.matrix.a'=gamete.trajecs[,2],
                       'aB.matrix.a'=gamete.trajecs[,3],'ab.matrix.a'=gamete.trajecs[,4],
                       'AB.matrix.b'=gamete.trajecs[,5],'Ab.matrix.b'=gamete.trajecs[,6],
                       'aB.matrix.b'=gamete.trajecs[,7],'ab.matrix.b'=gamete.trajecs[,8],
                       'A.matrix.a'=allele.trajecs[,1],'b.matrix.a'=allele.trajecs[,2],
                       'A.matrix.b'=allele.trajecs[,3],'b.matrix.b'=allele.trajecs[,4],
                       'zygote.w'=zygote.w)
    return(trajec.list)
}

simulate.trajecs.WF = function(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs) {
    pop.list = make.populations(source.freqs,pops)
    zygote.w = make.fitnessmatrix(alpha,alpha1,beta,beta1,s)
    gamete.trajecs = matrix(nrow=(t+1),ncol=pops*8)
    gamete.trajecs[1,1:4] = pop.list[[3]]
    gamete.trajecs[1,5:8] = pop.list[[4]]
    allele.trajecs = matrix(nrow=(t+1),ncol=pops*4)
    allele.trajecs[1,1:2] = c(sum(pop.list[[3]][c(1,2)]),sum(pop.list[[3]][c(2,4)]))
    allele.trajecs[1,3:4] = c(sum(pop.list[[4]][c(1,2)]),sum(pop.list[[4]][c(2,4)]))
    for (j in 2:(t+1)) {
        pop.list = recombination.and.selection(pop.list,zygote.w,pops,r)
        pop.list = migration(pop.list,source.freqs,m,msub,pops)
        pop.list = make.zygotes.WF(pop.list,pops,N)
        gamete.trajecs[j,(1:(pops*8/2))]  = pop.list[[3]]
        gamete.trajecs[j,5:8] = pop.list[[4]]
        allele.trajecs[j,1:2] = c(sum(pop.list[[3]][c(1,2)]),sum(pop.list[[3]][c(2,4)]))
        allele.trajecs[j,3:4] = c(sum(pop.list[[4]][c(1,2)]),sum(pop.list[[4]][c(2,4)]))
    }
    trajec.list = list('AB.matrix.a'=gamete.trajecs[,1],'Ab.matrix.a'=gamete.trajecs[,2],
                       'aB.matrix.a'=gamete.trajecs[,3],'ab.matrix.a'=gamete.trajecs[,4],
                       'AB.matrix.b'=gamete.trajecs[,5],'Ab.matrix.b'=gamete.trajecs[,6],
                       'aB.matrix.b'=gamete.trajecs[,7],'ab.matrix.b'=gamete.trajecs[,8],
                       'A.matrix.a'=allele.trajecs[,1],'b.matrix.a'=allele.trajecs[,2],
                       'A.matrix.b'=allele.trajecs[,3],'b.matrix.b'=allele.trajecs[,4],
                       'zygote.w'=zygote.w)
    return(trajec.list)
}

WF.trajec.reps = function(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs,reps) {
    AB.matrix.a = matrix(nrow=t+1,ncol=reps)
    Ab.matrix.a = matrix(nrow=t+1,ncol=reps)
    aB.matrix.a = matrix(nrow=t+1,ncol=reps)
    ab.matrix.a = matrix(nrow=t+1,ncol=reps)
    AB.matrix.b = matrix(nrow=t+1,ncol=reps)
    Ab.matrix.b = matrix(nrow=t+1,ncol=reps)
    aB.matrix.b = matrix(nrow=t+1,ncol=reps)
    ab.matrix.b = matrix(nrow=t+1,ncol=reps)
    A.matrix.a  = matrix(nrow=t+1,ncol=reps)
    b.matrix.a  = matrix(nrow=t+1,ncol=reps)
    A.matrix.b  = matrix(nrow=t+1,ncol=reps)
    b.matrix.b  = matrix(nrow=t+1,ncol=reps)
    if (reps>1) {
        pb = txtProgressBar(1,reps,1,style=3)
        for (b in 1:reps) {
            trajec.list.WF = simulate.trajecs.WF(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs)
            AB.matrix.a[,b] = trajec.list.WF[[1]]
            Ab.matrix.a[,b] = trajec.list.WF[[2]]
            aB.matrix.a[,b] = trajec.list.WF[[3]]
            ab.matrix.a[,b] = trajec.list.WF[[4]]
            AB.matrix.b[,b] = trajec.list.WF[[5]]
            Ab.matrix.b[,b] = trajec.list.WF[[6]]
            aB.matrix.b[,b] = trajec.list.WF[[7]]
            ab.matrix.b[,b] = trajec.list.WF[[8]]
            A.matrix.a[,b]  = trajec.list.WF[[9]]
            b.matrix.a[,b]  = trajec.list.WF[[10]]
            A.matrix.b[,b]  = trajec.list.WF[[11]]
            b.matrix.b[,b]  = trajec.list.WF[[12]]
            setTxtProgressBar(pb, b)
        }
        close(pb)
    }
    if (reps==1) {
        for (b in 1:reps) {
            trajec.list.WF = simulate.trajecs.WF(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs)
            gamete.trajecs.WF = trajec.list.WF[[1]]
            allele.trajecs.WF = trajec.list.WF[[2]]
            AB.matrix.a[,b] = gamete.trajecs.WF[,1]
            Ab.matrix.a[,b] = gamete.trajecs.WF[,2]
            aB.matrix.a[,b] = gamete.trajecs.WF[,3]
            ab.matrix.a[,b] = gamete.trajecs.WF[,4]
            AB.matrix.b[,b] = gamete.trajecs.WF[,5]
            Ab.matrix.b[,b] = gamete.trajecs.WF[,6]
            aB.matrix.b[,b] = gamete.trajecs.WF[,7]
            ab.matrix.b[,b] = gamete.trajecs.WF[,8]
            A.matrix.a[,b]  = allele.trajecs.WF[,1]
            b.matrix.a[,b]  = allele.trajecs.WF[,2]
            A.matrix.b[,b]  = allele.trajecs.WF[,3]
            b.matrix.b[,b]  = allele.trajecs.WF[,4]
        }
    }
    zygote.w = trajec.list.WF[[13]]
    WF.rep.list = list('AB.matrix.a'=AB.matrix.a,'Ab.matrix.a'=Ab.matrix.a,
                       'aB.matrix.a'=aB.matrix.a,'ab.matrix.a'=ab.matrix.a,
                       'AB.matrix.b'=AB.matrix.b,'Ab.matrix.b'=Ab.matrix.b,
                       'aB.matrix.b'=aB.matrix.b,'ab.matrix.b'=ab.matrix.b,
                       'A.matrix.a'=A.matrix.a,'b.matrix.a'=b.matrix.a,
                       'A.matrix.b'=A.matrix.b,'b.matrix.b'=b.matrix.b,'zygote.w'=zygote.w)
    return(WF.rep.list)
}

find.equilibrium.freqs = function(pops,N,m,msub,r,alpha,alpha1,beta,beta1,s,source.freqs) {
    pop.list = make.populations(source.freqs,pops)
    zygote.w = make.fitnessmatrix(alpha,alpha1,beta,beta1,s)
    gamete.trajecs = matrix(nrow=(1000+1),ncol=pops*8)
    gamete.trajecs[1,1:4] = pop.list[[3]]
    gamete.trajecs[1,5:8] = pop.list[[4]]
    allele.trajecs = matrix(nrow=(1000+1),ncol=pops*4)
    allele.trajecs[1,1:2] = c(sum(pop.list[[3]][c(1,2)]),sum(pop.list[[3]][c(2,4)]))
    allele.trajecs[1,3:4] = c(sum(pop.list[[4]][c(1,2)]),sum(pop.list[[4]][c(2,4)]))
    for (j in 2:(t+1)) {
        pop.list = recombination.and.selection(pop.list,zygote.w,pops,r)
        pop.list = migration(pop.list,source.freqs,m,msub,pops)
        pop.list = make.zygotes(pop.list,pops,N)
        gamete.trajecs[j,(1:(pops*8/2))]  = pop.list[[3]]
        gamete.trajecs[j,5:8] = pop.list[[4]]
        allele.trajecs[j,1:2] = c(sum(pop.list[[3]][c(1,2)]),sum(pop.list[[3]][c(2,4)]))
        allele.trajecs[j,3:4] = c(sum(pop.list[[4]][c(1,2)]),sum(pop.list[[4]][c(2,4)]))
        frequency.diffs = abs(allele.trajecs[j,]-allele.trajecs[(j-1),])
        if (length(which(frequency.diffs<1e-5))==4) {
            break
        }
    }
    eq.freq.list = list('AB.eq.a'=gamete.trajecs[j,1],'Ab.eq.a'=gamete.trajecs[j,2],
                        'aB.eq.a'=gamete.trajecs[j,3],'ab.eq.a'=gamete.trajecs[j,4],
                        'AB.eq.b'=gamete.trajecs[j,5],'Ab.eq.b'=gamete.trajecs[j,6],
                        'aB.eq.b'=gamete.trajecs[j,7],'ab.eq.b'=gamete.trajecs[j,8],
                        'A.eq.a'=allele.trajecs[j,1],'b.eq.a'=allele.trajecs[j,2],
                        'A.eq.b'=allele.trajecs[j,3],'b.eq.b'=allele.trajecs[j,4],
                        'zygote.w'=zygote.w)
    return(eq.freq.list)
}

plot.WF.versus.deterministic = function(WF.det.list,WF.sto.list) {
#    pdf("dmiWF_trajecs.pdf",width=10,height=7)
    library('RColorBrewer')
    colors = brewer.pal(4,'Paired')
    par(mfrow=c(1,2))
    plot(WF.sto.list[[1]][,1],col=colors[1],typ='l',main='Species 1',lwd=2,ylim=c(0,1),xlab='generation',ylab='frequency')
    lines(WF.sto.list[[2]][,1],col=colors[2],lwd=1)
    lines(WF.sto.list[[3]][,1],col=colors[3],lwd=1)
    lines(WF.sto.list[[4]][,1],col=colors[4],lwd=1)
    for (i in 2:reps) {
        lines(WF.sto.list[[1]][,i],col=colors[1],lwd=1)
        lines(WF.sto.list[[2]][,i],col=colors[2],lwd=1)
        lines(WF.sto.list[[3]][,i],col=colors[3],lwd=1)
        lines(WF.sto.list[[4]][,i],col=colors[4],lwd=1)
    }
    lines(WF.det.list[[1]],col='black',lwd=1)
    lines(WF.det.list[[2]],col='black',lwd=1)
    lines(WF.det.list[[3]],col='black',lwd=1)
    lines(WF.det.list[[4]],col='black',lwd=1)
    legend(10,.6,c("AB","Ab","aB","ab"),c(colors))
    plot(WF.sto.list[[5]][,1],col=colors[1],typ='l',main='Species 2',lwd=2,ylim=c(0,1),xlab='generation',ylab='frequency')
    lines(WF.sto.list[[6]][,1],col=colors[2],lwd=1)
    lines(WF.sto.list[[7]][,1],col=colors[3],lwd=1)
    lines(WF.sto.list[[8]][,1],col=colors[4],lwd=1)
    for (i in 2:reps) {
        lines(WF.sto.list[[5]][,i],col=colors[1],lwd=1)
        lines(WF.sto.list[[6]][,i],col=colors[2],lwd=1)
        lines(WF.sto.list[[7]][,i],col=colors[3],lwd=1)
        lines(WF.sto.list[[8]][,i],col=colors[4],lwd=1)
    }
    lines(WF.det.list[[5]],col='black',lwd=1)
    lines(WF.det.list[[6]],col='black',lwd=1)
    lines(WF.det.list[[7]],col='black',lwd=1)
    lines(WF.det.list[[8]],col='black',lwd=1)
#    dev.off()
}
