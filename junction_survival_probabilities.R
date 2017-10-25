
# gets the junction survival probabilities
transition.rates.eq = function(eq.freq.list,r,v,w,t,t.prime) {
    AB.eq = eq.freq.list$AB.eq.a
    Ab.eq = eq.freq.list$Ab.eq.a
    aB.eq = eq.freq.list$aB.eq.a
    ab.eq = eq.freq.list$ab.eq.a
    A.eq  = eq.freq.list$A.eq.a
    b.eq  = eq.freq.list$b.eq.a
    zygote.w   = eq.freq.list$zygote.w
    tau.vec = c(1:(min(t,t.prime)-1))
    for (tau in 1:(min(t,t.prime)-1)) {  
        gamete.trajecs = cbind(AB.eq,Ab.eq,aB.eq,ab.eq)
        allele.trajecs = cbind(A.eq,b.eq)
        trajec.list = list("gamete.trajecs"=gamete.trajecs,"allele.trajecs"=allele.trajecs,"zygote.w"=zygote.w)
        transition.rates.vec.WF = continuous.transition.rate.eq(trajec.list,r,v,w,t.prime,tau)                                              
        tau.vec[tau] = transition.rates.vec.WF 
    }   
    tau.result = sum(tau.vec)  
    return(tau.result)    
}                         

# gets the junction survival probabilities
transition.rates = function(WF.rep.list,r,v,w,t,t.prime) {
    AB.trajecs = WF.rep.list$AB.matrix.a
    Ab.trajecs = WF.rep.list$Ab.matrix.a
    aB.trajecs = WF.rep.list$aB.matrix.a
    ab.trajecs = WF.rep.list$ab.matrix.a
    A.trajecs  = WF.rep.list$A.matrix.a
    b.trajecs  = WF.rep.list$b.matrix.a
    zygote.w   = WF.rep.list$zygote.w
    if (is.vector(AB.trajecs)) {
        reps = 1
    } 
    if (!is.vector(AB.trajecs)) {
        reps = ncol(AB.trajecs)
        tau.matrix = matrix(nrow=(min(t,t.prime)-1),ncol=reps)
    }     
    for (tau in 1:(min(t,t.prime)-1)) {  
        if (reps>1) {
            transition.rates.vec.WF = c(1:reps) 
            for (i in 1:reps) {
                gamete.trajecs = cbind(AB.trajecs[,i],Ab.trajecs[,i],aB.trajecs[,i],ab.trajecs[,i])
                allele.trajecs = cbind(A.trajecs[,i],b.trajecs[,i])
                trajec.list = list("gamete.trajecs"=gamete.trajecs,"allele.trajecs"=allele.trajecs,"zygote.w"=zygote.w)
                transition.rates.vec.WF[i] = continuous.transition.rate(trajec.list,r,v,w,t.prime,tau)
            }
            tau.matrix[tau,] = transition.rates.vec.WF                                                                             
        }
        if (reps==1) {
            gamete.trajecs = cbind(AB.trajecs,Ab.trajecs,aB.trajecs,ab.trajecs)
            allele.trajecs = cbind(A.trajecs,b.trajecs)
            trajec.list = list("gamete.trajecs"=gamete.trajecs,"allele.trajecs"=allele.trajecs,"zygote.w"=zygote.w)
            transition.rates.vec.WF = continuous.transition.rate(trajec.list,r,v,w,t.prime,tau)                                              
            tau.results = sum(transition.rates.vec.WF)
        } 
    }        
    if (reps>1) {
        tau.results = colSums(tau.matrix)      
    }   
    return(tau.results)    
}                         


# in tau.matrix, columns are WF trajectory reps, and rows are the products for each tau in equation 2.                                                   
continuous.transition.rate = function(trajec.list,r,v,w,t.prime,tau) {  
    gamete.trajecs = trajec.list[[1]]  
    allele.trajecs = trajec.list[[2]]
    zygote.w       = trajec.list[[3]]
    AABB.w = zygote.w[1,1]; AABb.w = zygote.w[1,2]; AAbb.w = zygote.w[1,3]
    AaBB.w = zygote.w[2,1]; AaBb.w = zygote.w[2,2]; Aabb.w = zygote.w[2,3]
    aaBB.w = zygote.w[3,1]; aaBb.w = zygote.w[3,2]; aabb.w = zygote.w[3,3]    
    A.allele = allele.trajecs[,1]           
    B.allele = 1-allele.trajecs[,2]         
    a.allele = 1-allele.trajecs[,1]         
    b.allele = allele.trajecs[,2]  
    for (j in (tau+1):(t.prime)) {     
        x.1 = gamete.trajecs[j,1]                 
        x.2 = gamete.trajecs[j,2]               
        x.3 = gamete.trajecs[j,3]               
        x.4 = gamete.trajecs[j,4]                   
        initial.junction.prob = c(rep(0,4),rep(2*A.allele[j]*b.allele[j],4),rep(2*a.allele[j]*B.allele[j],4),rep(0,5))         
        initial.junction.prob = c(rep(0,4),rep(c(gamete.trajecs[j,]),2),rep(0,5))*initial.junction.prob
        initial.junction.prob[17] = 1-sum(initial.junction.prob[1:16])                                           
        P.transition.matrix = matrix(nrow=17,ncol=17)                  
        P.transition.matrix[1,c(1:16)]  = c(rep(AaBb.w,12),rep(0,4))*
                                          c(rep(.5,12),rep(0,4))*
                                          c(rep((1-r),4),rep(r,8),rep(1,4))*
                                          c(rep(1,4),rep((w/(v+w)),4),rep((v/(v+w)),4), rep(1,4))*
                                          c(rep(c(x.4,x.2,x.3,x.1),4))
        P.transition.matrix[2,c(1:16)]  = c(rep(AABb.w,8),rep(0,8))*
                                          c(rep(.5,8),rep(0,8))*
                                          c((((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,12))* 
                                          c(rep(1,4),rep(r,4),rep(1,8))*
                                          c(rep(1,4),rep((w/(v+w)),4),rep(1,8))*
                                          c(rep(1,4),c(x.4,x.2,x.3,x.1),rep(1,8))                                
        P.transition.matrix[3,c(1:16)]  = c(rep(AaBB.w,4),rep(0,4),rep(AaBB.w,4),rep(0,4))*
                                          c(rep(.5,4),rep(0,4),rep(.5,4),rep(0,4))*
                                          c((((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,12))* 
                                          c(rep(1,8),rep((r*(v/(v+w))),4),rep(1,4))*
                                          c(rep(1,8),c(x.4,x.2,x.3,x.1),rep(1,4))
        P.transition.matrix[4,c(1:16)]  = AABB.w*.5*c((((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(0,12)) 
        P.transition.matrix[5,c(1:16)]  = c(rep(0,4),rep(Aabb.w,4),rep(0,4),rep(Aabb.w,4))*
                                          c(rep(0,4),rep(.5,4),rep(0,4),rep(.5,4))*
                                          c(rep(1,4),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,8))* 
                                          c(rep(1,12),rep((r*(v/(v+w))),4))*
                                          c(rep(1,12),c(x.4,x.2,x.3,x.1))                                 
        P.transition.matrix[6,c(1:16)]  = c(rep(0,4),rep(AAbb.w,4),rep(0,8))*
                                          c(rep(0,4),rep(.5,4),rep(0,8))*
                                          c(rep(0,4),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(0,8)) 
        P.transition.matrix[7,c(1:16)]  = c(rep(AaBb.w,8),rep(0,4),rep(AaBb.w,4))*
                                          c(rep(.5,8),rep(0,4),rep(.5,4))*
                                          c(rep(r,4),rep((1-r),4),rep(1,4),rep(r,4))*
                                          c(rep((w/(v+w)),4),rep(1,8),rep((v/(v+w)),4))*
                                          c(rep(c(x.4,x.2,x.3,x.1),2),rep(1,4),c(x.4,x.2,x.3,x.1))
        P.transition.matrix[8,c(1:16)]  = c(rep(AABb.w,8),rep(0,8))*
                                          c(rep(.5,8),rep(0,8))*
                                          c(rep(r*(w/(v+w)),4),rep(1,12))*
                                          c(rep(1,4),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,8))*                                  
                                          c(c(x.4,x.2,x.3,x.1),rep(1,12))
        P.transition.matrix[9,c(1:16)]  = c(rep(0,8),rep(aaBb.w,8))*
                                          c(rep(0,8),rep(.5,8))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,4))*
                                          c(rep(1,12),rep(r*(w/(v+w)),4))*
                                          c(rep(1,12),c(x.4,x.2,x.3,x.1))
        P.transition.matrix[10,c(1:16)] = c(rep(AaBb.w,4),rep(0,4),rep(AaBb.w,8))*
                                          c(rep(.5,4),rep(0,4),rep(.5,8))*
                                          c(rep(r,4),rep(1,4),rep((1-r),4),rep(r,4))*
                                          c(rep(v/(v+w),4),rep(1,8),rep(w/(v+w),4))*
                                          c(c(x.4,x.2,x.3,x.1),rep(1,4),rep(c(x.4,x.2,x.3,x.1),2))
        P.transition.matrix[11,c(1:16)] = c(rep(0,8),rep(aaBB.w,4),rep(0,4))*
                                          c(rep(0,8),rep(.5,4),rep(0,4))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(1,4))                         
        P.transition.matrix[12,c(1:16)] = c(rep(AaBB.w,4),rep(0,4),rep(AaBB.w,4),rep(0,4))*
                                          c(rep(.5,4),rep(0,4),rep(.5,4),rep(0,4))*
                                          c(rep(r,4),rep(1,12))*
                                          c(rep(v/(v+w),4),rep(1,12))*
                                          c(c(x.4,x.2,x.3,x.1),rep(1,12))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,4))                                  
        P.transition.matrix[13,c(1:16)] = c(rep(0,12),rep(aabb.w,4))*
                                          c(rep(0,12),rep(.5,4))*
                                          c(rep(1,12),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)))                         
        P.transition.matrix[14,c(1:16)] = c(rep(0,4),rep(Aabb.w,4),rep(0,4),rep(Aabb.w,4))*
                                          c(rep(0,4),rep(.5,4),rep(0,4),rep(.5,4))*
                                          c(rep(1,4),rep(r*(v/(v+w)),4),rep(1,8))*
                                          c(rep(1,4),c(x.4,x.2,x.3,x.1),rep(1,8))*                                  	
                                          c(rep(1,12),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)))                                  
        P.transition.matrix[15,c(1:16)] = c(rep(0,8),rep(aaBb.w,8))*
                                          c(rep(0,8),rep(.5,8))*
                                          c(rep(1,8),rep(r*(w/(v+w)),4),rep(1,4))*
                                          c(rep(1,8),c(x.4,x.2,x.3,x.1),rep(1,4))*
                                          c(rep(1,12),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)))
        P.transition.matrix[16,c(1:16)] = c(rep(0,4),rep(AaBb.w,12))*
                                          c(rep(0,4),rep(.5,12))*
                                          c(rep(1,4),rep(r,8),rep((1-r),4))*
                                          c(rep(1,4),rep(v/(v+w),4),rep(w/(v+w),4),rep(1,4))*
                                          c(rep(1,4),rep(c(x.4,x.2,x.3,x.1),3))
        P.transition.matrix[1,17]  = 1-sum(P.transition.matrix[1,c(1:16)])
        P.transition.matrix[2,17]  = 1-sum(P.transition.matrix[2,c(1:16)])   
        P.transition.matrix[3,17]  = 1-sum(P.transition.matrix[3,c(1:16)]) 
        P.transition.matrix[4,17]  = 1-sum(P.transition.matrix[4,c(1:16)])         
        P.transition.matrix[5,17]  = 1-sum(P.transition.matrix[5,c(1:16)])         
        P.transition.matrix[6,17]  = 1-sum(P.transition.matrix[6,c(1:16)])     
        P.transition.matrix[7,17]  = 1-sum(P.transition.matrix[7,c(1:16)])     
        P.transition.matrix[8,17]  = 1-sum(P.transition.matrix[8,c(1:16)])     
        P.transition.matrix[9,17]  = 1-sum(P.transition.matrix[9,c(1:16)])     
        P.transition.matrix[10,17] = 1-sum(P.transition.matrix[10,c(1:16)])     
        P.transition.matrix[11,17] = 1-sum(P.transition.matrix[11,c(1:16)])     
        P.transition.matrix[12,17] = 1-sum(P.transition.matrix[12,c(1:16)])     
        P.transition.matrix[13,17] = 1-sum(P.transition.matrix[13,c(1:16)])   
        P.transition.matrix[14,17] = 1-sum(P.transition.matrix[14,c(1:16)]) 
        P.transition.matrix[15,17] = 1-sum(P.transition.matrix[15,c(1:16)]) 
        P.transition.matrix[16,17] = 1-sum(P.transition.matrix[16,c(1:16)]) 
        P.transition.matrix[17,]   = c(rep(0,16),1)            
        if (j==(tau+1)) {              
            product.vec = initial.junction.prob %*% P.transition.matrix       
        }                
        if (j!=(tau+1)) {      
            product.vec = product.vec %*% P.transition.matrix
        }                       
    }                                 
    junction.surv.prob = 1-product.vec[17]   
    return(junction.surv.prob)
}   


# in tau.matrix, columns are WF trajectory reps, and rows are the products for each tau in equation 2.                                                   
continuous.transition.rate.eq = function(trajec.list,r,v,w,t.prime,tau) {  
    gamete.trajecs = trajec.list[[1]]  
    allele.trajecs = trajec.list[[2]]
    zygote.w       = trajec.list[[3]]
    AABB.w = zygote.w[1,1]; AABb.w = zygote.w[1,2]; AAbb.w = zygote.w[1,3]
    AaBB.w = zygote.w[2,1]; AaBb.w = zygote.w[2,2]; Aabb.w = zygote.w[2,3]
    aaBB.w = zygote.w[3,1]; aaBb.w = zygote.w[3,2]; aabb.w = zygote.w[3,3]   
    A.allele = allele.trajecs[1]           
    B.allele = 1-allele.trajecs[2]         
    a.allele = 1-allele.trajecs[1]         
    b.allele = allele.trajecs[2]        
    for (j in (tau+1):(t.prime)) {     
        x.1 = gamete.trajecs[1]                 
        x.2 = gamete.trajecs[2]               
        x.3 = gamete.trajecs[3]               
        x.4 = gamete.trajecs[4]                   
        initial.junction.prob = c(rep(0,4),rep(2*A.allele*b.allele,4),rep(2*a.allele*B.allele,4),rep(0,5))         
        initial.junction.prob = c(rep(0,4),rep(c(gamete.trajecs),2),rep(0,5))*initial.junction.prob
        initial.junction.prob[17] = 1-sum(initial.junction.prob[1:16])                                           
        P.transition.matrix = matrix(nrow=17,ncol=17)                  
        P.transition.matrix[1,c(1:16)]  = c(rep(AaBb.w,12),rep(0,4))*
                                          c(rep(.5,12),rep(0,4))*
                                          c(rep((1-r),4),rep(r,8),rep(1,4))*
                                          c(rep(1,4),rep((w/(v+w)),4),rep((v/(v+w)),4), rep(1,4))*
                                          c(rep(c(x.4,x.2,x.3,x.1),4))
        P.transition.matrix[2,c(1:16)]  = c(rep(AABb.w,8),rep(0,8))*
                                          c(rep(.5,8),rep(0,8))*
                                          c((((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,12))* 
                                          c(rep(1,4),rep(r,4),rep(1,8))*
                                          c(rep(1,4),rep((w/(v+w)),4),rep(1,8))*
                                          c(rep(1,4),c(x.4,x.2,x.3,x.1),rep(1,8))                                
        P.transition.matrix[3,c(1:16)]  = c(rep(AaBB.w,4),rep(0,4),rep(AaBB.w,4),rep(0,4))*
                                          c(rep(.5,4),rep(0,4),rep(.5,4),rep(0,4))*
                                          c((((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,12))* 
                                          c(rep(1,8),rep((r*(v/(v+w))),4),rep(1,4))*
                                          c(rep(1,8),c(x.4,x.2,x.3,x.1),rep(1,4))
        P.transition.matrix[4,c(1:16)]  = AABB.w*.5*c((((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(0,12)) 
        P.transition.matrix[5,c(1:16)]  = c(rep(0,4),rep(Aabb.w,4),rep(0,4),rep(Aabb.w,4))*
                                          c(rep(0,4),rep(.5,4),rep(0,4),rep(.5,4))*
                                          c(rep(1,4),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,8))* 
                                          c(rep(1,12),rep((r*(v/(v+w))),4))*
                                          c(rep(1,12),c(x.4,x.2,x.3,x.1))                                 
        P.transition.matrix[6,c(1:16)]  = c(rep(0,4),rep(AAbb.w,4),rep(0,8))*
                                          c(rep(0,4),rep(.5,4),rep(0,8))*
                                          c(rep(0,4),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(0,8)) 
        P.transition.matrix[7,c(1:16)]  = c(rep(AaBb.w,8),rep(0,4),rep(AaBb.w,4))*
                                          c(rep(.5,8),rep(0,4),rep(.5,4))*
                                          c(rep(r,4),rep((1-r),4),rep(1,4),rep(r,4))*
                                          c(rep((w/(v+w)),4),rep(1,8),rep((v/(v+w)),4))*
                                          c(rep(c(x.4,x.2,x.3,x.1),2),rep(1,4),c(x.4,x.2,x.3,x.1))
        P.transition.matrix[8,c(1:16)]  = c(rep(AABb.w,8),rep(0,8))*
                                          c(rep(.5,8),rep(0,8))*
                                          c(rep(r*(w/(v+w)),4),rep(1,12))*
                                          c(rep(1,4),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,8))*                                  
                                          c(c(x.4,x.2,x.3,x.1),rep(1,12))
        P.transition.matrix[9,c(1:16)]  = c(rep(0,8),rep(aaBb.w,8))*
                                          c(rep(0,8),rep(.5,8))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,4))*
                                          c(rep(1,12),rep(r*(w/(v+w)),4))*
                                          c(rep(1,12),c(x.4,x.2,x.3,x.1))
        P.transition.matrix[10,c(1:16)] = c(rep(AaBb.w,4),rep(0,4),rep(AaBb.w,8))*
                                          c(rep(.5,4),rep(0,4),rep(.5,8))*
                                          c(rep(r,4),rep(1,4),rep((1-r),4),rep(r,4))*
                                          c(rep(v/(v+w),4),rep(1,8),rep(w/(v+w),4))*
                                          c(c(x.4,x.2,x.3,x.1),rep(1,4),rep(c(x.4,x.2,x.3,x.1),2))
        P.transition.matrix[11,c(1:16)] = c(rep(0,8),rep(aaBB.w,4),rep(0,4))*
                                          c(rep(0,8),rep(.5,4),rep(0,4))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(1,4))                         
        P.transition.matrix[12,c(1:16)] = c(rep(AaBB.w,4),rep(0,4),rep(AaBB.w,4),rep(0,4))*
                                          c(rep(.5,4),rep(0,4),rep(.5,4),rep(0,4))*
                                          c(rep(r,4),rep(1,12))*
                                          c(rep(v/(v+w),4),rep(1,12))*
                                          c(c(x.4,x.2,x.3,x.1),rep(1,12))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,4))                                  
        P.transition.matrix[13,c(1:16)] = c(rep(0,12),rep(aabb.w,4))*
                                          c(rep(0,12),rep(.5,4))*
                                          c(rep(1,12),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)))                         
        P.transition.matrix[14,c(1:16)] = c(rep(0,4),rep(Aabb.w,4),rep(0,4),rep(Aabb.w,4))*
                                          c(rep(0,4),rep(.5,4),rep(0,4),rep(.5,4))*
                                          c(rep(1,4),rep(r*(v/(v+w)),4),rep(1,8))*
                                          c(rep(1,4),c(x.4,x.2,x.3,x.1),rep(1,8))*                                  	
                                          c(rep(1,12),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)))                                  
        P.transition.matrix[15,c(1:16)] = c(rep(0,8),rep(aaBb.w,8))*
                                          c(rep(0,8),rep(.5,8))*
                                          c(rep(1,8),rep(r*(w/(v+w)),4),rep(1,4))*
                                          c(rep(1,8),c(x.4,x.2,x.3,x.1),rep(1,4))*
                                          c(rep(1,12),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)))
        P.transition.matrix[16,c(1:16)] = c(rep(0,4),rep(AaBb.w,12))*
                                          c(rep(0,4),rep(.5,12))*
                                          c(rep(1,4),rep(r,8),rep((1-r),4))*
                                          c(rep(1,4),rep(v/(v+w),4),rep(w/(v+w),4),rep(1,4))*
                                          c(rep(1,4),rep(c(x.4,x.2,x.3,x.1),3))
        P.transition.matrix[1,17]  = 1-sum(P.transition.matrix[1,c(1:16)])
        P.transition.matrix[2,17]  = 1-sum(P.transition.matrix[2,c(1:16)])   
        P.transition.matrix[3,17]  = 1-sum(P.transition.matrix[3,c(1:16)]) 
        P.transition.matrix[4,17]  = 1-sum(P.transition.matrix[4,c(1:16)])         
        P.transition.matrix[5,17]  = 1-sum(P.transition.matrix[5,c(1:16)])         
        P.transition.matrix[6,17]  = 1-sum(P.transition.matrix[6,c(1:16)])     
        P.transition.matrix[7,17]  = 1-sum(P.transition.matrix[7,c(1:16)])     
        P.transition.matrix[8,17]  = 1-sum(P.transition.matrix[8,c(1:16)])     
        P.transition.matrix[9,17]  = 1-sum(P.transition.matrix[9,c(1:16)])     
        P.transition.matrix[10,17] = 1-sum(P.transition.matrix[10,c(1:16)])     
        P.transition.matrix[11,17] = 1-sum(P.transition.matrix[11,c(1:16)])     
        P.transition.matrix[12,17] = 1-sum(P.transition.matrix[12,c(1:16)])     
        P.transition.matrix[13,17] = 1-sum(P.transition.matrix[13,c(1:16)])   
        P.transition.matrix[14,17] = 1-sum(P.transition.matrix[14,c(1:16)]) 
        P.transition.matrix[15,17] = 1-sum(P.transition.matrix[15,c(1:16)]) 
        P.transition.matrix[16,17] = 1-sum(P.transition.matrix[16,c(1:16)]) 
        P.transition.matrix[17,]   = c(rep(0,16),1)            
        if (j==(tau+1)) {              
            product.vec = initial.junction.prob %*% P.transition.matrix       
        }                
        if (j!=(tau+1)) {      
            product.vec = product.vec %*% P.transition.matrix
        }                       
    }                                 
    junction.surv.prob = 1-product.vec[17]   
    return(junction.surv.prob)
}   

# Equation 2
make.transition.matrix = function(freq.list,r,v,w,t.vec) {
    rownames = t.vec
    colnames = t.vec
    t.matrix = matrix(nrow=length(t.vec),ncol=length(t.vec))
    index.matrix.a = apply(t.matrix,1,function(X) X=t.vec)
    index.matrix.b = t(index.matrix.a)
    index.matrix = cbind(as.vector(index.matrix.a),as.vector(index.matrix.b))
    dimnames(t.matrix) = list(rownames,colnames)
    if (length(freq.list[[1]])==1) {
        pb = txtProgressBar(1,(length(t.vec)^2),1,style=3)                                        
        for (i in 2:(length(t.vec)^2)) {
            tau.vec = transition.rates.eq(freq.list,r,v,w,index.matrix[i,1],index.matrix[i,2])    
            t.matrix[index.matrix[i,1],index.matrix[i,2]] = mean(tau.vec) 
            setTxtProgressBar(pb, i)                                                            
        }  
        close(pb)                       
    }   
    if (length(freq.list[[1]])>1) {
        pb = txtProgressBar(1,(length(t.vec)^2),1,style=3)                                                         
        for (i in 2:(length(t.vec)^2)) {   
            tau.vec = transition.rates(freq.list,r,v,w,index.matrix[i,1],index.matrix[i,2])     
            t.matrix[index.matrix[i,1],index.matrix[i,2]] = mean(tau.vec)       
            setTxtProgressBar(pb, i)                                                            
        }   
        close(pb)                       
    }
    return(t.matrix)    
}    


junction.survival.func = function(trajec.list,r,v,w,gens.prob) { 
    gamete.trajecs = trajec.list[[1]]  
    allele.trajecs = trajec.list[[2]]
    zygote.w       = trajec.list[[3]]
    AABB.w = zygote.w[1,1]; AABb.w = zygote.w[1,2]; AAbb.w = zygote.w[1,3]
    AaBB.w = zygote.w[2,1]; AaBb.w = zygote.w[2,2]; Aabb.w = zygote.w[2,3]
    aaBB.w = zygote.w[3,1]; aaBb.w = zygote.w[3,2]; aabb.w = zygote.w[3,3]    
    A.allele = allele.trajecs[,1]           
    B.allele = 1-allele.trajecs[,2]         
    a.allele = 1-allele.trajecs[,1]         
    b.allele = allele.trajecs[,2]  
    starting.gen = (nrow(gamete.trajecs)-gens.prob)
    initial.junction.prob = c(rep(0,4),rep(2*A.allele[starting.gen]*b.allele[starting.gen],4),rep(2*a.allele[starting.gen]*B.allele[starting.gen],4),rep(0,5))         
    initial.junction.prob = c(rep(0,4),rep(c(gamete.trajecs[starting.gen,]),2),rep(0,5))*initial.junction.prob
    initial.junction.prob[17] = 1-sum(initial.junction.prob[1:16])
    for (j in 1:gens.prob) {                
        i = c(starting.gen:nrow(gamete.trajecs))[j]
        x.1 = gamete.trajecs[i,1]                 
        x.2 = gamete.trajecs[i,2]               
        x.3 = gamete.trajecs[i,3]               
        x.4 = gamete.trajecs[i,4]                                      
        P.transition.matrix = matrix(nrow=17,ncol=17)
        P.transition.matrix[1,c(1:16)]  = c(rep(AaBb.w,12),rep(0,4))*
                                          c(rep(.5,12),rep(0,4))*
                                          c(rep((1-r),4),rep(r,8),rep(1,4))*
                                          c(rep(1,4),rep((w/(v+w)),4),rep((v/(v+w)),4), rep(1,4))*
                                          c(rep(c(x.4,x.2,x.3,x.1),4))
        P.transition.matrix[2,c(1:16)]  = c(rep(AABb.w,8),rep(0,8))*
                                          c(rep(.5,8),rep(0,8))*
                                          c((((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,12))* 
                                          c(rep(1,4),rep(r,4),rep(1,8))*
                                          c(rep(1,4),rep((w/(v+w)),4),rep(1,8))*
                                          c(rep(1,4),c(x.4,x.2,x.3,x.1),rep(1,8))                                
        P.transition.matrix[3,c(1:16)]  = c(rep(AaBB.w,4),rep(0,4),rep(AaBB.w,4),rep(0,4))*
                                          c(rep(.5,4),rep(0,4),rep(.5,4),rep(0,4))*
                                          c((((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,12))* 
                                          c(rep(1,8),rep((r*(v/(v+w))),4),rep(1,4))*
                                          c(rep(1,8),c(x.4,x.2,x.3,x.1),rep(1,4))
        P.transition.matrix[4,c(1:16)]  = AABB.w*.5*c((((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(0,12)) 
        P.transition.matrix[5,c(1:16)]  = c(rep(0,4),rep(Aabb.w,4),rep(0,4),rep(Aabb.w,4))*
                                          c(rep(0,4),rep(.5,4),rep(0,4),rep(.5,4))*
                                          c(rep(1,4),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,8))* 
                                          c(rep(1,12),rep((r*(v/(v+w))),4))*
                                          c(rep(1,12),c(x.4,x.2,x.3,x.1))                                 
        P.transition.matrix[6,c(1:16)]  = c(rep(0,4),rep(AAbb.w,4),rep(0,8))*
                                          c(rep(0,4),rep(.5,4),rep(0,8))*
                                          c(rep(0,4),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(0,8)) 
        P.transition.matrix[7,c(1:16)]  = c(rep(AaBb.w,8),rep(0,4),rep(AaBb.w,4))*
                                          c(rep(.5,8),rep(0,4),rep(.5,4))*
                                          c(rep(r,4),rep((1-r),4),rep(1,4),rep(r,4))*
                                          c(rep((w/(v+w)),4),rep(1,8),rep((v/(v+w)),4))*
                                          c(rep(c(x.4,x.2,x.3,x.1),2),rep(1,4),c(x.4,x.2,x.3,x.1))
        P.transition.matrix[8,c(1:16)]  = c(rep(AABb.w,8),rep(0,8))*
                                          c(rep(.5,8),rep(0,8))*
                                          c(rep(r*(w/(v+w)),4),rep(1,12))*
                                          c(rep(1,4),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,8))*                                  
                                          c(c(x.4,x.2,x.3,x.1),rep(1,12))
        P.transition.matrix[9,c(1:16)]  = c(rep(0,8),rep(aaBb.w,8))*
                                          c(rep(0,8),rep(.5,8))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)),rep(1,4))*
                                          c(rep(1,12),rep(r*(w/(v+w)),4))*
                                          c(rep(1,12),c(x.4,x.2,x.3,x.1))
        P.transition.matrix[10,c(1:16)] = c(rep(AaBb.w,4),rep(0,4),rep(AaBb.w,8))*
                                          c(rep(.5,4),rep(0,4),rep(.5,8))*
                                          c(rep(r,4),rep(1,4),rep((1-r),4),rep(r,4))*
                                          c(rep(v/(v+w),4),rep(1,8),rep(w/(v+w),4))*
                                          c(c(x.4,x.2,x.3,x.1),rep(1,4),rep(c(x.4,x.2,x.3,x.1),2))
        P.transition.matrix[11,c(1:16)] = c(rep(0,8),rep(aaBB.w,4),rep(0,4))*
                                          c(rep(0,8),rep(.5,4),rep(0,4))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)), rep(1,4))                         
        P.transition.matrix[12,c(1:16)] = c(rep(AaBB.w,4),rep(0,4),rep(AaBB.w,4),rep(0,4))*
                                          c(rep(.5,4),rep(0,4),rep(.5,4),rep(0,4))*
                                          c(rep(r,4),rep(1,12))*
                                          c(rep(v/(v+w),4),rep(1,12))*
                                          c(c(x.4,x.2,x.3,x.1),rep(1,12))*
                                          c(rep(1,8),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)),rep(1,4))                                  
        P.transition.matrix[13,c(1:16)] = c(rep(0,12),rep(aabb.w,4))*
                                          c(rep(0,12),rep(.5,4))*
                                          c(rep(1,12),(((1-r)*x.4) + (r*x.4)), (((1-r)*x.2) + (r*x.2)), (((1-r)*x.3) + (r*x.3)), (((1-r)*x.1) + (r*x.1)))                         
        P.transition.matrix[14,c(1:16)] = c(rep(0,4),rep(Aabb.w,4),rep(0,4),rep(Aabb.w,4))*
                                          c(rep(0,4),rep(.5,4),rep(0,4),rep(.5,4))*
                                          c(rep(1,4),rep(r*(v/(v+w)),4),rep(1,8))*
                                          c(rep(1,4),c(x.4,x.2,x.3,x.1),rep(1,8))*                                  	
                                          c(rep(1,12),(((1-r)*x.4) + (r*(w/(v+w))*x.4)), (((1-r)*x.2) + (r*(w/(v+w))*x.2)), (((1-r)*x.3) + (r*(w/(v+w))*x.3)), (((1-r)*x.1) + (r*(w/(v+w))*x.1)))                                  
        P.transition.matrix[15,c(1:16)] = c(rep(0,8),rep(aaBb.w,8))*
                                          c(rep(0,8),rep(.5,8))*
                                          c(rep(1,8),rep(r*(w/(v+w)),4),rep(1,4))*
                                          c(rep(1,8),c(x.4,x.2,x.3,x.1),rep(1,4))*
                                          c(rep(1,12),(((1-r)*x.4) + (r*(v/(v+w))*x.4)), (((1-r)*x.2) + (r*(v/(v+w))*x.2)), (((1-r)*x.3) + (r*(v/(v+w))*x.3)), (((1-r)*x.1) + (r*(v/(v+w))*x.1)))
        P.transition.matrix[16,c(1:16)] = c(rep(0,4),rep(AaBb.w,12))*
                                          c(rep(0,4),rep(.5,12))*
                                          c(rep(1,4),rep(r,8),rep((1-r),4))*
                                          c(rep(1,4),rep(v/(v+w),4),rep(w/(v+w),4),rep(1,4))*
                                          c(rep(1,4),rep(c(x.4,x.2,x.3,x.1),3))
        P.transition.matrix[1,17]  = 1-sum(P.transition.matrix[1,c(1:16)])
        P.transition.matrix[2,17]  = 1-sum(P.transition.matrix[2,c(1:16)])   
        P.transition.matrix[3,17]  = 1-sum(P.transition.matrix[3,c(1:16)]) 
        P.transition.matrix[4,17]  = 1-sum(P.transition.matrix[4,c(1:16)])         
        P.transition.matrix[5,17]  = 1-sum(P.transition.matrix[5,c(1:16)])         
        P.transition.matrix[6,17]  = 1-sum(P.transition.matrix[6,c(1:16)])     
        P.transition.matrix[7,17]  = 1-sum(P.transition.matrix[7,c(1:16)])     
        P.transition.matrix[8,17]  = 1-sum(P.transition.matrix[8,c(1:16)])     
        P.transition.matrix[9,17]  = 1-sum(P.transition.matrix[9,c(1:16)])     
        P.transition.matrix[10,17] = 1-sum(P.transition.matrix[10,c(1:16)])     
        P.transition.matrix[11,17] = 1-sum(P.transition.matrix[11,c(1:16)])     
        P.transition.matrix[12,17] = 1-sum(P.transition.matrix[12,c(1:16)])     
        P.transition.matrix[13,17] = 1-sum(P.transition.matrix[13,c(1:16)])   
        P.transition.matrix[14,17] = 1-sum(P.transition.matrix[14,c(1:16)]) 
        P.transition.matrix[15,17] = 1-sum(P.transition.matrix[15,c(1:16)]) 
        P.transition.matrix[16,17] = 1-sum(P.transition.matrix[16,c(1:16)]) 
        P.transition.matrix[17,]   = c(rep(0,16),1)            
        if (j==1) {              
            product.vec = initial.junction.prob %*% P.transition.matrix       
        }                
        if (j!=1) {      
            product.vec = product.vec %*% P.transition.matrix
        }                       
    }                                 
    junction.surv.prob = 1-product.vec[17]   
    junction.surv.prob
    return(junction.surv.prob)
}        


survival.probs.WF = function(trajec.list,r,v,w,gens.prob) {
    A.trajecs  = trajec.list$A.matrix.a   
    if (!is.null(ncol(A.trajecs))) {
        AB.trajecs = trajec.list$AB.matrix.a
        Ab.trajecs = trajec.list$Ab.matrix.a
        aB.trajecs = trajec.list$aB.matrix.a
        ab.trajecs = trajec.list$ab.matrix.a
        A.trajecs  = trajec.list$A.matrix.a
        b.trajecs  = trajec.list$b.matrix.a
        zygote.w   = trajec.list$zygote.w
        survival.prob.vec.WF = c(1:ncol(A.trajecs))    
        for (i in 1:ncol(A.trajecs)) {
            gamete.trajecs = cbind(AB.trajecs[,i],Ab.trajecs[,i],aB.trajecs[,i],ab.trajecs[,i])
            allele.trajecs = cbind(A.trajecs[,i],b.trajecs[,i])
            WF.list = list("gamete.trajecs"=gamete.trajecs,"allele.trajecs"=allele.trajecs,"zygote.w"=zygote.w)
            survival.prob.vec.WF[i] = junction.survival.func(WF.list,r,v,w,gens.prob)
        }       
    }           
    if (is.null(ncol(A.trajecs))) {    
        AB.trajecs = trajec.list$AB.matrix.a   
        Ab.trajecs = trajec.list$Ab.matrix.a   
        aB.trajecs = trajec.list$aB.matrix.a   
        ab.trajecs = trajec.list$ab.matrix.a   
        A.trajecs  = trajec.list$A.matrix.a  
        b.trajecs  = trajec.list$b.matrix.a   
        zygote.w   = trajec.list$zygote.w   
        gamete.trajecs = cbind(AB.trajecs,Ab.trajecs,aB.trajecs,ab.trajecs)   
        allele.trajecs = cbind(A.trajecs,b.trajecs)   
        traj.list = list("gamete.trajecs"=gamete.trajecs,"allele.trajecs"=allele.trajecs,"zygote.w"=zygote.w) 
        survival.prob.vec.WF = junction.survival.func(traj.list,r,v,w,gens.prob)   
    }   
    return(survival.prob.vec.WF)
}
 

plot.survival.prob.versus.gens = function(WF.det.list,WF.sto.list,r,v,w) {  
    A.trajecs  = WF.sto.list$A.matrix.a   
    survival.time.matrix = matrix(ncol=nrow(A.trajecs),nrow=ncol(A.trajecs))  
    deterministic.vector = c(1:nrow(A.trajecs)) 
    pb = txtProgressBar(1,nrow(A.trajecs),1,style=3)    
    for (i in 1:100) {     
        gens.prob = i   
        survival.time.matrix[,i] = survival.probs.WF(WF.sto.list,r,v,w,gens.prob)    
        deterministic.vector[i] =  survival.probs.WF(WF.det.list,r,v,w,gens.prob)    
        setTxtProgressBar(pb, i)         
    }    
    close(pb)
    matplot(t(survival.time.matrix[,c(1:10)]),type='l',col="firebrick",lty=1,lwd=2)
    lines(deterministic.vector[c(1:10)],lwd=2,lty=2)
    return(survival.time.matrix)                                                       
}



     

