
#---------------Gamete frequency trajectory simulations-----------------

pops = 1        # number of subpopulations per species
N = 100000      # number of individuals per deme      
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
source("gamete_trajectories.r")
make.fitnessmatrix(alpha,alpha1,beta,beta1,s)
WF.det.list = simulate.trajecs.hs(N,m,r,t,alpha,alpha1,beta,beta1,s,source.freqs,starting.freqs) 
WF.sto.list = WF.trajec.reps.hs(N,m,r,t,alpha,alpha1,beta,beta1,s,source.freqs,starting.freqs,reps)
plot.WF.versus.deterministic.hs(WF.det.list,WF.sto.list)

#------------Frequency Trajectories--------------
source("gamete_trajectories.r")
WF.det.list = simulate.trajecs(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs) 
WF.sto.list = WF.trajec.reps(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs,reps)
plot.WF.versus.deterministic(WF.det.list,WF.sto.list)

#------------Survival Probabilities--------------
source("gamete_trajectories.r")
source("junction_survival_probabilities.r")
WF.det.list = simulate.trajecs(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs) 
WF.sto.list = WF.trajec.reps(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs,reps)
survival.time.matrix.WF = plot.survival.prob.versus.gens(WF.det.list,WF.sto.list,r,v,w)
        
#------------Transition Matrix (Eq. 2)-----------
source("gamete_trajectories.r")
source("junction_survival_probabilities.r")
WF.det.list = simulate.trajecs(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs) 
WF.sto.list = WF.trajec.reps(pops,N,m,msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs,reps)
eq.freq.list = find.equilibrium.freqs(pops,N,m,msub,r,alpha,alpha1,beta,beta1,s,source.freqs)
t.matrix.eq = make.transition.matrix(eq.freq.list,r,v,w,t.vec)

#---------------------dfuse Versus Tracts, No Selection------------------------

pars = list(
    "file.name" = "dmi",          # how to label the input and output files
    "demes.n" = 1,                # number of demes
    "demes.N" = 5000,             # adult and juvenile carrying capacity
    "mig.rate" = NULL,            # migration rate
    "reps.n" = 1,                 # number of replicates
    "gens.n" = 50,                # number of generations
    "sel.coef" = 0,               # selection coefficient
    "gen.print" = 50,             # print every x generations
    "markers.n" = 800,            # marker loci per chromosome
    "dmi.positions" = c(0.3,0.7), # positions of the dmi loci
    "bin.size" = .05,             # Morgan size to bin the tract lengths
    "final.gen" = TRUE,
    "color"     = NULL,
    "dfuse.dir" = "/home/joelsmith/Projects/dmis/code/dfuse_files/")      


plot.tracts.versus.dfuse = function(pars) {
    # Runs simulation with dfuse to get tract length distribution.
    source(paste(pars$dfuse.dir,"dfuse_functions.r",sep=''))
    results.list = run.dfuse(pars)
    bin.list = plot.tract.length.distr(results.list,pars)
    # Computes the expected tract length distribution under the Gravel model.
    source("/home/joelsmith/Projects/dmis/code/gravel_transition_matrix.r")     
    t.vec = c(1:pars$gens.n)  
    x.vec = seq(1,(1/pars$bin.size),1)
    m.vec = c(rep(pars$mig.rate,(length(t.vec)))) 
    m.vec[length(m.vec)] = 0.5
    initial.probs = c(rep((1/length(t.vec)),(length(t.vec))),rep(0,length(t.vec)))
    tracts = run.tract.length.functions(m.vec,t.vec,x.vec,initial.probs)
    tracts = cbind((pars$bin.size*x.vec)-(pars$bin.size/2),(tracts/sum(tracts)))
    lines(tracts,lty=2,lwd=2,col=pars$color)
}

#pdf("simulated_versus_predicted_variable_migration.pdf")
migration.vec = rev(c(.001,.1))
library("RColorBrewer")     
colors = brewer.pal(5,'Set1')[c(2,5,1,3)]   
for (i in c(1:length(migration.vec))) {
    pars$mig.rate = migration.vec[i]
    pars$color = colors[i]
    plot.tracts.versus.dfuse(pars)
    par(new=T)
}
legend("topright",c(as.character(rev(c(migration.vec))),"Predicted"),cex=1.4,col=c(colors[1:length(migration.vec)],'black'),title="Migration Rate",lwd=3,lty=c(1,1,2))
#dev.off()

#-----------------Tracts Pulse Migration Plots---------------------

source("gravel_transition_matrix.r")                     
t.vec = c(1:30)
x.vec = seq(1,20,1)   # tract length bins.
initial.probs = c(rep((1/length(t.vec)),(length(t.vec))),rep(0,length(t.vec)))

mig.a = c(rep(0,(length(t.vec)-1)),.6) 
mig.b = c(rep(0,(length(t.vec)-1)),.2) 
mig.c = c(rep(0,((length(t.vec)/2)-1)),.6,rep(0,((length(t.vec)/2)))) 
mig.d = c(rep(0,((length(t.vec)/2)-1)),.2,rep(0,((length(t.vec)/2)))) 
    
tracts.a = run.tract.length.functions(mig.a,t.vec,x.vec,initial.probs)
tracts.b = run.tract.length.functions(mig.b,t.vec,x.vec,initial.probs)
tracts.c = run.tract.length.functions(mig.c,t.vec,x.vec,initial.probs)
tracts.d = run.tract.length.functions(mig.d,t.vec,x.vec,initial.probs)
#tracts.a = tracts.a/sum(tracts.a)        
#tracts.b = tracts.b/sum(tracts.b)        
#tracts.c = tracts.c/sum(tracts.c)        
#tracts.d = tracts.d/sum(tracts.d)                   

#pdf("high_low_recent_old_comparisons.pdf")
library("RColorBrewer")     
colors = brewer.pal(5,'Set1')[c(2,5)]     
plot(cbind(((1/length(x.vec))*x.vec),tracts.a),typ='l',lwd=2.5,col=colors[1],ylim=c(0,.3),xlim=c(0,.4),cex.lab=1.2,cex.main=1.2,   
    xlab="Tract Length (Morgans)",ylab="Frequency",main="Tract Length Distributions")  
lines(cbind(((1/length(x.vec))*x.vec),tracts.b),lwd=2.5,col=colors[2])       
lines(cbind(((1/length(x.vec))*x.vec),tracts.c),lwd=2.5,col=colors[1],lty=2)   
lines(cbind(((1/length(x.vec))*x.vec),tracts.d),lwd=2.5,col=colors[2],lty=2)    
legend("topright",c("Old/Low","Old/High","Recent/Low","Recent/High"),   
               col=c(colors[1:2],colors[1:2]),    
               lty=c(1,1,2,2),lwd=2.3)          
#dev.off()
#--------------Tracts Continuous Migration Plots----------------
                  
source("tracts_js.r")                                         
t.vec = c(1:30)
x.vec = seq(1,20,1) 
initial.probs = c(rep(0,(length(t.vec)-1)),1,rep(0,length(t.vec)))
 
mig.a = rep(0.001,length(t.vec))
mig.a[length(mig.a)] = 1-mig.a[length(mig.a)]
mig.b = rep(0.03,length(t.vec))
mig.b[length(mig.b)] = 1-mig.b[length(mig.b)]
mig.c = rep(0.05,length(t.vec))
mig.c[length(mig.c)] = 1-mig.c[length(mig.c)]

tracts.a = run.tract.length.functions(mig.a,t.vec,x.vec,initial.probs)
tracts.b = run.tract.length.functions(mig.b,t.vec,x.vec,initial.probs)
tracts.c = run.tract.length.functions(mig.c,t.vec,x.vec,initial.probs)

tracts.a = tracts.a/sum(tracts.a)
tracts.b = tracts.b/sum(tracts.b)
tracts.c = tracts.c/sum(tracts.c)

tracts.a = log(tracts.a)
tracts.b = log(tracts.b)
tracts.c = log(tracts.c)

library("RColorBrewer")
colors = brewer.pal(6, "Set1")[c(1,2,5)]
plot(tracts.a,typ="l",lwd=3,col=colors[2])
lines(tracts.b,lwd=3,col=colors[1])
lines(tracts.c,lwd=3,col=colors[3])


#--------------Tracts Continuous Migration Plots with edges----------------
     
source("tract_length_distributions.r")   
t.vec = c(1:30)                                                                           
m.vec = rep(0.001,length(t.vec))   
#m.vec[length(m.vec)] = 1-m.vec[length(m.vec)]
x.vec = seq(1,40,1)    
initial.probs = c(rep(0,(length(t.vec)-1)),1,rep(0,length(t.vec)))   
    
tract.length.distr.list = run.tract.length.functions.with.edges(m.vec, t.vec, x.vec, initial.probs)
               
tract.length.distr.i = log(tract.length.distr.list[[1]]/sum(tract.length.distr.list[[1]]))
tract.length.distr.e = log(tract.length.distr.list[[2]]/sum(tract.length.distr.list[[2]]))
tract.length.distr.f = log(tract.length.distr.list[[3]]/sum(tract.length.distr.list[[3]]))
tract.length.distr.g = log(tract.length.distr.list[[4]]/sum(tract.length.distr.list[[4]]))

plot(tract.length.distr.g,typ='l',lwd=3,col="cadetblue",ylab='Relative Frequency',xlab='Morgans')  
lines(tract.length.distr.e,lwd=3,col='darkolivegreen')
legend('topright',legend = c("no edges","edges"),col=c("cadetblue","darkolivegreen"),lty=1,lwd=3)  
    
#lines(tract.length.distr.f,lwd=3,col='blue')   
#lines(tract.length.distr.i,lwd=3,col='orange') 

   
#------Plots Tract Length Distributions-------     
source("../dfuse/dmi_sim.R")
pars = list(
    "file.name" = "dmi",          # how to label the input and output files
    "demes.n" = 1,                # number of demes
    "demes.N" = 1000,             # adult and juvenile carrying capacity
    "mig.rate" = 0,               # migration rate
    "reps.n" = 1,                 # number of replicates
    "gens.n" = 75,                # number of generations
    "sel.coef" = 0,               # selection coefficient
    "gen.print" = 75,             # print every x generations
    "markers.n" = 500,            # marker loci per chromosome
    "dmi.positions" = c(0.3,0.7), # positions of the dmi loci
    "bin.size" = .05,             # Morgan size to bin the tract lengths
    "final.gen" = TRUE)            
                                   
results.list = run.dfuse(pars)
#pdf("tract_length_distr.pdf",height=8,width=10)
plot.tract.length.distr(results.list,pars)
par(new=TRUE)
pars$mig.rate=.1
results.list = run.dfuse(pars)
plot.tract.length.distr(results.list,pars)
#plot.individuals.junctions(results.list,pars)
#dev.off()

#--------------------forqs observed versus predicted tracts---------------------------  

source("/home/joelsmith/Projects/dmis/code/tracts_js.r")
source("/home/joelsmith/Projects/dmis/code/forqs_files/forqs_functions.r")                                         
params = list(
    "population.size" = 2500,
    "migration.rate" = .03,
    "generations" = 30,
    "chromosome.length" = 100000000,
    "output.file" = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_out/forqs_output_source_target",
    "bin.size" = 5000000)

mig.vec = c(.001,.03,.05)  
t.vec = c(1:30) 
x.vec = seq(1,params$chromosome.length/params$bin.size,1) 
initial.probs = c(rep(0,(length(t.vec)-1)),1,rep(0,length(t.vec))) 
   
mig.a = rep((1-mig.vec[1]),length(t.vec))
mig.b = rep((1-mig.vec[2]),length(t.vec))
mig.c = rep((1-mig.vec[3]),length(t.vec))

distr.a = run.tract.length.functions.with.edges(mig.a,t.vec,x.vec,initial.probs)
distr.b = run.tract.length.functions.with.edges(mig.b,t.vec,x.vec,initial.probs)
distr.c = run.tract.length.functions.with.edges(mig.c,t.vec,x.vec,initial.probs)
forqs.list = plot.multiple.simulations(params,mig.vec)    
total.edges.a = distr.a[[1]]+distr.a[[2]]+distr.a[[3]]
total.edges.b = distr.b[[1]]+distr.b[[2]]+distr.b[[3]]
total.edges.c = distr.c[[1]]+distr.c[[2]]+distr.c[[3]]
tracts.a = total.edges.a/sum(total.edges.a)
tracts.b = total.edges.b/sum(total.edges.b)
tracts.c = total.edges.c/sum(total.edges.c)

tracts.a.g = log(distr.a[[4]]/sum(distr.a[[4]]))
tracts.b.g = log(distr.b[[4]]/sum(distr.b[[4]]))
tracts.c.g = log(distr.c[[4]]/sum(distr.c[[4]]))
forqs.bins.a = log(forqs.list[[1]]/sum(forqs.list[[1]]))
forqs.bins.b = log(forqs.list[[2]]/sum(forqs.list[[2]]))
forqs.bins.c = log(forqs.list[[3]]/sum(forqs.list[[3]]))
plot(tracts.a.g,col='orange',typ='l',lty=1,lwd=2,ylim=c(-11,0),ylab='Relative Frequency')  
lines(tracts.b.g,col='red',lty=1,lwd=2)    
lines(tracts.c.g,col='blue',lty=1,lwd=2) 
points(forqs.bins.a,col='blue',pch=19,lwd=2)    
points(forqs.bins.b,col='red',pch=19,lwd=2) 
points(forqs.bins.c,col='orange',pch=19,lwd=2)    
legend("topright", title='Migration Rate',
        legend=c("0.001 Forqs","0.03   Forqs","0.05   Forqs",
                 paste(mig.vec[1]," Predicted"),
                 paste(mig.vec[2],"   Predicted"),
                 paste(mig.vec[3],"   Predicted")), 
       col=c(rep(c("blue","red","orange"),3)), 
       pch=c(20,20,20,NA,NA,NA),
       lty=c(NA,NA,NA,1,1,1),
       lwd=c(2,2,2,2,2,2),
       bty="n", border=F, ncol=2)  

#---------------------Forqs plots--------------------------------------------
source("forqs_functions.r")

params = list(  
    "population.size" = 500,
    "generations" = 30,
    "chromosome.length" = 1000000,
    "output.file" = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_out/forqs_output_source_target",
    "bin.size" = 50000)

mig.vec = c(.001,.03,.05)

run.multiple.simulations(params,mig.vec)
plot.multiple.simulations(params,mig.vec)

plot.chroms(chroms.matrix)

#---------------------dfuse tract lengths around DMIs------------------------

source("../dfuse/dmi_sim.R")
pars = list(
    "file.name" = "dmi",          # how to label the input and output files
    "demes.n" = 1,                # number of demes
    "demes.N" = 500,             # adult and juvenile carrying capacity
    "mig.rate" = .01,             # migration rate
    "reps.n" = 20,                 # number of replicates
    "gens.n" = 30,                # number of generations
    "sel.coef" = .5,              # selection coefficient
    "gen.print" = 30,             # print every x generations
    "markers.n" = 200,            # marker loci per chromosome
    "dmi.positions" = c(0.3,0.7), # positions of the dmi loci
    "bin.size" = .05,              # Morgan size to bin the tract lengths
    "final.gen" = TRUE,           # Only plots last generation
    "color"     = NULL)           # Specifies color for plot     
#pdf("tract_length_distr_dmi_neutral.pdf",height=8,width=10)
plot.dmi.tract.lengths(pars)
#dev.off()
     
#------Plots Haplotypes----- 
source("../dfuse/dmi_sim.R")                                       
pars = list(
    "file.name" = "dmi",          # how to label the input and output files
    "demes.n" = 1,                # number of demes
    "demes.N" = 100,              # adult and juvenile carrying capacity
    "mig.rate" = .01,               # migration rate
    "reps.n" = 1,                 # number of replicates
    "gens.n" = 30,                # number of generations
    "sel.coef" = 0,              # selection coefficient
    "gen.print" = 30,             # print every x generations
    "markers.n" = 100,            # marker loci per chromosome
    "dmi.positions" = c(0.3,0.7), # positions of the dmi loci
    "final.gen"= TRUE) 
              
results.list = run.dfuse(pars)
#pdf("haplotypes_neutral.pdf",height=8,width=10)
plot.individuals.junctions(results.list,pars)
#dev.off()
    
    
#------Plots Junction Counts------- 
source("../dfuse/dmi_sim.R")                                                               
pars = list(
    "file.name" = "dmi",          # how to label the input and output files
    "demes.n" = 1,                # number of demes
    "demes.N" = 100,              # adult and juvenile carrying capacity
    "mig.rate" = 0.01,            # migration rate
    "reps.n" = 20,                # number of replicates
    "gens.n" = 30,               # number of generations
    "sel.coef" = .5,              # selection coefficient
    "gen.print" = 10,             # print every x generations
    "markers.n" = 100,            # marker loci per chromosome
    "dmi.positions" = c(0.3,0.7)) # positions of the dmi loci

results.list = run.dfuse(pars)
junction.list = get.junction.densities(results.list)
#pdf("junctions.pdf",height=8,width=10)
plot.junction.densities(results.list,junction.list,pars)
#dev.off()

#plot.junction.sums(junction.list,results.list,pars)

#------Plots Ancestry Proportions-------   

source("../dfuse/dmi_sim.R")                                                             
pars = list(
    "file.name" = "dmi",          # how to label the input and output files
    "demes.n" = 1,                # number of demes
    "demes.N" = 100,              # adult and juvenile carrying capacity
    "mig.rate" = 0.01,            # migration rate
    "reps.n" = 20,                # number of replicates
    "gens.n" = 30,               # number of generations
    "sel.coef" = .5,              # selection coefficient
    "gen.print" = 10,             # print every x generations
    "markers.n" = 100,            # marker loci per chromosome
    "dmi.positions" = c(0.3,0.7)) # positions of the dmi loci

results.list = run.dfuse(pars)
ancestry.list = get.ancestry.proportions(results.list)
#pdf("ancestry.pdf",height=8,width=10)
plot.ancestry.proportions(results.list,ancestry.list,pars)
#dev.off()

#-------- Comparison with Gavrilets analytical results

gavrilets.solutions = function(m,msub,alpha,beta,r,s) {
    AB.star = m/(s+(1-s)*r)
    Ab.star = (((1-s)*r)/((s+(1-s)*r)))*(m/beta)
    aB.star.a = (((1-s)*r)/(s+((1-s)*r)))*(m/alpha)
    aB.star.b = (((1-s)*r)/(s+((1-s)*r)))*(m/msub)
    ab.star.a = 1-AB.star-Ab.star-aB.star.a
    ab.star.b = 1-AB.star-Ab.star-aB.star.b
    gamete.freqs = c(AB.star, Ab.star, aB.star.a, ab.star.a, aB.star.b, ab.star.b)
    allele.freqs = c(sum(gamete.freqs[c(1,2)]),sum(gamete.freqs[c(2,4)]),sum(gamete.freqs[c(2,6)]))
    star.list = list('gamete.freqs'=gamete.freqs,'allele.freqs'=allele.freqs)
    return(star.list)
}
  
m.vec = seq(0.0,0.1,0.001)
species.a.freqs = matrix(nrow=length(m.vec),ncol=2) 
species.b.freqs = matrix(nrow=length(m.vec),ncol=2)
gavrilets.freqs = matrix(nrow=length(m.vec),ncol=3)
for (i in 1:length(m.vec)) {
    trajec.list = simulate.trajecs(pops,N,m.vec[i],msub,r,t,alpha,alpha1,beta,beta1,s,source.freqs)  
    star.list =  gavrilets.solutions(m.vec[i],msub,alpha,beta,r,s) 
    gamete.trajecs = trajec.list[[1]]                                  
    allele.freqs = trajec.list[[2]]
    equilibrium.freqs = allele.freqs[nrow(allele.freqs),]
    species.a.freqs[i,] = equilibrium.freqs[c(1,2)]
    species.b.freqs[i,] = equilibrium.freqs[c(3,4)]
    gavrilets.freqs[i,] = star.list[[2]]
}

#pdf('equilibrium_freqs.pdf',height=7,width=7)
library('RColorBrewer')
colors = brewer.pal(4,'Dark2')   
plot(m.vec,species.b.freqs[,1],col=colors[2],typ='l',main='Equilibrium Frequencies',lwd=2,ylim=c(0,1),xlab='m',ylab='frequency')
lines(m.vec,species.b.freqs[,2],col=colors[4],lwd=2)
lines(m.vec,gavrilets.freqs[,1],col=colors[2],lwd=2,lty=3)
lines(m.vec,gavrilets.freqs[,3],col=colors[4],lwd=2,lty=3)
legend(.001,.6,c('A numerical','	b numerical','A analytical','b analytical'),lty=c(1,1,3,3),lwd=2,col=c(colors[c(2,4)],colors[c(2,4)])) 
#dev.off()


#----------Erlang and Gamma check for the tracts model------------

derlang = function(x,k,lambda) {
    return(((lambda^k)*(x^(k-1))*(exp(-lambda*x)))/factorial(k-1))
}

x = 2
k = 2
lambda = 10

dgamma(x,k,lambda)
derlang(x,k,lambda)

#-------------Power plots for tract lengths-----------------------

tractlength.powerplot = function(samplesize.vec,background.tract.length,alpha.vec,reps) {
    for (i in 1:length(samplesize.vec)) {
        print(samplesize.vec[i])
        rate.bg = 1/background.tract.length
        n = samplesize.vec[i]
        alpha.matrix = matrix(nrow=reps,ncol=length(alpha.vec))
        pb = txtProgressBar(1,length(alpha.vec),1,style=3) 
        for (j in 1:length(alpha.vec)) {    
            rate.al = rate.bg/alpha.vec[j]
            for (k in 1:reps) {
                tract.lengths.al = rexp(n,rate.al)
                tract.lengths.bg = rexp(n,rate.bg)                       
                mle.al = 1/mean(tract.lengths.al)
                mle.bg = 1/mean(tract.lengths.bg)
                likelihood.al = sum(dexp(tract.lengths.al,rate=mle.al,log=T))
                likelihood.bg = sum(dexp(tract.lengths.al,rate=mle.bg,log=T))
                likelihood.ratio = 2*(likelihood.al-likelihood.bg)    
                alpha.matrix[k,j] = pchisq(likelihood.ratio,df=1)           
            }
            setTxtProgressBar(pb, j)            
        }
        close(pb)  
        power.vec = apply(alpha.matrix,2,function(X) length(which(X>=0.95)))/reps  
        library('RColorBrewer')
        colors = brewer.pal(8,'Set1')                                                    
        if (i==1) {      
            plot(alpha.vec,power.vec,col=colors[i],ylab='power',xlab='alpha',ylim=c(0,1),pch=20)
        } else {
            points(alpha.vec,power.vec,col=colors[i],pch=20)
        }
    }   
    legend(1.3,.6,c("10","100","1000"),col=(colors[1:3]),pch=15)
}

samplesize.vec          = c(10,100,1000)     # the number of tracts
background.tract.length = 50000              # null tract length
alpha.vec               = seq(1.01,1.5,0.01) # a vector of values to adjust the alternative tract length. 
reps                    = 10000              # the number of datasets to simulate for each alpha.

tractlength.powerplot(samplesize.vec,background.tract.length,alpha.vec,reps)


#-------------Power plots for number of junctions-----------------------

junctions.powerplot = function(samplesize.vec,background.junction.rate,region.bp,alpha.vec,reps) {
    for (i in 1:length(samplesize.vec)) {
        print(samplesize.vec[i])
        n = samplesize.vec[i]
        rate.bg = background.junction.rate
        alpha.matrix = matrix(nrow=reps,ncol=length(alpha.vec))
        pb = txtProgressBar(1,length(alpha.vec),1,style=3) 
        for (j in 1:length(alpha.vec)) {    
            rate.al = rate.bg/alpha.vec[j]
            for (k in 1:reps) {
                junctions.al = rbinom(n,region.bp,prob=rate.al)
                junctions.bg = rbinom(n,region.bp,prob=rate.bg)                       
                mle.al = mean(junctions.al)/region.bp
                mle.bg = mean(junctions.bg)/region.bp
                likelihood.al = sum(dbinom(junctions.al,size=region.bp,prob=mle.al,log=T))
                likelihood.bg = sum(dbinom(junctions.al,size=region.bp,prob=mle.bg,log=T))
                likelihood.ratio = 2*(likelihood.al-likelihood.bg)    
                alpha.matrix[k,j] = pchisq(likelihood.ratio,df=1)           
            }
            setTxtProgressBar(pb, j)            
        }
        close(pb)  
        power.vec = apply(alpha.matrix,2,function(X) length(which(X>=0.95)))/reps  
        library('RColorBrewer')
        colors = brewer.pal(8,'Set1')                                                                                                                 
        if (i==1) {      
            plot(alpha.vec,power.vec,col=colors[i],ylab='power',xlab='alpha',ylim=c(0,1),pch=20)
        } else {
            points(alpha.vec,power.vec,col=colors[i],pch=20)
        }
    }
    legend(1.03,.6,c("10","100","1000"),col=(colors[1:3]),pch=15)
}


samplesize.vec = c(10,100,1000)   # The number of junctions.
region.bp = 10000                 # The size of the region. 
background.junction.rate = .1     # The null rate of junctions observed per basepair
alpha.vec = seq(1.001,1.05,0.001) # A vector of values to adjust the alternative junction rate.
reps = 10000                      # The number of datasets to simulate for each alpha.

junctions.powerplot(samplesize.vec,background.junction.rate,region.bp,alpha.vec,reps)



