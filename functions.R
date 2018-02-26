
#-----------------simulation parameters------------------

popsize = 500   
chrlen = 1000  
dmipos = c(.5,.5)  
windowsize = .2                    
chrnum = 2
reps = 1
bins = 10
migrate = .1  
nmarkers = 1000   
tstart = 100             
model = "pathway"
type = "boxplots"
demefile = "input/demefile"        
episfile = "input/episfile"        
outfile = "output/dmi"               
outfile.a = "output/neu"                      
outfile.b = "output/sel"     

# null model
seldmi = 0 
selanc = 0   
domdmi = 0  

# model (a)
seldmi = .9
selanc = .5  
domdmi = 0 

# model (b)
seldmi = .9  
selanc = 0   
domdmi = 0 

# model (c)
seldmi = .9  
selanc = 0   
domdmi = 3 

# models (d - f)
seldmi = .9  
selanc = 0  
domdmi = 0

                          
#---------------------run simulations--------------------
 
run_dfuse(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi,selanc,domdmi,reps,demefile,episfile,outfile,chrnum,model)    
                 
#--------------------plot haplotypes---------------------           

plot_haplotypes(outfile, dmipos, chrnum)      
                                                            
#----------------compute summary statistics--------------                                        
                                         
compute_summary_stats(outfile, dmipos, chrnum, windowsize)                                              
    
#-----------------plot summary statistics----------------   

plot_summary_stats(outfile, dmipos, windowsize, reps)
                  
#----------------run simulation comparison---------------                       
                      
run_comparisons(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi,selanc,domdmi,reps,demefile,episfile,outfile.a,outfile.b,chrnum,model)                     
                                                              
#-----plot neutral versus selection distributions--------             

plot_distributions(outfile.a,outfile.b,reps,bins,type)   

#---------------compare distributions--------------------

compare_distributions(outfile.a,outfile.b,reps) 
        
#---------------find time to equilibrium-----------------

find_equilibrium(popsize,chrlen,migrate,nmarkers,dmipos,seldmi,selanc,domdmi,demefile,episfile,outfile,chrnum,model)

#--------------------------------------------------------    

compare_distributions = function(outfile.a,outfile.b,reps) {
    recover()        
    library(entropy)
}


find_equilibrium = function(popsize,chrlen,migrate,nmarkers,dmipos,seldmi,selanc,domdmi,demefile,episfile,outfile,chrnum,model) {
    ancestry.a = 0
    ancestry.b = 0
    junctions.a = 0
    junctions.b = 0
    anc.diff.table = c(0,0,0,0)
    for (tstart in seq(1,150,1)) {
        print(tstart)
        run_dfuse(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi,selanc,domdmi,reps=1,demefile,episfile,outfile,chrnum,model)
        compute_summary_stats(outfile, dmipos, chrnum, windowsize)   
        load(paste(outfile,".stats.RDATA",sep=''))
        ancestry.a.prime = summarylist$ancestry$mean.a
        ancestry.b.prime = summarylist$ancestry$mean.b
        junctions.a.prime = summarylist$junctions$mean.a
        junctions.b.prime = summarylist$junctions$mean.b
        anc.a.diff = abs(ancestry.a.prime-ancestry.a)
        anc.b.diff = abs(ancestry.b.prime-ancestry.b)
        jun.a.diff = abs(junctions.a.prime-junctions.a)
        jun.b.diff = abs(junctions.b.prime-junctions.b)
        anc.diff.table = rbind(anc.diff.table,c(anc.a.diff,anc.b.diff,jun.a.diff,jun.b.diff))
    }
    diff.table =  anc.diff.table[-1,]
    par(mfrow=c(1,2))
    #ancestry plot
    matplot(diff.table[,-c(3,4)],type='l',lty=1,lwd=1.5,col=c("dodgerblue","darkorange"),main="Ancestry Proportion",xlab="Simulation Generations",ylab="Difference")
    #junction plot
    matplot(diff.table[,-c(1,2)],type='l',lty=1,lwd=1.5,col=c("dodgerblue","darkorange"),main="Junction Counts",xlab="Simulation Generations",ylab="Difference")
    legend("bottomright",c("Chromosome 1","Chromosome 2"),col=c("dodgerblue","darkorange"),lty=1,lwd=2)
}

plot_distributions = function(outfile.a,outfile.b,reps,bins,type) {
    load(paste(outfile.a,".stats.RDATA",sep=''))     
    ancestry.neu = summarylist$ancestry     
    junction.neu = summarylist$junctions     
    tracts.neu = summarylist$tracts  
    load(paste(outfile.b,".stats.RDATA",sep=''))                 
    ancestry.sel = summarylist$ancestry   
    junction.sel = summarylist$junctions   
    tracts.sel = summarylist$tracts        
    # ancestry and junction distributions are first
    anc.neu.a = cbind(ancestry.neu$mean.a,rep('Neutral',reps))  
    anc.neu.b = cbind(ancestry.neu$mean.b,rep('Neutral',reps))  
    anc.sel.a = cbind(ancestry.sel$mean.a,rep('Selection',reps))    
    anc.sel.b = cbind(ancestry.sel$mean.b,rep('Selection',reps))      
    jun.neu.a = cbind(junction.neu$mean.a,rep('Neutral',reps))    
    jun.neu.b = cbind(junction.neu$mean.b,rep('Neutral',reps))               
    jun.sel.a = cbind(junction.sel$mean.a,rep('Selection',reps))      
    jun.sel.b = cbind(junction.sel$mean.b,rep('Selection',reps))  
    anc.a = cbind(rbind(anc.neu.a,anc.sel.a),rep("Locus 1",reps*2))
    anc.b = cbind(rbind(anc.neu.b,anc.sel.b),rep("Locus 2",reps*2))  
    jun.a = cbind(rbind(jun.neu.a,jun.sel.a),rep("Locus 1",reps*2))
    jun.b = cbind(rbind(jun.neu.b,jun.sel.b),rep("Locus 2",reps*2))  
    anc = rbind(anc.a,anc.b)
    jun = rbind(jun.a,jun.b)
    anc.df = data.frame(Ancestry=as.numeric(anc[,1]),Locus=anc[,3],Simulation=anc[,2])
    jun.df = data.frame(Junctions=as.numeric(jun[,1]),Locus=jun[,3],Simulation=jun[,2])    
    # now the tracts
    tracts.df = make_tracts_df(tracts.sel,tracts.neu,reps,bins)
    tracts.df = tracts.df[-which(tracts.df$Length==1),]    
    # plotting
    library('RColorBrewer')
    library('ggplot2')
    library('gridExtra') 
    library('grid')   
    tracts.a.a = tracts.df[intersect(which(tracts.df$Locus=="Locus 1"),which(tracts.df$Ancestry=="A")),]
    tracts.a.b = tracts.df[intersect(which(tracts.df$Locus=="Locus 1"),which(tracts.df$Ancestry=="B")),]
    tracts.b.a = tracts.df[intersect(which(tracts.df$Locus=="Locus 2"),which(tracts.df$Ancestry=="A")),]
    tracts.b.b = tracts.df[intersect(which(tracts.df$Locus=="Locus 2"),which(tracts.df$Ancestry=="B")),]  
    if (type=="boxplots") {
        p1 = ggplot(aes(y = Ancestry, x = Locus, fill = Simulation), data = anc.df) + geom_boxplot() + theme(legend.position="none") + labs(x ="")
        p2 = ggplot(aes(y = Junctions, x = Locus, fill = Simulation), data = jun.df) + geom_boxplot() + guides(fill=FALSE) + theme(legend.position="right") + labs(x ="")
        p3 = ggplot(aes(y = Frequency, x = as.character(Length), fill = Simulation), data = tracts.a.a) + geom_boxplot() + theme(legend.position="none") + labs(x ="",title="Locus 1")
        p4 = ggplot(aes(y = Frequency, x = as.character(Length), fill = Simulation), data = tracts.b.b) + geom_boxplot() + theme(legend.position="none") + labs(x ="",y="",title="Locus 2") 
        p5 = ggplot(aes(y = Frequency, x = as.character(Length), fill = Simulation), data = tracts.a.b) + geom_boxplot() + theme(legend.position="none") + labs(x ="Tract Length")
        p6 = ggplot(aes(y = Frequency, x = as.character(Length), fill = Simulation), data = tracts.b.a) + geom_boxplot() + theme(legend.position="none") + labs(x ="Tract Length",y="")   
        grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,nrow=3,position="bottom")         
    } 
    if (type=="lineplots") {
        tract.means.a.a = compute_tract_means(tracts.a.a)
        tract.means.a.b = compute_tract_means(tracts.a.b)
        tract.means.b.a = compute_tract_means(tracts.b.a)
        tract.means.b.b = compute_tract_means(tracts.b.b)       
        p1 = ggplot(aes(y = Ancestry, x = Locus, fill = Simulation), data = anc.df) + geom_boxplot() + theme(legend.position="none") + labs(x ="")
        p2 = ggplot(aes(y = Junctions, x = Locus, fill = Simulation), data = jun.df) + geom_boxplot() + guides(fill=FALSE) + theme(legend.position="right") + labs(x ="")
        p3 = ggplot(data=tract.means.a.a, aes(x=as.character(Length), y=Mean, ymin=Lower, ymax=Upper, fill=Simulation, group=Simulation)) + geom_line() + geom_ribbon(alpha=0.5) + theme(legend.position="none") + labs(x ="",title="Locus 1")
        p4 = ggplot(data=tract.means.a.b, aes(x=as.character(Length), y=Mean, ymin=Lower, ymax=Upper, fill=Simulation, group=Simulation)) + geom_line() + geom_ribbon(alpha=0.5) + theme(legend.position="none") + labs(x ="",y="",title="Locus 2")    
        p5 = ggplot(data=tract.means.b.a, aes(x=as.character(Length), y=Mean, ymin=Lower, ymax=Upper, fill=Simulation, group=Simulation)) + geom_line() + geom_ribbon(alpha=0.5) + theme(legend.position="none") + labs(x ="Tract Length")
        p6 = ggplot(data=tract.means.b.b, aes(x=as.character(Length), y=Mean, ymin=Lower, ymax=Upper, fill=Simulation, group=Simulation)) + geom_line() + geom_ribbon(alpha=0.5) + theme(legend.position="none") + labs(x ="Tract Length",y="")
        grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,nrow=3,position="bottom")         
    }
}

compute_tract_means = function(tractsdata) {
    L.vec = unique(tractsdata$Length)    
    tract.matrix.sel = matrix(nrow=length(L.vec),ncol=6)
    tract.matrix.neu = matrix(nrow=length(unique(L.vec)),ncol=6)
    tracts.sel = tractsdata[which(tractsdata$Simulation=="Selection"),]
    tracts.neu = tractsdata[which(tractsdata$Simulation=="Neutral"),]   
    for (i in 1:length(L.vec)) {
        tract.matrix.sel[i,1] = L.vec[i]
        tract.matrix.neu[i,1] = L.vec[i]
        tract.matrix.sel[i,2] = mean(tracts.sel$Frequency[which(tracts.sel$Length==L.vec[i])])
        tract.matrix.neu[i,2] = mean(tracts.neu$Frequency[which(tracts.neu$Length==L.vec[i])])
        tract.matrix.sel[i,c(3,4)] = quantile(tracts.sel$Frequency[which(tracts.sel$Length==L.vec[i])],probs=c(.05,.95))
        tract.matrix.neu[i,c(3,4)] = quantile(tracts.neu$Frequency[which(tracts.neu$Length==L.vec[i])],probs=c(.05,.95))
    }
    tract.df.sel = data.frame(tract.matrix.sel)    
    tract.df.neu = data.frame(tract.matrix.neu)   
    tract.df.sel[,5] = tracts.sel[1,5] 
    tract.df.neu[,5] = tracts.neu[1,5] 
    tract.df.sel[,6] = tracts.sel[1,6] 
    tract.df.neu[,6] = tracts.neu[1,6] 
    tract.df = rbind(tract.df.sel,tract.df.neu)    
    colnames(tract.df) = c("Length","Mean","Lower","Upper","Simulation","Locus")
    return(tract.df)
} 
 
 
make_tracts_df = function(tracts.sel,tracts.neu,reps,bins) {
    tracts.df = data.frame(Replicate = integer(),Length = integer(),Counts = integer(),Ancestry = character(),Simulation = character(),Locus = character())
    for (i in 1:reps) {
        sel.rep.a = cbind(tracts.sel$a.lens[i,],tracts.sel$a.anc[i,])        
        sel.rep.b = cbind(tracts.sel$b.lens[i,],tracts.sel$b.anc[i,])
        neu.rep.a = cbind(tracts.neu$a.lens[i,],tracts.neu$a.anc[i,])
        neu.rep.b = cbind(tracts.neu$b.lens[i,],tracts.neu$b.anc[i,])
        breakpoints = seq(0,1,1/bins)
        sel.hist.a.a = hist(sel.rep.a[which(sel.rep.a[,2]==0),1],plot=F,breaks=breakpoints)
        sel.hist.a.b = hist(sel.rep.a[which(sel.rep.a[,2]==1),1],plot=F,breaks=breakpoints)
        sel.hist.b.a = hist(sel.rep.b[which(sel.rep.b[,2]==0),1],plot=F,breaks=breakpoints)
        sel.hist.b.b = hist(sel.rep.b[which(sel.rep.b[,2]==1),1],plot=F,breaks=breakpoints)        
        neu.hist.a.a = hist(neu.rep.a[which(neu.rep.a[,2]==0),1],plot=F,breaks=breakpoints)
        neu.hist.a.b = hist(neu.rep.a[which(neu.rep.a[,2]==1),1],plot=F,breaks=breakpoints)
        neu.hist.b.a = hist(neu.rep.b[which(neu.rep.b[,2]==0),1],plot=F,breaks=breakpoints)
        neu.hist.b.b = hist(neu.rep.b[which(neu.rep.b[,2]==1),1],plot=F,breaks=breakpoints)   
        sel.a.a = data.frame(Replicate = rep(i,length(sel.hist.a.a$counts)), Length = sel.hist.a.a$breaks[-1],
                             Frequency = sel.hist.a.a$counts/sum(sel.hist.a.a$counts), Ancestry = rep("A",length(sel.hist.a.a$counts)),
                             Simulation = rep("Selection",length(sel.hist.a.a$counts)), 
                             Locus = rep("Locus 1",length(sel.hist.a.a$counts)))
        sel.a.b = data.frame(Replicate = rep(i,length(sel.hist.a.b$counts)), Length = sel.hist.a.b$breaks[-1],
                             Frequency = sel.hist.a.b$counts/sum(sel.hist.a.b$counts), Ancestry = rep("B",length(sel.hist.a.b$counts)),
                             Simulation = rep("Selection",length(sel.hist.a.b$counts)), 
                             Locus = rep("Locus 1",length(sel.hist.a.b$counts)))                                                   
        sel.b.a = data.frame(Replicate = rep(i,length(sel.hist.b.a$counts)), Length = sel.hist.b.a$breaks[-1],
                             Frequency = sel.hist.b.a$counts/sum(sel.hist.b.a$counts), Ancestry = rep("A",length(sel.hist.b.a$counts)),
                             Simulation = rep("Selection",length(sel.hist.b.a$counts)), 
                             Locus = rep("Locus 2",length(sel.hist.b.a$counts)))
        sel.b.b = data.frame(Replicate = rep(i,length(sel.hist.b.b$counts)), Length = sel.hist.b.b$breaks[-1],
                             Frequency = sel.hist.b.b$counts/sum(sel.hist.b.b$counts), Ancestry = rep("B",length(sel.hist.b.b$counts)),
                             Simulation = rep("Selection",length(sel.hist.b.b$counts)), 
                             Locus = rep("Locus 2",length(sel.hist.b.b$counts)))   
        neu.a.a = data.frame(Replicate = rep(i,length(neu.hist.a.a$counts)), Length = neu.hist.a.a$breaks[-1],
                             Frequency = neu.hist.a.a$counts/sum(neu.hist.a.a$counts), Ancestry = rep("A",length(neu.hist.a.a$counts)),
                             Simulation = rep("Neutral",length(neu.hist.a.a$counts)), 
                             Locus = rep("Locus 1",length(neu.hist.a.a$counts)))
        neu.a.b = data.frame(Replicate = rep(i,length(neu.hist.a.b$counts)), Length = neu.hist.a.b$breaks[-1],
                             Frequency = neu.hist.a.b$counts/sum(neu.hist.a.b$counts), Ancestry = rep("B",length(neu.hist.a.b$counts)),
                             Simulation = rep("Neutral",length(neu.hist.a.b$counts)), 
                             Locus = rep("Locus 1",length(neu.hist.a.b$counts)))                                                   
        neu.b.a = data.frame(Replicate = rep(i,length(neu.hist.b.a$counts)), Length = neu.hist.b.a$breaks[-1],
                             Frequency = neu.hist.b.a$counts/sum(neu.hist.b.a$counts), Ancestry = rep("A",length(neu.hist.b.a$counts)),
                             Simulation = rep("Neutral",length(neu.hist.b.a$counts)), 
                             Locus = rep("Locus 2",length(neu.hist.b.a$counts)))
        neu.b.b = data.frame(Replicate = rep(i,length(neu.hist.b.b$counts)), Length = neu.hist.b.b$breaks[-1],
                             Frequency = neu.hist.b.b$counts/sum(neu.hist.b.b$counts), Ancestry = rep("B",length(neu.hist.b.b$counts)),
                             Simulation = rep("Neutral",length(neu.hist.b.b$counts)), 
                             Locus = rep("Locus 2",length(neu.hist.b.b$counts)))
        new.rep = rbind(sel.a.a,sel.a.b,sel.b.a,sel.b.b,neu.a.a,neu.a.b,neu.b.a,neu.b.b)
        tracts.df = rbind(tracts.df,new.rep)
    }
    return(tracts.df)
} 
                                          
run_dfuse = function(popsize, chrlen, migrate, nmarkers, tstart, dmipos, seldmi, selanc, domdmi, reps, demefile, episfile, outfile, chrnum, model) {
    # Write the epistasis and deme files for dfuse    
    # loci on same or different chromosomes
    if (chrnum==1) { 
        episinfo = matrix(nrow=8,ncol=3)  
        if (model=="dmi") { 
            episinfo[,1] = c(2,0,0,0,0,0,dmipos[1],0)   
            episinfo[,2] = c(1,1,0,1,0,0,dmipos[2],1)   
        }
        if (model=="pathway") {
            episinfo[,1] = c(2,0,0,0,1,0,dmipos[1],0)        
            episinfo[,2] = c(1,1,0,1,1,1,dmipos[2],1)            
        }  
        episinfo[1,3] = 2    
    }                 
    if (chrnum==2) {                                       
        episinfo = matrix(nrow=8,ncol=3)
        if (model=="dmi") {               
            episinfo[,1] = c(2,0,0,0,0,0,dmipos[1],0)        
            episinfo[,2] = c(1,1,0,1,0,1,dmipos[2],1)
        }    
        if (model=="pathway") {
            episinfo[,1] = c(2,0,0,0,1,0,dmipos[1],0)        
            episinfo[,2] = c(1,1,0,1,1,1,dmipos[2],1)        
        }      
        episinfo[1,3] = 2                             
    }                                                    
    demeinfo = c(1, popsize, popsize)                        
    write(demeinfo,file=demefile,ncolumns=1)                                       
    write.table(episinfo,file=episfile,row.names=F,col.names=F,quote=F,na="")           
    # Run dfuse                                                                   
    system(paste("dfuse/dfuse -d",demefile,"-g",tstart,"-m",migrate,                  
                 "-G",tstart,"-o",outfile,"-l",nmarkers,                           
                 "-e",episfile,"-c",seldmi,"-a",selanc,                         
                 "-R",domdmi,"-S 'fertility' -r",reps))                        
    # Get the stuff we want from the dfuse output                                    
    markers = scan(file=paste(outfile,".haplo",sep=''),nlines=1,what="character")  
    positions = as.numeric(markers[-c(1:2)])
    if (chrnum==1) { 
        listtable = read.table(paste(outfile,".haplo",sep=''),sep=',')[,c(c(1:4),c(9:(nmarkers+8)))]  
    }
    if (chrnum==2) {
        listtable = read.table(paste(outfile,".haplo",sep=''),sep=',')[,c(c(1:4),c(9:((nmarkers*2)+8)))]      
    }
    # Convert to matrix
    haptable = no_list_func(listtable,"char")             
    if (length(which(haptable[,5]=="NA"))>0) {            
        haptable = haptable[-which(haptable[,5]=="NA"),]    
    }
    haptable[which(haptable=="a")] = 0   
    haptable[which(haptable=="b")] = 1  
    haplotypes = apply(haptable[-1,],2,as.numeric)
    if (chrnum==1) {
        colnames(haplotypes) = c(haptable[1,c(1:4)],positions) 
    }
    if (chrnum==2) {
        colnames(haplotypes) = c(haptable[1,c(1:4)],positions,positions) 
    }    
    save(haplotypes,file=paste(outfile,".haplo.RDATA",sep=''))
}   
                
compute_summary_stats = function(outfile,dmipos,chrnum,windowsize) {
    load(paste(outfile,".haplo.RDATA",sep=''))   
    reps = length(unique(haplotypes[,1]))       
    if (reps==1) {  
        ancestry = ancestry_proportion(haplotypes, chrnum, windowsize, dmipos)  
        junctions = junction_count(haplotypes, chrnum, windowsize, dmipos)  
        tracts = tract_lengths(haplotypes, dmipos, chrnum) 
    } 
    if (reps>1) {
        if (chrnum==1) {
            ancestry.a = matrix(ncol=ncol(haplotypes)-4,nrow=reps)  
            junctions.a = matrix(ncol=ncol(haplotypes)-4,nrow=reps)  
            anc.window.a = c(1:reps)
            anc.window.b = c(1:reps)
            jun.window.a = c(1:reps)
            jun.window.b = c(1:reps)
            colnames(ancestry.a) = colnames(haplotypes)[-c(1:4)]
            colnames(junctions.a) = colnames(haplotypes)[-c(1:4)]
        }
        if (chrnum==2) {
            colnames = colnames(haplotypes)
            names = colnames[-c(1:min(which(colnames==1)))]
            ancestry.a = matrix(ncol=length(names),nrow=reps)  
            ancestry.b = matrix(ncol=length(names),nrow=reps) 
            junctions.a = matrix(ncol=length(names),nrow=reps) 
            junctions.b = matrix(ncol=length(names),nrow=reps) 
            anc.window.a = c(1:reps)
            anc.window.b = c(1:reps)
            jun.window.a = c(1:reps)
            jun.window.b = c(1:reps)
            colnames(ancestry.a) = names  
            colnames(junctions.a) = names  
            colnames(ancestry.b) = names  
            colnames(junctions.b) = names    
        } 
        a.lentable = matrix(ncol=popsize*2,nrow=reps)  
        b.lentable = matrix(ncol=popsize*2,nrow=reps)  
        a.anctable = matrix(ncol=popsize*2,nrow=reps)  
        b.anctable = matrix(ncol=popsize*2,nrow=reps) 
        for (i in 1:reps) { 
            haprep = haplotypes[which(haplotypes[,1]==i),]
            if (chrnum==1) {  
                ancestry.list = ancestry_proportion(haprep,chrnum,windowsize,dmipos)
                junctions.list = junction_count(haprep,chrnum,windowsize,dmipos)                
                ancestry.a[i,] = ancestry.list[[1]]
                junctions.a[i,] = junctions.list[[1]]
                anc.window.a[i] = ancestry.list[[2]]
                anc.window.b[i] = ancestry.list[[3]]
                jun.window.a[i] = junctions.list[[2]]
                jun.window.b[i] = junctions.list[[3]]
                ancestry = list("ancestry"=ancestry.a,"mean.a"=anc.window.a,"mean.b"=anc.window.b)  
                junctions = list("junctions"=junctions.a,"mean.a"=jun.window.a,"mean.b"=jun.window.b)                  
            }
            tracts = tract_lengths(haprep,dmipos,chrnum)  
            a.lentable[i,c(1:length(tracts$a.lens))] = tracts$a.lens  
            b.lentable[i,c(1:length(tracts$b.lens))] = tracts$b.lens   
            a.anctable[i,c(1:length(tracts$a.anc))] = tracts$a.anc    
            b.anctable[i,c(1:length(tracts$b.anc))] = tracts$b.anc   
            if (chrnum==2) {     
                ancestry.list = ancestry_proportion(haprep,chrnum,windowsize,dmipos)
                junctions.list = junction_count(haprep,chrnum,windowsize,dmipos)                            
                ancestry.a[i,] = ancestry.list[[1]]
                ancestry.b[i,] = ancestry.list[[2]] 
                junctions.a[i,] = junctions.list[[1]] 
                junctions.b[i,] = junctions.list[[2]]
                anc.window.a[i] = ancestry.list[[3]]
                anc.window.b[i] = ancestry.list[[4]]
                jun.window.a[i] = junctions.list[[3]]
                jun.window.b[i] = junctions.list[[4]]
                ancestry = list("ancestry.a"=ancestry.a, "ancestry.b"=ancestry.b,"mean.a"=anc.window.a,"mean.b"=anc.window.b)  
                junctions = list("junctions.a"=junctions.a,"junctions.b"=junctions.b,"mean.a"=jun.window.a,"mean.b"=jun.window.b)                  
            }                    
        }      
        tracts = list("a.lens"=a.lentable,"b.lens"=b.lentable,"a.anc"=a.anctable,"b.anc"=b.anctable)                                                                      
    }  
    summarylist = list("ancestry" = ancestry,"junctions"=junctions,"tracts"=tracts)
    save(summarylist,file=paste(outfile,".stats.RDATA",sep=''))       
}             
          
ancestry_proportion = function(haprep,chrnum,windowsize,dmipos) {   
    if (chrnum==1) {  
        colnames = colnames(haprep) 
        pos = colnames[-c(1:4)]    
        ancestry = apply(haprep[,-c(1:4)],2,mean) 
        left.a = max(which(pos<=(dmipos[1]-(windowsize/2))))+1
        left.b = max(which(pos<=(dmipos[2]-(windowsize/2))))+1
        right.a = min(which(pos>=(dmipos[1]+(windowsize/2))))-1
        right.b = min(which(pos>=(dmipos[2]+(windowsize/2))))-1
        mean.a = mean(ancestry[left.a:right.a])
        mean.b = mean(ancestry[left.b:right.b])   
        ancestry = list("ancestry"=ancestry,"mean.a"=mean.a,"mean.b"=mean.b)
    }   
    if (chrnum==2) {   
        colnames = colnames(haprep) 
        pos = colnames[-c(1:min(which(colnames==1)))]
        firstchrom = -c(c(1:4),c((min(which(colnames==1))+1):ncol(haprep)))  
        seconchrom = -c(1:min(which(colnames==1)))  
        ancestry.a = apply(haprep[,firstchrom],2,mean)    
        ancestry.b = apply(haprep[,seconchrom],2,mean)  
        left = max(which(pos<=(dmipos[1]-(windowsize/2))))+1
        right = min(which(pos>=(dmipos[1]+(windowsize/2))))-1
        mean.a = mean(ancestry.a[left:right])   
        mean.b = mean(ancestry.b[left:right])   
        ancestry = list("ancestry.a"=ancestry.a,"ancestry.b"=ancestry.b,"mean.a"=mean.a,"mean.b"=mean.b)  
    }    
    return(ancestry)  
}     
           
junction_count = function(haprep,chrnum,windowsize,dmipos) {
    if (chrnum==1) {
        colnames = colnames(haprep) 
        pos = colnames[-c(1:4)]    
        diffs = apply(haprep[,-c(1:4)],1,function(X) abs(diff(X)))
        junctions = c(0,apply(t(diffs),2,sum))
        left.a = max(which(pos<=(dmipos[1]-(windowsize/2))))+1
        left.b = max(which(pos<=(dmipos[2]-(windowsize/2))))+1
        right.a = min(which(pos>=(dmipos[1]+(windowsize/2))))-1
        right.b = min(which(pos>=(dmipos[2]+(windowsize/2))))-1
        mean.a = mean(junctions[left.a:right.a])
        mean.b = mean(junctions[left.b:right.b])   
        junctions = list("junctions"=junctions,"mean.a"=mean.a,"mean.b"=mean.b)
    }
    if (chrnum==2) {
        colnames = colnames(haprep)
        pos = colnames[-c(1:min(which(colnames==1)))]        
        firstchrom = -c(c(1:4),c((min(which(colnames==1))+1):ncol(haprep)))
        seconchrom = -c(1:min(which(colnames==1)))
        diffs.a = apply(haprep[,firstchrom],1,function(X) abs(diff(X)))
        diffs.b = apply(haprep[,seconchrom],1,function(X) abs(diff(X)))      
        junctions.a = c(0,apply(t(diffs.a),2,sum))
        junctions.b = c(0,apply(t(diffs.b),2,sum))
        left = max(which(pos<=(dmipos[1]-(windowsize/2))))+1
        right = min(which(pos>=(dmipos[1]+(windowsize/2))))-1
        mean.a = mean(junctions.a[left:right])   
        mean.b = mean(junctions.b[left:right])               
        junctions = list("junctions.a"=junctions.a,"junctions.b"=junctions.b,"mean.a"=mean.a,"mean.b"=mean.b)     
    }
    return(junctions)
} 
  
tract_lengths = function(haprep,dmipos,chrnum) { 
    # For one chromosome
    if (chrnum==1) {
        pos = as.numeric(colnames(haprep)[-c(1:4)])   
        haps = haprep[,-c(1:4)]               
        anclist = tract_ancestries(pos, dmipos, haps, chrnum)   
        dmiindex = anclist$dmiindex                                            
        diffs = t(apply(haps,1,function(X) abs(diff(X))))      
        diffs = cbind(rep(1,nrow(diffs)),diffs)          
        diffs[,ncol(diffs)] = 1                      
        junctpos = apply(diffs,1,function(X) which(X==1))         
        a.haps = matrix(nrow=nrow(haprep),ncol=2)         
        b.haps  = matrix(nrow=nrow(haprep),ncol=2)       
        a.haps [,c(1,2)] = cbind(rep(0,nrow(haps)),rep(1,nrow(haps)))             
        b.haps [,c(1,2)] = cbind(rep(0,nrow(haps)),rep(1,nrow(haps))) 
        if (typeof(junctpos)=="list") {        
            for (i in 1:nrow(diffs)) {               
                a.haps [i,1] = pos[junctpos[[i]][max(which(junctpos[[i]]<dmiindex[1]))]]  
                a.haps [i,2] = pos[junctpos[[i]][min(which(junctpos[[i]]>dmiindex[1]))]]    
                b.haps [i,1] = pos[junctpos[[i]][max(which(junctpos[[i]]<dmiindex[2]))]]  
                b.haps [i,2] = pos[junctpos[[i]][min(which(junctpos[[i]]>dmiindex[2]))]]            
            }   
        }                               
        a.lens = a.haps[,2]-a.haps[,1]  
        b.lens = b.haps[,2]-b.haps[,1]  
        a.anc = anclist$anctable[,1]   
        b.anc = anclist$anctable[,2]           
        tracts = list("a.lens"=a.lens,"b.lens"=b.lens,"a.anc"=a.anc,"b.anc"=b.anc)   
    }
    # For two chromosomes
    if (chrnum==2) {
        colnames = colnames(haprep)
        firstchrom = -c(c(1:4),c((min(which(colnames==1))+1):ncol(haprep)))
        seconchrom = -c(1:min(which(colnames==1)))
        haps.a = haprep[,firstchrom]
        haps.b = haprep[,seconchrom]
        pos = colnames(haps.a)
        haps = list("haps.a"=haps.a,"haps.b"=haps.b)
        anclist = tract_ancestries(pos, dmipos, haps, chrnum)    
        dmiindex = anclist$dmiindex                                                    
        diffs.a = t(apply(haps.a,1,function(X) abs(diff(X))))      
        diffs.b = t(apply(haps.b,1,function(X) abs(diff(X))))  
        diffs.a = cbind(rep(1,nrow(diffs.a)),diffs.a)      
        diffs.b = cbind(rep(1,nrow(diffs.b)),diffs.b)         
        diffs.a[,ncol(diffs.a)] = 1   
        diffs.b[,ncol(diffs.b)] = 1                                                 
        junctpos.a = apply(diffs.a,1,function(X) which(X==1))         
        junctpos.b = apply(diffs.b,1,function(X) which(X==1))         
        a.haps = matrix(nrow=nrow(haps.a),ncol=2)         
        b.haps  = matrix(nrow=nrow(haps.b),ncol=2)       
        a.haps [,c(1,2)] = cbind(rep(0,nrow(haps.a)),rep(1,nrow(haps.a)))             
        b.haps [,c(1,2)] = cbind(rep(0,nrow(haps.b)),rep(1,nrow(haps.b))) 
        if ((typeof(junctpos.a)=="list")&&(typeof(junctpos.b)=="list")) {        
            for (i in 1:nrow(diffs.a)) {  
                a.haps[i,1] = pos[junctpos.a[[i]][max(which(junctpos.a[[i]]<dmiindex))]]  
                a.haps[i,2] = pos[junctpos.a[[i]][min(which(junctpos.a[[i]]>dmiindex))]]    
                b.haps[i,1] = pos[junctpos.b[[i]][max(which(junctpos.b[[i]]<dmiindex))]]  
                b.haps[i,2] = pos[junctpos.b[[i]][min(which(junctpos.b[[i]]>dmiindex))]]            
            }   
        }                               
        a.lens = as.numeric(a.haps[,2])-as.numeric(a.haps[,1])  
        b.lens = as.numeric(b.haps[,2])-as.numeric(b.haps[,1])  
        a.anc = anclist$anctable[,1]   
        b.anc = anclist$anctable[,2]           
        tracts = list("a.lens"=a.lens,"b.lens"=b.lens,"a.anc"=a.anc,"b.anc"=b.anc)                 
    }        
    return(tracts)         
}                               
                 
tract_ancestries = function(pos,dmipos,haps,chrnum) {        
    if (chrnum==1) {    
        if (length(intersect(pos,dmipos))==0) {
            dmiindex = which(!is.na(match(sort(c(pos,dmipos)),dmipos)))      
            dmitab = cbind(rep(2,nrow(haps)),rep(2,nrow(haps)))       
            newhaps = cbind(haps,dmitab)                           
            colnames(newhaps) = c(colnames(haps),dmipos)        
            newhaps = newhaps[,order(as.numeric(colnames(newhaps)))]    
            a.loc = newhaps[,c(dmiindex[1]-1,dmiindex[1]+1)]    
            b.loc = newhaps[,c(dmiindex[2]-1,dmiindex[2]+1)]   
            a.diffs = which(c(a.loc[,1]-a.loc[,2])!=0)           
            b.diffs = which(c(b.loc[,1]-b.loc[,2])!=0)   
            a.anc = a.loc[,1]
            b.anc = b.loc[,1]   
            a.anc[a.diffs] = rbinom(length(a.diffs),1,.5)
            b.anc[b.diffs] = rbinom(length(b.diffs),1,.5)
            anctable = cbind(a.anc,b.anc)
        }                                                                          
        if (length(intersect(pos,dmipos))==2) {                    
            dmiindex = which(!is.na(match(pos,dmipos)))                
            anctable = haps[,dmiindex]                        
        }   
    }
    if (chrnum==2) {    
        dmipos = dmipos[1]
        haps.a = haps[[1]]
        haps.b = haps[[2]]
        if (length(intersect(pos,dmipos[1]))==0) {                            
            dmiindex = which(!is.na(match(sort(c(pos,dmipos)),dmipos)))   
            dmitab = cbind(rep(2,nrow(haps.a)))        
            newhaps.a = cbind(haps.a,dmitab)                 
            newhaps.b = cbind(haps.b,dmitab)                        
            colnames(newhaps.a) = c(colnames(haps.a),dmipos)     
            colnames(newhaps.b) = c(colnames(haps.b),dmipos)                       
            newhaps.a = newhaps.a[,order(as.numeric(colnames(newhaps.a)))]    
            newhaps.b = newhaps.b[,order(as.numeric(colnames(newhaps.b)))]    
            a.loc = newhaps.a[,c(dmiindex-1,dmiindex+1)]       
            b.loc = newhaps.b[,c(dmiindex-1,dmiindex+1)]       
            a.diffs = which(c(a.loc[,1]-a.loc[,2])!=0)            
            b.diffs = which(c(b.loc[,1]-b.loc[,2])!=0)     
            a.anc = a.loc[,1]                              
            b.anc = b.loc[,1]                                  
            a.anc[a.diffs] = rbinom(length(a.diffs),1,.5)   
            b.anc[b.diffs] = rbinom(length(b.diffs),1,.5)   
            anctable = cbind(a.anc,b.anc)                   
        }   
        if (length(intersect(pos,dmipos))==2) {                    
            dmiindex = which(!is.na(match(pos,dmipos)))                
            anctable = cbind(haps.a[,dmiindex],haps.b[,dmiindex])                        
        }   
    }    
    anclist = list("anctable"=anctable,"dmiindex"=dmiindex)
    return(anclist) 
}   
  
plot_haplotypes = function(outfile,dmipos,chrnum) {
    load(paste(outfile,".haplo.RDATA",sep=''))  
    reps = length(unique(haplotypes[,1]))
    if (chrnum==1) {
        if (reps==1) {
            image(t(haplotypes[,-c(1:4)]),col=c('dodgerblue','orange'))
            abline(v=dmipos,lty=2)
        } 
        if (reps>1) {
            par(mfrow=c((reps/5),5))
            for (i in 1:reps) {
                haprep = haplotypes[which(haplotypes[,1]==i),-c(1:4)]
                image(t(haprep),col=c('dodgerblue','orange'))
                abline(v=c(dmipos),lty=2)
            }
        } 
    }
    if (chrnum==2) {
        colnames = colnames(haplotypes)
        firstchrom = -c(c(1:4),c((min(which(colnames==1))+1):ncol(haplotypes)))
        seconchrom = -c(1:min(which(colnames==1)))
        haps.a = haplotypes[,firstchrom]
        haps.b = haplotypes[,seconchrom]    
        if (reps==1) {
            par(mfrow=c(1,2)) 
            par(mar=c(5,5,4,0))
            image(t(haps.a),col=c('dodgerblue','orange'),ylab="Haplotype",main="Chromosome 1",cex.lab=1.4,cex.main=1.2)
            abline(v=dmipos[[1]],lty=2)   
            par(mar=c(5,2,4,3))   
            image(t(haps.b),col=c('dodgerblue','orange'),main="Chromosome 2",yaxt='n',cex.main=1.2) 
            abline(v=dmipos[[2]],lty=2)   
            mtext("Position (Morgans)",side=1,line=2.5,at=-.1,cex=1.4)  
        }               
        if (reps>1) {   
            par(mfrow=c(reps,2))  
            for (i in 1:reps) {    
                haprep.a = haps.a[which(haplotypes[,1]==i),]  
                haprep.b = haps.b[which(haplotypes[,1]==i),]  
                par(mar=c(1,5,1,0))
                image(t(haprep.a),col=c('dodgerblue','orange'),ylab="",main="",xaxt='n',cex.lab=1.4,cex.main=1.2)
                abline(v=dmipos[[1]])   
                par(mar=c(1,2,1,3))   
                image(t(haprep.b),col=c('dodgerblue','orange'),main="",yaxt='n',xaxt='n',cex.main=1.2) 
                abline(v=dmipos[[2]])  
            }
        } 
    }
}  
 
plot_summary_stats = function(outfile, dmipos, windowsize, reps) {
    load(paste(outfile,".stats.RDATA",sep=''))     
    ancestry = summarylist$ancestry
    junctions = summarylist$junctions
    tracts = summarylist$tracts
    if (chrnum==1) {
        if (reps==1) {
            layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))
            par(mar=c(4,4.5,3,2))
            plot(names(ancestry),ancestry,type='l',lwd=1.5,col='darkolivegreen4',
                 ylab='Ancestry Proportion',xlab='Position (Morgans)',ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
            abline(v=dmipos,lty=2,lwd=1.5)
            par(mar=c(5,4.5,3,2))
            plot(names(junctions),junctions,type='l',lwd=1.5,col='darkolivegreen4',
                 ylab='Junction Count',xlab='Position (Morgans)',cex.axis=1.5,cex.lab=1.5)    
            abline(v=dmipos,lty=2,lwd=1.5)     
            a.tracts = cbind(tracts$a.lens,tracts$a.anc)
            b.tracts = cbind(tracts$b.lens,tracts$b.anc)
            par(mar=c(5,4.5,2,2))
            plot(density(a.tracts[which(a.tracts[,2]==0)]),col="dodgerblue",lwd=1.5,ylim=c(0,2),
                 xlab="Tract Length (Morgans)",main='Locus A',cex.main=1.3,cex.lab=1.5,cex.axis=1.5)
            lines(density(a.tracts[which(a.tracts[,2]==1)]),col="darkorange",lwd=1.5,
                  xlab="Tract Length (Morgans)")           
            par(mar=c(5,4.5,2,2))         
            plot(density(b.tracts[which(b.tracts[,2]==0)]),col="dodgerblue",lwd=1.5,ylim=c(0,2),   
                 xlab="Tract Length (Morgans)",main='Locus B',cex.main=1.3,cex.lab=1.5,cex.axis=1.5)   
            lines(density(b.tracts[which(b.tracts[,2]==1)]),col="darkorange",lwd=1.5,   
                  xlab="Tract Length (Morgans)")                          
        }   
        if (reps!=1) {    
            layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))
            par(mar=c(4,4.5,3,2))
            matplot(colnames(ancestry[[1]]),t(ancestry[[1]]),type='l',lty=1,col='darkolivegreen4',   
                    ylab='Ancestry Proportion',xlab='Position (Morgans)',ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
            abline(v=dmipos,lty=2,lwd=1.5)
            par(mar=c(5,4.5,3,2))
            matplot(colnames(junctions[[1]]),t(junctions[[1]]),type='l',lty=1,col='darkolivegreen4',
                 ylab='Junction Count',xlab='Position (Morgans)',cex.axis=1.5,cex.lab=1.5)    
            abline(v=dmipos,lty=2,lwd=1.5)    
            a.tracts = cbind(tracts$a.lens[1,],tracts$a.anc[1,])
            par(mar=c(5,4.5,2,2))
            plot(density(a.tracts[which(a.tracts[,2]==0)]),col="dodgerblue",lwd=1,
                 xlab="Tract Length (Morgans)",main='Locus A',cex.main=1.3,cex.lab=1.5,cex.axis=1.5,ylim=c(0,2))
            lines(density(a.tracts[which(a.tracts[,2]==1)]),col="darkorange",lwd=1,
                  xlab="Tract Length (Morgans)")                    
            for (i in 2:reps) {       
                a.tracts = cbind(tracts$a.lens[i,],tracts$a.anc[i,])                              
                lines(density(a.tracts[which(a.tracts[,2]==0)]),col="dodgerblue",lwd=1)                   
                lines(density(a.tracts[which(a.tracts[,2]==1)]),col="darkorange",lwd=1)   
            }
            b.tracts = cbind(tracts$b.lens[1,],tracts$b.anc[1,])
            par(mar=c(5,4.5,2,2))
            plot(density(b.tracts[which(b.tracts[,2]==0)]),col="dodgerblue",lwd=1,
                 xlab="Tract Length (Morgans)",main='Locus B',cex.main=1.3,cex.lab=1.5,cex.axis=1.5,ylim=c(0,2))
            lines(density(b.tracts[which(b.tracts[,2]==1)]),col="darkorange",lwd=1,
                  xlab="Tract Length (Morgans)")                    
            for (i in 2:reps) {       
                b.tracts = cbind(tracts$b.lens[i,],tracts$b.anc[i,])                              
                lines(density(b.tracts[which(b.tracts[,2]==0)]),col="dodgerblue",lwd=1)                   
                lines(density(b.tracts[which(b.tracts[,2]==1)]),col="darkorange",lwd=1)   
            }                      
        }   
    } 
    if (chrnum==2) {
        if (reps==1) {   
            layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
            par(mar=c(4,4.5,3,0))
            plot(names(ancestry[[1]]),ancestry[[1]],type='l',lwd=1.5,col='darkolivegreen4',main="Chromosome 1",
                 ylab='Ancestry Proportion',xlab='Position (Morgans)',ylim=c(0,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.3)
            abline(v=dmipos[1],lty=2,lwd=1.5)
            par(mar=c(4,2.5,3,2)) 
            plot(names(ancestry[[2]]),ancestry[[2]],type='l',lwd=1.5,col='darkolivegreen4',yaxt='n',main="Chromosome 2",
                 ylab='',xlab='Position (Morgans)',ylim=c(0,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.3)
            abline(v=dmipos[2],lty=2,lwd=1.5)                
            maxy = max(c(max(junctions[[1]]),max(junctions[[2]])))           
            par(mar=c(4,4.5,3,0))
            plot(names(junctions[[1]]),junctions[[1]],type='l',lwd=1.5,col='darkolivegreen4',
                 ylab='Junction Count',xlab='Position (Morgans)',cex.axis=1.5,cex.lab=1.5,ylim=c(1,maxy))
            abline(v=dmipos[1],lty=2,lwd=1.5)
            par(mar=c(4,2.5,3,2)) 
            plot(names(junctions[[2]]),junctions[[2]],type='l',lwd=1.5,col='darkolivegreen4',yaxt='n',
                 ylab='',xlab='Position (Morgans)',cex.axis=1.5,cex.lab=1.5,ylim=c(1,maxy))
            abline(v=dmipos[2],lty=2,lwd=1.5)                
            a.tracts = cbind(tracts$a.lens,tracts$a.anc)
            b.tracts = cbind(tracts$b.lens,tracts$b.anc)
            par(mar=c(5,4.5,2,2))
            plot(density(a.tracts[which(a.tracts[,2]==0)]),col="dodgerblue",lwd=1.5,
                 xlab="Tract Length (Morgans)",main='',cex.main=1.3,cex.lab=1.5,cex.axis=1.5,ylim=c(0,2))
            lines(density(a.tracts[which(a.tracts[,2]==1)]),col="darkorange",lwd=1.5,
                  xlab="Tract Length (Morgans)")           
            par(mar=c(5,4.5,2,2))         
            plot(density(b.tracts[which(b.tracts[,2]==0)]),col="dodgerblue",lwd=1.5,   
                 xlab="Tract Length (Morgans)",main='',cex.main=1.3,cex.lab=1.5,cex.axis=1.5,ylim=c(0,2))   
            lines(density(b.tracts[which(b.tracts[,2]==1)]),col="darkorange",lwd=1.5,   
                  xlab="Tract Length (Morgans)")                                        
        }       
        if (reps!=1) {    
            layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
            par(mar=c(5,4.5,2,0))
            matplot(colnames(ancestry[[1]]),t(ancestry[[1]]),type='l',lty=1,col='darkolivegreen4',,main="Chromosome 1",   
                    ylab='Ancestry Proportion',xlab='Position (Morgans)',ylim=c(0,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.3)
            abline(v=dmipos,lty=2,lwd=1.5)
            par(mar=c(5,2.5,2,2))
            matplot(colnames(ancestry[[1]]),t(ancestry[[2]]),type='l',lty=1,col='darkolivegreen4',yaxt='n',cex.main=1.3,   
                    ylab='',xlab='Position (Morgans)',ylim=c(0,1),cex.axis=1.5,cex.lab=1.5,main="Chromosome 2")
            abline(v=dmipos,lty=2,lwd=1.5)            
            par(mar=c(5,4.5,2,0))
            matplot(colnames(junctions[[1]]),t(junctions[[1]]),type='l',lty=1,col='darkolivegreen4',
                 ylab='Junction Count',xlab='Position (Morgans)',cex.axis=1.5,cex.lab=1.5)    
            abline(v=dmipos,lty=2,lwd=1.5)
            par(mar=c(5,2.5,2,2))
            matplot(colnames(junctions[[1]]),t(junctions[[2]]),type='l',lty=1,col='darkolivegreen4',yaxt='n',
                 ylab='',xlab='Position (Morgans)',cex.axis=1.5,cex.lab=1.5)    
            abline(v=dmipos,lty=2,lwd=1.5)                 
            a.tracts = cbind(tracts$a.lens[1,],tracts$a.anc[1,])
            par(mar=c(5,4.5,2,0))
            plot(density(a.tracts[which(a.tracts[,2]==0)]),col="dodgerblue",lwd=1,
                 xlab="Tract Length (Morgans)",main='Locus A',cex.main=1.3,cex.lab=1.5,cex.axis=1.5,ylim=c(0,2))
            lines(density(a.tracts[which(a.tracts[,2]==1)]),col="darkorange",lwd=1,
                  xlab="Tract Length (Morgans)")                    
            for (i in 2:reps) {       
                a.tracts = cbind(tracts$a.lens[i,],tracts$a.anc[i,])                              
                lines(density(a.tracts[which(a.tracts[,2]==0)]),col="dodgerblue",lwd=1)                   
                lines(density(a.tracts[which(a.tracts[,2]==1)]),col="darkorange",lwd=1)   
            }
            b.tracts = cbind(tracts$b.lens[1,],tracts$b.anc[1,])
            par(mar=c(5,2.5,2,2))
            plot(density(b.tracts[which(b.tracts[,2]==0)]),col="dodgerblue",lwd=1,yaxt='n',
                 xlab="Tract Length (Morgans)",main='Locus B',cex.main=1.3,cex.lab=1.5,cex.axis=1.5,ylim=c(0,2))
            lines(density(b.tracts[which(b.tracts[,2]==1)]),col="darkorange",lwd=1,
                  xlab="Tract Length (Morgans)")                    
            for (i in 2:reps) {       
                b.tracts = cbind(tracts$b.lens[i,],tracts$b.anc[i,])                              
                lines(density(b.tracts[which(b.tracts[,2]==0)]),col="dodgerblue",lwd=1)                   
                lines(density(b.tracts[which(b.tracts[,2]==1)]),col="darkorange",lwd=1)   
            }               
        }           
    }        
}
 
run_comparisons = function(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi,selanc,domdmi,reps,demefile,episfile,outfile.a,outfile.b,chrnum,model) {
    # The neutral simulation:     
    run_dfuse(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi=0,selanc=0,domdmi=0,reps,demefile,episfile,outfile.a,chrnum,model)
    compute_summary_stats(outfile.a, dmipos, chrnum, windowsize)            
    # The selection simulation:         
    run_dfuse(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi,selanc,domdmi,reps,demefile,episfile,outfile.b,chrnum,model)
    compute_summary_stats(outfile.b, dmipos, chrnum, windowsize)  
}      
                   
no_list_func = function(data,type) {    
	new.matrix = matrix(nrow=nrow(data),ncol=ncol(data))   
	if (type=="num") {   
	    for (i in 1:ncol(data)) {   
		    new.matrix[,i] = as.numeric(paste(unlist(data[,i])))   
	    }   
	}   
	if (type=="char") {   
		for (i in 1:ncol(data)) {   
		    new.matrix[,i] = as.character(paste(unlist(data[,i])))   
	    }   
	}   
	return(new.matrix)   
}    
  
grid_arrange_shared_legend = function(..., ncol = 2, nrow = 3, position = c("bottom", "right")) {
	plots <- list(...)
	position <- match.arg(position)
	g <- ggplotGrob(plots[[1]] + 
	theme(legend.position = position))$grobs
	legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
	lheight <- sum(legend$height)
	lwidth <- sum(legend$width)
	gl <- lapply(plots, function(x) x +
	theme(legend.position = "none"))
	gl <- c(gl, ncol = ncol, nrow = nrow)
	combined <- switch(position,
	                   "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
	                   legend,ncol = 1,
					   heights = unit.c(unit(1, "npc") - lheight, lheight)),
					   "right" = arrangeGrob(do.call(arrangeGrob, gl),
				       legend, ncol = 1,
					   widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
	grid.newpage()
	grid.draw(combined)
	# return gtable invisibly
	invisible(combined)
}    


add.transparency = function(colors,alpha) {
    for (i in 1:length(colors)) {
        colors[i] = paste(colors[i],alpha,sep='')
    }
    return(colors)
}
 
g_legend = function(a.gplot){
    tmp = ggplot_gtable(ggplot_build(a.gplot))
    leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend = tmp$grobs[[leg]]
    return(legend)
}


