
#-----simulation parameters-------

popsize = 100   
chrlen = 1000  
dmipos = c(.5,.5)   
chrnum = 2
reps = 10
migrate = .01  
nmarkers = 100   
tstart = 5             
seldmi = .1   
selanc = 0   
domdmi = 2  
windowsize = 20                   
demefile ="input/demefile"        
episfile ="input/episfile"       
outfile = "output/dmi"               
                             
#---------run simulations----------
 
run_dfuse(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi,selanc,domdmi,reps,demefile,episfile,outfile,chrnum)    
                                                            
#----compute summary statistics----                                        
                                         
compute_summary_stats(outfile, dmipos, chrnum)                                              
    
#-----plot summary statistics------   

plot_summary_stats(outfile, dmipos, windowsize, reps)
             
#----------plot haplotypes---------           

plot_haplotypes(outfile, dmipos, chrnum)             
                                          
#----------------------------------    

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
            matplot(colnames(ancestry),t(ancestry),type='l',lty=1,col='darkolivegreen4',   
                    ylab='Ancestry Proportion',xlab='Position (Morgans)',ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
            abline(v=dmipos,lty=2,lwd=1.5)
            par(mar=c(5,4.5,3,2))
            matplot(colnames(junctions),t(junctions),type='l',lty=1,col='darkolivegreen4',
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

                              
run_dfuse = function(popsize, chrlen, migrate, nmarkers, tstart, dmipos, seldmi, selanc, domdmi, reps, demefile, episfile, outfile, chrnum) {
    # Write the epistasis and deme files for dfuse    
    # DMIs on same or different chromosomes
    if (chrnum==1) { 
        episinfo = matrix(nrow=8,ncol=3)   
        episinfo[,1] = c(2,0,0,0,0,0,dmipos[1],0)   
        episinfo[,2] = c(1,1,0,1,0,0,dmipos[2],1)   
        episinfo[1,3] = 2    
    }  
    if (chrnum==2) {                                       
        episinfo = matrix(nrow=8,ncol=3)               
        episinfo[,1] = c(2,0,0,0,0,0,dmipos[1],0)        
        episinfo[,2] = c(1,1,0,1,0,1,dmipos[2],1)    
        episinfo[1,3] = 2                             
    }                                                   
    demeinfo = c(1, popsize, popsize)                        
    write(demeinfo,file=demefile,ncolumns=1)                                       
    write.table(episinfo,file=episfile,row.names=F,col.names=F,quote=F,na="")           
    # Run dfuse                                                                   
    system(paste("dfuse -d",demefile,"-g",tstart,"-m",migrate,                  
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
               
compute_summary_stats = function(outfile,dmipos,chrnum) {
    load(paste(outfile,".haplo.RDATA",sep=''))   
    reps = length(unique(haplotypes[,1]))       
    if (reps==1) {  
        ancestry = ancestry_proportion(haplotypes, chrnum)  
        junctions = junction_count(haplotypes, chrnum)  
        tracts = tract_lengths(haplotypes, dmipos, chrnum) 
    } 
    if (reps>1) {
        if (chrnum==1) {
            ancestry = matrix(ncol=ncol(haplotypes)-4,nrow=reps)  
            junctions = matrix(ncol=ncol(haplotypes)-4,nrow=reps)  
            colnames(ancestry) = colnames(haplotypes)[-c(1:4)]
            colnames(junctions) = colnames(haplotypes)[-c(1:4)]
        }
        if (chrnum==2) {
            colnames = colnames(haplotypes)
            names = colnames[-c(1:min(which(colnames==1)))]
            ancestry.a = matrix(ncol=length(names),nrow=reps)  
            ancestry.b = matrix(ncol=length(names),nrow=reps) 
            junctions.a = matrix(ncol=length(names),nrow=reps) 
            junctions.b = matrix(ncol=length(names),nrow=reps) 
            colnames(ancestry.a) = names
            colnames(junctions.a) = names    
        } 
        a.lentable = matrix(ncol=popsize*2,nrow=reps)  
        b.lentable = matrix(ncol=popsize*2,nrow=reps)  
        a.anctable = matrix(ncol=popsize*2,nrow=reps)  
        b.anctable = matrix(ncol=popsize*2,nrow=reps) 
        for (i in 1:reps) { 
            haprep = haplotypes[which(haplotypes[,1]==i),]
            if (chrnum==1) {  
                ancestry[i,] = ancestry_proportion(haprep,chrnum)
                junctions[i,] = junction_count(haprep, chrnum)
            }
            tracts = tract_lengths(haprep,dmipos,chrnum)  
            a.lentable[i,c(1:length(tracts$a.lens))] = tracts$a.lens  
            b.lentable[i,c(1:length(tracts$b.lens))] = tracts$b.lens   
            a.anctable[i,c(1:length(tracts$a.anc))] = tracts$a.anc    
            b.anctable[i,c(1:length(tracts$b.anc))] = tracts$b.anc   
            if (chrnum==2) {     
                ancestry = ancestry_proportion(haprep,chrnum)
                junctions = junction_count(haprep, chrnum)             
                ancestry.a[i,] = ancestry[[1]]
                ancestry.b[i,] = ancestry[[2]] 
                junctions.a[i,] = junctions[[1]] 
                junctions.b[i,] = junctions[[2]]
                ancestry = list("ancestry.a"=ancestry.a, "ancestry.b"=ancestry.b)  
                junctions = list("junctions.a"=junctions.a,"junctions.b"=junctions.b)                  
            }                    
        }      
        tracts = list("a.lens"=a.lentable,"b.lens"=b.lentable,"a.anc"=a.anctable,"b.anc"=b.anctable)                                                                      
    }  
    summarylist = list("ancestry" = ancestry,"junctions"=junctions,"tracts"=tracts)
    save(summarylist,file=paste(outfile,".stats.RDATA",sep=''))       
}             
          
ancestry_proportion = function(haprep,chrnum) {
    if (chrnum==1) {
        ancestry = apply(haprep[,-c(1:4)],2,mean)
    }
    if (chrnum==2) {
        colnames = colnames(haprep)
        firstchrom = -c(c(1:4),c((min(which(colnames==1))+1):ncol(haprep)))
        seconchrom = -c(1:min(which(colnames==1)))
        ancestry.a = apply(haprep[,firstchrom],2,mean)  
        ancestry.b = apply(haprep[,seconchrom],2,mean)
        ancestry = list("ancestry.a"=ancestry.a, "ancestry.b"=ancestry.b)
    }
    return(ancestry)
}
           
junction_count = function(haprep,chrnum) {
    if (chrnum==1) {
        diffs = apply(haprep[,-c(1:4)],1,function(X) abs(diff(X)))
        junctions = c(0,apply(t(diffs),2,sum))
    }
    if (chrnum==2) {
        colnames = colnames(haprep)
        firstchrom = -c(c(1:4),c((min(which(colnames==1))+1):ncol(haprep)))
        seconchrom = -c(1:min(which(colnames==1)))
        diffs.a = apply(haprep[,firstchrom],1,function(X) abs(diff(X)))
        diffs.b = apply(haprep[,seconchrom],1,function(X) abs(diff(X)))      
        junctions.a = c(0,apply(t(diffs.a),2,sum))
        junctions.b = c(0,apply(t(diffs.b),2,sum))    
        junctions = list("junctions.a"=junctions.a,"junctions.b"=junctions.b)     
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
  




