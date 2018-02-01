
#-----simulation parameters-------

run_dfuse(
popsize = 10,
chrlen = 1000,
migrate = .01,
markers = 100,
tstart = 5,            
dmipos = c(.3,.7),
seldmi = 0,
selanc = 0,
domdmi = 2,
reps = 1,
demefile ="input/demefile",
episfile ="input/episfile",
outfile = "output/dmi") 

#----compute summary statistics----

compute_summary_stats(outfile)

#----------------------------------

run_dfuse = function(popsize, chrlen, migrate, markers, tstart, dmipos, seldmi, selanc, domdmi, reps, demefile, episfile, outfile) {
    recover()
    # Write the deme and epistasis files for dfuse
    demeinfo = c(1, popsize, popsize)   
    episinfo = matrix(nrow=8,ncol=3)   
    episinfo[,1] = c(2,0,0,0,0,0,dmipos[1],0)   
    episinfo[,2] = c(1,1,0,1,0,0,dmipos[2],1)   
    episinfo[1,3] = c(2)    
    write(demeinfo,file=demefile,ncolumns=1)     
    write.table(episinfo,file=episfile,row.names=F,col.names=F,quote=F,na="")   
    # Run dfuse
    system(paste("dfuse -d",demefile,"-r 1 -g",tstart,"-m",migrate,   
                 "-G",tstart,"-o",outfile,"-l",markers,   
                 "-e",episfile,"-c",seldmi,"-a",selanc,   
                 "-R",domdmi,"-S 'fertility' -r",reps))             
    # Get the stuff we want from the dfuse output   
    positions = scan(file=paste(outfile,".haplo",sep=''),nlines=1,what="character")
    listtable = read.table(paste(outfile,".haplo",sep=''),sep=',')  
    # Convert to matrix
    mattable = no_list_func(listtable,"char")             
    haptable = mattable[,c(c(1:4),c(9:ncol(mattable)))]   
    if (length(which(haptable[,5]=="NA"))>0) {            
        haptable = haptable[-which(haptable[,5]=="NA"),]    
    }
    haptable[which(haptable=="a")] = 0   
    haptable[which(haptable=="b")] = 1  
    haplotypes = apply(haptable[-1,],2,as.numeric)
    colnames(haplotypes) = haptable[1,] 
    save(haplotypes,file=paste(outfile,"haplotypes.RDATA",sep=''))
}   
               
compute_summary_stats = function(outfile,dmipos) { 
    load(paste(outfile,"haplotypes.RDATA",sep=''))  
    reps = length(unique(haplotypes[,1]))      
    if (reps==1) {
        ancestry = ancestry_proportion(haplotypes)
        junctions = junction_count(haplotypes)
        tracts = tract_lengths(haplotypes)
    }
    if (reps>1) {
        ancestable = matrix(ncol=ncol(haplotypes)-4,nrow=reps)
        juncttable = matrix(ncol=ncol(haplotypes)-4,nrow=reps)
        for (i in 1:reps) {
            haprep = haplotypes[which(haplotypes[,1]==i),]
            ancestable[i,] = ancestry_proportion(haprep)
            juncttable[i,] = junction_count(haprep)
        }                                                   
    }         
}             
          
ancestry_proportion = function(haprep) {
    ancestry = apply(haprep[,-c(1:4)],2,mean)
    return(ancestry)
}
           
junction_count = function(haprep) {
    diffs = apply(haprep[,-c(1:4)],1,function(X) abs(diff(X)))
    junctions = c(0,apply(t(diffs),2,sum))
    return(junctions)
}
 
tract_lengths = function(haprep,dmipos) {
    recover()
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
  




