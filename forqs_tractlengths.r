
run.multiple.simulations = function(params, mig.vec) {
    temploc = '/home/joelsmith/Projects/dmis/code/forqs_files/forqs_temp/'
    for (i in 1:length(mig.vec)) {
        params$migration.rate = mig.vec[i]   
        run.forqs.simulation(params)
        chroms.matrix = get.chroms(params)
        tract.table.list = get.forqs.tract.length.distr(chroms.matrix)
        save(tract.table.list,file=paste(temploc,'tracts_',mig.vec[i],'.RDATA',sep=''))
        rm(chroms.matrix)
        rm(tract.table.list)
    }
}

plot.multiple.simulations = function(params, mig.vec) {
    forqs.list = list()
    temploc = '/home/joelsmith/Projects/dmis/code/forqs_files/forqs_temp/'
    for (i in 1:length(mig.vec)) {
        load(paste(temploc,'tracts_',mig.vec[i],".RDATA",sep=''))
        tract.bin.list = get.tract.length.bins(tract.table.list,params)
        tract.bin.plot = tract.bin.list[[1]][,2]/sum(tract.bin.list[[1]][,2])
        forqs.list[[i]] = tract.bin.plot
    }
    return(forqs.list)
}
   
run.forqs.simulation = function(params) {
    options(scipen=10)    
    outdir = paste(params$output.file,"mig_",params$migration.rate,sep='')
    forqs.input = vector()
    forqs.input[1]   = "Trajectory_Constant id_population_size_trajectory"
    forqs.input[2]  = paste("value = ",params$population.size,sep='')
    forqs.input[4]   = "Trajectory_Constant id_migration_rate_trajectory"
    forqs.input[5]  = paste("value = ",params$migration.rate,sep='')
    forqs.input[7]   = "PopulationConfigGenerator_LinearSteppingStone popconfig_generator"
    forqs.input[8]  = paste("generation_count = ",params$generations,sep='')   
    forqs.input[9]  = "population_count = 2"
    forqs.input[10] = "chromosome_pair_count = 1"
    forqs.input[11] = paste("chromosome_lengths = ",params$chromosome.length,sep='')   
    forqs.input[12] = "population_size = id_population_size_trajectory"   
    forqs.input[13] = "migration_rate:from:to = id_migration_rate_trajectory 1 2"          
    forqs.input[15] = "RecombinationPositionGenerator_SingleCrossover rpg"
    forqs.input[17] = "Reporter_Population reporter_population"
    forqs.input[19] = "Reporter_Timer reporter_timer"
    forqs.input[21] = "SimulatorConfig"  
    forqs.input[22] = paste("output_directory = ",outdir,sep='')    
    forqs.input[23] = "population_config_generator = popconfig_generator"              
    forqs.input[24] = "recombination_position_generator = rpg"    
    forqs.input[25] = "reporter = reporter_timer"    
    forqs.input[26] = "reporter = reporter_population"     
    forqs.input[which(is.na(forqs.input))] = " "
    input.file = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_in/forqs_source_target.txt"
    write(forqs.input,file=input.file,ncolumns=1)      
    system("rm forqs.seed") 
    system(paste("rm -r -f ",params$output.file,sep=''))
    system(paste("forqs ",input.file,sep=''))
}

get.forqs.tract.length.distr = function(chroms.matrix) {
    pos = chroms.matrix[1,]
    haps.matrix = chroms.matrix[-1,]    
    tract.table = c(0,0)
    for (i in 1:nrow(haps.matrix)) {
        switches = which(diff(haps.matrix[i,])!=0)
        switches.sites = c(1,switches,ncol(haps.matrix))
        hap.ids = haps.matrix[i,switches.sites[-1]]
        switches.inner.pos = pos[switches] + round((pos[switches+1] - pos[switches])/2) 
        hap.lengths = diff(c(0,switches.inner.pos,pos[length(pos)]))
        tract.table = rbind(tract.table,cbind(hap.lengths,hap.ids))
    }
    tract.table = tract.table[-1,]
    tracts.a = tract.table[which(tract.table[,2]==1),]
    tracts.b = tract.table[which(tract.table[,2]==2),]
    tract.table.list = list("tracts.a" = tracts.a,"tracts.b" = tracts.b)
    for (i in 1:2) {
        tract.table = tract.table.list[[i]]
        tract.bins = unique(tract.table[,1])
        tract.length.table = cbind(tract.bins,rep(0,length(tract.bins)))
        for (j in 1:length(tract.bins)) {
            tract.length.table[j,2] = length(which(tract.table[,1]%in%tract.bins[j]==TRUE))
        }
        tract.table.list[[i]] = tract.length.table[order(tract.length.table[,1]),]
   }
   return(tract.table.list)
}
   
get.chroms = function(params) {
    outdir = paste(params$output.file,"mig_",params$migration.rate,sep='')
    temp.loc = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_temp"
    system(paste("tail -n+3 ",outdir,"/population_final_pop2.txt > ",temp.loc,"/pop_chunks.txt",sep=''))            
    system(paste("sed 's/[^0-9]/ /g' ",temp.loc,"/pop_chunks.txt > ",temp.loc,"/pop_data.txt",sep=''))
    system(paste("cat ",temp.loc,"/pop_data.txt | awk '{ print NF }' > ",temp.loc,"/pop_nums.txt",sep=''))        
    cn <- paste("V",1:max(scan(paste(temp.loc,"/pop_nums.txt",sep=''))),sep='')
    chunks = read.table(file = paste(temp.loc,"/pop_data.txt",sep=''), col.names = cn,fill=TRUE)
    chunks = no.list.func(chunks)
    chroms.list = make.chroms(chunks)
    chroms.matrix = chroms.list[[1]]  
    pos = c(chroms.list[[2]],params$chromosome.length)
    n = nrow(chroms.matrix)
    chroms.matrix[is.element(chroms.matrix,c(0:(n-1)))] = 1
    chroms.matrix[is.element(chroms.matrix,c(n:(n*2)))] = 2
    chroms.matrix = rbind(pos,chroms.matrix)
    return(chroms.matrix)    
}
         
plot.chroms = function(chroms.matrix) {
    colors = c("slateblue","sienna1")
    image(t(chroms.matrix[-1,]),col=colors) 
}    
     
get.tract.length.bins = function(tract.table.list,params) {
    options(scipen=10)      
    for (j in 1:2) {
        tract.length.distr = tract.table.list[[j]]
        bins = seq(1,params$chromosome.length,params$bin.size)
        bin.table = cbind(bins,c((bins[-1]-1),params$chromosome.length))
        new.bins = matrix(nrow=nrow(bin.table),ncol=2)
        new.bins[,1] = bin.table[,2]#-(params$bin.size/2)
        for (i in 1:nrow(bin.table)) { 
            uppers = which(tract.length.distr[,1] >= bin.table[i,1])
            lowers = which(tract.length.distr[,1] <= bin.table[i,2])
            new.bins[i,2] = sum(tract.length.distr[intersect(uppers,lowers),2])      
        }  
        tract.table.list[[j]] = new.bins
    }
    return(tract.table.list)
}     
                                                               
make.chroms = function(data) {
    ind.pos = data[,seq(3,ncol(data),by=2)]
    pos = c(0,unique(sort(as.numeric(ind.pos[which(ind.pos != "NA")]))))
    ids = data[,seq(2,ncol(data),by=2)]
    chrom.matrix = matrix(nrow = nrow(data),ncol = (length(pos)+1))
    for (i in 1:nrow(data)) {
        chunk.index = length(which(ind.pos[i,]!="NA"))
        if (chunk.index==0) {
            chrom.matrix[i,] = as.numeric(ids[i,1])
            next
        }
        chunk.start = 1
        for (j in 1:chunk.index) {
            chunk.end = which(pos == ind.pos[i,j])
            chrom.matrix[i,chunk.start:chunk.end] = as.numeric(ids[i,j])
            chunk.start = chunk.end   
        }
        chrom.matrix[i,c(max(which(!is.na(chrom.matrix[i,]))):ncol(chrom.matrix))] = ids[i,max(which(ids[i,]!="NA"))]
    }
    chrom.matrix = apply(chrom.matrix,2,as.numeric)
    chrom.list = list("chrom.matrix"=chrom.matrix,"pos"=pos)
    return(chrom.list)
}

no.list.func = function(data) {
    new.matrix = matrix(nrow = nrow(data), ncol = ncol(data))
	for (i in 1:ncol(data)) {
        new.matrix[,i] = (paste(unlist(data[,i])))
    }
	return(new.matrix)
}

#---------------------Forqs plots--------------------------------------------
source("forqs_tractlengths.r")

params = list(  
    "population.size" = 1000,
    "generations" = 30,
    "chromosome.length" = 100000000,
    "output.file" = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_out/",
    "bin.size" = 500000)

mig.vec = c(.001,.03,.05)

run.multiple.simulations(params,mig.vec)
binlist = plot.multiple.simulations(params,mig.vec)
binmatrix = matrix(unlist(binlist),ncol=3,byrow=TRUE)
temploc = '/Projects/dmis/code/forqs_files/forqs_temp/'
write.csv(binmatrix,file=paste(temploc,'binmatrix.csv',sep=''),col.names=FALSE)    




