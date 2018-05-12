

source("../../functions.R")
#-----------------simulation parameters------------------

popsize = 500   
chrlen = 1000 
dmipos = c(.5,.5)  
windowsize = .2                    
chrnum = 2
reps = 100
bins = 10
migrate = .1  
nmarkers = 1000   
tstart = 100             
model = "pathway"
type = "boxplots"
ones = "no"
demefile = "input/demefile"        
episfile = "input/episfile"        
outfile = "output/dmi"         
outfile.a = "neu"                      
outfile.b = "moda"     

# neutral
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
                       
#----------------compute summary statistics--------------     

compute_summary_stats(outfile, dmipos, chrnum, windowsize)                                              
    
#--------------------plot haplotypes---------------------           

plot_haplotypes(outfile, dmipos, chrnum)      
                                                      
#-----------------plot summary statistics----------------   

plot_summary_stats(outfile, dmipos, windowsize, reps)
                  
#----------------run simulation comparison---------------                       
                      
run_comparisons(popsize,chrlen,migrate,nmarkers,tstart,dmipos,seldmi,selanc,domdmi,reps,demefile,episfile,outfile.a,outfile.b,chrnum,model)                     
                                                              
#-----plot neutral versus selection distributions--------             

plot_distributions(outfile.a,outfile.b,reps,bins,type)   

#---------------compare distributions--------------------      
 
compare_distributions(outfile.a,outfile.b,reps,title)

#---------------compare tract lengths--------------------      
 
outfile.a = "neu"                      
outfile.b = "moda"
title = "Model A"
pdf("tracts_model_a",height=9,width=9)  
compare_tracts(outfile.a,outfile.b,reps,title,ones)
dev.off()
#---------------find time to equilibrium-----------------        

find_equilibrium(popsize,chrlen,migrate,nmarkers,dmipos,seldmi,selanc,domdmi,demefile,episfile,outfile,chrnum,model)

#--------------------------------------------------------        

