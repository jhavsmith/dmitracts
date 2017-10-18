
%reset
%matplotlib
#----------parameters----------
popsize = 5000   
chrlen = 1000000
ratevec = [.001,0.03,0.05]
tstart = 30 
npts = 20 
maxlen = 1 
pop = 1 
Ls = [1]
#-------------------------------

exec(open("/home/joelsmith/Projects/dmis/code/tractlengths/forqs_tractlengths.py").read())

run_forqs(ratevec, tstart, popsize, chrlen)
get_tracts(ratevec, chrlen, maxlen, npts)
plot_figure(ratevec, tstart, Ls, pop, maxlen, npts)


