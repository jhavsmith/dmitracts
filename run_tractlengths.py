
%reset
%matplotlib

#-----simulation parameters-------
popsize = 2000
chrlen = 1000000
rate = .01
tstart = 50
dmipos = (.3,.7)
seldmi = .1
selanc = 0
domdmi = 2
#---------------------------------
ratevec = [.001,.05,.1]
markers = 1000
npts = 20
maxlen = 1
pop = 1
Ls = [1]
migscheme = 'hybrid'



exec(open("/home/joelsmith/Projects/dmis/code/dmitracts/forqs_tractlengths.py").read())
exec(open("/home/joelsmith/Projects/dmis/code/dmitracts/dfuse_tractlengths.py").read())

print_fitness_matrix(seldmi,selanc,domdmi)

#----------------DMI tract length distributions-------------------
run_dfuse_dmis(rate, tstart, popsize, markers, dmipos, seldmi, selanc, domdmi)
get_dfuse_tracts_dmis(rate, npts, maxlen, chrlen, dmipos)
plot_figure_dfuse_dmis(rate, tstart, Ls, pop, maxlen, npts, chrlen)
plot_image_dfuse_dmis(rate)

#-----------neutral tract length distributions-------------------
run_dfuse(ratevec, tstart, popsize, markers)
get_dfuse_tracts(ratevec,npts,maxlen,chrlen)
plot_figure_dfuse(ratevec, tstart, Ls, pop, maxlen, npts, chrlen)
plot_image_dfuse(ratevec)

run_forqs(ratevec, tstart, popsize, chrlen, migscheme)
get_forqs_tracts(ratevec, chrlen, maxlen, npts, migscheme)
plot_figure_forqs(ratevec, tstart, Ls, pop, maxlen, npts, migscheme)
#----------------------------------------------------------------
