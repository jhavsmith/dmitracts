
%reset
%matplotlib

#-----neutral parameters-------
popsize = 1000
chrlen = 1000000
markers = 1000
ratevec = [.001,.03,.05]
tstart = 30
npts = 20
maxlen = 1
pop = 1
Ls = [1]
migscheme = 'hybrid'
#------selection parameters-----
dmipos = (.3,.7)
seldmi = .9
selanc = 0
domdmi = 0
#-------------------------------

exec(open("/home/joelsmith/Projects/dmis/code/tractlengths/forqs_tractlengths.py").read())
exec(open("/home/joelsmith/Projects/dmis/code/tractlengths/dfuse_tractlengths.py").read())

#----------------dmi tract length distributions-------------------
run_dfuse_dmis(ratevec, tstart, popsize, markers, dmipos, seldmi, selanc, domdmi)
get_dfuse_tracts(ratevec,npts,maxlen,chrlen)
plot_figure_dfuse(ratevec, tstart, Ls, pop, maxlen, npts, migscheme)

#-----------neutral tract length distributions-------------------
run_dfuse(ratevec, tstart, popsize, markers)
get_dfuse_tracts(ratevec,npts,maxlen,chrlen)
plot_figure_dfuse(ratevec, tstart, Ls, pop, maxlen, npts, migscheme)

run_forqs(ratevec, tstart, popsize, chrlen, migscheme)
get_forqs_tracts(ratevec, chrlen, maxlen, npts, migscheme)
plot_figure_forqs(ratevec, tstart, Ls, pop, maxlen, npts, migscheme)
#----------------------------------------------------------------
