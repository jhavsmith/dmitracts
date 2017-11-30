import matplotlib.pyplot as plt
import subprocess as sp
import pickle as pkl
import pandas as pd
import numpy as np
import math
import pdb

#------------------Neutral Simulations----------------------

def run_dfuse(ratevec, tstart, popsize, markers):
    codedir = '../dfuse_files/'
    make_deme_file(popsize)
    for rate in ratevec:
        runstr = "dfuse -d "+codedir+"inputfiles/demefile_dmi -r 1 -g "+ \
                 str(tstart)+" -m "+str(rate*2)+" -G "+str(tstart)+" -o "+ \
                 codedir+"outputfiles/dmi"+" -l "+str(markers)
        sp.run(runstr,shell=True)
        haptable = pd.read_csv(codedir+"outputfiles/dmi.haplo",skiprows=[0])
        haptable.drop(haptable.columns[0:8], axis=1, inplace=True)
        haptable = haptable[haptable.isin(['a','b']).all(1)]
        mainfile = codedir+"outputfiles/dmi.main"
        pos = open(mainfile).readline().rstrip()
        pos = np.asarray(pos.split(' ')[2:])
        haptable.columns = pos
        haptable[haptable.isin(['a'])] = 1
        haptable[haptable.isin(['b'])] = 2
        hapfile = codedir+"tempfiles/haptable_mig_"+str(rate)+".csv"
        haptable.to_csv(hapfile,index=False)

def make_deme_file(popsize):
    inputdir = '../dfuse_files/inputfiles/'
    thefile = open(inputdir+"demefile_dmi", 'w+')
    input_list = [''] * 4
    input_list[0] = '1'
    input_list[1] = str(int(popsize))
    input_list[2] = str(int(popsize))
    for item in input_list:
        thefile.write("%s\n" % item)
    thefile.close()

def get_dfuse_tracts(ratevec,npts,maxlen,chrlen):
    codedir = "../dfuse_files/tempfiles/"
    bins = np.array([np.arange(0,maxlen*(1+.5/npts),float(maxlen)/npts)])[0]
    bintable = pd.DataFrame(index=bins[:-1],columns=ratevec)
    for rate in ratevec:
        chromdict = make_dfuse_chroms(rate,chrlen)
        tractdict = make_tracts(chromdict,chrlen)
        tractbins = make_tractbins(tractdict,chrlen,maxlen,npts)
        popahist = tractbins.get('popahist')
        popbhist = tractbins.get('popbhist')
        relpopa = popahist/np.sum(popahist)
        relpopb = popbhist/np.sum(popbhist)
        bintable[rate] = np.mean(np.vstack((relpopa,relpopb)),axis=0)
    binfile = codedir+"dfuse_bintable.pkl"
    pkl.dump(bintable,open(binfile, "wb" ))

def make_dfuse_chroms(rate, chrlen):
    codedir = "../dfuse_files/tempfiles/"
    haptable = pd.read_csv(codedir+"haptable_mig_"+str(rate)+".csv")
    chroms = np.asarray(haptable)
    pos = np.asarray(haptable.columns,dtype='float64')*chrlen
    chromdict = {"chroms":chroms,"pos":pos}
    return(chromdict)

def plot_figure_dfuse(ratevec,tstart,Ls,pop,maxlen,npts,chrlen):
    """Plots the predicted tract length distributions against the
    simulated tract length distributions from dfuse."""
    bindir = '/home/joelsmith/Projects/dmis/code/dfuse_files/tempfiles/'
    bintable = pkl.load(open(bindir+"dfuse_bintable.pkl", "rb"))
    bins = np.append(bintable.index.values,1)
    rows = bintable.shape[0]
    xAxis = np.array(list(range(1,rows+1)))[:rows-1] / np.array(rows)
    for rate in ratevec:
        mig = continuous_mig_hybridzone(rate,tstart)
        model = demographic_model(mig)
        tractsOut = model.expectperbin(Ls,pop,bins)
        tractDist = np.array(tractsOut)
        relTractDist = tractDist / np.sum(tractDist)
        obs = np.asarray(np.where(np.logical_not(bintable[rate]==0)))[0][:-1]
        obsvalues = bintable.iloc[obs][rate]
        plt.plot(bins[:npts],np.log(relTractDist[:npts:]),label='_nolegend_')
        plt.scatter(x=bintable.index[obs],y=np.log(obsvalues))
        plt.xlabel('Tract Length (Morgans)')
        plt.ylabel('log(Relative Frequency)')
        plt.legend(loc='upper right',title='Migration Rate')

def plot_image_dfuse(ratevec):
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
    codedir = "../dfuse_files/tempfiles/"
    chroms1 = np.asarray(pd.read_csv(codedir+"haptable_mig_"+str(ratevec[0])+".csv"))
    chroms2 = np.asarray(pd.read_csv(codedir+"haptable_mig_"+str(ratevec[1])+".csv"))
    chroms3 = np.asarray(pd.read_csv(codedir+"haptable_mig_"+str(ratevec[2])+".csv"))
    mycmap = plt.cm.get_cmap('tab10')
    ax1.imshow(chroms1,extent=[0,1,0,chroms1.shape[0]],aspect='auto',cmap=mycmap,vmin=1,vmax=7)
    ax2.imshow(chroms2,extent=[0,1,0,chroms2.shape[0]],aspect='auto',cmap=mycmap,vmin=1,vmax=7)
    ax3.imshow(chroms3,extent=[0,1,0,chroms3.shape[0]],aspect='auto',cmap=mycmap,vmin=1,vmax=7)
    ax1.set_title(ratevec[0])
    ax2.set_title(ratevec[1])
    ax3.set_title(ratevec[2])
    ax3.set_xlabel('Position (Morgans)')


#------------------Simulations with DMIs----------------------

def print_fitness_matrix(seldmi,selanc,domdmi):
    if domdmi==0:
        h0=1
        h1=1
    if domdmi==1:
        h0=1/4
        h1=1/2
    if domdmi==2:
        h0=0
        h1=1/2
    if domdmi==3:
        h0=0
        h1=0
    fitness_matrix = np.zeros((3,3))
    fitness_matrix[[0,2],[0,2]] = 1
    fitness_matrix[[1,2],[0,1]] = 1-(selanc/2)
    fitness_matrix[[0,1],[1,2]] = 1-seldmi*h1
    fitness_matrix[1,1] = 1-seldmi*h0
    fitness_matrix[2,0] = 1-selanc
    fitness_matrix[0,2] = 1-seldmi
    print('\nFitness Matrix:\n')
    print(fitness_matrix)

def run_dfuse_dmis(rate, tstart, popsize, markers, dmipos, seldmi, selanc, domdmi):
    print_fitness_matrix(seldmi, selanc, domdmi)
    codedir = '../dfuse_files/'
    make_deme_file(popsize)
    make_epis_file(dmipos)
    runstr = "dfuse -d "+codedir+"inputfiles/demefile_dmi -r 1 -g "+ \
             str(tstart)+" -m "+str(rate*2)+" -G "+str(tstart)+" -o "+ \
             codedir+"outputfiles/dmi"+" -l "+str(markers)+" -e "+ \
             codedir+"inputfiles/episfile_dmi -c "+str(seldmi)+ \
             " -a "+str(selanc)+" -R "+str(domdmi)+" -S 'fertility'"
    sp.run(runstr,shell=True)
    haptable = pd.read_csv(codedir+"outputfiles/dmi.haplo",skiprows=[0])
    haptable.drop(haptable.columns[0:8], axis=1, inplace=True)
    haptable = haptable[haptable.isin(['a','b']).all(1)]
    mainfile = codedir+"outputfiles/dmi.main"
    pos = open(mainfile).readline().rstrip()
    pos = np.asarray(pos.split(' ')[2:])
    haptable.columns = pos
    haptable[haptable.isin(['a'])] = 1
    haptable[haptable.isin(['b'])] = 2
    hapfile = codedir+"tempfiles/haptable_mig_"+str(rate)+".csv"
    haptable.to_csv(hapfile,index=False)

def make_epis_file(dmipos):
    inputdir = '../dfuse_files/inputfiles/'
    thefile = open(inputdir+"episfile_dmi", 'w+')
    input_list = [''] * 8
    input_list[0] = '2 1 2'
    input_list[1] = '0 1'
    input_list[2] = '0 0'
    input_list[3] = '0 1'
    input_list[4] = '0 0'
    input_list[5] = '0 0'
    input_list[6] = str(dmipos[0])+' '+str(dmipos[1])
    input_list[7] = '0 1'
    for item in input_list:
        thefile.write("%s\n" % item)
    thefile.close()

def get_dfuse_tracts_dmis(rate,npts,maxlen,chrlen,dmipos):
    codedir = "../dfuse_files/tempfiles/"
    bins = np.array([np.arange(0,maxlen*(1+.5/npts),float(maxlen)/npts)])[0]
    chromdict = make_dfuse_chroms(rate,chrlen)
    tractdict = make_tracts_dmis(chromdict,chrlen,dmipos)
    tractbinsa = make_tractbins_dmis(tractdict['loca'],chrlen,maxlen,npts)
    tractbinsb = make_tractbins_dmis(tractdict['locb'],chrlen,maxlen,npts)
    apopahist = tractbinsa.get('popahist')
    apopbhist = tractbinsa.get('popbhist')
    bpopahist = tractbinsb.get('popahist')
    bpopbhist = tractbinsb.get('popbhist')
    bintablea = pd.DataFrame(data={'col1':apopahist,'col2':apopbhist})
    bintableb = pd.DataFrame(data={'col1':bpopahist,'col2':bpopbhist})
    binfilea = codedir+"dfuse_bintable_a.pkl"
    binfileb = codedir+"dfuse_bintable_b.pkl"
    pkl.dump(bintablea,open(binfilea, "wb" ))
    pkl.dump(bintableb,open(binfileb, "wb" ))

def make_tractbins_dmis(tractdict,chrlen,maxlen,npts):
    bins = np.array([np.arange(0,maxlen*(1+.5/npts),float(maxlen)/npts)])[0]
    popa = np.reshape(tractdict['popa'],(tractdict['popa'].shape[1],2))
    popb = np.reshape(tractdict['popb'],(tractdict['popb'].shape[1],2))
    popaMorgans = popa[:,0]/chrlen
    popbMorgans = popb[:,0]/chrlen
    countsa = np.histogram(popaMorgans,bins)[0]
    countsb = np.histogram(popbMorgans,bins)[0]
    popahist = countsa/(np.sum(countsa)+np.sum(countsb))
    popbhist = countsb/(np.sum(countsa)+np.sum(countsb))
    tractbins = {'popahist':popahist,'popbhist':popbhist}
    return(tractbins)

def plot_figure_dfuse_dmis(ratevec,tstart,Ls,pop,maxlen,npts,chrlen):
    """Plots the predicted tract length distributions against the
    simulated tract length distributions from dfuse."""
    bindir = '/home/joelsmith/Projects/dmis/code/dfuse_files/tempfiles/'
    bintablea = pkl.load(open(bindir+"dfuse_bintable_a.pkl", "rb"))
    bintableb = pkl.load(open(bindir+"dfuse_bintable_b.pkl", "rb"))
    bins = np.append(bintablea.index.values,1)
    rows = bintablea.shape[0]
    xAxis = np.array(list(range(1,rows+1)))[:rows-1] / np.array(rows)
    # tracts
    mig = continuous_mig_hybridzone(rate,tstart)
    model = demographic_model(mig)
    tractsOut = model.expectperbin(Ls,pop,bins)
    tractDist = np.array(tractsOut)
    relTractDist = tractDist / np.sum(tractDist)
    #
    f, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,sharey=True)
    ax1.plot(bins[:npts],np.log(relTractDist[:npts:]),label='_nolegend_')
    ax2.plot(bins[:npts],np.log(relTractDist[:npts:]),label='_nolegend_')
    obsa = np.asarray(np.where(np.logical_not(bintablea['col1']==0)))[0][:-1]
    obsb = np.asarray(np.where(np.logical_not(bintablea['col2']==0)))[0][:-1]
    obsvaluesa = bintablea.iloc[obsa]['col1']
    obsvaluesb = bintablea.iloc[obsb]['col2']
    ax1.scatter(x=bintablea.index[obsa],y=np.log(obsvaluesa))
    ax1.scatter(x=bintablea.index[obsb],y=np.log(obsvaluesb))
    ax1.set_xlabel('Tract Length (Morgans)')
    ax1.set_ylabel('log(Relative Frequency)')
    obsa = np.asarray(np.where(np.logical_not(bintableb['col1']==0)))[0][:-1]
    obsb = np.asarray(np.where(np.logical_not(bintableb['col2']==0)))[0][:-1]
    obsvaluesa = bintableb.iloc[obsa]['col1']
    obsvaluesb = bintableb.iloc[obsb]['col2']
    ax2.scatter(x=bintableb.index[obsa],y=np.log(obsvaluesa))
    ax2.scatter(x=bintableb.index[obsb],y=np.log(obsvaluesb))
    ax2.set_xlabel('Tract Length (Morgans)')





def make_tracts_dmis(chromdict,chrlen,dmipos):
    chroms = chromdict['chroms']
    pos = np.int_(chromdict['pos'])
    loca = np.where(pos==dmipos[0]*chrlen)[0]
    locb = np.where(pos==dmipos[1]*chrlen)[0]
    switches = np.apply_along_axis(np.diff,1,chroms)
    switches[:,0] = 1
    index = 0
    lentablea = np.zeros((switches.shape[0],2),dtype='int')
    lentableb = np.zeros((switches.shape[0],2),dtype='int')
    for row in range(chroms.shape[0]):
        rowswitches = np.append(np.nonzero(switches[row,:])[0],len(pos)-1)
        lefta = rowswitches[np.max(np.where(rowswitches<loca)[0])]
        righta = rowswitches[np.min(np.where(rowswitches>loca)[0])]
        leftb = rowswitches[np.max(np.where(rowswitches<locb)[0])]
        rightb = rowswitches[np.min(np.where(rowswitches>locb)[0])]
        lena = pos[righta]-pos[lefta]
        lenb = pos[rightb]-pos[leftb]
        anca = chroms[row,loca]
        ancb = chroms[row,locb]
        lentablea[index:index+1,] = np.array([lena,anca])
        lentableb[index:index+1,] = np.array([lenb,ancb])
        index = index+1
    locapopa = lentablea[np.where(lentablea[:,1]==1),]
    locapopb = lentablea[np.where(lentablea[:,1]==2),]
    locbpopa = lentableb[np.where(lentableb[:,1]==1),]
    locbpopb = lentableb[np.where(lentableb[:,1]==2),]
    loca = {"popa":locapopa,"popb":locapopb}
    locb = {"popa":locbpopa,"popb":locbpopb}
    tractdict = {"loca":loca,"locb":locb}
    return(tractdict)

def plot_image_dfuse_dmis(rate):
    codedir = "../dfuse_files/tempfiles/"
    chroms = np.asarray(pd.read_csv(codedir+"haptable_mig_"+str(rate)+".csv"))
    mycmap = plt.cm.get_cmap('tab10')
    plt.imshow(chroms,extent=[0,1,0,chroms.shape[0]],aspect='auto',cmap=mycmap,vmin=1,vmax=7)
    plt.xlabel('Position (Morgans)')
    plt.ylabel('Chromosome')
