import matplotlib.pyplot as plt
import subprocess as sp
import pickle as pkl
import pandas as pd
import numpy as np
import math
import pdb

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
        hapfile = codedir+"tempfiles/haptable_mig_"+str(rate)+".pkl"
        pkl.dump(haptable,open(hapfile, "wb" ))

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
    haptable = pkl.load(open(codedir+"haptable_mig_"+str(rate)+".pkl", "rb"))
    chroms = np.asarray(haptable)
    pos = np.asarray(haptable.columns,dtype='float64')*chrlen
    chromdict = {"chroms":chroms,"pos":pos}
    return(chromdict)

def plot_figure_dfuse(ratevec,tstart,Ls,pop,maxlen,npts,migscheme):
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

#------------------Simulations with DMIs----------------------

def run_dfuse_dmis(ratevec, tstart, popsize, markers, dmipos, seldmi, selanc, domdmi):
    print_fitness_matrix(seldmi, selanc, domdmi)
    codedir = '../dfuse_files/'
    make_deme_file(popsize)
    make_epis_file(dmipos)
    for rate in ratevec:
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
        hapfile = codedir+"tempfiles/haptable_mig_"+str(rate)+".pkl"
        pkl.dump(haptable,open(hapfile, "wb" ))

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
    fitness_matrix[[1,2],[0,1]] = 1-seldmi*h1
    fitness_matrix[[0,1],[1,2]] = 1-(selanc/2)
    fitness_matrix[1,1] = 1-seldmi*h0
    fitness_matrix[2,0] = 1-seldmi
    fitness_matrix[0,2] = 1-selanc
    print('\nFitness Matrix:\n')
    print(fitness_matrix)
