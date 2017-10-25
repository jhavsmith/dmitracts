import matplotlib.pyplot as plt
import subprocess as sp
import pickle as pkl
import pandas as pd
import numpy as np
import math

exec(open("/home/joelsmith/Projects/software/tracts/tracts.py").read())

def run_forqs(ratevec, tstart, popsize, chrlen, migscheme):
    """Runs 1 of 2 different migration scenarios.

    Scenario 1:
    A target population initially composed of individuals
    of one ancestry receives migrants from the source population at
    some rate each generation to the present

    Scenario 2:
    In the first generation, 2 populations form an admixed population with
    equal ancestry contributions (.5). Both source populations then contributions
    ancestry to this admixed population at a constant rate to the present.
    """
    codedir = "/home/joelsmith/Projects/dmis/code/"
    sp.run("rm -r -f "+codedir+"forqs_files/forqs_out",shell=True)
    if migscheme == 'source':
        for rate in ratevec:
            input_list = [''] * 27
            input_list[0] = "Trajectory_Constant id_population_size_trajectory"
            input_list[1] = "value = " + str(popsize)
            input_list[3] = "Trajectory_Constant id_migration_rate_trajectory"
            input_list[4] = "value = " + str(rate)
            input_list[6] = "PopulationConfigGenerator_LinearSteppingStone popconfig_generator"
            input_list[7] = "generation_count = " + str(tstart)
            input_list[8] = "population_count = 2"
            input_list[9] = "chromosome_pair_count = 1"
            input_list[10] = "chromosome_lengths = " + str(chrlen)
            input_list[11] = "population_size = id_population_size_trajectory"
            input_list[12] = "migration_rate:from:to = id_migration_rate_trajectory 1 2"
            input_list[14] = "RecombinationPositionGenerator_Uniform rpg"
            input_list[15] = "rate = 1"
            input_list[17] = "Reporter_Population reporter_population"
            input_list[19] = "Reporter_Timer reporter_timer"
            input_list[21] = "SimulatorConfig"
            input_list[22] = "output_directory = "+codedir+"forqs_files/forqs_out/mig_"+str(rate)+"/"
            input_list[23] = "population_config_generator = popconfig_generator"
            input_list[24] = "recombination_position_generator = rpg"
            input_list[25] = "reporter = reporter_timer"
            input_list[26] = "reporter = reporter_population"
            thefile = open(codedir+"forqs_files/forqs_in/forqs_in_mig_"+str(rate)+".txt", 'w+')
            for item in input_list:
                thefile.write("%s\n" % item)
            thefile.close()
            sp.run("rm forqs.seed",shell=True)
            sp.run("forqs "+codedir+"forqs_files/forqs_in/forqs_in_mig_"+str(rate)+".txt",shell=True)
    if migscheme == 'hybrid':
        for rate in ratevec:
            input_list = [''] * 35
            input_list[0] = "Trajectory_Constant id_population_size_trajectory"
            input_list[1] = "value = " + str(popsize)
            input_list[3] = "Trajectory_Constant id_migration_rate_0"
            input_list[4] = "value = .5"
            input_list[6] = "Trajectory_Constant id_migration_rate_1"
            input_list[7] = "value = " + str(rate)
            input_list[9] = "Trajectory_GenerationComposite id_migration_rate"
            input_list[10] = "generation:trajectory = 0 id_migration_rate_0"
            input_list[11] = "generation:trajectory = 2 id_migration_rate_1"
            input_list[13] = "PopulationConfigGenerator_LinearSteppingStone popconfig_generator"
            input_list[14] = "generation_count = " + str(tstart)
            input_list[15] = "population_count = 3"
            input_list[16] = "chromosome_pair_count = 1"
            input_list[17] = "chromosome_lengths = " + str(chrlen)
            input_list[18] = "population_size = id_population_size_trajectory"
            input_list[19] = "migration_rate:from:to = id_migration_rate 1 3"
            input_list[20] = "migration_rate:from:to = id_migration_rate 2 3"
            input_list[22] = "RecombinationPositionGenerator_Uniform rpg"
            input_list[23] = "rate = 1"
            input_list[25] = "Reporter_Population reporter_population"
            input_list[27] = "Reporter_Timer reporter_timer"
            input_list[29] = "SimulatorConfig"
            input_list[30] = "output_directory = "+codedir+"forqs_files/forqs_out/mig_"+str(rate)+"/"
            input_list[31] = "population_config_generator = popconfig_generator"
            input_list[32] = "recombination_position_generator = rpg"
            input_list[33] = "reporter = reporter_timer"
            input_list[34] = "reporter = reporter_population"
            thefile = open(codedir+"forqs_files/forqs_in/forqs_in_mig_"+str(rate)+".txt", 'w+')
            for item in input_list:
                thefile.write("%s\n" % item)
            thefile.close()
            sp.run("rm forqs.seed",shell=True)
            sp.run("forqs "+codedir+"forqs_files/forqs_in/forqs_in_mig_"+str(rate)+".txt",shell=True)

def get_forqs_tracts(ratevec, chrlen, maxlen, npts, migscheme):
    bins = np.array([np.arange(0,maxlen*(1+.5/npts),float(maxlen)/npts)])[0]
    bintable = pd.DataFrame(index=bins[:-1],columns=ratevec)
    for rate in ratevec:
        data = get_data(rate, migscheme)
        chromdict = make_chroms(data['chunks'],chrlen,rate)
        tractdict = make_tracts(chromdict,chrlen)
        tractbins = make_tractbins(tractdict,chrlen,maxlen,npts)
        popahist = tractbins.get('popahist')
        bintable[rate] = popahist/np.sum(popahist)
    bindir = '/home/joelsmith/Projects/dmis/code/forqs_files/forqs_temp/'
    pkl.dump(bintable,open(bindir+"forqs_bintable.pkl", "wb" ))

def get_data(rate, migscheme):
    tempdir = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_temp/"
    outdir = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_out/mig_"+str(rate)+"/"
    if migscheme == 'hybrid':
        sp.run("tail -n+3 "+outdir+"population_final_pop3.txt > "+tempdir+"pop_chunks.txt",shell=True)
    if migscheme == 'source':
        sp.run("tail -n+3 "+outdir+"population_final_pop2.txt > "+tempdir+"pop_chunks.txt",shell=True)
    sp.run("sed 's/[^0-9]/ /g' "+tempdir+"pop_chunks.txt > "+tempdir+"pop_data.txt",shell=True)
    sp.run("cat "+tempdir+"pop_data.txt | awk '{ print NF }' > "+tempdir+"pop_nums.txt",shell=True)
    longrow = max(open(tempdir+"pop_data.txt", 'r'), key=len)
    cn = list(range(1,len(longrow.split())+1))
    chunks = pd.read_table(tempdir+"pop_data.txt", names=cn, delim_whitespace=True)
    chunks = chunks.as_matrix()
    poptext = np.loadtxt(tempdir+"pop_nums.txt")
    popnums = poptext[np.nonzero(poptext)]
    sp.run("rm "+tempdir+"pop_nums.txt",shell=True)
    sp.run("rm "+tempdir+"pop_data.txt",shell=True)
    sp.run("rm "+tempdir+"pop_chunks.txt",shell=True)
    data = {'chunks':chunks,'popnums':popnums}
    return(data)

def make_chroms(chunks, chrlen, rate):
    indpos = chunks[:,2::2]
    pos =  np.int_([value for value in np.unique(indpos) if not math.isnan(value)])
    ids = chunks[:,1::2]
    chrommatrix = np.zeros(shape=(len(ids),len(pos)))
    for row in range(chrommatrix.shape[0]):
        chunkindex = len(indpos[row,np.logical_not(np.isnan(indpos[row,]))])
        if chunkindex==0:
            chrommatrix[row,:] = ids[row,0]
        chunkstart = 0
        for chunk in range(chunkindex):
            chunkend = list(np.where(pos==indpos[row,chunk])[0])[0]
            chrommatrix[row,chunkstart:chunkend] = ids[row,chunk]
            chunkstart = chunkend
        truevalues = ids[row,np.logical_not(np.isnan(ids[row,]))]
        chrommatrix[row,np.where(chrommatrix[row,]==0)] = truevalues[-1]
    chroms = np.copy(chrommatrix)
    n = chrommatrix.shape[0]
    chroms[np.where(chroms<(n-1))] = 1
    chroms[np.where(chroms>(n-1))] = 2
    haptable = pd.DataFrame(data = chroms,columns=pos,dtype='i')
    tempdir = "/home/joelsmith/Projects/dmis/code/forqs_files/forqs_temp/"
    hapfile = tempdir+"haptable_mig_"+str(rate)+".csv"
    haptable.to_csv(hapfile,index=False)
    pos = np.append(pos,chrlen)
    chromdict = {"chroms":chroms,"pos":pos}
    return(chromdict)

def make_tracts(chromdict,chrlen):
    chroms = chromdict['chroms']
    pos = chromdict['pos']
    switches = np.apply_along_axis(np.diff,1,chroms)
    switchvec = np.apply_over_axes(np.count_nonzero,switches,1)
    tractstot = np.sum(switchvec)+len(np.nonzero(switchvec)[0])
    switchmatrix = switches[np.nonzero(switchvec)[0],]
    hapsmatrix = chroms[np.nonzero(switchvec)[0],]
    tractsmatrix = np.zeros((tractstot,2),dtype='int')
    noswitches = chroms[np.where(switchvec==0),0][0]
    fulltracts = np.transpose(np.vstack((np.repeat(chrlen,len(noswitches)),noswitches)))
    index = 0
    for row in range(switchmatrix.shape[0]):
        switchsites = np.logical_not(switchmatrix[row,]==0)
        switchpos = pos[np.append(np.where(switchsites),len(pos)-1)]
        switchpos = np.append(0,switchpos)
        hapids = hapsmatrix[row,np.append(np.where(switchsites),len(pos)-2)]
        haplengths = np.diff(switchpos)
        hapstack = np.transpose(np.vstack((haplengths,hapids)))
        tractsmatrix[index:(index+len(hapids)),] = hapstack
        index = index + (len(hapids))
    tractsmatrix = np.vstack((tractsmatrix,fulltracts))
    tractsmatrix = tractsmatrix.astype(int)
    popa = tractsmatrix[np.where(tractsmatrix[:,1]==1),]
    popb = tractsmatrix[np.where(tractsmatrix[:,1]==2),]
    tractdict = {"popa":popa,"popb":popb}
    return(tractdict)

def make_tractbins(tractdict,chrlen,maxlen,npts):
    bins = np.array([np.arange(0,maxlen*(1+.5/npts),float(maxlen)/npts)])[0]
    popa = np.reshape(tractdict['popa'],(tractdict['popa'].shape[1],2))
    popb = np.reshape(tractdict['popb'],(tractdict['popb'].shape[1],2))
    popaMorgans = popa[:,0]/chrlen
    popbMorgans = popb[:,0]/chrlen
    popahist = np.histogram(popaMorgans,bins,density=True)[0]
    popbhist = np.histogram(popbMorgans,bins,density=True)[0]
    tractbins = {'popahist':popahist,'popbhist':popbhist}
    return(tractbins)

def plot_figure_forqs(ratevec,tstart,Ls,pop,maxlen,npts,migscheme):
    """Plots the predicted tract length distributions against the
    simulated tract length distributions from forqs."""
    bindir = '/home/joelsmith/Projects/dmis/code/forqs_files/forqs_temp/'
    bintable = pkl.load(open(bindir+"forqs_bintable.pkl", "rb"))
    bins = np.append(bintable.index.values,1)
    rows = bintable.shape[0]
    xAxis = np.array(list(range(1,rows+1)))[:rows-1] / np.array(rows)
    for rate in ratevec:
        if migscheme=='hybrid':
            mig = continuous_mig_hybridzone(rate,tstart)
        if migscheme=='source':
            mig = continuous_mig(rate,tstart)
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
