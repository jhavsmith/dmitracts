
import numpy as np
import itertools as it
import os
import pylab
import tkinter as Tk
import tkinter.filedialog
import pdb

try:
    from scipy.misc.common import factorial
except ImportError:
    from scipy.misc import factorial

from scipy.special import gammainc, gammaln
import sys


class demographic_model(object):
    def __init__(self, mig, psivec, max_remaining_tracts=1e-5, max_morgans=100):
        """ Migratory model takes as an input a vector containing the migration
            proportions over the last generations. Each row is a time, each
            column is a population. row zero corresponds to the current
            generation. The migration rate at the last generation (time $T$) is
            the "founding generation" and should sum up to 1. Assume that
            non-admixed individuals have been removed.

            max_remaining_tracts is the proportion of tracts that are allowed
            to be incomplete after cutoff Lambda
            (Appendix 2 in Gravel: doi: 10.1534/genetics.112.139808)
            cutoff=1-\sum(b_i)

            max_morgans is used to impose a cutoff to the number of Markov transitions.
            If the simulated morgan lengths of tracts in an infinite genome is more than
            max_morgans, issue a warning and stop generating new transitions
        """
        small=1e-10
        self.mig = mig
        (self.ngen, self.npop) = mig.shape
        self.max_remaining_tracts=max_remaining_tracts #tolerance for incomplete
        # convergence
        self.max_morgans=max_morgans
        # the total migration per generation
        self.totmig = mig.sum(axis=1)

        # test for reasonableness of migration matrix

        if abs(self.totmig[-1] - 1) > small:
            eprint("founding migration should sum up to 1. Now:", mig[-1, :],
                    "sum up to ", self.totmig[-1])
            raise ValueError("founding migration sum is not 1")

        if self.totmig[0] > small:
            eprint("migrants at last generation should be removed from",
                    "sample!")
            eprint("currently", self.totmig[0])
            raise ValueError("migrants from last generation are not removed")

        self.totmig[0] = 0

        if self.totmig[1] > small:
            eprint("migrants at penultimate generation should be removed from",
                    "sample!")
            eprint("currently", self.totmig[1])
            raise ValueError(
                    "migrants from penultimate generation are not removed")

        if ((self.totmig > 1).any() or (mig < 0).any()):
            eprint("migration rates should be between 0 and 1")
            eprint("currently", mig)
            raise ValueError("mig")
        if (mig[:-1] == 1).any():

            eprint("warning: population was completely replaced after",
                    "founding event")

        # Tracts represents ancestry tracts as a markov model
        # states are defined by source population and generation of arrival.
        # i.e. if one population contributes migrants in two generations, in the MM
        # (see fig 3 in MOLA paper)

        # identify states where migration occurred as these are the relevant
        # states in our Markov model. Each state is a tuple of the form:
        # (generation, population)

        self.states = list(map(tuple, np.array(mig.nonzero()).transpose()))
        self.nstates = len(self.states)
        self.npops = mig.shape[1]

        # get the equilibrium distribution in each state
        self.equil = np.zeros(self.nstates)
        self.stateINpop = [[] for pop in range(self.npops)]
        self.stateOUTpop = [[] for pop in range(self.npops)]

        for i, state in enumerate(self.states):
            self.stateINpop[state[1]].append(i)
            for other in range(1, self.npops+1):
                self.stateOUTpop[(state[1]+other)%self.npops].append(i)
            self.equil[i] = mig[state]*(1-self.totmig)[1:state[0]].prod()

        self.equil /= self.equil.sum()

        # calculate the ancestry proportions as a function of time
        self.proportions = np.zeros(mig.shape)

        # could be optimized using array operations and precomputing survivals
        for pop in range(self.npop):
            for time in range(self.ngen):
                for g in range(time, self.ngen):
                    self.proportions[time, pop] += \
                        mig[g, pop]*(1-self.totmig)[time:g].prod()

        # calculate the transition matrix
        self.dicTpopTau = {}

        # we could precompute prod
        for (t, pop) in self.states:
            for tau in range(t):
                prod = (1-self.totmig)[tau+1:t].prod()
                self.dicTpopTau[(t, pop, tau)] = mig[t, pop]*prod

        self.mat = np.zeros((len(self.states), len(self.states)))
        for nump, (tp, popp) in enumerate(self.states):
            for num, (t, pop) in enumerate(self.states):
                tot = 0
                for tau in range(1, min(t, tp)):
                    tot += self.dicTpopTau[(tp, popp, tau)] * psivec[(tau-1)]
                self.mat[num, nump] = tot

        self.uniformizemat()
        self.ndists = []
        for i in range(self.npops):
            self.ndists.append(self.popNdist(i))
        self.switchdensity()

    def gen_variance(self, popnum):
        """ 1. Calculate the expected genealogy variance in the model.
            2. Calculate the e(d) (equation 3 in MOLA (Models of Local ancestry) paper)
            3. Generations go from 0 to self.ngen-1.
            """
        legterm = [self.proportions[self.ngen - d, popnum]**2 * \
                np.prod(1 - self.totmig[:(self.ngen - d)])
                for d in range(1, self.ngen)]
        trunkterm = [np.sum([self.mig[u, popnum] * \
                np.prod(1 - self.totmig[:u])
                for u in range(self.ngen-d)])
                for d in range(1, self.ngen)]

        # Now calculate the actual variance.
        return np.sum([2**(d-self.ngen) * (legterm[d-1]+trunkterm[d-1])\
            for d in range(1, self.ngen)]) +\
                    self.proportions[0, popnum] *\
                    (1/2.**(self.ngen-1)-self.proportions[0, popnum])

    def uniformizemat(self):
        """ Uniformize the transition matrix so that each state has the same
            total transition rate. """
        self.unifmat = self.mat.copy()
        lmat = len(self.mat)
        # identify the highest non-self total transition rate
        outgoing = (self.mat - np.diag(self.mat.diagonal())).sum(axis=1)

        self.maxrate = outgoing.max() #max outgoing rate
        # reset the self-transition rate. We forget about the prior self-transition rates,
        # as they do not matter for the trajectories.
        for i in range(lmat):
            self.unifmat[i, i] = self.maxrate-outgoing[i]
        self.unifmat /= self.maxrate

    def popNdist(self, pop):
        """ Calculate the distribution of number of steps before exiting
            population. """
        if len(self.stateINpop[pop]) == 0:
            return []
        # get the equilibrium distribution in tracts OUTSIDE pop.
        tempequil = self.equil.copy()
        tempequil[self.stateINpop[pop]] = 0
        # Apply one evolution step
        new = np.dot(tempequil, self.unifmat)
        # select states in relevant population

        newrest = new[self.stateINpop[pop]]
        newrest /= newrest.sum() #normalize probabilities
        # reduce the matrix to apply only to states of current population
        shortmat = self.unifmat[np.meshgrid(self.stateINpop[pop],self.stateINpop[pop])].transpose()
        # calculate the amount that fall out of the state
        escapes = 1 - shortmat.sum(axis=1)
        # decide on the number of iterations

        # if the expected length of tracts becomes more than maxmorgans,
        # bail our and issue warning.
        nit = int(self.max_morgans* self.maxrate)

        nDistribution = []
        for i in range(nit): # nit is the max number of iterations.
            # will exit loop earlier if tracts are all complete.
            nDistribution.append(np.dot(escapes, newrest))
            newrest = np.dot(newrest, shortmat)
            if newrest.sum()<self.max_remaining_tracts:#we stop when there are at most
            # we request that the proportions of tracts that are still
            # incomplete be at most self.cutoff tracts left
                break
        if newrest.sum()>self.max_remaining_tracts:
            print("Warning: After %d time steps, %f of tracts are incomplete" %
                    (nit,newrest.sum()))
            print("This can happen when one population has really long",
                    "tracts.")

        nDistribution.append(newrest.sum())
        return nDistribution

    def Erlang(self, i, x, T):
        if i > 10:
            lg = i*np.log(T)+(i-1)*np.log(x)-T*x-gammaln(i)
            return np.exp(lg)
        return T**i*x**(i - 1)*np.exp(- T*x)/factorial(i - 1)

    def inners(self, L, x, pop):
        """ Calculate the length distribution of tract lengths not hitting a
            chromosome edge. """
        if(x > L):
            return 0
        else:
            return np.sum(
                    self.ndists[pop][i] *
                        (L-x) *
                        self.Erlang(i+1, x, self.maxrate)
                        for i in range(len(self.ndists[pop])))

    def outers(self, L, x, pop):
        """ Calculate the length distribution of tract lengths hitting a single
            chromosome edge. """
        if(x > L):
            return 0
        else:
            nd = self.ndists[pop]
            mx = self.maxrate * x
            return 2 * np.sum(
                    nd[i] * (1 - gammainc(i+1, mx))
                    for i in range(
                        len(nd))
            ) + 2 * (1-np.sum(nd))

    def full(self, L, pop):
        """ The expected fraction of full-chromosome tracts, p. 63 May 24,
            2011. """
        return np.sum(
                self.ndists[pop][i] * (((i+1) / float(self.maxrate) - L) +
                    L * gammainc(i + 1, self.maxrate * L) -
                    float(i+1) / self.maxrate * gammainc(i+2, self.maxrate*L))
                for i in range(len(self.ndists[pop]))
        ) + (1 - np.sum(self.ndists[pop])) * \
                (len(self.ndists[pop])/self.maxrate - L)

    def Z(self, L, pop):
        """the normalizing factor, to ensure that the tract density is 1."""
        return L + np.sum(
                self.ndists[pop][i]*(i+1)/self.maxrate
                for i in range(len(self.ndists[pop]))
        ) + (1 - np.sum([self.ndists[pop]])) * \
                len(self.ndists[pop])/self.maxrate

    def switchdensity(self):
        """ Calculate the density of ancestry switchpoints per morgan in our
            model. """
        self.switchdensities = np.zeros((self.npops, self.npops))
        # could optimize by precomputing survivals earlier
        self.survivals = [(1 - self.totmig[:i]).prod()
                for i in range(self.ngen)]
        for pop1 in range(self.npops):
            for pop2 in range(pop1):
                self.switchdensities[pop1, pop2] = \
                        np.sum(
                                [2 * self.proportions[i+1, pop1] * \
                                        self.proportions[i+1, pop2]* \
                                        self.survivals[i+1]
                                    for i in range(1, self.ngen-1)])
                self.switchdensities[pop2, pop1] = \
                        self.switchdensities[pop1, pop2]

        self.totSwitchDens = self.switchdensities.sum(axis=1)

    def expectperbin(self, Ls, pop, bins):
        """ The expected number of tracts per bin for a diploid individual with
            distribution of chromosome lengths given by Ls. The bin should be a
            list with n+1 breakpoints for n bins. We will always add an extra
            value for the full chromosomes as an extra bin at the end. The last
            bin should not go beyond the end of the longest chromosome. For
            now, perform poor man's integral by using the bin midpoint value
            times width. """
        self.totalPerInd = \
                [L*self.totSwitchDens[pop]+2.*self.proportions[0, pop]
                        for L in Ls]
        self.totalfull = \
                np.sum(
                        [(L*self.totSwitchDens[pop]+2. * \
                                self.proportions[0, pop]) * \
                                self.full(L, pop)/self.Z(L, pop)
                            for L in Ls])
        lsval = []
        for binNum in range(len(bins) - 1):
            mid = (bins[binNum] + bins[binNum+1]) / 2.
            val = np.sum(
                    [(L*self.totSwitchDens[pop] + \
                            2. * self.proportions[0, pop]) * \
                            (self.inners(L, mid, pop) + \
                                self.outers(L, mid, pop)) / self.Z(L, pop)
                        for L in Ls]) \
                    * (bins[binNum+1] - bins[binNum])
            lsval.append(max(val, 1e-17))

        lsval.append(max(self.totalfull, 1e-17))
        return lsval

def continuous_mig_hybridzone(rate, tstart):
    """ A simple model in which a population initially composed of a single
    ancestry receives migrants at a constant rate each generation to the
    present.
    """
    gen = int(np.ceil(tstart))+1
    mig = np.zeros((gen+1, 2))
    mig[-1,:] = np.array([.5, .5])
    mig[2:-1,:] = np.array([rate, rate])

    return mig
