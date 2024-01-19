# Script for algorithms to use
import numpy as np
from hyperopt import tpe, hp, fmin, rand, STATUS_OK, Trials, atpe
from multiprocessing import Pool
import emcee
import os
import loss_model_global_mpi_sp_eo
import pickle

def a_tpe(prior, f, args, outfile):
    """
    Performs the optimisation of the search space using the TPE algorithm
    :param prior: dictionary of prior distribution of the parameters
    :param f: loss function (takes param dictionary as parameter)
    :param args: [number of iterations desired]
    :return: test: optimized parameter vector, best: dictionary
    """

    space = {}          # empty dictionary to append

    for i in range(0, len(prior['parameter_name'])):   # conditional check for which distribution to use

        distribution = prior['distribution_type'][i]
        if distribution == 'uniform':
            space[prior['parameter_name'][i]] = hp.uniform(prior['parameter_name'][i], prior['lower_bound'][i], prior['higher_bound'][i])
        if distribution == 'random':
            value = hp.randint(prior[0][i], prior[2][i])

    it = args           # number of search iterations
    if it <= 0 or type(it) != int:
        raise Exception('Invalid iteration value, must be positive integer')
    if outfile is None:
        trials = Trials()
    else:
        trials = pickle.load(open(outfile, "rb"))
    best = fmin(
        fn=loss_model_global_mpi_sp_eo.loss,           # Objective Function to optimize
        space=space,    # hyperparameter's Search space
        algo=tpe.suggest,  # Optimization algorithm (representative TPE)
        max_evals=it,   # Number of optimization attempts
        trials=trials   # Handles of diagnosis quantities
    )

    return trials


def a_atpe(prior, f, args):
    """
    Performs the optimisation of the search space using the  adaptive TPE algorithm
    :param prior: array containing prior knowledge on parameters
    :param f: loss function (takes param dictionary as parameter)
    :param args: [number of iterations desired]
    :return: test: optimized parameter vector, best: dictionary
    """

    space = {}              # empty dictionary to append

    for i in range(0, len(prior[0])):   # conditional check for which distribution to use
        if prior[3][i] == 'uniform':
            key = hp.uniform(prior[0][i], prior[1][i], prior[2][i])
        if prior[3][i] == 'random':
            key = hp.randint(prior[0][i], prior[2][i])

        space.update({prior[0][i]: key})

    it = args               # number of search iterations
    trials = Trials()
    best = fmin(
        fn=f,               # Objective Function to optimize
        space=space,        # hyperparameter's Search space
        algo=atpe.suggest,  # Optimization algorithm (representative TPE)
        max_evals=it,       # Number of optimization attempts
        trials=trials       # Handles of diagnosis quantities
    )

    test = np.array(list(best.values()))
    return test, best


def a_rand(prior, f, args):
    """
    Performs the optimisation of the search space using the random search algorithm
    :param prior: dictionary of prior distribution of the parameters
    :param f: loss function (takes param dictionary as parameter)
    :param args: [number of iterations desired]
    :return: test: optimized parameter vector, best: dictionary
    """

    space = {}              # empty dictionary to append

    for i in range(0, len(prior[0])):   # conditional check for which distribution to use
        if prior[3][i] == 'uniform':
            key = hp.uniform(prior[0][i], prior[1][i], prior[2][i])
        if prior[3][i] == 'random':
            key = hp.randint(prior[0][i], prior[2][i])

        space.update({prior[0][i]: key})

    it = args               # number of search iterations
    trials = Trials()
    best = fmin(
        fn=f,               # Objective Function to optimize
        space=space,        # hyperparameter's Search space
        algo=rand.suggest,  # Optimization algorithm (representative TPE)
        max_evals=it,       # Number of optimization attempts
        trials=trials       # Handles of diagnosis quantities
    )

    test = np.array(list(best.values()))
    return test, best


def loss(param, prior, f):
    """
    Converts input from MCMC algo into dict that is compatible with the input needed for the loss function
    :param param: parameters to pass to pass to loss func
    :param prior: parameters name for dict input
    :param f: loss function
    :return: 1/loss since MCMC maximises likelihood instead of minimizing
    """
    dic = {}    # Create empty dictionnary
    for i in range(0,len(param)):    # fill dict with param name and value
        key = param[i]
        dic.update({prior.parameter_name.values[i]: key})
    chi = loss_model_global.loss(dic)   # evaluate loss for param vector
    #chi = chi.get('loss')
    return 1/chi


def priore(param, prior):
    """
    Computes prior distribution of each parameter
    :param param: parameter vector
    :param prior: array containing info on prior distribution of params
    :return: likelihood contribution of prior distribution
    """
    p = np.zeros(len(param))
    for i in range(0,len(param)):   # checking if each parameter is in prior distribution range
        if prior.distribution_type.values[i] == 'uniform':
            if prior.lower_bound.values[i] < param[i] < prior.higher_bound.values[i]:
                p = 0
            else:
                return - np.inf
        return p


def combine(param,prior,f):
    """
    Combines prior likelihood and loss likelihood
    :param param: parameter vector
    :param prior: array of info on prior
    :param f: loss function
    :return: value of the likelihood function at parameter vector coordinates
    """
    p = priore(param,prior)
    if not np.isfinite(p):
        return -np.inf

    return p + loss(param, prior, f)


def a_mcmc(prior, f, args):  # need to add args
    """
    Performs optimisation of search space using the Monte Carlo Markov Chain algorithme
    :param prior: prior knowledge on the parameter
    :param f: loss function
    :param args: number of walkers, number of iterations
    :return: vector of optimized params and other stats
    """
    nparam = prior.shape[0]         # number of parameter to optimize
    lb = prior.lower_bound.values   # lower bound of parameters
    hb = prior.higher_bound.values  # upper bound of parameters
    distribution = prior.distribution_type.values
    nwalkers = args[0]       # number of walkers
    it = args[1]             # number of iterations

    # initializing walkers given prior parameter distributions
    p0 = np.zeros([nwalkers, nparam])
    for i in range(0, nparam):
        if distribution[i] == 'uniform':
            p0[:, i] = np.random.uniform(lb[i], hb[i], size=nwalkers)  # array of random startpoint of walkers
        if distribution[i] == 'normal':
            p0[:, i] = np.random.normal(lb[i], hb[i], size=nwalkers)  # array of random startpoint of walkers

    # setting algorithm run
    with Pool(processes=os.environ.get('SLURM_CPUS_PER_TASK', default=4)) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, nparam, combine,pool=pool, args=[prior, f])  # Specify algo params
        sampler.run_mcmc(p0, it, progress=True,skip_initial_state_check=True)  # Runs algorithm

    # Setting outputs
    samples = sampler.get_chain()  # value of each param of each walker at each iteration

    flat_samples = sampler.get_chain(flat=True)
    param1 = np.mean(flat_samples[:, 0])
    param2 = np.mean(flat_samples[:, 1])
    param3 = np.mean(flat_samples[:, 2])
    param4 = np.mean(flat_samples[:, 3])
    param5 = np.mean(flat_samples[:, 4])

    param = np.array([param1, param2, param3, param4, param5])

    stats = {'Samples': samples}

    return param, stats

# END OF TPE ALGORITHM FUNCTION SCRIPT
