# Main script for Bayesian Optimisation Framework
import algo_scores_mo
import argparse
import numpy as np
import pandas as pd
import pickle

def main(prior, loss, algorithm):
    # loading prior distribution of parameters
    prior = pd.read_csv(prior, sep='\t', header=2)

    # Launching optimiaztion
    if algorithm == 'tpe':
        opt_out = algo_scores_mo.a_tpe(prior, loss, 3000)  # Optimisation with tpe algorithm
        #best = opt_out.best_trial               # List of best parameters
        #stats = opt_out.trials                  # Info on each step in optimisation
        #loss = opt_out.losses()                 # List of loss score at each step in optimization

        # Saving outputs to dictionnary
        dict_out = open("s3_mo_opt_out_3000.pkl", "wb")
        pickle.dump(opt_out, dict_out)
        dict_out.close()
        # Saving output files
        #with open('opt_out.txt', 'wb') as fp:
        #    pickle.dump(opt_out, fp)
        #with open('opt_out_stats.txt', 'wb') as fp:
        #    pickle.dump(stats, fp)
        #with open('opt_out_loss.txt', 'wb') as fp:
        #    pickle.dump(loss, fp)

    if algorithm == 'mcmc':
        vector, stats = algo.a_mcmc(prior, loss, args=[68,150]) #Optimisation with MCMC algorithm

    #print(vector)
    #np.savetxt('vect.txt', vector)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
      description='Arguments needed for optimisation')

    parser.add_argument('--prior',
                        type=str,
                        help='file name containing parameter priors',
                        required=True)
    parser.add_argument('--loss',
                        type=str,
                        help='python file of loss function',
                        required=True)
    parser.add_argument('--algo',
                        type=str,
                        help='python file to desired algorithm',
                        required=True)
    args = parser.parse_args()
    main(prior=args.prior, loss=args.loss, algorithm=args.algo)


