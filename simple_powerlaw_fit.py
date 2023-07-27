import numpy as np

def simple_powerlaw_fit(data, xmin=None, xmax=None, Ngrid=1024,alpha_min=-4,alpha_max=4,quantiles=[.16,.5,.84]):
    """Computes the specified quantiles of the posterior distribution of the power law slope alpha, given a dataset and desired quantiles
    
    Example usage to get intervals on the power-law slope of the IMF over the interval [1,10]:
    `
    masses = np.loadtxt("my_IMF_data.dat")
    quantiles = simple_powerlaw_fit(masses, mmin=1, mmax=10)
    `
    """
    if xmin is None: xmin = data.min()
    if xmax is None: xmax = data.max()
    data = data[(data>xmin)*(data<xmax)] # prune data to specified limits

    alpha_grid = np.linspace(alpha_min,alpha_max,Ngrid) # initialize the 1D grid in parameter space
    lnprob = np.zeros_like(alpha_grid) # grid that stores the values of the posterior distribution

    normgrid =  (1 + alpha_grid) / (xmax**(1+alpha_grid)-xmin**(1+alpha_grid))  # grid of normalization of the distribution x^alpha over the limits [mmin,mmax]

    for d in data: # sum the posterior log-likelihood distribution
        lnprob += np.log(d**alpha_grid * normgrid)

    # convert log likelihood to likelihood
    lnprob -= lnprob.max()
    prob = np.exp(lnprob) # watch for overflow errors here
    prob /= np.trapz(prob,alpha_grid) # normalize

    q = np.interp(quantiles,np.cumsum(prob)/prob.sum(), alpha_grid)
    return q # returns quantiles of the posterior distribution
