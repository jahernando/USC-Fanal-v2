import numpy as np
import scipy.optimize as optimize

import matplotlib.pyplot as plt


"""
    Module for histogram fitting

    predefined functions: gaus, lin, exp, gausline, gausexp

"""


functions = ['gaus', 'line', 'exp', 'gausline', 'gausexp']

current_module = __import__(__name__)


def _predefined_function(fun, guess, x):

    fnames = None
    if (type(fun) == str):
        assert fun in functions
        if (guess is None):
            guess = eval('g'+fun)(x)
        fnames = getattr(current_module.hfit, 'n' + fun)
        fun    = getattr(current_module.hfit, 'f' + fun)

    return fun, guess, fnames


def hfit(x, bins, fun, guess = None, range = None):
    """ fit a histogram to a function with guess parameters
    inputs:
    x    : np.array, values to build the histogram
    bins : int (100), tuple, bins of the histogram
    fun  : callable(x, *parameters) or string, function to fit
           str  = ['gaus', 'line', 'exp' ], for gaussian, line fit
    guess: tuple (None), values of the guess/initial parameters for the fit
           if fun is a predefined function, no need to add initial gess parameters
    range: tuple (None), range of the values to histogram

    TODO : add a mask argument to mask the parameters to fix in the fit!
    TODO : consider the errors in the histogram!
    TODO : check that the returned errors are ok!
    """


    fun, guess, _ = _predefined_function(fun, guess, x)

    range = range if range is not None else (np.min(x), np.max(x))
    yc, xe = np.histogram(x, bins, range)
    xc = 0.5 * (xe[1:] + xe[:-1])
    #yerr = np.sqrt(np.maximum(yc, 2.4))
    #fpar, fcov = optimize.curve_fit(fun, xc, yc, guess, sigma = yerr, absolute_sigma = True)
    fpar, fcov = optimize.curve_fit(fun, xc, yc, guess)

    return fpar, np.sqrt(np.diag(fcov))


def hresiduals(x, bins, fun = 'gaus', pars = None, **kargs):

    fun, guess, _ = _predefined_function(fun, pars, x)

    ys, edges = np.histogram(x, bins, **kargs)
    xcs       = 0.5 * (edges[:-1] + edges[1:])
    yerr      = np.maximum(np.sqrt(ys), 1.)
    res       = (fun(xcs, *guess) - ys)/yerr
    chi2      = np.sum(res * res)
    ndf       = len(ys[ys > 0]) - len(pars)

    return res, edges, chi2, ndf



def hfitres(x, bins, fun, guess = None, **kargs):

    pars, epars     = hfit      (x, bins, fun, guess = guess, **kargs)
    _, _, chi2, ndf = hresiduals(x, bins, fun, pars, **kargs )

    return pars, epars, chi2, ndf


def str_parameters(pars, covpars, parnames = None, formate = '6.2f'):
    s = ''
    for i, par in enumerate(pars):
        namepar = r'$a_'+str(i)+'$' if parnames is None else parnames[i]
        covpar   = covpars[i]
        s += namepar + ' = '
        s += (('{0:'+formate+'}').format(par))   + r'$\pm$'
        s += (('{0:'+formate+'}').format(covpar))+ '\n'
    return s


def fgaus(x, a, b, c):
    """ return a gausian function
    """
    if c <= 0.:
        return np.inf
    return a * np.exp(- (x-b)**2 / (2* c**2) )


def ggaus(x):
    """ return guess parameters for a guassian function
    """
    return (len(x), np.mean(x), np.std(x))


ngaus = [r'$N_\mu$', r'$\mu$', r'$\sigma$']


def fline(x, a, b):
    """ return a line a* x + b
    """
    return a * x + b


def gline(x):
    """ return guess parameters for a line function
    """
    ys, xe = np.histogram(x, 2)
    xc = 0.5 * (xe[1:] + xe[:-1])
    a = (ys[1] - ys[0])/ (xc[1] - xc[0])
    b = ys[0] - a * xc[0]
    return a, b


nline = ['a', 'b']


def fexp(x, a, b):
    """ an exponential function a * exp(-b * x)
    """
    return a * np.exp( - b * x)


def gexp(x):
    """ guess parameters for an exponential
    """
    ys, xs = np.histogram(x, 2)
    xcs = 0.5 * (xs[1:] + xs[:-1])
    dx = xcs[1] - xcs[0]
    b = - (np.log(ys[1]) - np.log(ys[0]))/dx
    a = ys[0] * np.exp(b * xcs [0])
    return (a, b)


nexp = [r'$N_\tau$', r'$\tau$']

def fgausline(x, na, mu, sig, a, b):
    return fgaus(x, na, mu, sig) + fline(x, a, b)


def ggausline(x):
    return list(ggaus(x)) + list(gline(x))


ngausline = ngaus + nline


def fgausexp(x, na, mu, sig, nb, tau):
    return fgaus(x, na, mu, sig) + fexp(x, nb, tau)


def ggausexp(x):
    return list(ggaus(x)) + list(gexp(x))


ngausexp = ngaus + nexp

#------------
