import numpy             as np
import pandas as pd

from functools import reduce
#from collections.abc import Iterable

#--- with lists

def list_transpose(ll):
    """
    
    transpose a list m-length with n-length each item
    returns a list of n-items each of m-lenght

    Parameters
    ----------
    ll : list of lists

    Returns
    -------
    lt : list of lists

    """
    m = len(ll[0])
    lt = [[x[i] for x in ll] for i in range(m)]
    return lt


def list_to_df(ll, names):
    """
    
    Converts a list of list into a pandas DataFrame

    Parameters
    ----------
    ll    : list(list), contents of the DF colums
    names : list(str), name of the DF columnes

    Returns
    -------
    df    : DataFrame
    """

    assert len(ll) == len(names), 'required same number of lists and names'    

    df     = {}
    for i, name in enumerate(names): 
        df[name] = ll[i]
        
    df = pd.DataFrame(df)
    return df



#--- general utilies

def remove_nan(vals : np.array) -> np.array:
    """ returns the np.array without nan
    """
    return vals[~np.isnan(vals)]


def in_range(vals : np.array,
             range : tuple = None,
             upper_limit_in = False) -> np.array(bool):
    """ returns a np.array(bool) with the elements of val that are in range
    inputs:
        vals : np.array
        range: None, (x0, x1) or x0: None, all values; (c1, c1): x >= x0, x <xf;
               x0: x == x0                 
        upper_limit_int, if True x <= xf
    returns
        sel : np.array(bool) where True/False indicates if the elements of vals are
        in range
    """

    if (range is None): 
        return vals >= np.min(vals)

    if (isinstance(range, list) or (isinstance(range, tuple))):
        sel1 = (vals >= range[0])
        sel2 = (vals <= range[1]) if upper_limit_in else (vals < range[1])
        return sel1 & sel2
    
    return vals == range


def centers(xs : np.array) -> np.array:
    """ returns the center between the participn
    inputs:
        xs: np.array
    returns:
        np.array with the centers of xs (dimension len(xs)-1)
    """
    return 0.5* ( xs[1: ] + xs[: -1])




def stats(vals : np.array, range : tuple = None):
    vals = np.array(vals)
    vals = remove_nan(vals)
    sel  = in_range(vals, range)
    vv = vals[sel]
    mean, std, evts, oevts = np.mean(vv), np.std(vv), len(vv), len(vals) - len(vv)
    return evts, mean, std, oevts


def str_stats(vals, range = None, formate = '6.2f'):
    evts, mean, std, ovts = stats(vals, range)
    s  = 'entries '+str(evts)+'\n'
    s += (('mean {0:'+formate+'}').format(mean))+'\n'
    s += (('std  {0:'+formate+'}').format(std))
    return s


#-----

def efficiency(sel, n = None):
    """ compute the efficiency and uncertantie of a selection
    inputs:
        sel: np.array(bool), bool array with True/False
        n  : int, denominator, if n is None, use len(sel)
    """
    n    = n if n is not None else len(sel)
    eff  = np.sum(sel)/n
    ueff = np.sqrt(eff * (1- eff) / n)
    return eff, ueff

#-------

def selection(df, varname, varrange, oper = np.logical_and):
    """ apply the selection on a DataFrame requirend that the variable(s) are in the range
    
    inputs:
        df       : dataFrame
        varname  : str or list(str), name of list of names of the variable s 
        varrange : tuple(float, float), range of the selection, all (-np.inf, np.inf)
        oper     : bool operation, default and
    returns:
        sel      : np.array(bool) same same of DF with True/False
                   if the item fulfull variable value (varname) inside the range (varrange)
    """
    
    _isiter = lambda x: isinstance(x, list) or isinstance(x, tuple)
    if _isiter(varname):
        assert len(varname) == len(varrange), \
            'required same length of variables and ranges'
        sels = [selection(df, ivar, ivarran) for ivar, ivarran \
                      in zip(varname, varrange)]
        return reduce(oper, sels)
    
    return in_range(df[varname], varrange)
    

def selection_efficiency(df, varname, varrange):
    
    """
    Computes the efficiency of a selection in a DataFrame

    Parameters
    ----------
    df       : DataFrame,
    varname  : str or list(str), name of list of names of the variables of the selection
    varrange : (float, float) or list(float, float), range or ranges of the variables of the selection

    Returns
    -------
    eff      : (float, float), efficiency and uncertanty on the efficiency

    """
    return efficiency(selection(df, varname, varrange))


def selection_sample(df, varname, varrange):
    """
    
    return a sample of a DataFrame with the events that pass the selection
    

    Parameters
    ----------
    df       : DataFrame 
    varname  : str or list(str), name of list of names of the variable o f the selection
        DESCRIPTION.
    varrange : (float, float) or list(float, float), range or list of ranges of the variable sof the selection
        DESCRIPTION.

    Returns
    -------
    df      ; DataFrame

    """
    
    return df[selection(df, varname, varrange)]

