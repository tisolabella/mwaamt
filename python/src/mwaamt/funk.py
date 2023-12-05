#####################################################################
# IMPORTS
#####################################################################

import numpy as np


#####################################################################
# CONSTANT
#####################################################################




#####################################################################
# FUNCTION DEFINITIONS
#####################################################################

def singlefit(x, scale, aae):
    """The function for the overall power-law fit"""
    return scale * x ** ( - aae )


def doublefit(x, scale_1, alpha_1, scale_2, alpha_2):
    """The function for the double power_law fit"""
    return scale_1 * x ** ( - alpha_1) + scale_2 * x ** ( - alpha_2 )


def line(x, a, b):
    return a * x + b


def set_resolution(*best_parameters, iteration=2):
    """Function to increase the parameter scan resolution"""
    return_list = []
    resolution = 0.1 / (iteration + 1)
    for best_param in best_parameters:
        param_list = [best_param - 2 * resolution, 
                best_param - resolution, best_param,
                best_param + resolution, 
                best_param + 2 * resolution]
        return_list.append(param_list)
    return return_list

def get_chisq(data, expected, sigma, ndf):
    """Calculates the chi squared"""
    chi = 0
    for d, e, s in zip(data, expected, sigma):
        chi += (((d - e) ** 2) / (s ** 2))
    # Returns the chisquared and the reduced chisq
    return chi, (chi / ndf)


def average(l):
    """Returns the average value of a list"""
    s = 0
    for x in l:
        s += x
    return s / len(l)

def weighted_average(l, el):
    """Returns the average of a list weigthed by el"""
    s = 0
    ws = 0
    for x, ex in zip(l, el):
        s += x / ex ** 2
        ws += 1 / ex ** 2
    return s / ws, 1 / np.sqrt(ws)

def stddev(l):
    """Returns the standard deviation of a list"""
    ll = [x ** 2 for x in l]
    avg_sq, sq_avg = average(ll), average(l) ** 2
    return np.sqrt(avg_sq - sq_avg)

