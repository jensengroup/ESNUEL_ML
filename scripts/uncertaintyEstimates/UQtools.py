import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss

from sklearn.linear_model import LinearRegression

from scipy.stats import pearsonr
from scipy.stats import norm
from scipy.stats import bootstrap
from scipy.integrate import quad

from bisect import bisect_left


def order_sig_and_errors(sigmas, errors):
    ordered_df = pd.DataFrame()
    ordered_df["uq"] = sigmas
    ordered_df["errors"] = errors
    ordered_df["abs_z"] = abs(ordered_df.errors)/ordered_df.uq
    ordered_df = ordered_df.sort_values(by="uq")
    return ordered_df


def spearman_rank_corr(v1, v2):
    v1_ranked = ss.rankdata(v1)
    v2_ranked = ss.rankdata(v2)
    return pearsonr(v1_ranked, v2_ranked)


def rmse(x, axis=None):
    return np.sqrt((x**2).mean())


def get_bootstrap_intervals(errors_ordered, Nbins=10):
    """
    calculate the confidence intervals at a given p-level. 
    """
    ci_low = []
    ci_high = []
    N_total = len(errors_ordered)
    N_entries = math.ceil(N_total/Nbins)

    for i in range(0, N_total, N_entries):
        data = errors_ordered[i:i+N_entries]
        res = bootstrap((data,), rmse, vectorized=False)
        ci_low.append(res.confidence_interval[0])
        ci_high.append(res.confidence_interval[1])
    return ci_low, ci_high


def expected_rho(uncertainties):
    """
    for each uncertainty we draw a random Gaussian error to simulate the expected errors
    the spearman rank coeff. is then calculated between uncertainties and errors. 
    """
    sim_errors = []
    for sigma in uncertainties:
        error = np.abs(np.random.normal(0, sigma))
        sim_errors.append(error)
    
    rho, _ = spearman_rank_corr(uncertainties, sim_errors)
    return rho, sim_errors


def NLL(uncertainties, errors):
    NLL = 0
    for uncertainty, error in zip(uncertainties, errors):
        temp = math.log(2*np.pi*uncertainty**2)+(error)**2/uncertainty**2
        NLL += temp
    
    NLL = NLL/(2*len(uncertainties))
    return NLL


def calibration_curve(errors_sigma):
    N_errors = len(errors_sigma)
    gaus_pred = []
    errors_observed = []
    for i in np.arange(-10, 0+0.01, 0.01):
        gaus_int = 2*norm(loc=0, scale=1).cdf(i)
        gaus_pred.append(gaus_int)
        observed_errors = (errors_sigma > abs(i)).sum()
        errors_frac = observed_errors/N_errors
        errors_observed.append(errors_frac)
    
    return gaus_pred, errors_observed


def plot_calibration_curve(gaus_pred, errors_observed, mis_cal):
    fig, ax = plt.subplots(figsize=(5,5))
    ax.fill_between(gaus_pred, gaus_pred, errors_observed, color="purple", alpha=0.4, label="Miscalibration area = {:0.3f}".format(mis_cal))
    ax.plot(gaus_pred, errors_observed, color="purple", alpha=1)
    ax.plot(np.arange(0,1,0.01),np.arange(0,1,0.01), linestyle='-', color='dimgrey')
    ax.set_xlabel("Expected fraction of errors", fontsize=16)
    ax.set_ylabel("Observed fraction of errors", fontsize=16)
    ax.legend(fontsize=14, loc='upper left', frameon=False)
    plt.xticks(np.arange(0,1.1,0.2))
    plt.yticks(np.arange(0,1.1,0.2))
    plt.xlim(0,1)
    plt.ylim(0,1)
    return fig


def plot_Z_scores(errors, uncertainties):
    Z_scores = errors/uncertainties
    N_bins = 29
    xmin, xmax = -7,7
    y, bin_edges = np.histogram(Z_scores, bins=N_bins, range=(xmin, xmax))
    bin_width = bin_edges[1] - bin_edges[0]
    x = 0.5*(bin_edges[1:] + bin_edges[:-1])
    sy = np.sqrt(y)
    target_values = np.array([len(errors)*bin_width*norm.pdf(x_value) for x_value in x])
    fig, ax = plt.subplots(figsize=(7,5))
    ax.hist(Z_scores, bins=N_bins, range=(xmin, xmax), color='purple', alpha=0.3)
    ax.errorbar(x, y, sy, fmt='.', color='k')
    ax.plot(np.arange(-7, 7, 0.1), len(errors)*bin_width*norm.pdf(np.arange(-7, 7, 0.1), 0, 1), color='k')
    
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    ax.set_xlabel("Error (Z)", fontsize=16)
    ax.set_ylabel("Count", fontsize=16)
    ax.set_xlim([-7, 7])
    return fig, ax

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0], myList[1], 0, 1
    if pos == len(myList):
        return myList[-2], myList[-1], -2, -1
    before = myList[pos - 1]
    after = myList[pos]
    if myNumber<before or myNumber>after:
        print("problem")
    else:
        return before, after, pos-1, pos


def f_linear_segment(x, point_list=None, x_list=None):
    
    x1, x2, x1_idx, x2_idx = take_closest(x_list, x) 
    f = point_list[x1_idx]+(x-x1)/(x2-x1)*(point_list[x2_idx]-point_list[x1_idx])
    
    return f


def area_function(x, observed_list, predicted_list):
    h = abs((f_linear_segment(x, observed_list, predicted_list)-x))
    return h


def calibration_area(observed, predicted):
    area = 0
    x = min(predicted)
    while x < max(predicted):
        temp, _ = quad(area_function, x, x+0.001, args=(observed, predicted))
        area += temp
        x += 0.001
    return area


def chi_squared(x_values, x_sigmas, target_values):
    mask = x_values > 0
    chi_value = ((x_values[mask]-target_values[mask])/x_sigmas[mask])**2
    chi_value = np.sum(chi_value)
    

    N_free_cs = len(x_values[mask])
    print(N_free_cs)
    chi_prob =  ss.chi2.sf(chi_value, N_free_cs)
    return chi_value, chi_prob


def get_slope_metric(uq_ordered, errors_ordered, Nbins=10, include_bootstrap=True):
    """
    Calculates the error-based calibration metrices

    uq_ordered: list of uncertainties in increasing order
    error_ordered: list of observed errors corresponding to the uncertainties in uq_ordered
    NBins: integer deciding how many bins to use for the error-based calibration metric
    include_bootstrap: boolean deciding wiether to include 95% confidence intervals on RMSE values from bootstrapping
    """
    rmvs, rmses, ci_low, ci_high = get_rmvs_and_rmses(uq_ordered, errors_ordered, Nbins=Nbins, include_bootstrap=include_bootstrap)
    
    x = np.array(rmvs)
    y = np.array(rmses)

    # Obtain the coefficient of determination by calling the model with the score() function, then print the coefficient:
    slope, intercept, r_value, p_value, std_err = ss.linregress(x, y)
    r_sq = r_value**2
    print('R-squared:', r_sq)

    # Print the Intercept:
    print('Intercept:', intercept)

    # Print the Slope:
    print('Slope:', slope) 

    # Predict a Response and print it:
    # y_pred = model.predict(x)
    y_pred = slope * x + intercept
    
    fig, ax = plt.subplots(figsize=(5, 5))
    assymetric_errors = [np.array(rmses)-np.array(ci_low), np.array(ci_high)-np.array(rmses)]
    ax.errorbar(x, y, yerr = assymetric_errors, fmt="o", linewidth=2)

    if ax.get_xlim()[-1] > ax.get_ylim()[-1]:
        plt.xlim(0,ax.get_xlim()[-1])
        plt.ylim(0,ax.get_xlim()[-1])
        ax.plot(np.arange(0,ax.get_xlim()[-1],0.01),np.arange(0,ax.get_xlim()[-1],0.01), linestyle='-', color='dimgrey')
        plt.xticks(np.arange(0,ax.get_xlim()[-1]+1,10))
        plt.yticks(np.arange(0,ax.get_xlim()[-1]+1,10))
    else:
        plt.xlim(0,ax.get_ylim()[-1])
        plt.ylim(0,ax.get_ylim()[-1])
        ax.plot(np.arange(0,ax.get_ylim()[-1],0.01),np.arange(0,ax.get_ylim()[-1],0.01), linestyle='-', color='dimgrey')
        plt.xticks(np.arange(0,ax.get_ylim()[-1]+1,10))
        plt.yticks(np.arange(0,ax.get_ylim()[-1]+1,10))

    if intercept >= 0:
        ax.plot(rmvs, y_pred, linestyle="dashed", color='k', label="y = {:0.2f}".format(slope)+"x + {:0.2f}\n".format(intercept) + r'$R^2$ = '"{:0.2f}".format(r_sq))
    else:
        ax.plot(rmvs, y_pred, linestyle="dashed", color='k', label="y = {:0.2f}".format(slope)+"x - {:0.2f}\n".format(intercept) + r'$R^2$ = '"{:0.2f}".format(r_sq))

    ax.set_xlabel("RMV", fontsize=16)
    ax.set_ylabel("RMSE", fontsize=16)
    ax.legend(fontsize=16, frameon=False)
    return fig, slope, r_sq, intercept


def get_rmvs_and_rmses(uq_ordered, errors_ordered, Nbins=10, include_bootstrap=True):
    """
    uq orderes should be the list of uncertainties in increasing order and errors should be the corresponding errors
    Nbins determine how many bins the data should be divided into
    """

    N_total = len(uq_ordered)
    N_entries = math.ceil(N_total/Nbins)
    #print(N_entries)
    rmvs = [np.sqrt((uq_ordered[i:i+N_entries]**2).mean()) for i in range(0, N_total, N_entries)]
    #print(rmvs)
    rmses = [np.sqrt((errors_ordered[i:i+N_entries]**2).mean()) for i in range(0, N_total, N_entries)]
    if include_bootstrap:
        ci_low, ci_high = get_bootstrap_intervals(errors_ordered, Nbins=Nbins)
    else:
        ci_low, ci_high = None, None
    return rmvs, rmses, ci_low, ci_high





