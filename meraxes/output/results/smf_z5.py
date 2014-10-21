#!/usr/bin/env python

"""Plot the z=5 SMF.

Usage: smf_z5.py <fname> [Hubble_h]

Arguments:
    Hubble_h   Hubble constant [default: 0.702]
"""

import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt
from ssimpl import plotutils, meraxes, munge
import pandas as pd
from astropy import log

import os
import sys
__script_dir__ = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(__script_dir__)
from smf_common import smf_uncert

__author__ = "Simon Mutch"
__date__   = "2014-05-01"


def plot(gals, simprops, ax, h):

    print "Plotting z=5 SMF..."
    markers = ['o', 's', 'D', 'p', 'v', '^'][::-1]

    # N.B. NOTE SFR CUT TO MATCH KATSIANIS 2014 (via SMIT 2012) DATA (log(SFR)
    # > -0.44)
    sfr_limit = -0.44  # lower log10(sfr [Msol]) limit for galaxies
    sfr_limit *= (0.7**2)/(h**2)

    # generate the model smf
    # sm = stellar mass
    sm = np.log10(gals.StellarMass[(gals.StellarMass > 0)
                                   & (np.log10(gals.Sfr) > sfr_limit)] * 1.0e10)

    # n_dropped = gals.shape[0] - sm.shape[0]
    # if n_dropped > 0:
    #     log.warn("Dropped %d galaxies (%.1f%% of total) with stellar mass <= 0" %
    #             (n_dropped, float(n_dropped)/gals.shape[0]*100))

    smf, edges = munge.mass_function(sm, simprops["Volume"], "knuth",
                                     return_edges=True)

    # plot the model
    l, = ax.plot(smf[:,0], np.log10(smf[:,1]), ls="--", lw=4)

    # plot the total SMF (i.e. without SFR cut)
    sm = np.log10(gals.StellarMass[gals.StellarMass > 0] * 1.0e10)
    smf = munge.mass_function(sm, simprops["Volume"], edges,
                              poisson_uncert=True)
    ax.plot(smf[:,0], np.log10(smf[:,1]), ls='-', lw=4, color=l.get_color(),
            label="Meraxes")
    # yerr = [np.log10(smf[:,1]-smf[:,2]), np.log10(smf[:,1]+smf[:,2])]
    yerr = np.log10(smf_uncert(sm, simprops["Volume"], edges, 5.0))
    yerr[0][np.isinf(yerr[0])] = -7
    p = ax.fill_between(smf[:,0], yerr[0], yerr[1], color=l.get_color(),
                    alpha=0.3)
    r = plt.Rectangle([0,10], 0.5, 0.1, color=p.get_facecolor()[0],
                      alpha=p.get_alpha(), label="Meraxes - convolved+Poisson")
    ax.add_patch(r)

    # read the observed smf
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/smf/Gonzalez11_z5_smf.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 3,
        names = ["sm", "log_phi", "m_err", "p_err"])

    # convert obs to same hubble value and IMF
    obs.sm += 2.0*np.log10(0.7/h)
    obs.sm += 0.25  # IMF correction Chabrier -> Salpeter
    for col in ["log_phi", "m_err", "p_err"]:
        obs[col] -= 3.0*np.log10(0.7/h)

    # plot the observations
    ax.errorbar(obs.sm, obs.log_phi, yerr=[obs.m_err, obs.p_err],
                label="Gonzalez et al. 2011", ls="None",
                lw=2, capsize=2.5, marker=markers.pop(), mec='None')

    # do it all again for the next set of observations
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/smf/Katsianis_z5_smf.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 3,
        names = ["sm", "log_phi", "err"])

    # convert obs to same hubble value and IMF
    obs.sm += 2.0*np.log10(0.702/h)
    for col in ["log_phi", "err"]:
        obs[col] -= 3.0*np.log10(0.702/h)

    # plot the observations
    ax.errorbar(obs.sm, obs.log_phi, yerr=obs.err,
                label="Katsianis et al. 2014", ls="None",
                lw=2, capsize=2.5, marker=markers.pop(), mec='None')

    # and again...
    obs = pd.read_table(os.path.join(__script_dir__,
        "../../utils/obs_datasets/smf/Duncan14_MF_z5.txt"),
        delim_whitespace=True,
        header = None,
        skiprows = 8,
        names = ["sm", "phi", "merr", "perr"])
    obs.merr[obs.merr >= obs.phi] = obs.phi - 1e-10

    # convert obs to same hubble value and IMF
    obs.sm += 2.0*np.log10(0.702/h)
    obs.sm += 0.25  # IMF correction Chabrier -> Salpeter
    for col in ["phi", "merr", "perr"]:
        obs[col] /= (0.7**3/h**3)

    # plot the observations
    ax.errorbar(obs.sm, np.log10(obs.phi),
                yerr=[np.log10(obs.phi / (obs.phi - obs.merr)),
                      np.log10(1.0 + (obs.perr / obs.phi))],
                label="Duncan et al. 2014", ls="None",
                lw=2, capsize=2.5, marker=markers.pop(), mec='None')

    # add some text
    ax.text(0.95,0.95, "z=5\nh={:.2f}\nSalpeter IMF\n".format(h)+
            r"log$_{10}$(SFR)"+" > {:.2f}".format(sfr_limit)+r"M$_{\odot}$/yr",
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes)

    ax.set_xlim([7.5,11])
    ax.set_ylim([-6,-1])

    ax.set_xlabel(r"$\log_{10}(M_* / {\rm M_{\odot}})$")
    ax.set_ylabel(r"$\log_{10}(\phi / ({\rm dex^{-1}\,Mpc^{-3}}))$")


if __name__ == '__main__':

    args = docopt(__doc__)
    fname = args['<fname>']
    if args['Hubble_h'] is False:
        h = 0.702
    else:
        h = float(args['Hubble_h'])

    snap, redshift = meraxes.io.check_for_redshift(fname, 5.0, tol=0.1)

    props = ("StellarMass", "Mvir", "Type", "Sfr")
    gals, simprops = meraxes.io.read_gals(fname, snapshot=snap, props=props,
            sim_props=True, h=h)
    gals = gals.view(np.recarray)

    fig, ax = plt.subplots(1,1)
    plot(gals, simprops, ax, h)
    ax.yaxis.set_tick_params(which='both', color='w')
    ax.legend(loc="lower left")
    fig.tight_layout()
    output_fname = os.path.join(os.path.dirname(fname), "plots/smf-z5.png")
    plt.savefig(output_fname)
