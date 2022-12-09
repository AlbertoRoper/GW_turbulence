"""
cosmology_plots.py is a Python routine that contains functions relevant
to generate the plots of the cosmology Jupyter notebook where Friedmann
equations are solved (cosmology.ipynb)

Author: Alberto Roper Pol
Date: 27/11/2022
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# get working directory, where the runs and routines should be stored
dir0 = os.getcwd() + '/'
HOME = dir0 + '/..'

# import GW_turbulence routines that are relevant for the project
os.chdir(HOME)
import plot_sets
os.chdir(dir0)

def plot_Omega_vs_a(a, Om_tot, OmR0, OmM0, aEQ, aL, save=True):
    
    """
    Function that plots the ratio of energy density to present-time
    critical energy density as a function of the scale factor a.
    
    Arguments:
        a -- array of scale factors
        Om_tot -- ratio of total energy density to present-time critical
                  energy density
        OmR0 -- ratio of present-time radiation to critical energy density
        OmM0 -- ratio of present-time matter to critical energy density
        aEQ -- scale factor at the time of equipartition
        aL -- scale factor at the time of dark-energy domination
        save -- option to save the plot under 'friedmann/plots/Omega_vs_a.pdf'
    """

    plt.figure(figsize=(8, 5))
    plt.plot(a, Om_tot, lw=2, color='blue')
    plt.loglog()
    plt.xlim(1e-16, 1e2)
    plt.xticks(np.logspace(-15, 0, 6))
    plt.ylim(1e-5, 1e65)
    plt.yticks(np.logspace(0, 60, 5))
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel(r'$\Omega(a)$', fontsize=36)
    plot_sets.axes_lines()

    plt.text(6e-13, 5e23, r'$\sim\!a^{-4}$', fontsize=28)
    plt.text(8e-4, 1e12, r'$\sim\!a^{-3}$', fontsize=28)
    xx = np.logspace(-15.6, -4.3)
    plt.plot(xx, .01*OmR0*(xx)**(-4), color='black', lw=.8)
    xx = np.logspace(-3.2, -.5)
    plt.plot(xx, 1e2*OmM0*xx**(-3), color='black', lw=.8)

    plt.vlines(aEQ, 1e-10, 1e70, color='black', lw=.8)
    plt.vlines(aL, 1e-10, 1e70, color='black', lw=.8)
    plt.text(5e-4, 1e35, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(1.5, 1e3, r'$a_{\Lambda}$', fontsize=28)

    if save:
        plt.savefig('friedmann/plots/Omega_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_w_vs_a(a, w, OmL0, aEQ, aL, save=True):
    
    """
    Function that plots the equation of state w such that p = w rho
    as a function of the scale factor a.
    
    Arguments:
        a -- array of scale factors
        w -- ratio of pressure to energy density
        OmL0 -- ratio of present-time dark to critical energy density
        aEQ -- scale factor at the time of equipartition
        aL -- scale factor at the time of dark-energy domination
        save -- option to save the plot under 'friedmann/plots/w_vs_a.pdf'
    """
    
    
    plt.figure(figsize=(8, 5))
    plt.plot(a, w, color='blue')
    plt.xscale('log')
    plt.ylim(-1.1, .5)
    plt.xlim(1e-16, 1e3)
    plt.xticks(np.logspace(-15, 3, 7))
    plt.yticks(np.linspace(-1, .5, 4))
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel('$w$', fontsize=36)
    plot_sets.axes_lines()

    plt.hlines(1/3, 1e-8, 1e-3, color='black', lw=.8)
    plt.hlines(0, 1e-5, 1e7, color='black', lw=.8)
    plt.hlines(-1, 1e0, 2e7, color='black', lw=.8)
    plt.hlines(-OmL0, 1e0, 2e7, color='black', lw=.8)

    plt.text(6e0, -.6, r'$\Omega_{\Lambda}^0$', fontsize=28)
    plt.text(1e-12, .15, r'$w = {1\over3}$', fontsize=28)
    plt.text(2e-6, -.5, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(2, .26, r'$a_{\Lambda}$', fontsize=28)

    plt.vlines(aEQ, -2, 3, color='black', lw=.8)
    plt.vlines(aL, -2, 3, color='black', lw=.8)

    if save:
        plt.savefig('friedmann/plots/w_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_Hubble_vs_a(a, ad, ap, H0, aEQ, aL, save=True):
    
    """
    Function that plots the Hubble rates \dot a/a and a'/a divided by H0 as
    a function of the scale factor a.
    
    Arguments:
        a -- array of scale factors
        ad -- cosmic time derivative of a
        ap -- conformal time derivative of a
        H0 -- current Hubble rate
        aEQ -- scale factor at the time of equipartition
        aL -- scale factor at the time of dark-energy domination
        save -- option to save the plot under 'friedmann/plots/Hubble_vs_a.pdf'
                and 'friedmann/plots/confHubble_vs_a.pdf'
    """
    
    ### plot the Hubble rate \dot{a}/a vs a
    plt.figure(figsize=(8, 5))
    plt.plot(a, ad/a/H0, color='blue')
    plt.loglog()
    plt.xlim(1e-18, 1e0)
    plt.xticks(np.logspace(-18, 0, 7))
    plt.ylim(1e-2, 1e35)
    plt.yticks(np.logspace(0, 32, 5))
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel('$H/H_0$', fontsize=36)
    plot_sets.axes_lines()
    
    plt.vlines(aL, 1e-20, 1e40, color='black', lw=.8)
    plt.vlines(aEQ, 1e-20, 1e40, color='black', lw=.8)

    xx = np.logspace(-14, -5)
    plt.plot(xx, 3e-1/xx**2, color='black', lw=.8)
    xx = np.logspace(-3, -.5)
    plt.plot(xx, 20/xx**1.5, color='black', lw=.8)

    plt.text(1e-10, 1e21, r'$a^{-2}$', fontsize=28)
    plt.text(7e-3, 1e6, r'$a^{-{3\over2}}$', fontsize=28)
    plt.text(6e-6, 1e20, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(6e-2, 1e18, r'$a_{\Lambda}$', fontsize=28)
    
    if save:
        plt.savefig('friedmann/plots/Hubble_vs_a.pdf',
                    bbox_inches='tight')
    
    ### plot the conformal Hubble rate a'/a vs a

    plt.figure(figsize=(8, 5))
    plt.plot(a, ap/a/H0, color='blue')
    plt.loglog()
    plt.xlim(1e-18, 1e0)
    plt.xticks(np.logspace(-18, 0, 7))
    plt.ylim(1e-1, 1e16)
    plt.yticks(np.logspace(0, 16, 5))
    plt.vlines(aL, 1e-20, 1e30, color='black', lw=.8)
    plt.vlines(aEQ, 1e-20, 1e30, color='black', lw=.8)
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel('$a_0^{-1} \, {\cal H}/H_0$', fontsize=36)
    plot_sets.axes_lines()

    xx = np.logspace(-14, -5)
    plt.plot(xx, .05/xx, color='black', lw=.8)
    xx = np.logspace(-3, -.5)
    plt.plot(xx, 4/xx**.5, color='black', lw=.8)

    plt.text(1e-10, 3e9, r'$a^{-1}$', fontsize=28)
    plt.text(5e-3, 2e2, r'$a^{-{1\over2}}$', fontsize=28)
    plt.text(6e-6, 1e7, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(6e-2, 1e9, r'$a_{\Lambda}$', fontsize=28)
    
    if save:
        plt.savefig('friedmann/plots/confHubble_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_derHubble_vs_a(a, add, app, H0, aEQ, aL, a_acc, save=True):
    
    """
    Function that plots the acceleration of the universe divided by H0^2 as
    a function of the scale factor a.
    
    Arguments:
        a -- array of scale factors
        add -- second cosmic time derivative of a
        app -- second conformal time derivative of a
        H0 -- current Hubble rate
        aEQ -- scale factor at the time of equipartition
        aL -- scale factor at the time of dark-energy domination
        a_acc -- scale factor at which the Universe acceleration becomes positive
        save -- option to save the plot under 'friedmann/plots/add_vs_a.pdf'
                and 'friedmann/plots/app_vs_a.pdf'
    """
    
    #### plot the accelerated rate of Universe's expansion \ddot{a}/a vs a
    plt.figure(figsize=(8, 5))
    plt.plot(a, abs(add/a/H0**2), color='blue')
    plt.loglog()
    plt.xlim(1e-18, 1e0)
    plt.xticks(np.logspace(-18, 0, 7))
    plt.ylim(1e-2, 1e70)
    plt.yticks(np.logspace(0, 64, 5))
    plot_sets.axes_lines()
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel(r'$|\ddot{a}|/(aH_0^2)$', fontsize=36)
    
    plt.vlines(a_acc, 1e-20, 1e70, color='black', lw=.8)
    plt.vlines(aEQ, 1e-20, 1e70, color='black', lw=.8)

    xx = np.logspace(-14, -5)
    plt.plot(xx, 1e-2/xx**4, color='black', lw=.8)
    xx = np.logspace(-3, -.5)
    plt.plot(xx, 50/xx**3, color='black', lw=.8)

    plt.text(1e-10, 1e40, r'$a^{-4}$', fontsize=28)
    plt.text(1e-2, 1e10, r'$a^{-3}$', fontsize=28)
    plt.text(6e-6, 1e48, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(1e-2, 1e36, r'$a_{\rm acc}$', fontsize=28)
    
    if save:
        plt.savefig('friedmann/plots/add_vs_a.pdf',
                    bbox_inches='tight')
    
    #### plot the conformal accelerated rate of Universe's
    #### expansion a''/a vs a
    plt.figure(figsize=(8, 5))
    plt.plot(a, app/a/H0**2, color='blue')
    plt.loglog()
    plt.xlim(1e-18, 1e0)
    plt.xticks(np.logspace(-18, 0, 7))
    plt.ylim(1e-1, 1e18)
    plt.yticks(np.logspace(0, 16, 5))
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel("$a_0^{-2} \, a''/(a H_0^2)$", fontsize=36)
    plot_sets.axes_lines()
    
    plt.vlines(aL, 1e-20, 1e30, color='black', lw=.8)
    plt.vlines(aEQ, 1e-20, 1e30, color='black', lw=.8)

    xx = np.logspace(-14, -2)
    plt.plot(xx, 1/xx, color='black', lw=.8)

    plt.text(1e-10, 6e10, r'$a^{-1}$', fontsize=28)
    plt.text(6e-6, 1e13, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(6e-2, 1e9, r'$a_{\Lambda}$', fontsize=28)
    
    if save:
        plt.savefig('friedmann/plots/app_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_eta_vs_t(t, eta, t0, eta0, tEQ, etaEQ, save=True):
    
    """
    Function that plots confomal vs cosmic time.
    
    Arguments:
        t -- array of cosmic times
        eta -- array of conformal times
        t0 -- value of cosmic present time
        eta0 -- value of conformal present time
        tEQ -- cosmic time at equipartition
        etaEQ -- conformal time at equipartition
        save -- option to save the plot under 'friedmann/plots/eta_vs_t.pdf'
    """
    
    import astropy.units as u
    
    plt.figure(figsize=(8, 5))
    plt.plot(t, eta, color='blue')
    plt.loglog()
    plt.xlim(1e-14, t0.value*1e1)
    plt.ylim(1e3, 1e19)
    plt.xticks(np.logspace(-14, 18, 5))
    plt.yticks(np.logspace(4, 18, 8))
    plt.xlabel('$t$ [s]', fontsize=36)
    plt.ylabel(r'$a_0 \, \eta$ [s]', fontsize=36)
    plot_sets.axes_lines()
    
    xx = np.logspace(-10, 6)
    plt.plot(xx, 4e10*xx**.5, color='black', lw=.8)
    xx = np.logspace(14, 17)
    plt.plot(xx, 5e11*xx**(1/3), color='black', lw=.8)
    plt.vlines(t0.value, 1e0, 1e30, color='black', lw=.8)
    plt.vlines(tEQ.to(u.s).value, 1e0, 1e30, color='black', lw=.8)
    
    plt.text(1e-3, 7e10, r'$\sim\!t^{1\over2}$', fontsize=28)
    plt.text(1e14, 1e15, r'$\sim\!t^{1\over3}$', fontsize=28)
    plt.text(2e9, 1e8, r'$t_{\rm EQ}$', fontsize=28)
    plt.text(1e16, 1e10, r'$t_0$', fontsize=28)
    
    if save:
        plt.savefig('friedmann/plots/eta_vs_t.pdf',
                    bbox_inches='tight')
        
def plot_a_vs_eta(eta, a, eta0, save=True):
    
    """
    Function that plots scale factor vs conformal time.
    
    Arguments:
        eta -- array of conformal times
        a -- array of scale factors
        eta0 -- value of conformal present time
        save -- option to save the plot under 'friedmann/plots/a_vs_eta.pdf'
    """
    
    plt.figure(figsize=(8, 5))
    plt.plot(eta, np.log10(a), color='blue')
    plt.xscale('log')
    plt.xlim(1e4, 4e18)
    plt.ylim(-17, 2)
    plt.xticks(np.logspace(4, 18, 8))
    plt.yticks(np.linspace(-16, 2, 10))
    plt.xlabel(r'$a_0 \, \eta \ [$s$]$', fontsize=36)
    plt.ylabel(r'$\log_{10}(a/a_0)$', fontsize=36)
    plot_sets.axes_lines()

    xx = np.logspace(5, 15)
    plt.plot(xx, np.log10(4e-21*xx), color='black', lw=.8)
    xx = np.logspace(17.4, 18.1)
    plt.plot(xx, np.log10(1.5e-37*xx**2), color='black', lw=.8)
    
    plt.hlines(1, 1e3, 1e20, color='black', ls='dashed', lw=.8)
    plt.vlines(eta0.value, -18, 2, color='black', ls='dashed', lw=.8)

    plt.text(3e9, np.log10(4e-13), r'$\sim\!\eta$', fontsize=28)
    plt.text(1e16, np.log10(1e-1), r'$\sim\!\eta^2$', fontsize=28)
    plt.text(2e17, -4, r'$\eta_0$', fontsize=28)

    if save:
        plt.savefig('friedmann/plots/a_vs_eta.pdf',
                    bbox_inches='tight')
        
def plot_a_vs_t(t, a, H0, tEQ, tL, t0, save=True):
    
    """
    Function that plots scale factor vs cosmic time.
    
    Arguments:
        t -- array of cosmic times
        a -- array of scale factors
        H0 -- present-time Hubble rate
        tL -- value of cosmic time at dark energy domination
        t0 -- value of cosmic present time
        save -- option to save the plot under 'friedmann/plots/a_vs_t.pdf'
    """
    
    import astropy.units as u
    
    plt.figure(figsize=(8, 5))
    plt.plot(t, np.log10(a), color='blue')
    plt.xscale('log')
    plt.xlim(1e-9, 2e19)
    plt.ylim(-17, 2)
    plt.xticks(np.logspace(-10, 20, 7))
    plt.yticks(np.linspace(-16, 2, 10))
    plt.xlabel(r'$t \ [$s$]$', fontsize=36)
    plt.ylabel(r'$\log_{10}(a/a_0)$', fontsize=36)
    plot_sets.axes_lines()

    xx = np.logspace(-23, -7)
    plt.plot(xx/H0, np.log10(.04*np.sqrt(xx)), color='black', lw=.8)
    xx = np.logspace(-4, -1)
    plt.plot(xx/H0, np.log10(.1*xx**(2/3)), color='black', lw=.8)
    xx = np.logspace(-1, 1)
    plt.plot(.6*xx/H0, np.log10(np.exp(xx)), color='black', lw=.8)

    plt.vlines(tEQ.to(u.s).value, -17, 2, color='black', ls='dashed', lw=.8)
    plt.vlines(tL.to(u.s).value, -17, 2, color='black', ls='dashed', lw=.8)
    plt.vlines(t0.to(u.s).value, -18, 2, color='black', ls='dashed', lw=.8)
    
    plt.text(3e9, -10, r'$t_{\rm EQ}$', fontsize=28)
    plt.text(8e15, -8, r'$t_{\Lambda}$', fontsize=28)
    plt.text(1e18, -11, r'$t_0$', fontsize=28)
    plt.text(2e0, np.log10(6e-12), r'$\sim\!t^{1/2}$', fontsize=28)
    plt.text(1e13, np.log10(1e-5), r'$\sim\!t^{2/3}$', fontsize=28)
    plt.text(4e12, -1, r'$\sim\!e^{H_0\, t}$', fontsize=22)

    if save:
        plt.savefig('friedmann/plots/a_vs_t.pdf',
                    bbox_inches='tight')
        
def plot_gs_gS_vs_T(T, a, g, gS, save=True):
    
    """
    Function that plots the evolution of relativistic and adiabatic
    degrees of freedom vs the scale factor during the RD era.
    
    Arguments:
        T -- array of temperatures
        a -- array of scale factors
        g -- array of relativistic dofs
        gS -- array of adiabatic dofs
        save -- option to save the plot under 'friedmann/plots/g_vs_a.pdf'
    """
    
    plt.figure(figsize=(8, 5))
    plt.plot(a, g, color='blue', label=r'$g_*$')
    plt.plot(a, gS, color='blue', ls='dashed', label=r'$g_{\rm S}$')
    plt.loglog()
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel(r'$g$', fontsize=36)
    plt.xlim(1e-18, 1e-6)
    plt.xticks(np.logspace(-18, -6, 7))
    plt.ylim(2, 200)
    plt.legend(frameon=False, loc='upper right', fontsize=28)
    plot_sets.axes_lines()

    if save:
        plt.savefig('friedmann/plots/g_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_Omrad_rat_vs_a(a, Om_rat, save=True):
    
    """
    Function that plots the factor of the radiation energy density due
    to the variation of degrees of freedom vs the scale factor during
    the RD era.
    
    Arguments:
        a -- array of scale factors
        Om_rat -- factor of Omega radiation
        save -- option to save the plot under
                'friedmann/plots/rho_vs_a_rat.pdf'
    """
    
    plt.figure(figsize=(8, 5))
    plt.plot(a, Om_rat, '.', color='blue')
    plt.xscale('log')
    plt.ylim(0.3, 1.05)
    plt.xlim(1e-16, 1e-5)
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel(r'$\Biggl(\frac{g_*}{g_*^0}\Biggr)' + \
               r'\Biggl(\frac{g_{\rm S}}{g_{\rm S}^0}\Biggr)^{-{4\over3}}$',
               fontsize=36)

    if save:
        plt.savefig('friedmann/plots/rho_vs_a_rat.pdf',
                    bbox_inches='tight')
        
def plot_dofs_rad(a, w, Om_mat, Om_rad, save=True):
    
    """
    Function that plots the function 1 - 3w vs a
    
    Arguments:
        a -- array of scale factors
        w -- equation of state
        Om_mat -- matter energy density (normalized)
        Om_rad -- radiation energy density (normalized)
        save -- option to save the plot under
                'friedmann/plots/Dw_vs_a_rad.pdf'
    """
    
    plt.figure(figsize=(8, 5))

    plt.plot(a, 1 - 3*w, color='blue')
    plt.plot(a, Om_mat/Om_rad, color='blue', ls='dashed', lw=.8)
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel(r'$1 - 3w(a)$', fontsize=36)
    plt.loglog()
    plt.ylim(1e-13, 1e1)
    plt.xlim(1e-16, 1e3)
    plt.xticks(np.logspace(-15, 3, 7))
    plt.yticks(np.logspace(-12, 0, 7))
    plot_sets.axes_lines()
    
    xx = np.logspace(-15, -4)
    plt.plot(xx, 1e3*xx, color='black', lw=.8)
    plt.text(2e-9, 2e-7, r'$\sim\!a$', fontsize=28)

    if save:
        plt.savefig('friedmann/plots/Dw_vs_a_rad.pdf',
                    bbox_inches='tight')
    
def plot_etaH_vs_a(a, ap, eta, aEQ, aL, save=True):
    
    """
    Function that plots the product of conformal Hubble rate and
    conformal time, which becomes ~ 1 during RD era and ~ 2 during
    MD era.
    
    Arguments:
        a -- array of scale factors
        ap -- array of conformal time derivative of scale factors
        eta -- array of conformal time
        aEQ -- scale factor at equipartition
        aL -- scale factor at dark energy domination
        save -- option to save the plot under
                'friedmann/plots/etaH_vs_a_rad.pdf'
    """
    
    plt.figure(figsize=(8, 5))
    plt.plot(a, ap/a*eta - 1, color='blue')
    plt.xlabel(r'${a\over a_0}$', fontsize=36)
    plt.ylabel(r'$\eta {\cal H} - 1$', fontsize=36)
    plt.loglog()
    plt.ylim(1e-5, 3e1)
    plt.xlim(1e-16, 1e1)
    plt.xticks(np.logspace(-15, 0, 6))
    plt.yticks(np.logspace(-5, 1, 7))
    plot_sets.axes_lines()
    
    plt.hlines(1, 1e-18, 1e4, color='black', lw=.8, ls='dashed')
    plt.vlines(aEQ, 1e-6, 1e2, color='black', lw=.8)
    plt.vlines(aL, 1e-6, 1e2, color='black', lw=.8)
    plt.text(5e-4, 1e-3, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(5e-2, 5e-2, r'$a_{\Lambda}$', fontsize=28)
    
    if save:
        plt.savefig('friedmann/plots/etaH_vs_a_rad.pdf',
                    bbox_inches='tight')
        
def plot_an_vs_etan(eta_n, a_n, eta_n_0, eta_n_EQ, ast, epoch='EWPT', save=True):
    
    """
    Function that plots normalized scale factor vs normalized conformal time.
    """

    plt.figure(figsize=(8, 5))
    plt.plot(eta_n, a_n, color='blue')
    plt.xlim(1, 1e14)
    plt.ylim(1, 1e18)
    plt.loglog()
    plot_sets.axes_lines()
    plt.xlabel(r'$\eta/\eta_*$')
    plt.ylabel(r'$a$')
    plt.text(3e12, 1e12, '$\eta_0$', fontsize=28)
    plt.vlines(eta_n_0, 1e0, 1e20, color='black', lw=.8)
    plt.vlines(eta_n_EQ, 1e0, 1e20, color='black', lw=.8)
    plt.text(5e9, 1e5, r'$\eta_{\rm EQ}$', fontsize=28)
    plt.title(r'$a_* = %.3f \times 10^{-16}$'%(ast*1e16), pad=10)
    if save:
        print("saved figure 'friedmann/plots/HHn_vs_etan_%s.pdf'"%epoch)
        plt.savefig('friedmann/plots/an_vs_etan_' + epoch +'.pdf',
                    bbox_inches='tight')
        
def plot_HHn_vs_etan(eta_n, HH_n, eta_n_0, eta_n_EQ, ast, epoch='EWPT', save=True):

    plt.figure(figsize=(8, 5))
    plt.plot(eta_n, HH_n, color='blue')
    plt.xlim(1, 1e14)
    plt.ylim(1e-14, 1)
    plt.loglog()
    plot_sets.axes_lines()
    plt.xlabel(r'$\eta/\eta_*$')
    plt.ylabel(r'${\cal H} (\eta)$')
    plt.vlines(eta_n_0, 1e-15, 1e0, color='black', lw=.8)
    plt.vlines(eta_n_EQ, 1e-15, 1e20, color='black', lw=.8)
    plt.text(3e12, 1e-3, '$\eta_0$', fontsize=28)
    plt.text(5e9, 1e-12, r'$\eta_{\rm EQ}$', fontsize=28)
    plt.title(r'$a_* = %.3f \times 10^{-16}$'%(ast*1e16), pad=10)
    if save:
        print("saved figure 'friedmann/plots/HHn_vs_etan_%s.pdf'"%epoch)
        plt.savefig('friedmann/plots/HHn_vs_etan_' + epoch +'.pdf',
                    bbox_inches='tight')
        
def plot_app_a_n_vs_etan(eta_n, app_a_n, eta_n_0, eta_n_EQ, ast, epoch='EWPT', save=True):

    plt.figure(figsize=(8, 5))
    plt.plot(eta_n, app_a_n, color='blue')
    plt.xlim(1, 1e14)
    plt.ylim(1e-27, 1e-10)
    plt.loglog()
    plot_sets.axes_lines()
    plt.xlabel(r'$\eta/\eta_*$')
    plt.ylabel(r"$a''/a$")
    plt.vlines(eta_n_0, 1e-30, 1e0, color='black', lw=.8)
    plt.text(3e12, 1e-20, '$\eta_0$', fontsize=28)
    plt.vlines(eta_n_EQ, 1e-30, 1e20, color='black', lw=.8)
    plt.text(5e9, 1e-17, r'$\eta_{\rm EQ}$', fontsize=28)
    plt.title(r'$a_* = %.3f \times 10^{-16}$'%(ast*1e16), pad=10)
    if save:
        print("saved figure 'friedmann/plots/app_n_vs_etan_%s.pdf'"%epoch)
        plt.savefig('friedmann/plots/app_n_vs_etan_' + epoch +'.pdf',
                    bbox_inches='tight')
        
def plot_factor_app_a_normalized(a, app_a_n, fact_app, epoch='EWPT', save=True):
    
    plt.figure(figsize=(12, 8))
    plt.plot(a, app_a_n/a, color='blue')
    plt.plot(a, fact_app, color='blue', ls='dashed')
    plt.loglog()
    plt.xlim(1e-18, 1e0)
    plt.ylim(.5, 1e4)
    plot_sets.axes_lines()
    plt.xlabel(r'$\frac{a}{a_0}$', fontsize=34)
    plt.ylabel(r"$\frac{a''}{a} \, \frac{a_0}{a_*} \, \frac{1}{{\cal H}}$", fontsize=34)
    plt.text(1e-9, 6e2, r'$\frac{\Omega_{{\rm mat}, 0}}' + \
             r'{\Omega_{{\rm rad}, 0}} \frac{g_*}{g_*^0}' + \
             r' \biggl(\frac{g_{\rm S}}{g_{\rm S}}^0\biggr)^{-{4\over3}}$', fontsize=30)
    plt.hlines(4.5e3, 1e-18, 1, color='black', ls='dashed', lw=.8)
    plt.text(5e-9, 5.5e3, r'$4.5 \times 10^3$')
    if save:
        print("saved figure 'friedmann/plots/factor_app_normalized_%s.pdf'"%epoch)
        plt.savefig('friedmann/plots/factor_app_normalized_' + epoch +'.pdf',
                    bbox_inches='tight')