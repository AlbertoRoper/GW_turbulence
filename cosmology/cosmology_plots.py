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

os.chdir('..')
import plot_sets
os.chdir('cosmology')

def plot_Omega_vs_a(a, Om_tot, OmR0, OmM0, aEQ, aL, save=True):

    plt.figure(figsize=(12, 8))

    plt.plot(a, Om_tot, lw=2, color='blue')

    xx = np.logspace(-15.6, -4.3)
    plt.plot(xx, .04*OmR0*(xx)**(-4), color='black', lw=.8)

    xx = np.logspace(-3.2, -.5)
    plt.plot(xx, .07*OmM0*xx**(-3), color='black', lw=.8)

    plt.loglog()

    plt.xticks(np.logspace(-16, 2, 10))
    plt.yticks(np.logspace(0, 60, 7))
    plt.xlim(1e-16, 1e2)
    plt.ylim(1e-5, 1e60)

    plt.text(6e-13, 5e29, r'$\sim\!a^{-4}$', fontsize=32)
    plt.text(4e-4, 1e-2, r'$\sim\!a^{-3}$', fontsize=32)

    plt.xlabel(r'${a\over a_0}$', fontsize=40)
    plt.ylabel(r'$\Omega(a)$', fontsize=36)
    
    plt.vlines(aEQ, 1e-10, 1e70, color='black', lw=.8)
    plt.vlines(aL, 1e-10, 1e70, color='black', lw=.8)

    plt.text(5e-4, 1e15, r'$a_{\rm EQ}$', fontsize=32)
    plt.text(1.5, 1e3, r'$a_{\Lambda}$', fontsize=32)

    plot_sets.axes_lines()

    if save:
        plt.savefig('friedmann/plots/Omega_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_w_vs_a(a, w, OmL0, aEQ, aL, save=True):
    
    plt.figure(figsize=(8, 5))
    plt.plot(a, w, color='blue')

    plt.hlines(1/3, 1e-8, 1e-3, color='black', lw=.8)
    plt.hlines(0, 1e-5, 1e7, color='black', lw=.8)
    plt.hlines(-1, 1e0, 2e7, color='black', lw=.8)
    plt.hlines(-OmL0, 1e0, 2e7, color='black', lw=.8)

    plt.text(6e0, -.6, r'$\Omega_{\Lambda}^0$', fontsize=28)
    plt.text(1e-12, .15, r'$w = {1\over3}$', fontsize=28)

    plt.xscale('log')
    plt.ylim(-1.1, .5)
    plt.xlim(1e-16, 1e3)

    plt.xticks(np.logspace(-15, 3, 7))
    plt.yticks(np.linspace(-1, .5, 4))

    plot_sets.axes_lines()

    plt.text(2e-6, -.5, r'$a_{\rm EQ}$', fontsize=28)
    plt.text(2, .26, r'$a_{\Lambda}$', fontsize=28)

    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel('$w$')

    plt.vlines(aEQ, -2, 3, color='black', lw=.8)
    plt.vlines(aL, -2, 3, color='black', lw=.8)

    if save:
        plt.savefig('friedmann/plots/w_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_Hubble_vs_a(a, ad, H0, aL, aEQ, ap, save=True):
    
    ### plot the Hubble rate \dot{a}/a vs a
    plt.figure(figsize=(8, 5))
    plt.plot(a, ad/a/H0, color='blue')
    plt.loglog()
    plt.xlim(1e-18, 1e0)
    plt.xticks(np.logspace(-18, 0, 7))
    plt.ylim(1e-1, 1e34)
    plt.yticks(np.logspace(0, 32, 5))
    plot_sets.axes_lines()
    plt.vlines(aL, 1e-20, 1e40, color='black', lw=.8)
    plt.vlines(aEQ, 1e-20, 1e40, color='black', lw=.8)
    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel('$H/H_0$')

    xx = np.logspace(-14, -5)
    plt.plot(xx, 1/xx**2, color='black', lw=.8)

    xx = np.logspace(-3, -.5)
    plt.plot(xx, 50/xx**1.5, color='black', lw=.8)

    plt.text(1e-10, 1e21, r'$a^{-2}$')
    plt.text(1e-2, 1e6, r'$a^{-{3\over2}}$')
    plt.text(1e-5, 1e20, r'$a_{\rm EQ}$')
    plt.text(1e-1, 1e18, r'$a_{\Lambda}$')
    
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
    plot_sets.axes_lines()
    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel('$a_0^{-1} \, {\cal H}/H_0$')

    xx = np.logspace(-14, -5)
    plt.plot(xx, .1/xx, color='black', lw=.8)

    xx = np.logspace(-3, -.5)
    plt.plot(xx, 4/xx**.5, color='black', lw=.8)

    plt.text(1e-10, 5e9, r'$a^{-1}$')
    plt.text(1e-2, 3e2, r'$a^{-{1\over2}}$')

    plt.text(1e-5, 1e13, r'$a_{\rm EQ}$')
    plt.text(1e-1, 1e9, r'$a_{\Lambda}$')
    
    if save:
        plt.savefig('friedmann/plots/confHubble_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_derHubble_vs_a(a, add, H0, aL, aEQ, app, a_acc, save=True):
    
    #### plot the accelerated rate of Universe's expansion \ddot{a}/a vs a
        
    plt.figure(figsize=(8, 5))
    plt.plot(a, abs(add/a/H0**2), color='blue')
    plt.loglog()
    plt.xlim(1e-18, 1e0)
    plt.xticks(np.logspace(-18, 0, 7))
    plt.ylim(1e-2, 1e70)
    plt.yticks(np.logspace(0, 64, 5))
    plot_sets.axes_lines()
    plt.vlines(a_acc, 1e-20, 1e70, color='black', lw=.8)
    plt.vlines(aEQ, 1e-20, 1e70, color='black', lw=.8)
    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel(r'$|\ddot{a}|/(aH_0^2)$')

    xx = np.logspace(-14, -5)
    plt.plot(xx, 1/xx**4, color='black', lw=.8)

    xx = np.logspace(-3, -.5)
    plt.plot(xx, 50/xx**3, color='black', lw=.8)

    plt.text(1e-10, 1e42, r'$a^{-4}$')
    plt.text(1e-2, 1e10, r'$a^{-3}$')
    plt.text(1e-5, 1e48, r'$a_{\rm EQ}$')
    plt.text(2e-2, 1e36, r'$a_{\rm acc}$')
    
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
    plt.vlines(aL, 1e-20, 1e30, color='black', lw=.8)
    plt.vlines(aEQ, 1e-20, 1e30, color='black', lw=.8)
    plot_sets.axes_lines()
    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel("$a_0^{-2} \, a''/(a H_0^2)$")

    xx = np.logspace(-14, -2)
    plt.plot(xx, 4/xx, color='black', lw=.8)

    plt.text(1e-10, 1e11, r'$a^{-1}$')
    plt.text(1e-5, 1e13, r'$a_{\rm EQ}$')
    plt.text(1e-1, 1e9, r'$a_{\Lambda}$')
    
    if save:
        plt.savefig('friedmann/plots/app_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_eta_vs_t(t, eta, t0, eta0, tEQ, etaEQ, save=True):
    
    import astropy.units as u
    
    plt.figure(figsize=(8, 5))
    plt.plot(t, eta, color='blue')
    xx = np.logspace(-10, 6)
    plt.plot(xx, 4e10*xx**.5, color='black', lw=.8)
    xx = np.logspace(14, 18)
    plt.plot(xx, 4e12*xx**(1/3), color='black', lw=.8)
    plt.loglog()
    plt.ylim(1e3, 1e19)
    plt.xlim(1e-14, t0.value*1e1)
    plt.xlabel('$t$')
    plt.ylabel(r'$a_0 \, \eta$')
    plt.vlines(t0.value, 1e0, 1e30, color='black', lw=.8)
    plt.hlines(eta0.value, 1e-20, 1e30, color='black', lw=.8)
    plt.vlines(tEQ.to(u.s).value, 1e0, 1e30, color='black', lw=.8)
    plt.hlines(etaEQ.to(u.s).value, 1e-20, 1e30, color='black', lw=.8)
    plt.xticks(np.logspace(-14, 18, 5))
    plt.yticks(np.logspace(4, 18, 8))
    plot_sets.axes_lines()
    plt.text(1e-3, 3e10, r'$\sim\!t^{1\over2}$')
    plt.text(4e14, 2e16, r'$\sim\!t^{1\over3}$')
    plt.text(4e9, 1e8, 'EQ')
    plt.text(1e16, 1e10, r'$t_0$')
    
    if save:
        plt.savefig('friedmann/plots/eta_vs_t.pdf',
                    bbox_inches='tight')
        
def plot_a_vs_eta(eta, a, eta_0, save=True):
    
    ##### plot the conformal time evolution of the scale factor

    plt.figure(figsize=(8, 5))

    plt.plot(eta, np.log10(a), color='blue')

    plt.xscale('log')
    plot_sets.axes_lines()
    plt.xlabel(r'$a_0 \, \eta \ [$s$]$')
    plt.ylabel(r'$\log_{10}(a/a_0)$')
    plt.xlim(1e4, 4e18)
    plt.ylim(-17, 2)
    plt.xticks(np.logspace(4, 18, 8))
    plt.yticks(np.linspace(-16, 2, 10))

    xx = np.logspace(5, 15)
    plt.plot(xx, np.log10(4e-21*xx), color='black', lw=.8)

    xx = np.logspace(17.4, 18.2)
    plt.plot(xx, np.log10(1.5e-37*xx**2), color='black', lw=.8)

    plt.text(3e9, np.log10(1e-13), r'$\sim\!\eta$', fontsize=28)
    plt.text(1e16, np.log10(1e-1), r'$\sim\!\eta^2$', fontsize=28)

    plt.hlines(1, 1e3, 1e20, color='black', ls='dashed', lw=.8)
    plt.vlines(eta_0.value, -18, 2, color='black', ls='dashed', lw=.8)

    plt.text(3e17, -4, r'$\eta_0$')

    if save:
        plt.savefig('friedmann/plots/a_vs_eta.pdf',
                    bbox_inches='tight')
        
def plot_a_vs_t(t, a, H0, tEQ, tL, t0, save=True):
        
    ##### plot the cosmic time evolution of the scale factor

    plt.figure(figsize=(8, 5))

    plt.plot(t, np.log10(a), color='blue')

    plt.xscale('log')
    plot_sets.axes_lines()
    plt.xlabel(r'$t \ [$s$]$')
    plt.ylabel(r'$\log_{10}(a/a_0)$')
    plt.xlim(1e-9, 2e19)
    plt.ylim(-17, 2)
    plt.xticks(np.logspace(-10, 20, 7))
    plt.yticks(np.linspace(-16, 2, 10))

    plt.text(7e0, np.log10(4e-12), r'$\sim\!t^{1/2}$', fontsize=28)
    plt.text(1e13, np.log10(1e-5), r'$\sim\!t^{2/3}$', fontsize=28)
    plt.text(4e12, -1, r'$\sim\!e^{H_0\, t}$', fontsize=22)

    xx = np.logspace(-23, -7)
    plt.plot(xx/H0, np.log10(.02*np.sqrt(xx)), color='black', lw=.8)

    xx = np.logspace(-4, -1)
    plt.plot(xx/H0, np.log10(.1*xx**(2/3)), color='black', lw=.8)

    xx = np.logspace(-1, 1)
    plt.plot(.6*xx/H0, np.log10(np.exp(xx)), color='black', lw=.8)

    plt.vlines(tEQ.value, -17, 2, color='black', ls='dashed', lw=.8)
    plt.vlines(tL.value, -17, 2, color='black', ls='dashed', lw=.8)

    plt.vlines(t0.value, -18, 2, color='black', ls='dashed', lw=.8)
    plt.text(1e10, -10, r'$t_{\rm EQ}$')
    plt.text(1e16, -8, r'$t_{\Lambda}$')
    plt.text(1e18, -11, r'$t_0$')

    if save:
        plt.savefig('friedmann/plots/a_vs_t.pdf',
                    bbox_inches='tight')
        
def plot_gs_gS_vs_T(T, a, g, gS, save=True):
    
    plt.figure(figsize=(8, 5))

    plt.plot(a, g, color='blue', label=r'$g_*$')
    plt.plot(a, gS, color='blue', ls='dashed', label=r'$g_{\rm S}$')

    plot_sets.axes_lines()
    plt.legend(frameon=False, loc='upper right', fontsize=28)

    plt.xlim(1e-18, 1e-6)
    plt.ylim(2, 200)

    plt.xticks(np.logspace(-18, -6, 7))

    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel(r'$g$')

    plt.loglog()
    if save:
        plt.savefig('friedmann/plots/g_vs_a.pdf',
                    bbox_inches='tight')
        
def plot_Omrad_rat_vs_a(a, Om_rat, save=True):
    
    plt.figure(figsize=(8, 5))

    plt.xscale('log')
    plt.ylim(0.3, 1.05)
    plt.xlim(1e-16, 1e-5)

    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel(r'$\Biggl(\frac{g_*}{g_*^0}\Biggr)' + \
               r'\Biggl(\frac{g_{\rm S}}{g_{\rm S}^0}\Biggr)^{-{4\over3}}$')

    plt.plot(a, Om_rat, '.', color='blue')

    if save:
        plt.savefig('friedmann/plots/rho_vs_a_rat.pdf',
                    bbox_inches='tight')
        
def plot_dofs_rad(a, w, Om_mat, Om_rad, save=True):
    
    plt.figure(figsize=(8, 5))

    plt.plot(a, 1 - 3*w, color='blue')
    plt.plot(a, Om_mat/Om_rad, color='blue', ls='dashed', lw=.8)

    xx = np.logspace(-15, -4)
    plt.plot(xx, 1e3*xx, color='black', lw=.8)

    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel(r'$1 - 3w(a)$')
    plt.loglog()
    plt.ylim(1e-13, 1e1)
    plt.xlim(1e-16, 1e3)
    plt.xticks(np.logspace(-15, 3, 7))
    plt.yticks(np.logspace(-12, 0, 7))
    plot_sets.axes_lines()
    plt.text(2e-9, 2e-7, r'$\sim\!a$', fontsize=28)

    if save:
        plt.savefig('friedmann/plots/Dw_vs_a_rad.pdf',
                    bbox_inches='tight')
    
def plot_etaH_vs_a(a, ap, eta, aEQ, aL, save=True):
    
    plt.figure(figsize=(8, 5))
    plt.plot(a, ap/a*eta - 1, color='blue')

    plt.xlabel(r'${a\over a_0}$')
    plt.ylabel(r'$\eta {\cal H} - 1$')
    plt.hlines(1, 1e-18, 1e4, color='black', lw=.8, ls='dashed')
    plt.loglog()
    plt.ylim(1e-5, 3e1)
    plt.xlim(1e-16, 1e1)
    plt.xticks(np.logspace(-15, 0, 6))
    plt.yticks(np.logspace(-5, 1, 7))
    plot_sets.axes_lines()
    plt.vlines(aEQ, 1e-6, 1e2, color='black', lw=.8)
    plt.vlines(aL, 1e-6, 1e2, color='black', lw=.8)
    plt.text(5e-4, 1e-3, r'$a_{\rm EQ}$')
    plt.text(1e-1, 1e-4, r'$a_{\Lambda}$')
    
    if save:
        plt.savefig('friedmann/plots/etaH_vs_a_rad.pdf',
                    bbox_inches='tight')