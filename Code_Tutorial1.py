#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 18:01:31 2018

@author: Matthias
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.special as sp
a_0 = 1.



def R_nl(r,n,l):
    rho = 2*r/(n*a_0)
    return np.sqrt( (2./(n*a_0))**2 * sp.factorial(n-l-1)/(2*n*sp.factorial(n+l))) * np.exp(-rho/2.)*rho**l*sp.eval_genlaguerre(n-l-1,2*l+1,rho)

def characteristic_radius(n,l,rmax,stepsize):
    '''
    computes the characteristic radius of the (n,l)-orbital
    n: Hauptquantenzahl
    l: Nebenquantenzahl
    rmax: upper boundary of 'integration'
    stepsize: stepsize of integration
    '''
    integration_range = np.linspace(0,rmax,float(rmax)/stepsize)
    sum_r = 0
    weight_sum = 0
    for r in integration_range:
        sum_r+=r*(R_nl(r,n,l))**2*4*np.pi*r**2
        weight_sum += (R_nl(r,n,l))**2*4*np.pi*r**2

    characteristic_r = sum_r/weight_sum
#    characteristic_r = sum_r

    return characteristic_r


def plot_R_nl(n,l,rmax,stepsize):
    '''
    plots R_nl in(n,l) in the specified range
    n: Hauptquantenzahl
    l: Nebenquantenzahl
    rmax: upper boundary of 'integration'
    stepsize: stepsize of integration
    '''
    R_nl_vals = []
    R_nl_sqrd  = []
    r_range = np.linspace(0,rmax, float(rmax)/stepsize)
    for r in r_range:
#        R_nl_vals.append(R_nl(r,n,l))
        R_nl_sqrd.append(R_nl(r,n,l)**2*4*np.pi*r**2)
#    plt.plot(r_range,R_nl_vals, 'r')
    plt.plot(r_range, R_nl_sqrd, 'b')
    
    return 0

def localization_cmap(nmin,nmax,lmin,lmax,stepsize):
    '''
    creates colorplot of localization depending on n (x-axis) and l (y-axis)
    nmin: minimal n-value
    nmax: maximum n value
    lmin: minimum l-velua
    lmax: maximum l-value
    stepsize: stepsize of r-values for the integration in computing the expectation value of r
    '''
    rmax = 200.
    lrange = range(lmin,lmax+1)
    nrange = range(nmin,nmax+1)
    
#    Z: array of characteristic r's
    Z = []
    
    for l in lrange:
        L = []
        for n in nrange:
            if l < n:
                L.append(characteristic_radius(n,l,rmax*n/4,stepsize))
            else:
                L.append(0)
        Z.append(L)
    
    p = plt.pcolor(Z, cmap='autumn')

    return Z
    
hello = localization_cmap(1,10 ,0,9,0.2)
#print(hello)
#print(np.random.rand(6,10))

#characteristic_r_6 = []

#for l in range(0,20):
#    characteristic_r_6.append(characteristic_radius(20,l,50.,0.01))
#
#plt.plot(range(0,20),characteristic_r_6)
#print(characteristic_r_6)

#plot_R_nl(10,0,500,0.01)
#print(characteristic_radius(10,0,200.,0.01))

#for n in range(1,8):
#    characteristic_r_6.append(characteristic_radius(n,0,200.,0.01))
#
#plt.plot(range(1,8),characteristic_r_6)

#
#plot_R_nl(3,0,50,0.01)
