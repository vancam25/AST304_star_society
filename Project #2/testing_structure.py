#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 16:33:36 2022

@author: steven
"""

import structure as stc
import astro_const as ac
import numpy as np

from observations import MassRadiusObservations

def test_routine(M_test, delta_m, eta, xi):
    
    Pc = stc.pressure_guess(M_test, ac.mue)
    
    m_res, r_res, p_res = stc.integrate(Pc, delta_m, eta, xi, ac.mue)
    
    print(m_res[-1], r_res[-1], p_res[-1])
    
    return(m_res[-1], r_res[-1], p_res[-1])

delta_m = ac.Msun
eta = 10**(-10)
xi = 0.1

test_routine(ac.Msun, delta_m, eta, xi)