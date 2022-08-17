# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 07:40:45 2019

@author: cimat
"""

import numpy as np

##p number = 0 , default example
##p number Healthy:  4=<p_number=<54
##p number unHealthy:  64=<p_number=<92

for i in np.arange(0,1):
    p_number = i
    exec(open("twalk_ogtt_incretins.py").read())
