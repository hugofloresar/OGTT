# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 07:40:45 2019

@author: cimat
"""

import numpy as np

##p number = 0 , default example
##p number NGT:  1=<p_number=<32
##p number with some anomaly:  39=<p_number=<60

for i in np.arange(0,1):
    p_number = i
    exec(open("twalk_ogtt_incretins.py").read())
