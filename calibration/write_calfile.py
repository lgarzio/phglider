#!/usr/bin/env python

"""
Author: Lori Garzio on 3/3/2021
Last modified: 3/3/2021
"""
import json

data = {'f6': 6.3113E-22,
        'f5': -3.5667E-18,
        'f4': 6.0166E-15,
        'f3': 6.1183E-13,
        'f2': -1.1932E-08,
        'f1': 1.3118E-05,
        'k0': -1.421745,
        'k2': -1.0438E-03,
        }

with open('sbe10528_20200821.txt', 'w') as outfile:
    json.dump(data, outfile)
