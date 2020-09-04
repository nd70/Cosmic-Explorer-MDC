from __future__ import division
import numpy as np


def count_notches(datafile, cutoff=100.0):
    data = np.loadtxt(datafile, delimiter='\t', skiprows=3)
    cols = ['t0', 'sigma', 'freq', 'cc_real', 'cc_imag']
    data_dict = dict(zip(cols, data.T))
    reals = data_dict['cc_real']
    freqs = data_dict['freq']
    df = freqs[2] - freqs[1]
    segs = np.array(list(set(data_dict['t0'])))
    notched = 0
    for ix, num in enumerate(reals):
        if float(num) == 0.0 and float(freqs[ix]) <= cutoff:
            notched += 1
    seg_len = len(reals) / len(segs)
    lost = notched / seg_len
    print('{0} bins notched below {1}Hz (equiv {2} segs lost)'.format(notched, cutoff, lost))
    print('  -> {0:.2f}% of all data notched'.format(notched * 100.0 / len(reals)))
    print('  -> {0:.2f}% of data below {1}Hz notched'.format(notched *\
            100.0 / (cutoff  * len(segs)/ df), cutoff))
