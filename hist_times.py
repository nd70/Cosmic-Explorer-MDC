from __future__ import division
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np


def get_cbc_data(datafile):
    data = np.genfromtxt(datafile, delimiter=' ', filling_values=0)
    cols = ['number', 't0', 'tc', 'm1', 'm2', 'z', 'dist', 'ra', 'decl',
               'psi', 'inc', 'params.startPhase', 'type', 'grb', 'rho']

    data_dict = dict(zip(cols, data.T))
    return data_dict


def chirp_times(event_params,
                  Msun  = 1,
                  f0    = 137.0,
                  f_min = 7.0,
                  f_max = 500.0):

    # Unpack params
    t0 = event_params['t0']
    tc = event_params['tc']
    z  = event_params['z']
    m1 = event_params['m1']
    m2 = event_params['m2']

    # Calculate times at f_min and f_max
    Mc   = np.power(m1 * m2, 3/5) * (1 + z) / np.power(m1 + m2, 1/5)
    tmax = tc - np.power(f0/f_max, 8/3) * np.power(1.21 * Msun/Mc, 5/3)
    tmin = tc - np.power(f0/f_min, 8/3) * np.power(1.21 * Msun/Mc, 5/3)

    # Chirp duration
    dt = np.abs(tmax - tmin)

    return dt


# BBH Only
bbh_data = get_cbc_data('Lists/all_bbh.txt')
times = chirp_times(bbh_data)
plt.hist(times, bins=np.arange(0, 150, 1))
plt.title('{} BBH Events: $t_0$ to $t_c$'.format(len(times)))
plt.savefig('all_bbh.png')
plt.close()

times = chirp_times(bbh_data, f_max=10)
plt.hist(times, bins=np.arange(0, 150, 1))
plt.title('{} BBH Events: Time to 10Hz'.format(len(times)))
plt.savefig('all_bbh_10Hz.png')
plt.close()

times = chirp_times(bbh_data, f_max=15)
plt.hist(times, bins=np.arange(0, 150, 1))
plt.title('{} BBH Events: Time to 15Hz'.format(len(times)))
plt.savefig('all_bbh_15Hz.png')
plt.close()

times = chirp_times(bbh_data, f_max=100)
plt.hist(times, bins=np.arange(0, 150, 1))
plt.title('{} BBH Events: Time to 100Hz'.format(len(times)))
plt.savefig('all_bbh_100Hz.png')
plt.close()

# BNS + BBH
cbc_data = get_cbc_data('Lists/all_bns.txt')
times = chirp_times(cbc_data)
plt.hist(times, bins=np.arange(0, 2000, 10), label='BNS')
plt.hist(chirp_times(bbh_data),bins=np.arange(0, 150, 10), label='BBH')
plt.title('{} BBH + BNS Events: $t_0$ to $t_c$'.format(len(times)))
plt.savefig('all_cbc.png')
plt.close()

times = chirp_times(cbc_data, f_max=10)
plt.hist(times, bins=np.arange(0, 2000, 10), label='BNS')
plt.hist(chirp_times(bbh_data),bins=np.arange(0, 150, 10), label='BBH')
plt.title('BBH + BNS Time to 10Hz'.format(len(times)))
plt.savefig('all_cbc_10Hz.png')
plt.close()

times = chirp_times(cbc_data, f_max=15)
plt.hist(times, bins=np.arange(0, 2000, 10), label='BNS')
plt.hist(chirp_times(bbh_data),bins=np.arange(0, 150, 10), label='BBH')
plt.title('BBH + BNS Time to 15Hz'.format(len(times)))
plt.savefig('all_cbc_15Hz.png')
plt.close()

times = chirp_times(cbc_data, f_max=100)
plt.hist(times, bins=np.arange(0, 2000, 10), label='BNS')
plt.hist(chirp_times(bbh_data),bins=np.arange(0, 150, 10), label='BBH')
plt.title('BBH + BNS Time to 100Hz'.format(len(times)))
plt.savefig('all_cbc_100Hz.png')
plt.close()

# BNS only
bns_data = get_cbc_data('Lists/all_bns.txt')
times = chirp_times(bns_data)
plt.hist(times, bins=np.arange(0, 2000, 10))
plt.title('{} BNS Events: $t_0$ to $t_c$'.format(len(times)))
plt.savefig('all_bns.png')
plt.close()

times = chirp_times(bns_data, f_max=10)
plt.hist(times, bins=np.arange(0, 2000, 10))
plt.title('{} BNS Events: Time to 10Hz'.format(len(times)))
plt.savefig('all_bns_10Hz.png')
plt.close()

times = chirp_times(bns_data, f_max=15)
plt.hist(times, bins=np.arange(0, 2000, 10))
plt.title('{} BNS Events: Time to 15Hz'.format(len(times)))
plt.savefig('all_bns_15Hz.png')
plt.close()

times = chirp_times(bns_data, f_max=100)
plt.hist(times, bins=np.arange(0, 2000, 10))
plt.title('{} BNS Events: Time to 100Hz'.format(len(times)))
plt.savefig('all_bns_100Hz.png')
plt.close()
