from __future__ import division
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
import scipy.io as sio
from scipy.ndimage.filters import gaussian_filter


def get_cc_data(datafile):
    if datafile.split('.')[-1] == 'dat':
        data      = np.loadtxt(datafile, delimiter='\t', skiprows=3)
        cols      = ['t0', 'sigma', 'freq', 'cc_real', 'cc_imag']
        data_dict = dict(zip(cols, data.T))

    elif datafile.split('.')[-1] == 'mat':
        data       = sio.loadmat(datafile)
        temp       = np.zeros((data['cc_real'].shape[1], 3))
        temp[:, 0] = data['cc_real'][0, :]
        temp[:, 1] = data['sigma'][0, :]
        temp[:, 2] = data['t0'][0, :]
        data_dict  = dict(zip(['cc_real', 'sigma', 't0'], temp.T))

    return data_dict


def plot_specgram(datafile, cleaned=True, normalize=True, df=1, dt=4, t_start=0, t_end=2000):
    cc_data = get_cc_data(datafile)
    times   = np.arange(t_start * dt//2, t_end * dt//2, dt//2)
    freqs   = np.arange(7, 500 + df, df)
    SNR     = np.zeros(shape=(freqs.shape[0], len(times)))

    for i in range(len(times)):
        st = len(freqs) * i
        et = st + len(freqs)
        SNR[:, i] = np.sqrt(np.square(cc_data['cc_real'][st:et])) / cc_data['sigma'][st:et]

    if normalize:
        for i in range(SNR.shape[0]):
            # SNR[i, :] /= np.median(SNR[i, :])
            SNR[i, :] /= np.median(np.array([x for x in SNR[i, :] if x != 0]))
        vmax = 10
    else:
        vmax = 1

    plt.pcolormesh(times, freqs, SNR, cmap='magma', vmin=0, vmax=vmax)
    plt.colorbar(label='SNR/Median SNR')
    plt.ylabel('Frequency (Hz)')
    plt.ylim([7, 200])
    plt.xlabel('Time (s)')

    if cleaned:
        name = 'chirp_subtracted.png'
        plt.title('CC Spectra Subtracted ({})'.format(int(cc_data['t0'][0])))
    else:
        name = 'chirp.png'
        plt.title('CC Spectra ({})'.format(int(cc_data['t0'][0])))

    plt.savefig(name)
    plt.close()



def find_edges(datafile, cleaned=True, normalize=True, df=1, dt=4, t_start=0, t_end=2000):
    cc_data = get_cc_data(datafile)
    times   = np.arange(t_start * dt//2, t_end * dt//2, dt//2)
    freqs   = np.arange(7, 500 + df, df)
    SNR     = np.zeros(shape=(freqs.shape[0], len(times)))

    for i in range(len(times)):
        st = len(freqs) * i
        et = st + len(freqs)
        SNR[:, i] = np.sqrt(np.square(cc_data['cc_real'][st:et])) / cc_data['sigma'][st:et]


    if normalize:
        for i in range(SNR.shape[0]):
            SNR[i, :] /= np.median(SNR[i, :])
        vmax = 10
    else:
        vmax = 1

    SNR = gaussian_filter(SNR, 1.0)
    max_snr = np.max(SNR)
    min_snr = np.min(SNR)

    for _ in range(2):
        for i in range(SNR.shape[0]):
            SNR[i, :] -= np.array(list(SNR[i, 1:]) + [0])

        for i in range(SNR.shape[0]):
            SNR[i, :] -= np.array([0] + list(SNR[i, :-1]))

    vmax = 1
    for i in range(SNR.shape[0]):
        for j in range(SNR.shape[1]):
            if SNR[i, j] >= max_snr / 20.0:
                SNR[i, j] == 1
            else:
                SNR[i, j] == 0

    SNR = gaussian_filter(SNR, 2.0)
    for i in range(SNR.shape[0]):
        for j in range(SNR.shape[1]):
            if SNR[i, j] >= 0.2:
                SNR[i, j] == 1
            else:
                SNR[i, j] == 0

    plt.pcolormesh(times, freqs, SNR, cmap='magma', vmin=0, vmax=vmax)
    plt.colorbar(label='SNR/Median SNR')
    plt.ylabel('Frequency (Hz)')
    plt.ylim([7, 200])
    # plt.xlim([350, 410])
    plt.xlabel('Time (s)')

    if cleaned:
        name = 'chirp_subtracted_edge.png'
        plt.title('CC Spectra Subtracted ({})'.format(int(cc_data['t0'][0])))
    else:
        name = 'chirp_edge.png'
        plt.title('CC Spectra ({})'.format(int(cc_data['t0'][0])))

    plt.savefig(name)
    plt.close()


# BBH
bbh_cleaned  = 'noise_bbh/H1L1_ccspectra.job1.trial1.mat'
bbh_original = ('/home/andrew.matas/MDC_frames/CE_single_job/cc_runs/noise_bbh/'
                'output/HL_a0/H1L1_ccspectra.job1.trial1.dat')
# BNS
bns_cleaned  = 'noise_bns/H1L1_ccspectra.job1.trial1.mat'
bns_original = ('/home/andrew.matas/MDC_frames/CE_single_job/cc_runs/noise_bns/'
                'output/HL_a0/H1L1_ccspectra.job1.trial1.dat')

print('Making cleaned spectrogram')
plot_specgram(bns_cleaned, normalize=True, df=0.25, t_start=0, t_end=1200)
find_edges(bns_original, normalize=True, df=0.25, t_start=0, t_end=1200)
print('Making original spectrogram')
plot_specgram(bns_original, cleaned=False, normalize=True, df=0.25, t_start=0, t_end=1200)
