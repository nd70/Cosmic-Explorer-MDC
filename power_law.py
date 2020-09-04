from __future__ import division
import numpy as np
import re
import os
import sys
import argparse
from operator import itemgetter
from ConfigParser import ConfigParser
import scipy.io as sio
import time
start_time = time.time()


# Change dir so this can be run from anywhere
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)


def parse_command_line():
    """
    parse command line
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--section", "-s",
                        help    = "section from config file",
                        default = "4s",
                        dest    = "section",
                        type    = str)

    parser.add_argument("--ini-file", "-i",
        				help    = "config file to read params from",
                        default = "configs.ini",
        				dest    = "ini_file",
                        type    = str)

    parser.add_argument("--verbose", "-v",
        				help    = "verbose printing",
                        default = False,
        				dest    = "verbose",
                        action  = 'store_true')

    params = parser.parse_args()
    return params


# Set Params
params   = parse_command_line()
section  = params.section
ini_file = params.ini_file
verbose  = params.verbose
settings = ConfigParser()
settings.read(ini_file)
cbc_PATH = settings.get(section, 'cbc_PATH')
cc_PATH  = settings.get(section, 'cc_PATH')
df       = settings.getfloat(section, 'df')
f0       = settings.getfloat(section, 'f0')
f_min    = settings.getfloat(section, 'f_min')
f_max    = settings.getfloat(section, 'f_max')
Msun     = settings.getfloat(section, 'Msun')
outdir   = settings.get(section, 'outdir')
seg_dur  = settings.getint(section, 'seg_dur')
forward  = settings.getint(section, 'forward')
backward = settings.getint(section, 'backward')

if not os.path.isdir(outdir):
    os.system('mkdir -p {}'.format(outdir))


def sort_files(file_list, keyword='ccspectra'):
    """
    collects the ccspectra files and sorts them in numerical order

    Parameters
    ----------
    file_list : `str`
        path to files

    Returns
    -------
    ordered : `list`
        list of ordered ccspectra files
    """

    reg  = re.compile(r'[0-9]+')
    temp = []
    for f in file_list:
        if keyword in f:
            num = int(reg.findall(f)[-2])
        else:
            num = int(reg.findall(f)[-1])
        temp.append((f, num))

    sort_out = sorted(temp, key=itemgetter(1))
    ordered  = [s[0] for s in sort_out]

    return ordered


def find_event_file(ccspec_files, cbc_data, seg_dur):
    """
    find all of the CBCs in a given job file

    Parameters
    ----------
    ccspec_files : `list`
        list containing the paths to all of the ccspectra job files

    cbc_data : `ndarray`
        numpy ndarray of shape (event, parameters) containing the
        CBC event data

    seg_dur : `int`
        segment duration used in the job file

    Returns
    -------
    event_file_list : `list`
        list of dicts containing the parameters of each event for each job file
        Ex. [{'Job1':{'times':[list-of-times], 'data':[cbc-data]}},
             {'Job2':{'times':[list-of-times], 'data':[cbc-data]}}
             ...
             ]

    header : `str`
        information from the first 3 lines of the jobfile
    """

    event_file_list = []
    for ccspec_file in ccspec_files:
        with open(ccspec_file) as cc:
            lines = cc.readlines()
            header = lines[:3]
            tmin, tmax = float(lines[3].split('\t')[0]), float(lines[-1].split('\t')[0])

        events = []
        for event in range(cbc_data.shape[0]):
            if (cbc_data[event][1] < tmax) and (float(cbc_data[event][2]) >= tmin):
                events.append(cbc_data[event])

        event_file_list.append({ccspec_file: dict(times=np.arange(tmin, tmax +\
                                                  seg_dur // 2,  seg_dur // 2),
                                                  data=np.array(events))})
    return event_file_list, header


def power_law(event_params, time,
              f_min = 7.0,
              f_max = 500.0,
              Msun  = 1.0,
              f0    = 134.0):

    """
    Assume we have the power law:

    f(t) = alpha / (t_c - t + delta)^(3/8)

    where delta keeps us from exploding at t_c = t_0.
    We require f(t_c) = f_max, which sets delta:

    delta = (alpha / f_max)^(8/3)

    Plugging this in and solving for alpha we find:

    alpha = [ (t_c - t_0) / (f_min^(-8/3) - f_max^(-8/3)) ]^(3/8)
    """

    # collect params
    event_params = [float(x) for x in event_params]
    _, t0, tc, _, _, _, _, _, _, _, _, _, _, _, _ = event_params

    # find alpha
    freq_diff = np.power(np.power(f_min, -8.0 / 3.0) - np.power(f_max, -8.0 / 3.0), -1.0)
    alpha = np.power((tc - t0) * freq_diff, 3.0 / 8.0)

    # assemble the chirp
    f_gw = 0.75 * alpha * np.power(tc - time + np.power(alpha / f_max, 8.0 / 3.0), -3.0 / 8.0)

    return f_gw


def smart_notching(mask, notch_st, freq_st, freq_end,
                   seg_dur  = 4,
                   df       = 0.25,
                   f_min    = 7,
                   f_max    = 500,
                   st_time  = 0,
                   forward  = 0,
                   backward = 0):
    """
    notch a frequency mask by calculating the positions of the bins to notch
    instead of looping over the (enormous) jobfile

    Parameters
    ----------
    mask : `ndarray`
        numpy array (initially of all 1s) frequency mask. "Notching" means
        setting bins to 0.

    notch_st : `int`
        segment where notching begins

    freq_st : `float`
        frequency bin where notching begins

    freq_end : `float`
        frequency bin where notching ends

    seg_dur : `int`
        segment duration used in the job file

    df : `float`
        frequency resolution

    f_min : `float`
        minimum frequency corresponding to tmin

    f_max : `float`
        maximum frequency corresponding to tmax

    st_time : `int`
        time at the start of the jobfile

    forward : `int`
        selects how many segments to notch at the current frequency in
        the future ("to the right" if you're looking at the spectrogram)

    backward : `int`
        selects how many segments to notch at the current frequency in
        the past ("to the left" if you're looking at the spectrogram)

    Returns
    -------
    mask : `ndarray`
        notched numpy ndarray frequency mask
    """

    # Calculate start and stop index for notching
    jumps = (notch_st - st_time) / (seg_dur // 2)  # how many segments to jump forward
    seg_start = (jumps * (f_max - f_min) / df) + jumps  # start of segment to notch
    index_st = int(seg_start + ((freq_st - f_min) / df))  # bin to start notching
    index_end = index_st + int((freq_end - freq_st) / df) + 1  # final bin to notch
    mask[index_st:index_end] = 0.0

    seg_size = (f_max - f_min) / df
    if forward > 0:
        for i in range(forward):
            try:
                forward_st = int(index_st + seg_size * (i + 1))
                forward_et = int(index_end + seg_size * (i + 1) + 1)
                mask[forward_st:forward_et] = 0.0
            except:
                pass

    if backward > 0:
        for i in range(backward):
            try:
                backward_st = int(index_st - seg_size * (i + 1))
                backward_et = int(index_end - seg_size * (i + 1) + 1)
                mask[backward_st:backward_et] = 0.0
            except:
                pass

    return mask


# get data files
ccspec_files = sort_files([cc_PATH + f for f in os.listdir(cc_PATH)
                           if "ccspec" in f])
sensint_files = sort_files([cc_PATH + f for f in os.listdir(cc_PATH)
                            if "sensints" in f])
# Get cbc event data
# ['number', 't0', 'tc', 'm1', 'm2', 'z', 'dist', 'ra', 'decl',
# 'psi', 'inc', 'params.startPhase', 'type', 'grb', 'rho']
event_file_list = []
cbc_data = np.genfromtxt(cbc_PATH, delimiter=' ', filling_values=0)

# Collect the events that belong in each file
event_file_list, header = find_event_file(ccspec_files, cbc_data, seg_dur)

# Start looping over chirps
for ii, jobfile in enumerate(event_file_list):
    if verbose:
        print('Running job {0}/{1}'.format(ii + 1, len(ccspec_files)))

    # unpack the event file list
    event_file = jobfile.keys()[0]
    jobfile = jobfile[event_file]
    data = jobfile['data']
    times = jobfile['times']
    job_start_time = times[0]

    # Make the mask for this job
    steps = int((f_max - f_min) / df) + 1
    ft_data = np.ones(shape=(len(times) * steps))

    # loop over each event within the job
    for ix in range(data.shape[0]):
        tmin, tmax = data[ix, 1], data[ix, 2]

        cc_min_diff = np.abs(times - tmin)
        cc_max_diff = np.abs(times - tmax)
        min_index = np.where(cc_min_diff == np.min(cc_min_diff))[0][0]
        max_index = np.where(cc_max_diff == np.min(cc_max_diff))[0][0]

        # Find how many segments need to be notched
        seg_st = times[min_index]
        segs_to_notch = int(((times[max_index] - times[min_index]) // (seg_dur // 2)))

        # Notch them
        for i in range(segs_to_notch):
            if verbose:
                sys.stdout.write('\r  * Notching Event {0}/{1}'.format(ix + 1, data.shape[0]))
                sys.stdout.flush()

            if i + 1 == segs_to_notch:
                seg_et = tmax
            else:
                seg_et = seg_st + (seg_dur // 2)

            f_low = power_law(data[ix, :], seg_st, f0=f0)
            f_high = power_law(data[ix, :], seg_et, f0=f0)

            if f_low < f_min:
                f_low = f_min
            if f_high > f_max:
                f_high = f_max

            ft_data = smart_notching(ft_data, seg_st, f_low, f_high,
                                     seg_dur  = seg_dur,
                                     df       = df,
                                     f_min    = f_min,
                                     f_max    = f_max,
                                     st_time  = job_start_time,
                                     forward  = forward,
                                     backward = backward)

            seg_st += (seg_dur // 2)

    # Write data to file
    if verbose:
        sys.stdout.write('\r\n  * Writing and saving ccspectra file... ')
        sys.stdout.flush()

    ccspec_jobfile = np.loadtxt(event_file, skiprows=3, delimiter='\t')
    cc_real = ccspec_jobfile[:, -2] * ft_data
    cc_imag = ccspec_jobfile[:, -1] * ft_data
    out_file = outdir + '/' + event_file.split('/')[-1]
    out_file = out_file[:-3] + 'mat'
    sio.savemat(out_file, {'t0'     : ccspec_jobfile[:, 0],
                           'sigma'  : ccspec_jobfile[:, 1],
                           'freq'   : ccspec_jobfile[:, 2],
                           'cc_real': cc_real,
                           'cc_imag': cc_imag})

    if verbose:
        sys.stdout.write('Done\n')
        sys.stdout.write('\r  * Writing and saving sensint file... ')
        sys.stdout.flush()

    sensint_jobfile = np.loadtxt(sensint_files[ii], skiprows=3, delimiter='\t')
    sens = sensint_jobfile[:, -1] * ft_data
    out_file = outdir + '/' + sensint_files[ii].split('/')[-1]
    out_file = out_file[:-3] + 'mat'
    sio.savemat(out_file, {'t0'     : sensint_jobfile[:, 0],
                          'sigma'   : sensint_jobfile[:, 1],
                          'freq'    : sensint_jobfile[:, 2],
                          'sens_int': sens})

if verbose:
    sys.stdout.write('Done\n')
    print('Total run time: {:.2f}s'.format(time.time() - start_time))
