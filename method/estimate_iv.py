"""
Estimation of the time constant and steady-state I-V relationship for the
staircase protocol.
"""
from __future__ import print_function
import sys
sys.path.append('../lib/')
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

Dt = 500e-3  # s
t_start = 3.4  # s
t_skip = 5e-3  # s
t_skip_2 = 5e-3  # s
n_steps = 19
fit_method = fit_single_exp


def fit_single_exp(current, times,
                       t_start, t_end,
                       t_trim, t_fit_until,
                       debug=None):
    """
    Fitting a 3-parameter single exponential to the tail current.

    Input
    =====
    - current: (array) Current time series
    - times: (array) Time points of `current`
    - t_start: (float) Starting time of the tail current
    - t_end: (float) Ending time of the tail current
    - t_trim: (float) Duration of the beginning of the current not to be fitted
    - t_fit_until: (float) Duration of the time series to be fitted
    - debug: (optional, (str, int)) Output debug plots if not `None`, expect
             the output folder name and the number of step of the tail current

    Output
    =====
    - (i, tau): The steady state current of tail current, and its time constant
    """
    from scipy.optimize import curve_fit
    def exp_func(t, a, b, c):
        # do a "proper exponential" decay fit
        # i.e. shift the t to t' where t' has zero at the start of the 
        # voltage step
        return a * np.exp( -b * (t - t_start)) + c
    time_window = np.where(np.logical_and(times > t_start, times <= t_end))[0]
    i_trim = np.argmin(np.abs(times - (t_start + t_trim))) - time_window[0]
    i_fit_until = np.argmin(np.abs(times - (t_start + t_fit_until))) \
                  - time_window[0]
    # trim off the first i_trim (100ms) in case it is still shooting up...
    x = times[time_window[0] + i_trim:time_window[0] + i_fit_until]
    y = current[time_window[0] + i_trim:
                time_window[0] + i_fit_until]
    # Compute init guess
    current_i, current_f = np.mean(y[0:50]), np.mean(y[-50:-1])
    a = current_i - current_f
    b = 1. / 500e-3
    c = current_f
    p0 = [a, b, c]
    if np.abs(a) < 40:
        if debug is not None:
            debugdir, i_step = debug
            plt.plot(times[time_window[0] - 500:time_window[-1] + 500],
                     current[time_window[0] - 500:time_window[-1] + 500],
                     c='#d62728')
            plt.axhline(y=c, ls='--', c='C2')
            plt.savefig('%s/s%s.png' % (debugdir, i_step))
            plt.close()
        return c, np.NaN
    popt, pcov = curve_fit(exp_func, x, y, p0=p0)
    tau = 1e3 / popt[1]
    out = popt[2]
    if debug is not None:
        debugdir, i_step = debug
        fig = plt.figure()
        plt.plot(times[time_window[0] - 500:time_window[-1] + 500],
                 current[time_window[0] - 500:time_window[-1] + 500],
                 c='#d62728')
        fitted = exp_func(times[time_window[0] + i_trim:
                                time_window[-1] + 500], *popt)
        plt.plot(times[time_window[0] + i_trim:time_window[-1] + 500], fitted,
                 '--', c='#1f77b4')
        plt.axhline(y=out, ls='--', c='C2')
        plt.axvline(x=times[time_window[0] + i_trim])
        plt.axvline(x=times[time_window[0] + i_fit_until])
        plt.savefig('%s/s%s.png' % (debugdir, i_step))
        plt.close()
    return out, tau


def get_iv(i, v, t, out, debug=True):
    """
    Estimate the I-V relationship using the staircase protocol recording.

    Input
    =====
    - i: (array) Recorded current
    - v: (array) Voltage protocol
    - t: (array) Time points of the recording
    - out: (str) Output folder name which will be created in `../out`
    - debug: (optional, bool) Create debug plots if True

    Output
    =====
    - (i_s, v_s, tau_s): current of the IV curve, voltage of the IV curve, and
                         time constant of the tau-V curve
    All results are saved in `../out/[out]` directory.
    """
    outdir = './out/' + out
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if debug == True:
        debugdir = outdir + '/debug'
        if not os.path.isdir(debugdir):
            os.makedirs(debugdir)
    else:
        debugdir = None

    print('Step, tau [ms]')

    plt.plot(t, i)
    plt.xlabel('t [s]')
    plt.ylabel('I [pA]')
    plt.savefig('%s/it' % outdir)
    plt.close()


    i_s = []
    v_s = []
    tau_s = []
    # 3 big steps before staircase
    t_i = 0.9  # s
    for i_step, t_step in enumerate([1.0, 0.5, 1.0]):
        i_step = i_step - 3
        t_f = t_i + t_step
        if i_step == -3:
            i_i, tau_i = fit_method(i, t, t_i, t_f, t_skip, 1. - t_skip_2,
                    debug=(debugdir, i_step))
        else:
            i_i, tau_i = fit_method(i, t, t_i, t_f, t_skip, Dt - t_skip_2,
                    debug=(debugdir, i_step))
        if np.isfinite(tau_i):
            print(i_step, tau_i)
        time_window = np.where(np.logical_and(t > t_i, t <= t_f))[0]
        i_s.append(i_i)
        v_s.append(np.mean(v[time_window]))
        tau_s.append(tau_i)
        t_i = t_f
    # staircase
    for i_step in range(n_steps):
        t_i, t_f = t_start + i_step * Dt, t_start + (i_step + 1) * Dt
        i_i, tau_i = fit_method(i, t, t_i, t_f, t_skip, Dt - t_skip_2,
                debug=(debugdir, i_step))
        if np.isfinite(tau_i):
            print(i_step, tau_i)
        time_window = np.where(np.logical_and(t > t_i, t <= t_f))[0]
        i_s.append(i_i)
        v_s.append(np.mean(v[time_window]))
        tau_s.append(tau_i)
    # last 2 steps after staircase
    t_i = 12.9  # s
    for i_step, t_step in enumerate([1.0, 0.5]):
        i_step = i_step + n_steps
        t_f = t_i + t_step
        i_i, tau_i = fit_method(i, t, t_i, t_f, t_skip, Dt - t_skip_2,
                debug=(debugdir, i_step))
        if np.isfinite(tau_i):
            print(i_step, tau_i)
        time_window = np.where(np.logical_and(t > t_i, t <= t_f))[0]
        i_s.append(i_i)
        v_s.append(np.mean(v[time_window]))
        tau_s.append(tau_i)
        t_i = t_f

    plt.plot(v_s, i_s, 'o-')
    plt.xlabel('V [mV]')
    plt.ylabel('I [pA]')
    plt.savefig('%s/iv' % outdir)
    plt.close()

    np.savetxt('%s/v_s.txt' % outdir, v_s)
    np.savetxt('%s/i_s.txt' % outdir, i_s)
    np.savetxt('%s/tau_s.txt' % outdir, tau_s)

    return i_s, v_s, tau_s

