#!/usr/bin/env python3
import sys
sys.path.append('method')
import os
import numpy as np
import matplotlib.pyplot as plt
import heka_reader as hr
import estimate_iv

cached = '--cached' in sys.argv

savedir = 'fig'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

f1 = 'data-rev/cho-herg/hergcho-staircaseramp-A13-after-sweep1.csv'
ft = 'data/cho-cell/herg25oc1-staircaseramp-times.csv'

t = np.loadtxt(ft, delimiter=',', skiprows=1)
i1 = np.loadtxt(f1, delimiter=',', skiprows=1)

v = np.loadtxt('data/protocol-staircaseramp.csv', delimiter=',',
               skiprows=1)[::2, 1]

iv_i1s, iv_v1s, iv_tau1s = [], [], []
liv_i1s = []
f1s = 'data-rev/cho-herg/selected-hergcho.txt'
selected = []
with open(f1s, 'r') as f:
    for l in f:
        if not l.startswith('#'):
            selected.append(l.split()[0])
for s in selected:
    f = 'data-rev/cho-herg/hergcho-staircaseramp-%s-after-sweep1.csv' % s
    i = np.loadtxt(f, delimiter=',', skiprows=1)
    if cached:
        iv_i1s.append(np.loadtxt('out/auto-1-%s/i_s.txt' % s))
        iv_v1s.append(np.loadtxt('out/auto-1-%s/v_s.txt' % s))
        iv_tau1s.append(np.loadtxt('out/auto-1-%s/tau_s.txt' % s))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t, out='auto-1-%s' % s)
        iv_i1s.append(ivi)
        iv_v1s.append(ivv)
        iv_tau1s.append(ivt)

iv_i2s, iv_v2s, iv_tau2s = [], [], []
liv_i2s = []
f2s = 'data-rev/cho-herg-2/selected-hergcho2.txt'
selected = []
with open(f2s, 'r') as f:
    for l in f:
        if not l.startswith('#'):
            selected.append(l.split()[0])
for s in selected:
    f = 'data-rev/cho-herg-2/hergcho2-staircaseramp-%s-after-sweep1.csv' % s
    i = np.loadtxt(f, delimiter=',', skiprows=1)
    if cached:
        iv_i2s.append(np.loadtxt('out/auto-2-%s/i_s.txt' % s))
        iv_v2s.append(np.loadtxt('out/auto-2-%s/v_s.txt' % s))
        iv_tau2s.append(np.loadtxt('out/auto-2-%s/tau_s.txt' % s))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t, out='auto-2-%s' % s)
        iv_i2s.append(ivi)
        iv_v2s.append(ivv)
        iv_tau2s.append(ivt)

iv_i3s, iv_v3s, iv_tau3s = [], [], []
liv_i3s = []
selected = [
    [0, 1, 0],
    [1, 1, 0],
    [3, 3, 0],
    [4, 3, 0],
]
for s in selected:
    i = hr.load('data-rev/cho-empty/CHO_I_20201113.dat', s)[0]
    if cached:
        iv_i3s.append(np.loadtxt('out/auto-3-%s/i_s.txt' % s[0]))
        iv_v3s.append(np.loadtxt('out/auto-3-%s/v_s.txt' % s[0]))
        iv_tau3s.append(np.loadtxt('out/auto-3-%s/tau_s.txt' % s[0]))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t, out='auto-3-%s' % s[0])
        iv_i3s.append(ivi)
        iv_v3s.append(ivv)
        iv_tau3s.append(ivt)

iv_i4s, iv_v4s, iv_tau4s = [], [], []
liv_i4s = []
f4s = 'data-rev/cho-empty-auto/selected-emptycho.txt'
selected = []
with open(f4s, 'r') as f:
    for l in f:
        if not l.startswith('#'):
            selected.append(l.split()[0])
for s in selected:
    f = 'data-rev/cho-empty-auto/emptycho-staircaseramp-%s-before-sweep1.csv' % s
    i = np.loadtxt(f, delimiter=',', skiprows=1)
    if cached:
        iv_i4s.append(np.loadtxt('out/auto-4-%s/i_s.txt' % s))
        iv_v4s.append(np.loadtxt('out/auto-4-%s/v_s.txt' % s))
        iv_tau4s.append(np.loadtxt('out/auto-4-%s/tau_s.txt' % s))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t, out='auto-4-%s' % s)
        iv_i4s.append(ivi)
        iv_v4s.append(ivv)
        iv_tau4s.append(ivt)


# Plot
fig, axes = plt.subplots(1, 2, figsize=(8, 2.75))

def plot_example(current, times, ax):
    from scipy.optimize import curve_fit

    def exp_func(t, a, b, c):
        # do a 'proper exponential' decay fit
        # i.e. shift the t to t' where t' has zero at the start of the 
        # voltage step
        return a * np.exp( -b * (t - t_start)) + c

    t_start = 0.9  # s
    t_step = 1.0  # s
    t_end = t_start + t_step
    t_trim = 5e-3  # s
    t_fit_until = 1. - 5e-3  # s

    time_window = np.where(np.logical_and(times > t_start, times <= t_end))[0]
    i_trim = np.argmin(np.abs(times - (t_start + t_trim))) - time_window[0]
    i_fit_until = np.argmin(np.abs(times - (t_start + t_fit_until))) \
                  - time_window[0]
    # Trim off the first i_trim (100ms) in case it is still shooting up...
    x = times[time_window[0] + i_trim:time_window[0] + i_fit_until]
    y = current[time_window[0] + i_trim:
                time_window[0] + i_fit_until]
    # Compute init guess
    current_i, current_f = np.mean(y[0:50]), np.mean(y[-50:-1])
    a = current_i - current_f
    b = 1. / 500e-3
    c = current_f
    p0 = [a, b, c]

    # Fit
    popt, pcov = curve_fit(exp_func, x, y, p0=p0)

    # Plot
    ax.plot(times[time_window[0] - 500:time_window[-1] + 500],
            current[time_window[0] - 500:time_window[-1] + 500],
            c='C0')
    fitted = exp_func(times[time_window[0] + i_trim:
                            time_window[-1] + 500], *popt)
    ax.plot(times[time_window[0] + i_trim:time_window[-1] + 500], fitted,
            '--', c='C3')
    ax.axhline(y=popt[2], ls='--', c='C2')
    ax.axvline(x=times[time_window[0] + i_trim], ls=':', c='#7f7f7f')
    ax.axvline(x=times[time_window[0] + i_fit_until], ls=':', c='#7f7f7f')
    ax.set_xlim([times[time_window[0] - 500], times[time_window[-1] + 500 - 1]])
    ax.set_ylim([-450, 400])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Current (pA)')

    return ax

axes[0] = plot_example(i1, t, axes[0])

bins = np.linspace(50, 1250, 20)
def hist(values, label, c, ax):
    import matplotlib
    alpha = 0.3
    rbg = matplotlib.colors.to_rgb(c)
    kwargs = dict(bins=bins,
                  histtype='stepfilled',
                  edgecolor=c, ls='--', lw=1.2,
                  fc=rbg + (alpha,))
    ax.hist(values, label=label, **kwargs)
    return ax

hist(np.array(iv_tau1s)[:, 0], label='(A)', c='C0', ax=axes[1])
hist(np.array(iv_tau2s)[:, 0], label='(B)', c='C1', ax=axes[1])
hist(np.array(iv_tau3s)[:, 0], label='(C)', c='C2', ax=axes[1])
hist(np.array(iv_tau4s)[:, 0], label='(D)', c='C3', ax=axes[1])
axes[1].set_xlim([50, 1250])
axes[1].legend()
axes[1].set_xlabel('Time constant (ms)')
axes[1].set_ylabel('Frequency')

# plt.subplots_adjust(top=0.975, bottom=0.1, right=0.975, left=0.075)
plt.savefig('%s/figure-3.png' % savedir, bbox_inches='tight', pad_inches=0,
        dpi=200)
plt.savefig('%s/figure-3.pdf' % savedir, bbox_inches='tight', pad_inches=0,
        format='pdf')

print('Done')