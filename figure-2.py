#!/usr/bin/env python3
import sys
sys.path.append('method')
import os
import numpy as np
import matplotlib.pyplot as plt
import estimate_iv
  
cached = '--cached' in sys.argv

def compute_leak(traces, i=[15000, 16250, 18000, 19250]):
  ya = np.mean(traces[0][i[0]:i[1]]); yb = np.mean(traces[0][i[2]:i[3]])
  xa = np.mean(traces[1][i[0]:i[1]]); xb = np.mean(traces[1][i[2]:i[3]])
  m = (ya - yb) / (xa - xb)
  c = ya - m * xa
  return m * traces[1] + c

def compute_leak_iv(traces, i=[15000, 16250, 18000, 19250]):
  ya = np.mean(traces[0][i[0]:i[1]]); yb = np.mean(traces[0][i[2]:i[3]])
  xa = np.mean(traces[1][i[0]:i[1]]); xb = np.mean(traces[1][i[2]:i[3]])
  m = (ya - yb) / (xa - xb)
  c = ya - m * xa
  v = np.linspace(-160, 80, 100)
  return m * v + c, v

def load(f, i, ds=2):
    # f: file path
    # i: index list, [group_ind, series_ind, sweep_ind]
    # ds: downsample
    #
    # return
    # ====
    # out: Recorded trace
    # times: Recorded trace's times
    import heka_reader
    b = heka_reader.Bundle(f)
    nt = 0  # Usually 0 is current
    out = b.data[i + [nt]] * 1e12  # A -> pA
    info = b.pul[i[0]][i[1]][i[2]][0]
    times = np.linspace(info.XStart,
            info.XStart + info.XInterval * (len(out)-1),
            len(out))
    return out[::ds], times[::ds]


savedir = 'fig'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

f1 = 'data-rev/cho-herg/hergcho-staircaseramp-A13-after-sweep1.csv'
f2 = 'data-rev/cho-herg-2/hergcho2-staircaseramp-A02-after-sweep1.csv'
f3 = 'data-rev/cho-empty/CHO_I_20201113.dat'
f4 = 'data-rev/cho-empty-auto/emptycho-staircaseramp-A01-before-sweep1.csv'
f5 = 'data/no-cell/nocell-staircaseramp-B01-after.csv'
ft = 'data/cho-cell/herg25oc1-staircaseramp-times.csv'

t = np.loadtxt(ft, delimiter=',', skiprows=1)
i1 = np.loadtxt(f1, delimiter=',', skiprows=1)
i2 = np.loadtxt(f2, delimiter=',', skiprows=1)
i3 = load(f3, [4, 3, 0])[0]
i4 = np.loadtxt(f4, delimiter=',', skiprows=1)
i5 = np.loadtxt(f5, delimiter=',', skiprows=1)

v = np.loadtxt('data/protocol-staircaseramp.csv', delimiter=',',
        skiprows=1)[::2, 1]

traces1 = [i1, v]
traces2 = [i2, v]
traces3 = [i3, v]
traces4 = [i4, v]
traces5 = [i5, v]

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
    liv_i1s.append(compute_leak([i, v]))

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
    liv_i2s.append(compute_leak([i, v]))

iv_i3s, iv_v3s, iv_tau3s = [], [], []
liv_i3s = []
selected = [
    [0, 1, 0],
    [1, 1, 0],
    [3, 3, 0],
    [4, 3, 0],
]
for s in selected:
    i = load(f3, s)[0]
    if cached:
        iv_i3s.append(np.loadtxt('out/auto-3-%s/i_s.txt' % s[0]))
        iv_v3s.append(np.loadtxt('out/auto-3-%s/v_s.txt' % s[0]))
        iv_tau3s.append(np.loadtxt('out/auto-3-%s/tau_s.txt' % s[0]))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t, out='auto-3-%s' % s[0])
        iv_i3s.append(ivi)
        iv_v3s.append(ivv)
        iv_tau3s.append(ivt)
    liv_i3s.append(compute_leak([i, v]))

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
    liv_i4s.append(compute_leak([i, v]))

if cached:
    iv_i5 = np.loadtxt('out/auto-5/i_s.txt')
    iv_v5 = np.loadtxt('out/auto-5/v_s.txt')
    iv_tau5 = np.loadtxt('out/auto-5/tau_s.txt')
else:
    iv_i5, iv_v5, iv_tau5 = estimate_iv.get_iv(i5, v, t, out='auto-5')

# Plot
fig = plt.figure(figsize=(8, 6))
grid = plt.GridSpec(11, 5, hspace=0.125, wspace=0.2)
axes = np.empty([6, 2], dtype=object)

# Voltage
axes[0, 0] = fig.add_subplot(grid[:1, :3])
axes[0, 0].plot(t, v, c='#7f7f7f')
axes[0, 0].set_xlim([t[0], t[-1]])
axes[0, 0].set_xticks([])

# Current
axes[1, 0] = fig.add_subplot(grid[1:3, :3])
axes[1, 0].plot(t, i1, c='C0')
axes[1, 0].plot(t, compute_leak(traces1), c='C1', ls='--')
axes[1, 0].set_ylim([-340, 440])
axes[1, 0].set_xlim([t[0], t[-1]])
axes[1, 0].set_xticks([])
axes[1, 0].text(-0.2, 0.85, '(A)', transform=axes[1, 0].transAxes, size=12,
        weight='bold')

axes[2, 0] = fig.add_subplot(grid[3:5, :3])
axes[2, 0].plot(t, i2, c='C0')
axes[2, 0].plot(t, compute_leak(traces2), c='C1', ls='--')
axes[2, 0].set_ylim([-420, 460])
axes[2, 0].set_xlim([t[0], t[-1]])
axes[2, 0].set_xticks([])
axes[2, 0].text(-0.2, 0.85, '(B)', transform=axes[2, 0].transAxes, size=12,
        weight='bold')

axes[3, 0] = fig.add_subplot(grid[5:7, :3])
axes[3, 0].plot(t, i3, c='C0')
axes[3, 0].plot(t, compute_leak(traces3), c='C1', ls='--')
axes[3, 0].set_ylim([-550, 600])
axes[3, 0].set_xlim([t[0], t[-1]])
axes[3, 0].set_xticks([])
axes[3, 0].text(-0.2, 0.85, '(C)', transform=axes[3, 0].transAxes, size=12,
        weight='bold')

axes[4, 0] = fig.add_subplot(grid[7:9, :3])
axes[4, 0].plot(t, i4, c='C0')
axes[4, 0].plot(t, compute_leak(traces4), c='C1', ls='--')
axes[4, 0].set_ylim([-420, 460])
axes[4, 0].set_xlim([t[0], t[-1]])
axes[4, 0].set_xticks([])
axes[4, 0].text(-0.2, 0.85, '(D)', transform=axes[4, 0].transAxes, size=12,
        weight='bold')

axes[5, 0] = fig.add_subplot(grid[9:11, :3])
axes[5, 0].plot(t, i5, c='C0')
axes[5, 0].plot(t, compute_leak(traces5), c='C1', ls='--')
axes[5, 0].set_ylim([-6000, 4500])
axes[5, 0].set_xlim([t[0], t[-1]])
axes[5, 0].text(-0.2, 0.85, '(E)', transform=axes[5, 0].transAxes, size=12,
        weight='bold')


axes[0, 0].set_ylabel('Voltage\n(mV)', rotation=90)
axes[1, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[2, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[3, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[4, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[5, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[-1, 0].set_xlabel('Time (s)', fontsize=12)

# IV
#axes[0, 1].axis('off')

def boxplot(ivis, ivvs, ax):
    x = [-120, -80, -60, -40, -20, 0, 20, 40]
    for xx in x:
        i = np.where(ivvs == xx)[0]
        y = ivis[:, i]
        bp = ax.boxplot(y.reshape(-1, 1), positions=[xx], showfliers=False,
                        widths=6)
        for e in ['whiskers', 'caps']:
            plt.setp(bp[e], color='#7f7f7f')
        for e in ['boxes', 'fliers', 'means', 'medians']:
            plt.setp(bp[e], color='C0')

def normalise(x):
    xr = np.max(x) - np.min(x)
    return (x - np.min(x)) / xr

def normalise_by(x, y):
    yr = np.max(y) - np.min(y)
    return (x - np.min(y)) / yr

axes[1, 1] = fig.add_subplot(grid[1:3, 3:])
axes[1, 1].set_title('Normalised current')
niv_i1s = []
for iv_v, iv_i, li in zip(iv_v1s, iv_i1s, liv_i1s):
    niv_i1s.append(normalise_by(iv_i, li))
boxplot(np.array(niv_i1s), np.round(iv_v1s[0]), axes[1, 1])
iv_i, iv_v = compute_leak_iv([liv_i1s[0], v])
iv_i = normalise_by(iv_i, liv_i1s[0])
axes[1, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[1, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[1, 1].set_xticks([])

axes[2, 1] = fig.add_subplot(grid[3:5, 3:])
niv_i2s = []
for iv_v, iv_i, li in zip(iv_v2s, iv_i2s, liv_i2s):
    niv_i2s.append(normalise_by(iv_i, li))
boxplot(np.array(niv_i2s), np.round(iv_v2s[0]), axes[2, 1])
iv_i, iv_v = compute_leak_iv([liv_i2s[0], v])
iv_i = normalise_by(iv_i, liv_i2s[0])
axes[2, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[2, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[2, 1].set_xticks([])

axes[3, 1] = fig.add_subplot(grid[5:7, 3:])
niv_i3s = []
for iv_v, iv_i, li in zip(iv_v3s, iv_i3s, liv_i3s):
    niv_i3s.append(normalise_by(iv_i, li))
boxplot(np.array(niv_i3s), np.round(iv_v3s[0]), axes[3, 1])
iv_i, iv_v = compute_leak_iv([liv_i3s[0], v])
iv_i = normalise_by(iv_i, liv_i3s[0])
axes[3, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[3, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[3, 1].set_xticks([])

axes[4, 1] = fig.add_subplot(grid[7:9, 3:])
niv_i4s = []
for iv_v, iv_i, li in zip(iv_v4s, iv_i4s, liv_i4s):
    niv_i4s.append(normalise_by(iv_i, li))
boxplot(np.array(niv_i4s), np.round(iv_v4s[0]), axes[4, 1])
iv_i, iv_v = compute_leak_iv([liv_i4s[0], v])
iv_i = normalise_by(iv_i, liv_i4s[0])
axes[4, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[4, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[4, 1].set_xticks([])

axes[5, 1] = fig.add_subplot(grid[9:11, 3:])
niv_i5 = normalise_by(iv_i5, compute_leak(traces5))
axes[5, 1].plot(np.round(iv_v5), niv_i5, 'x', c='C0')
iv_i, iv_v = compute_leak_iv(traces5)
iv_i = normalise_by(iv_i, compute_leak(traces5))
axes[5, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[5, 1].set_xlim([iv_v[0], iv_v[-1]])

axes[-1, 1].set_xlabel('Voltage (mV)', fontsize=12)

fig.align_ylabels(axes[:, 0])

plt.subplots_adjust(top=0.975, bottom=0.1, right=0.975, left=0.075)
plt.savefig('%s/figure-2.png' % savedir, bbox_inches='tight', pad_inches=0,
        dpi=200)
plt.savefig('%s/figure-2.pdf' % savedir, bbox_inches='tight', pad_inches=0,
        format='pdf')

print('Done')
