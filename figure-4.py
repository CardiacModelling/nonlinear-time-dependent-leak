#!/usr/bin/env python3
import sys
sys.path.append('method')
import os
import numpy as np
import matplotlib.pyplot as plt
import heka_reader as hr
import estimate_iv
  
cached = '--cached' in sys.argv
if '--seed' in sys.argv:
    seed = int(sys.argv[sys.argv.index('--seed') + 1])
else:
    seed = 0
np.random.seed(seed)

ELEAK = []

def compute_leak(traces, i=[15000, 16250, 18000, 19250]):
  ya = np.mean(traces[0][i[0]:i[1]]); yb = np.mean(traces[0][i[2]:i[3]])
  xa = np.mean(traces[1][i[0]:i[1]]); xb = np.mean(traces[1][i[2]:i[3]])
  m = (ya - yb) / (xa - xb)
  c = ya - m * xa
  ELEAK.append(- c / m)
  return m * traces[1] + c

def compute_leak_iv(traces, i=[15000, 16250, 18000, 19250]):
  ya = np.mean(traces[0][i[0]:i[1]]); yb = np.mean(traces[0][i[2]:i[3]])
  xa = np.mean(traces[1][i[0]:i[1]]); xb = np.mean(traces[1][i[2]:i[3]])
  m = (ya - yb) / (xa - xb)
  c = ya - m * xa
  v = np.linspace(-160, 80, 100)
  return m * v + c, v

i_m2 = [38250, 39000, 40000, 40750]

savedir = 'fig'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

f1 = 'data/silicone/20191010_sg_rsol5_ls_nocomp.dat'
f2 = 'data/silicone/20191010_sg_rrsol_ls_noc3.dat'
f3 = 'data/silicone/20191010_sg_rrsolm_ls_noc.dat'
#f3 = 'data/silicone/20191009_sg_msol2_ls_nocomp.dat'

idx1 = [0, 8, 0]
idx2 = [0, 0, 1]
idx3 = [0, 0, 0]

i1, t1 = hr.load(f1, idx1)
i2, t2 = hr.load(f2, idx2)
i3, t3 = hr.load(f3, idx3)

v = np.loadtxt('data/protocol-staircaseramp.csv', delimiter=',',
        skiprows=1)[::2, 1]

traces1 = [i1, v]
traces2 = [i2, v]
traces3 = [i3, v]


iv_i1s, iv_v1s, iv_tau1s = [], [], []
liv_i1s = []
if cached:
    iv_i1s.append(np.loadtxt('out/manual-1/i_s.txt'))
    iv_v1s.append(np.loadtxt('out/manual-1/v_s.txt'))
    iv_tau1s.append(np.loadtxt('out/manual-1/tau_s.txt'))
else:
    ivi, ivv, ivt = estimate_iv.get_iv(i1, v, t1, out='manual-1')
    iv_i1s.append(ivi)
    iv_v1s.append(ivv)
    iv_tau1s.append(ivt)
liv_i1s.append(compute_leak([i1, v]))
selected = [
    [1, 7, 1],
    [2, 1, 0],
]
for s in selected:
    i = hr.load('data-rev/silicone/Sylguard_20201020.dat', s)[0]
    if cached:
        iv_i1s.append(np.loadtxt('out/manual-1-1-%s/i_s.txt' % s[0]))
        iv_v1s.append(np.loadtxt('out/manual-1-1-%s/v_s.txt' % s[0]))
        iv_tau1s.append(np.loadtxt('out/manual-1-1-%s/tau_s.txt' % s[0]))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t, out='manual-1-1-%s' % s[0])
        iv_i1s.append(ivi)
        iv_v1s.append(ivv)
        iv_tau1s.append(ivt)
    liv_i1s.append(compute_leak([i, v]))
selected = [
    [0, 3, 0],
    [1, 1, 0],
    [2, 1, 0],
    [3, 2, 0],
    [4, 4, 0],
    [5, 1, 0],
    [6, 5, 0],
    [7, 1, 0],
    [8, 5, 0],
    [10, 3, 1],
]
for s in selected:
    i = hr.load('data-rev/silicone/Sylguard_20201021.dat', s)[0]
    if cached:
        iv_i1s.append(np.loadtxt('out/manual-1-2-%s/i_s.txt' % s[0]))
        iv_v1s.append(np.loadtxt('out/manual-1-2-%s/v_s.txt' % s[0]))
        iv_tau1s.append(np.loadtxt('out/manual-1-2-%s/tau_s.txt' % s[0]))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t, out='manual-1-2-%s' % s[0])
        iv_i1s.append(ivi)
        iv_v1s.append(ivv)
        iv_tau1s.append(ivt)
    liv_i1s.append(compute_leak([i, v]))


iv_i2s, iv_v2s, iv_tau2s = [], [], []
liv_i2s = []
if cached:
    iv_i2s.append(np.loadtxt('out/manual-2/i_s.txt'))
    iv_v2s.append(np.loadtxt('out/manual-2/v_s.txt'))
    iv_tau2s.append(np.loadtxt('out/manual-2/tau_s.txt'))
else:
    ivi, ivv, ivt = estimate_iv.get_iv(i2, v, t2, out='manual-2')
    iv_i2s.append(ivi)
    iv_v2s.append(ivv)
    iv_tau2s.append(ivt)
liv_i2s.append(compute_leak([i2, v], i=i_m2))
selected = [
    [4, 12, 1],
    [6, 9, 1],
    [7, 4, 1],
]
for s in selected:
    i = hr.load('data-rev/silicone/Sylguard_II_III_20201109.dat', s)[0]
    if cached:
        iv_i2s.append(np.loadtxt('out/manual-2-%s/i_s.txt' % s[0]))
        iv_v2s.append(np.loadtxt('out/manual-2-%s/v_s.txt' % s[0]))
        iv_tau2s.append(np.loadtxt('out/manual-2-%s/tau_s.txt' % s[0]))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t2, out='manual-2-%s' % s[0])
        iv_i2s.append(ivi)
        iv_v2s.append(ivv)
        iv_tau2s.append(ivt)
    liv_i2s.append(compute_leak([i, v], i=i_m2))


iv_i3s, iv_v3s, iv_tau3s = [], [], []
liv_i3s = []
if cached:
    iv_i3s.append(np.loadtxt('out/manual-3/i_s.txt'))
    iv_v3s.append(np.loadtxt('out/manual-3/v_s.txt'))
    iv_tau3s.append(np.loadtxt('out/manual-3/tau_s.txt'))
else:
    ivi, ivv, ivt = estimate_iv.get_iv(i3, v, t3, out='manual-3')
    iv_i3s.append(ivi)
    iv_v3s.append(ivv)
    iv_tau3s.append(ivt)
liv_i3s.append(compute_leak([i3, v]))
selected = [
    [4, 17, 0],  # or 1 TODO!
    [6, 14, 0],  # Nice!
    [7, 6, 0],
]
for s in selected:
    i = hr.load('data-rev/silicone/Sylguard_II_III_20201109.dat', s)[0]
    if cached:
        iv_i3s.append(np.loadtxt('out/manual-3-%s/i_s.txt' % s[0]))
        iv_v3s.append(np.loadtxt('out/manual-3-%s/v_s.txt' % s[0]))
        iv_tau3s.append(np.loadtxt('out/manual-3-%s/tau_s.txt' % s[0]))
    else:
        ivi, ivv, ivt = estimate_iv.get_iv(i, v, t3, out='manual-3-%s' % s[0])
        iv_i3s.append(ivi)
        iv_v3s.append(ivv)
        iv_tau3s.append(ivt)
    liv_i3s.append(compute_leak([i, v]))


# Plot
fig = plt.figure(figsize=(8, 4))
grid = plt.GridSpec(7, 5, hspace=0.125, wspace=0.275)
axes = np.empty([4, 2], dtype=object)

# Voltage
axes[0, 0] = fig.add_subplot(grid[:1, :3])
axes[0, 0].plot(t1, v, c='#7f7f7f')
axes[0, 0].set_xlim([t1[0], t1[-1]])
axes[0, 0].set_xticks([])

# Current
axes[1, 0] = fig.add_subplot(grid[1:3, :3])
axes[1, 0].plot(t1, i1, c='C0')
axes[1, 0].plot(t1, compute_leak(traces1), c='C1', ls='--')
axes[1, 0].set_ylim([-300, 400])
axes[1, 0].set_xlim([t1[0], t1[-1]])
axes[1, 0].set_xticks([])
axes[1, 0].text(-0.2, 0.85, '(I)', transform=axes[1, 0].transAxes, size=12,
        weight='bold')

axes[2, 0] = fig.add_subplot(grid[3:5, :3])
axes[2, 0].plot(t2, i2, c='C0')
axes[2, 0].plot(t2, compute_leak(traces2, i=i_m2),
                c='C1', ls='--')
axes[2, 0].set_ylim([-400, 210])
axes[2, 0].set_xlim([t2[0], t2[-1]])
axes[2, 0].set_xticks([])
axes[2, 0].text(-0.2, 0.85, '(II)', transform=axes[2, 0].transAxes, size=12,
        weight='bold')

if True:
    i3 = hr.load('data-rev/silicone/Sylguard_II_III_20201109.dat', [6, 14, 0])[0]
    traces3 = [i3, v]

axes[3, 0] = fig.add_subplot(grid[5:7, :3])
axes[3, 0].plot(t3, i3, c='C0')
axes[3, 0].plot(t3, compute_leak(traces3), c='C1', ls='--')
#axes[3, 0].set_ylim([-1700, 800])
axes[3, 0].set_ylim([-400, 210])
axes[3, 0].set_xlim([t3[0], t3[-1]])
axes[3, 0].text(-0.2, 0.85, '(III)', transform=axes[3, 0].transAxes, size=12,
        weight='bold')

axes[0, 0].set_ylabel('Voltage\n(mV)', rotation=90)
axes[1, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[2, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[3, 0].set_ylabel('Current\n(pA)', rotation=90)
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

axes[2, 1] = fig.add_subplot(grid[3:5, 3:])
niv_i2s = []
for iv_v, iv_i, li in zip(iv_v2s, iv_i2s, liv_i2s):
    niv_i2s.append(normalise_by(iv_i, li))
boxplot(np.array(niv_i2s), np.round(iv_v2s[0]), axes[2, 1])
iv_i, iv_v = compute_leak_iv([liv_i2s[0], v], i=i_m2)
iv_i = normalise_by(iv_i, liv_i2s[0])
axes[2, 1].plot(iv_v, iv_i, c='C1', ls='--')

axes[3, 1] = fig.add_subplot(grid[5:7, 3:])
niv_i3s = []
for iv_v, iv_i, li in zip(iv_v3s, iv_i3s, liv_i3s):
    niv_i3s.append(normalise_by(iv_i, li))
boxplot(np.array(niv_i3s), np.round(iv_v3s[0]), axes[3, 1])
iv_i, iv_v = compute_leak_iv([liv_i3s[0], v])
iv_i = normalise_by(iv_i, liv_i3s[0])
axes[3, 1].plot(iv_v, iv_i, c='C1', ls='--')

for i in range(1, 4):
    axes[i, 1].set_ylim([-2.9, 3.7])
    axes[i, 1].set_yticks([-2, 0, 2])
    axes[i, 1].set_yticklabels([r'$-2$', r'$0$', r'$2$'])
    axes[i, 1].grid(lw=1, alpha=0.4)
    axes[i, 1].set_xlim([-160, 80])
    axes[i, 1].set_xticks([-120, -80, -40, 0, 40])

axes[1, 1].set_xticklabels([])
axes[2, 1].set_xticklabels([])
axes[3, 1].set_xticklabels([r'$-120$', r'$-80$', r'$-40$', r'$0$', r'$40$'])

axes[-1, 1].set_xlabel('Voltage (mV)', fontsize=12)

fig.align_ylabels(axes[:, 0])

plt.subplots_adjust(top=0.975, bottom=0.1, right=0.975, left=0.075)
plt.savefig('%s/figure-4.png' % savedir, bbox_inches='tight', pad_inches=0,
        dpi=200)
plt.savefig('%s/figure-4.pdf' % savedir, bbox_inches='tight', pad_inches=0,
        format='pdf')

# print('ELeak =', np.mean(ELEAK), '+/-', np.std(ELEAK) / np.sqrt(len(ELEAK)))
np.savetxt('out/eleak-fig4.txt', ELEAK)

print('Done')
