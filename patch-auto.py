#!/usr/bin/env python3
import sys
sys.path.append('method')
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import estimate_iv
  
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


savedir = 'fig'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

f1 = 'data/cho-cell/herg25oc1-staircaseramp-A13-after.csv'
f2 = 'data/cho-cell/mutantsrun2-staircaseramp-A01-after.csv'
f3 = 'data/no-cell/nocell-staircaseramp-B01-after.csv'
ft = 'data/cho-cell/herg25oc1-staircaseramp-times.csv'

t = np.loadtxt(ft, delimiter=',', skiprows=1)
i1 = np.loadtxt(f1, delimiter=',', skiprows=1)
i2 = np.loadtxt(f2, delimiter=',', skiprows=1)
i3 = np.loadtxt(f3, delimiter=',', skiprows=1)

v = np.loadtxt('data/protocol-staircaseramp.csv', delimiter=',',
        skiprows=1)[::2, 1]

traces1 = [i1, v]
traces2 = [i2, v]
traces3 = [i3, v]

iv_i1, iv_v1, iv_tau1 = estimate_iv.get_iv(i1, v, t, out='auto-1')
iv_i2, iv_v2, iv_tau2 = estimate_iv.get_iv(i2, v, t, out='auto-2')
iv_i3, iv_v3, iv_tau3 = estimate_iv.get_iv(i3, v, t, out='auto-3')

# Plot
fig = plt.figure(figsize=(8, 4))
grid = plt.GridSpec(7, 4, hspace=0.125, wspace=0.5)
axes = np.empty([4, 2], dtype=object)

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
axes[1, 0].text(-0.125, 0.85, '(A)', transform=axes[1, 0].transAxes, size=12,
        weight='bold')

axes[2, 0] = fig.add_subplot(grid[3:5, :3])
axes[2, 0].plot(t, i2, c='C0')
axes[2, 0].plot(t, compute_leak(traces2), c='C1', ls='--')
axes[2, 0].set_ylim([-420, 460])
axes[2, 0].set_xlim([t[0], t[-1]])
axes[2, 0].set_xticks([])
axes[2, 0].text(-0.125, 0.85, '(B)', transform=axes[2, 0].transAxes, size=12,
        weight='bold')

axes[3, 0] = fig.add_subplot(grid[5:7, :3])
axes[3, 0].plot(t, i3, c='C0')
axes[3, 0].plot(t, compute_leak(traces3), c='C1', ls='--')
axes[3, 0].set_ylim([-6000, 4500])
axes[3, 0].set_xlim([t[0], t[-1]])
axes[3, 0].text(-0.125, 0.85, '(C)', transform=axes[3, 0].transAxes, size=12,
        weight='bold')

axes[0, 0].set_ylabel('Voltage\n(mV)', rotation=90)
axes[1, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[2, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[3, 0].set_ylabel('Current\n(pA)', rotation=90)
axes[-1, 0].set_xlabel('Time (s)', fontsize=12)

# IV
#axes[0, 1].axis('off')

axes[1, 1] = fig.add_subplot(grid[1:3, 3])
axes[1, 1].plot(iv_v1, iv_i1, 'x', c='C0')
iv_i, iv_v = compute_leak_iv(traces1)
axes[1, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[1, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[1, 1].set_xticks([])

axes[2, 1] = fig.add_subplot(grid[3:5, 3])
axes[2, 1].plot(iv_v2, iv_i2, 'x', c='C0')
iv_i, iv_v = compute_leak_iv(traces2)
axes[2, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[2, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[2, 1].set_xticks([])

axes[3, 1] = fig.add_subplot(grid[5:7, 3])
axes[3, 1].plot(iv_v3, iv_i3, 'x', c='C0')
iv_i, iv_v = compute_leak_iv(traces3)
axes[3, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[3, 1].set_xlim([iv_v[0], iv_v[-1]])

axes[-1, 1].set_xlabel('Voltage (mV)', fontsize=12)

plt.subplots_adjust(top=0.975, bottom=0.1, right=0.975, left=0.075)
plt.savefig('%s/patch-auto.png' % savedir, bbox_inches='tight', pad_inches=0,
        dpi=200)

print('Done')
