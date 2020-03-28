#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import heka_reader
import estimate_iv

def load(f, i, vccc=False):
    # f: file path
    # i: index list, [group_ind, series_ind, sweep_ind]
    # vccc: bool, if True, read 4 traces, else 2.
    #
    # return
    # ====
    # out: Recorded traces
    # times: Recorded traces' times
    b = heka_reader.Bundle(f)
    if vccc:
        nt = 4
    else:
        nt = 2
    out = []
    for t in range(nt):
        out.append(b.data[i + [t]])
    info = b.pul[i[0]][i[1]][i[2]][0]
    times = np.linspace(info.XStart,
            info.XStart + info.XInterval * (len(out[0])-1),
            len(out[0]))
    return out, times

  
def compute_leak(traces, i=[30000, 32500, 36000, 38500]):
  ya = np.mean(traces[0][i[0]:i[1]]); yb = np.mean(traces[0][i[2]:i[3]])
  xa = np.mean(traces[1][i[0]:i[1]]); xb = np.mean(traces[1][i[2]:i[3]])
  m = (ya - yb) / (xa - xb)
  c = ya - m * xa
  return m * traces[1] + c

def compute_leak_iv(traces, i=[30000, 32500, 36000, 38500]):
  ya = np.mean(traces[0][i[0]:i[1]]); yb = np.mean(traces[0][i[2]:i[3]])
  xa = np.mean(traces[1][i[0]:i[1]]); xb = np.mean(traces[1][i[2]:i[3]])
  m = (ya - yb) / (xa - xb)
  c = ya - m * xa
  v = np.linspace(-160e-3, 80e-3, 100)
  return m * v + c, v


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

traces1, times1 = load(f1, idx1)
traces2, times2 = load(f2, idx2)
traces3, times3 = load(f3, idx3)

i1 = traces1[0] * 1e12
i2 = traces2[0] * 1e12
i3 = traces3[0] * 1e12

v = np.loadtxt('protocol-staircaseramp.csv', delimiter=',', skiprows=1)[:, 1]

i1 = i1[::2]
t1 = times1[::2]
i2 = i2[::2]
t2 = times2[::2]
i3 = i3[::2]
t3 = times3[::2]
v = v[::2]

iv_i1, iv_v1, iv_tau1 = estimate_iv.get_iv(i1, v, t1, out='manual-1')
iv_i2, iv_v2, iv_tau2 = estimate_iv.get_iv(i2, v, t2, out='manual-2')
iv_i3, iv_v3, iv_tau3 = estimate_iv.get_iv(i3, v, t3, out='manual-3')

# Plot
fig = plt.figure(figsize=(8, 4))
grid = plt.GridSpec(7, 4, hspace=0.125, wspace=0.5)
axes = np.empty([4, 2], dtype=object)

# Voltage
axes[0, 0] = fig.add_subplot(grid[:1, :3])
axes[0, 0].plot(t1, v, c='#7f7f7f')
axes[0, 0].set_xlim([t1[0], t1[-1]])
axes[0, 0].set_xticks([])

# Current
axes[1, 0] = fig.add_subplot(grid[1:3, :3])
axes[1, 0].plot(times1, traces1[0] * 1e12, c='C0')
axes[1, 0].plot(times1, compute_leak(traces1) * 1e12, c='C1', ls='--')
axes[1, 0].set_ylim([-300, 400])
axes[1, 0].set_xlim([t1[0], t1[-1]])
axes[1, 0].set_xticks([])

axes[2, 0] = fig.add_subplot(grid[3:5, :3])
axes[2, 0].plot(times2, traces2[0] * 1e12, c='C0')
axes[2, 0].plot(times2, compute_leak(traces2,
    i=[76500, 78000, 80000, 81500]) * 1e12, c='C1', ls='--')
axes[2, 0].set_ylim([-400, 210])
axes[2, 0].set_xlim([t2[0], t2[-1]])
axes[2, 0].set_xticks([])

axes[3, 0] = fig.add_subplot(grid[5:7, :3])
axes[3, 0].plot(times3, traces3[0] * 1e12, c='C0')
axes[3, 0].plot(times3, compute_leak(traces3) * 1e12, c='C1', ls='--')
axes[3, 0].set_ylim([-1700, 800])
#axes[3, 0].set_ylim([-1500, 600])
axes[3, 0].set_xlim([t3[0], t3[-1]])

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
iv_i *= 1e12
iv_v *= 1e3
axes[1, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[1, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[1, 1].set_xticks([])

axes[2, 1] = fig.add_subplot(grid[3:5, 3])
axes[2, 1].plot(iv_v2, iv_i2, 'x', c='C0')
iv_i, iv_v = compute_leak_iv(traces2, i=[76500, 78000, 80000, 81500])
iv_i *= 1e12
iv_v *= 1e3
axes[2, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[2, 1].set_xlim([iv_v[0], iv_v[-1]])
axes[2, 1].set_xticks([])

axes[3, 1] = fig.add_subplot(grid[5:7, 3])
axes[3, 1].plot(iv_v3, iv_i3, 'x', c='C0')
iv_i, iv_v = compute_leak_iv(traces3)
iv_i *= 1e12
iv_v *= 1e3
axes[3, 1].plot(iv_v, iv_i, c='C1', ls='--')
axes[3, 1].set_xlim([iv_v[0], iv_v[-1]])

axes[-1, 1].set_xlabel('Voltage (mV)', fontsize=12)

plt.subplots_adjust(top=0.975, bottom=0.1, right=0.975, left=0.075)
plt.savefig('%s/patch-manual.png' % savedir, bbox_inches='tight', pad_inches=0,
        dpi=200)

print('Done')
