#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import heka_reader

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


savedir = 'fig/test'
if not os.path.isdir(savedir):
    os.makedirs(savedir)

f1 = 'data/silicone/20191010_sg_rsol5_ls_nocomp.dat'
f2 = 'data/silicone/20191010_sg_rrsol_ls_noc3.dat'
#f3 = 'data/silicone/20191010_sg_rrsolm_ls_noc.dat'
f3 = 'data/silicone/20191009_sg_msol2_ls_nocomp.dat'

i1 = [0, 8, 0]
i2 = [0, 0, 1]
i3 = [0, 0, 0]

traces1, times1 = load(f1, i1)
traces2, times2 = load(f2, i2)
traces3, times3 = load(f3, i3)

plt.figure(figsize=(8, 4))
plt.plot(times1, traces1[0] * 1e12)
plt.plot(times1, compute_leak(traces1) * 1e12)
plt.ylim([-500, 500])
plt.ylabel('Current (pA)')
plt.xlabel('Time (s)')
plt.savefig('%s/Roche_solution.png' % savedir, dpi=200)
plt.close()

plt.figure(figsize=(8, 4))
plt.plot(times2, traces2[0] * 1e12)
plt.plot(times2, compute_leak(traces2, i=[76500, 78000, 80000, 81500]) * 1e12)
plt.ylim([-600, 400])
plt.ylabel('Current (pA)')
plt.xlabel('Time (s)')
plt.savefig('%s/Roche_solution_insideout.png' % savedir, dpi=200)
plt.close()

plt.figure(figsize=(8, 4))
plt.plot(times3, traces3[0] * 1e12)
plt.plot(times3, compute_leak(traces3) * 1e12)
plt.ylim([-1500, 600])
plt.ylabel('Current (pA)')
plt.xlabel('Time (s)')
plt.savefig('%s/Manual_solution.png' % savedir, dpi=200)
plt.close()
