import numpy as np
from matplotlib import rc
import math
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True) #This is included for latex type setting
rc('axes', linewidth=1)#Change the width of axis
import matplotlib.pyplot as plt
#plt.gcf().subplots_adjust(bottom=0.15,left=0.15)
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
a1 = np.loadtxt('DensityPhaseSeg1b_18_Sig0_DenIso04_Rad04_Grid05.dat')
b1 = np.loadtxt('DensityPhaseSeg2b_18_Sig20_DenIso04_Rad04_Grid05.dat')
c1 = np.loadtxt('DensityPhaseSeg3b_18_Sig30_DenIso04_Rad04_Grid05.dat')
d1 = np.loadtxt('DensityPhaseSeg4b_18_Sig50_DenIso04_Rad04_Grid05.dat')
a2 = np.loadtxt('DensityPhaseSeg1b_18a_Sig0_DenIso04_Rad04_Grid05.dat')
b2 = np.loadtxt('DensityPhaseSeg2b_18a_Sig20_DenIso04_Rad04_Grid05.dat')
c2 = np.loadtxt('DensityPhaseSeg3b_18a_Sig30_DenIso04_Rad04_Grid05.dat')
d2 = np.loadtxt('DensityPhaseSeg4b_18a_Sig50_DenIso04_Rad04_Grid05.dat')
a3 = np.loadtxt('DensityPhaseSeg1b_18b_Sig0_DenIso04_Rad04_Grid05.dat')
b3 = np.loadtxt('DensityPhaseSeg2b_18b_Sig20_DenIso04_Rad04_Grid05.dat')
c3 = np.loadtxt('DensityPhaseSeg3b_18b_Sig30_DenIso04_Rad04_Grid05.dat')
d3 = np.loadtxt('DensityPhaseSeg4b_18b_Sig50_DenIso04_Rad04_Grid05.dat')
a4 = np.loadtxt('DensityPhaseSeg1b_18c_Sig0_DenIso04_Rad04_Grid05.dat')
b4 = np.loadtxt('DensityPhaseSeg2b_18c_Sig20_DenIso04_Rad04_Grid05.dat')
c4 = np.loadtxt('DensityPhaseSeg3b_18c_Sig30_DenIso04_Rad04_Grid05.dat')
d4 = np.loadtxt('DensityPhaseSeg4b_18c_Sig50_DenIso04_Rad04_Grid05.dat')
d5 = np.loadtxt('DensityPhaseSeg4b_18e_Sig50_DenIso04_Rad04_Grid05.dat')
d6 = np.loadtxt('DensityPhaseSeg4b_18f_Sig50_DenIso04_Rad04_Grid05.dat')
d7 = np.loadtxt('DensityPhaseSeg4b_18g_Sig50_DenIso04_Rad04_Grid05.dat')

#a1 = np.loadtxt('RadialDensPhaseSeg1b_14_Sig0NoCm.dat');
#a2 = np.loadtxt('RadialDensPhaseSeg1b_14a_Sig0NoCm.dat');
#a3 = np.loadtxt('RadialDensPhaseSeg1b_14b_Sig0NoCm.dat');
#a4 = np.loadtxt('RadialDensPhaseSeg1b_14c_Sig0NoCm.dat');
#a1_1 = np.loadtxt('RadialDensPhaseSeg1b_14_Sig0NoCm_NCr.dat');
#a2_1 = np.loadtxt('RadialDensPhaseSeg1b_14a_Sig0NoCm_NCr.dat');
#a3_1 = np.loadtxt('RadialDensPhaseSeg1b_14b_Sig0NoCm_NCr.dat');
#a4_1 = np.loadtxt('RadialDensPhaseSeg1b_14c_Sig0NoCm_NCr.dat');
#b1 = np.loadtxt('RadialDensPhaseSeg2b_14_Sig20NoCm.dat');
#b2 = np.loadtxt('RadialDensPhaseSeg2b_14a_Sig20NoCm.dat');
#b3 = np.loadtxt('RadialDensPhaseSeg2b_14b_Sig20NoCm.dat');
#b4 = np.loadtxt('RadialDensPhaseSeg2b_14c_Sig20NoCm.dat');
#b1_1 = np.loadtxt('RadialDensPhaseSeg2b_14_Sig20NoCm_NCr.dat');
#b2_1 = np.loadtxt('RadialDensPhaseSeg2b_14a_Sig20NoCm_NCr.dat');
#b3_1 = np.loadtxt('RadialDensPhaseSeg2b_14b_Sig20NoCm_NCr.dat');
#b4_1 = np.loadtxt('RadialDensPhaseSeg2b_14c_Sig20NoCm_NCr.dat');
#c1 = np.loadtxt('RadialDensPhaseSeg3b_14_Sig30NoCm.dat');
#c2 = np.loadtxt('RadialDensPhaseSeg3b_14a_Sig30NoCm.dat');
#c3 = np.loadtxt('RadialDensPhaseSeg3b_14b_Sig30NoCm.dat');
#c4 = np.loadtxt('RadialDensPhaseSeg3b_14c_Sig30NoCm.dat');
#c1_1 = np.loadtxt('RadialDensPhaseSeg3b_14_Sig30NoCm_NCr.dat');
#c2_1 = np.loadtxt('RadialDensPhaseSeg3b_14a_Sig30NoCm_NCr.dat');
#c3_1 = np.loadtxt('RadialDensPhaseSeg3b_14b_Sig30NoCm_NCr.dat');
#c4_1 = np.loadtxt('RadialDensPhaseSeg3b_14c_Sig30NoCm_NCr.dat');
#d1 = np.loadtxt('RadialDensPhaseSeg4b_14_Sig50NoCm.dat');
#d2 = np.loadtxt('RadialDensPhaseSeg4b_14a_Sig50NoCm.dat');
#d3 = np.loadtxt('RadialDensPhaseSeg4b_14b_Sig50NoCm.dat');
#d4 = np.loadtxt('RadialDensPhaseSeg4b_14c_Sig50NoCm.dat');
#d1_1 = np.loadtxt('RadialDensPhaseSeg4b_14_Sig50NoCm_NCr.dat');
#d2_1 = np.loadtxt('RadialDensPhaseSeg4b_14a_Sig50NoCm_NCr.dat');
#d3_1 = np.loadtxt('RadialDensPhaseSeg4b_14b_Sig50NoCm_NCr.dat');
#d4_1 = np.loadtxt('RadialDensPhaseSeg4b_14c_Sig50NoCm_NCr.dat');
#e1 = np.loadtxt('RadialDensPhaseSeg3b_20_Sig30N_NoCm.dat');
#e1_1 = np.loadtxt('RadialDensPhaseSeg3b_20_Sig30N_NoCm_NCr.dat');

# pi = math.atan2(1, 1);
Rs = 28.0;
ts = 12;
Vs = (4*math.pi*(28**3))/3;
print(math.pi)
Va = 4*math.pi*(28**2);
Nb = 1305*23;
a4_temp = a4[:,0:10]
b4_temp = b4[:,0:10]
c3_temp = c3[:,0:10]
c4_temp = c4[:,0:10]
d4_temp = d4[:,0:10]
d5_temp = d5[:,0:10]
d6_temp = d6[:,0:10]
d7_temp = d7[:,0:10]
print(a4_temp.size)
w, h = 10, a1[:,1].size+a2[:,1].size+a3[:,1].size+a4[:,1].size;
a = [[0 for x in range(w)] for y in range(h)]
w, h = 10, b1[:,1].size+b2[:,1].size+b3[:,1].size+b4[:,1].size;
b = [[0 for x in range(w)] for y in range(h)]
#w, h = 10, c1[:,1].size+c2[:,1].size+c3_temp[:,1].size+c4_temp[:,1].size;
w, h = 10, c1[:,1].size+c2[:,1].size+c3[:,1].size+c4[:,1].size;
c = [[0 for x in range(w)] for y in range(h)]
w, h = 10, d1[:,1].size+d2[:,1].size+d3[:,1].size+d4[:,1].size+d5[:,1].size;#+d6[:,1].size+d7[:,1].size;
d = [[0 for x in range(w)] for y in range(h)]

#w, h = 10, e1[:,1].size;
#e = [[0 for x in range(w)] for y in range(h)]
a = np.concatenate((a1, a2, a3, a4_temp))
b = np.concatenate((b1, b2, b3, b4_temp))
c = np.concatenate((c1, c2, c3_temp, c4_temp))
d = np.concatenate((d1, d2, d3, d4_temp, d5_temp))#, d6_temp, d7_temp))
#e = e1

#w, h = 10, a1_1[:,1].size+a2_1[:,1].size+a3_1[:,1].size+a4_1[:,1].size;
#a_1 = [[0 for x in range(w)] for y in range(h)]
#w, h = 10, b1_1[:,1].size+b2_1[:,1].size+b3_1[:,1].size+b4_1[:,1].size;
#b_1 = [[0 for x in range(w)] for y in range(h)]
#w, h = 10, c1_1[:,1].size+c2_1[:,1].size+c3_1[:,1].size+c4_1[:,1].size;
#c_1 = [[0 for x in range(w)] for y in range(h)]
#w, h = 10, d1_1[:,1].size+d2_1[:,1].size+d3_1[:,1].size+d4_1[:,1].size;
#d_1 = [[0 for x in range(w)] for y in range(h)]
#w, h = 10, e1_1[:,1].size;
#e_1 = [[0 for x in range(w)] for y in range(h)]
#a_1 = np.concatenate((a1_1, a2_1, a3_1, a4_1))
#b_1 = np.concatenate((b1_1, b2_1, b3_1, b4_1))
#c_1 = np.concatenate((c1_1, c2_1, c3_1, c4_1))
#d_1 = np.concatenate((d1_1, d2_1, d3_1, d4_1))
#e_1 = e1_1
#anet = a[a[:,1].size-1,1]-a_1[a_1[:,1].size-1,1];
#bnet = b[b[:,1].size-1,1]-b_1[b_1[:,1].size-1,1];
#cnet = c[c[:,1].size-1,1]-c_1[c_1[:,1].size-1,1];
#dnet = d[d[:,1].size-1,1]-d_1[d_1[:,1].size-1,1];
#enet = e[e[:,1].size-1,1]-e_1[e_1[:,1].size-1,1];
#net = [abs(enet), abs(anet), abs(bnet), abs(cnet), abs(dnet)];
#sig = [-30, 0, 20, 30, 50];
#print net
#print sig
#DenSa = (Nb-a[1:a.size,5])/(Vs-a[1:a.size,4]);
#DenSb = (Nb-b[1:b.size,5])/(Vs-b[1:b.size,4]);
#DenSc = (Nb-c[1:c.size,5])/(Vs-c[1:c.size,4]);
#DenSd = (Nb-d[1:d.size,5])/(Vs-d[1:d.size,4]);

## Combining Particle number density and volume fraction on a single plot
#fig, ax1 = plt.subplots()
#ax1.plot(a[ts:a.size,0],a[ts:a.size,4]/a[ts:a.size,3],'b-', label=r"$\sigma_0 = 0$", linewidth=2)
#ax1.plot(b[ts:b.size,0],b[ts:b.size,4]/b[ts:b.size,3],'r-', label=r"$\sigma_0 = 20$", linewidth=2)
#ax1.plot(c[ts:c.size,0],c[ts:c.size,4]/c[ts:c.size,3],'k-', label=r"$\sigma_0 = 30$", linewidth=2)
#ax1.plot(d[ts:d.size,0],d[ts:d.size,4]/d[ts:d.size,3],'g-', label=r"$\sigma_0 = 50$", linewidth=2)
#ax1.locator_params(axis='y', tight=True, nbins=4)
#ax1.locator_params(axis='x', tight=True, nbins=4)
#ax1.set_ylabel(r'$\rho_c$',fontsize=24)
#ax1.set_xlabel(r'$t$',fontsize=24)

#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# Volume fraction plot
#ax2.plot(a[ts:a.size,0],a[ts:a.size,3]/Vs,'b-', label=r"$\sigma_0 = 0$", linewidth=2)
#ax2.plot(b[ts:b.size,0],b[ts:b.size,3]/Vs,'r-', label=r"$\sigma_0 = 20$", linewidth=2)
#ax2.plot(c[ts:c.size,0],c[ts:c.size,3]/Vs,'k-', label=r"$\sigma_0 = 30$", linewidth=2)
#ax2.plot(d[ts:d.size,0],d[ts:d.size,3]/Vs,'g-', label=r"$\sigma_0 = 50$", linewidth=2)
#ax2.locator_params(axis='y', tight=True, nbins=4)
#ax2.locator_params(axis='x', tight=True, nbins=4)
#ax2.set_ylabel(r'$V_{cl}/V$',fontsize=24)
#ax2.setp(ax2.get_yticklabels(), fontsize=24)
#ax2.legend(loc='upper left',fontsize=20)
#legend = plt.legend(framealpha=1, frameon=True, loc='upper left',fontsize=16)
#legend.get_frame().set_edgecolor('k')
#plt.xlim(0, 520)
#plt.savefig('Fig2.eps')
#plt.show()

## Adding inset of volume fraction and particle fraction in the density plot


fig, ax1 = plt.subplots()
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)
ax1.plot(a[ts:a.size,0],a[ts:a.size,4]/a[ts:a.size,3],'k-', label=r"$\sigma_0 = 0$", linewidth=2)
ax1.plot(b[ts:b.size,0],b[ts:b.size,4]/b[ts:b.size,3],'b-', label=r"$\sigma_0 = 20$", linewidth=2)
ax1.plot(c[ts:c.size,0],c[ts:c.size,4]/c[ts:c.size,3],'g-', label=r"$\sigma_0 = 30$", linewidth=2)
ax1.plot(d[ts:d.size,0],d[ts:d.size,4]/d[ts:d.size,3],'r-', label=r"$\sigma_0 = 50$", linewidth=2)
ax1.set_ylabel(r'$\rho_{hcr}$',fontsize=22)
ax1.set_xlabel(r'$t$',fontsize=22)
ax1.set_ylim(0.38, 0.9)
ax1.set_xlim(0, 520)
ax1.locator_params(axis='y', tight=True, nbins=4)
ax1.locator_params(axis='x', tight=True, nbins=4)
legend = ax1.legend(framealpha=1, frameon=False,loc='best', bbox_to_anchor=(0., 0.01, 0.55, 0.55),fontsize=8)
#legend.get_frame().set_edgecolor('k')
plt.setp(ax1.get_yticklabels(), fontsize=20)
plt.setp(ax1.get_xticklabels(), fontsize=20)

ax2 = plt.axes([0,0,0,0])
#plt.gcf().subplots_adjust(bottom=0.15,left=0.15)
ip = InsetPosition(ax1, [0.62,0.6,0.35,0.35])
ax2.set_axes_locator(ip)

ax2.plot(a[ts:a.size,0],a[ts:a.size,3]/Vs,'k-', linewidth=1)
ax2.plot(b[ts:b.size,0],b[ts:b.size,3]/Vs,'b-', linewidth=1)
ax2.plot(c[ts:c.size,0],c[ts:c.size,3]/Vs,'g-', linewidth=1)
ax2.plot(d[ts:d.size,0],d[ts:d.size,3]/Vs,'r-', linewidth=1)

ax2.locator_params(axis='y', tight=True, nbins=4)
# ax2.locator_params(axis='x', tight=True, nbins=4)
ax2.set_ylabel(r'$V_{hcr}/V$',fontsize=12)
#ax2.set_xlabel(r'$t$',fontsize=12)
plt.setp(ax2.get_yticklabels(), fontsize=12)
#plt.setp(ax2.get_xticklabels(), fontsize=12)
ax2.set_xticklabels([]) ## TO KEEP ONLY THE TICKS
#ax2.axes.get_xaxis().set_visible(False) ## TO REMOVE THE TICKS AND LABELS
ax2.set_ylim(0, 0.401)
ax2.set_xlim(0, 520)

ax3 = plt.axes([0,0,0,1])
#plt.gcf().subplots_adjust(bottom=0.15,left=0.15)
ip = InsetPosition(ax1, [0.14,0.6,0.35,0.35])
ax3.set_axes_locator(ip)

ax3.plot(a[ts:a.size,0],a[ts:a.size,4]/Nb,'k-', label=r"$\sigma_0 = 0$", linewidth=1)
ax3.plot(b[ts:b.size,0],b[ts:b.size,4]/Nb,'b-', label=r"$\sigma_0 = 20$", linewidth=1)
ax3.plot(c[ts:c.size,0],c[ts:c.size,4]/Nb,'g-', label=r"$\sigma_0 = 30$", linewidth=1)
ax3.plot(d[ts:d.size,0],d[ts:d.size,4]/Nb,'r-', label=r"$\sigma_0 = 50$", linewidth=1)

plt.locator_params(axis='y', tight=True, nbins=4)
#plt.locator_params(axis='x', tight=True, nbins=4)
plt.ylabel(r'$N_{hcr}/N$',fontsize=12)
#plt.xlabel(r'$t$',fontsize=12)
plt.setp(ax3.get_yticklabels(), fontsize=12)
#plt.setp(ax3.get_xticklabels(), fontsize=12)
#ax3.legend(framealpha=1, frameon=True, loc='lower left',fontsize=8)
ax3.set_xticklabels([]) ## TO KEEP ONLY THE TICKS
#ax3.axes.get_xaxis().set_visible(False)
ax3.set_ylim(0, 0.6)
ax3.set_xlim(0, 520)


plt.savefig('Fig3b.eps')
plt.show()
