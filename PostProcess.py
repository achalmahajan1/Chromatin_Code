import numpy as np
from matplotlib import rc
import math
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True) #This is included for latex type setting
rc('axes', linewidth=2)#Change the width of axis
import matplotlib.pyplot as plt
plt.gcf().subplots_adjust(bottom=0.16,left=0.17)
rad1 = [3, 6, 8, 12, 16, 20, 24];
print, (rad1)
Rs = 28.0;
rad = [x/Rs for x in rad1];
## Nematic alignment data
avg_kc003_sig0_2 = [0.0774, 0.0268, 0.0149, 0.0101, 0.0073, 0.0057, 0.0044];
avg_kc003_sig20_2 = [0.0960, 0.0428, 0.0257, 0.0170, 0.0119, 0.0090, 0.0074];
avg_kc003_sig30_2 = [0.1191, 0.0642, 0.0425, 0.0320, 0.0239, 0.0180, 0.0133];
avg_kc003_sig50_2 = [0.1564, 0.0988, 0.0762, 0.0618, 0.0448, 0.0302, 0.0216];

avg_kc003_sig0_2Cr = [0.0833, 0.0431, 0.0272, 0.0184, 0.0126, 0.0094, 0.0073];
avg_kc003_sig20_2Cr = [0.0857, 0.0439, 0.0264, 0.0183, 0.0132, 0.0106, 0.0088];
avg_kc003_sig30_2Cr = [0.0852, 0.0444, 0.0275, 0.0195, 0.0136, 0.0104, 0.0085];
avg_kc003_sig50_2Cr = [0.0916, 0.0513, 0.0367, 0.0268, 0.0185, 0.0137, 0.0113];

avg_kc003_sig0_2NCr = [0.0984, 0.0323, 0.0174, 0.0116, 0.0083, 0.0064, 0.0051];
avg_kc003_sig20_2NCr = [0.1300, 0.0554, 0.0329, 0.0211, 0.0145, 0.0109, 0.0089];
avg_kc003_sig30_2NCr = [0.1692, 0.0898, 0.0588, 0.0441, 0.0334, 0.0255, 0.0191];
avg_kc003_sig50_2NCr = [0.2168, 0.1354, 0.1040, 0.0846, 0.0617, 0.0416, 0.0295];

stdev_kc003_sig0_2 = [0.0023, 0.0011, 9.7567e-04, 9.2790e-04, 8.7281e-04, 8.1348e-04, 8.0596e-04];
stdev_kc003_sig20_2 = [0.0030, 0.0023, 0.0018, 0.0015, 0.0011, 9.7519e-04, 0.0011];
stdev_kc003_sig30_2 = [0.0040, 0.0023, 0.0021, 0.0016, 0.0011, 9.0009e-04, 9.0009e-04];
stdev_kc003_sig50_2 = [0.0132, 0.0100, 0.0084, 0.0062, 0.0046, 0.0036, 0.0026];

stdev_kc003_sig0_2Cr = [0.0062, 0.0022, 0.0018, 0.0016, 0.0015, 0.0014, 0.0014];
stdev_kc003_sig20_2Cr = [0.0055, 0.0020, 0.0017, 0.0018, 0.0018, 0.0020, 0.0021];
stdev_kc003_sig30_2Cr = [0.0045, 0.0025, 0.0020, 0.0018, 0.0019, 0.0019, 0.0020];
stdev_kc003_sig50_2Cr = [0.0042, 0.0044, 0.0050, 0.0038, 0.0031, 0.0028, 0.0029];

stdev_kc003_sig0_2NCr = [0.0030, 0.0013, 0.0011, 9.9404e-04, 9.7304e-04, 9.1753e-04, 9.2736e-04];
stdev_kc003_sig20_2NCr = [0.0038, 0.0029, 0.0024, 0.0018, 0.0014, 0.0012, 0.0013];
stdev_kc003_sig30_2NCr = [0.0044, 0.0025, 0.0018, 0.0017, 0.0014, 0.0012, 0.0014];
stdev_kc003_sig50_2NCr = [0.0118, 0.0100, 0.0095, 0.0073, 0.0055, 0.0045, 0.0037];

## Peristence length data
sig = [0, 20, 30, 50];
Lp_av = [0.8490, 0.9973, 1.2972, 1.7998];
Lp_std = [0.005, 0.0060, 0.0204, 0.0204];
Lp_av_Cr = [0.666, 0.666, 0.6836, 0.7437];
Lp_std_Cr = [0.0052, 0.0048, 0.0076, 0.0063];
Lp_av_NCr = [0.982, 1.1655, 1.7077, 2.8171];
Lp_std_NCr = [0.0051, 0.0163, 0.0644, 0.1924];

## Mean Square displacement data
#a1 = np.loadtxt('MSD_N1KNc23_PhaseSeg1b_18.dat');
#a2 = np.loadtxt('MSD_N1KNc23_PhaseSeg3b_18_NCr.dat');
#a3 = np.loadtxt('MSD_N1KNc23_PhaseSeg3b_18.dat');
#a4 = np.loadtxt('MSD_N1KNc23_PhaseSeg3b_18_Cr.dat');
a1 = np.loadtxt('MSD_N1KNc23_PhaseSeg1b_18.dat');
a2 = np.loadtxt('MSD_N1KNc23_PhaseSeg2b_18.dat');
a3 = np.loadtxt('MSD_N1KNc23_PhaseSeg3b_18.dat');
a4 = np.loadtxt('MSD_N1KNc23_PhaseSeg4b_18.dat');
dt = 20;
tend = 36440;
tstart = 0;
MSD0 = a1[tstart:tend:dt,1]
t0 = a1[tstart:tend:dt,0]-10;
tend = 35830;
tstart = 0;
MSD20 = a2[tstart:tend:dt,1]
t20 = a2[tstart:tend:dt,0]-10;
tend = 34540;
tstart = 0;
MSD30 = a3[tstart:tend:dt,1]
t30 = a3[tstart:tend:dt,0]-10;
tend = 31710;
tstart = 0;
MSD50 = a4[tstart:tend:dt,1]
t50 = a4[tstart:tend:dt,0]-10;
tend = 12000;
tstart = 800;

t0a = a1[tstart:tend:dt,0]-10;
t20a = a2[tstart:tend:dt,0]-10;
t30a = a3[tstart:tend:dt,0]-10;
t50a = a4[tstart:tend:dt,0]-10;
## For comparison of with, without crosslinks and full chain
#a = 3.206;
#b = 0.4423;
#c = -3.124;
#MSD0a = a*t0a**b+c;
#a = 1.573;
#b = 1.035;
#c = -7.59;
#MSD20a = a*t20a**b+c;
a = 1.912;
b = 1.041;
c = -5.712;
MSD30a = a*t30a**b+c;
#a = 0.6793;
#b = 1.207;
#c = 1.956;
#MSD50a = a*t50a**b+c;
## For comparison of activity
#a = 3.206;
#b = 0.4423;
#c = -3.124;
#MSD0a = a*t0a**b+c;
#a = 0.4781;
#b = 0.9208;
#c = 1.381;
#MSD20a = a*t20a**b+c;
#a = 1.912;
#b = 1.041;
#c = -5.712;
#MSD30a = a*t30a**b+c;

ax2 = plt.subplot(111)
## Plot for the nematic alignment variation with the levels of activity
#ax2.errorbar(rad,avg_kc003_sig0_2,yerr=stdev_kc003_sig0_2,fmt='.b',ms=14, solid_capstyle='projecting', capsize=2,label=r"$\sigma_0 = 0$")
#ax2.errorbar(rad,avg_kc003_sig20_2,yerr=stdev_kc003_sig20_2,fmt='.r',ms=14, solid_capstyle='projecting', capsize=2,label=r"$\sigma_0 = 20$")
#ax2.errorbar(rad,avg_kc003_sig30_2,yerr=stdev_kc003_sig30_2,fmt='.k',ms=14, solid_capstyle='projecting', capsize=2,label=r"$\sigma_0 = 30$")
#ax2.errorbar(rad,avg_kc003_sig50_2,yerr=stdev_kc003_sig50_2, fmt='.g',ms=14, solid_capstyle='projecting', capsize=2,label=r"$\sigma_0 = 50$")

## Plot for the nematic alignment variation with the different domains (Full chain; heterochromatin; euchromatin)
#ax2.errorbar(rad,avg_kc003_sig30_2NCr,yerr=stdev_kc003_sig0_2,fmt='.b',ms=14, solid_capstyle='projecting', capsize=2,label=r"Non-crosslinks")
#ax2.errorbar(rad,avg_kc003_sig30_2,yerr=stdev_kc003_sig20_2,fmt='.k',ms=14, solid_capstyle='projecting', capsize=2,label=r"Entire chain")
#ax2.errorbar(rad,avg_kc003_sig30_2Cr,yerr=stdev_kc003_sig30_2,fmt='.r',ms=14, solid_capstyle='projecting', capsize=2,label=r"Crosslinks")

## Plot for the persistence length
#ax2.errorbar(sig,Lp_av_NCr,yerr=Lp_std_NCr,fmt='.b',ms=14, solid_capstyle='projecting', capsize=2,label=r"Non-crosslinks")
#ax2.errorbar(sig,Lp_av,yerr=Lp_std,fmt='.k',ms=14, solid_capstyle='projecting', capsize=2,label=r"Entire chain")
#ax2.errorbar(sig,Lp_av_Cr,yerr=Lp_std_Cr,fmt='.r',ms=14, solid_capstyle='projecting', capsize=2,label=r"Crosslinks")
## Plot for mean-square displacement
#ax2.loglog(t0,MSD0,'b-', label=r"$\sigma_0 = 0$", linewidth=2)
#ax2.loglog(t20,MSD20,'r-', label=r"Non-crosslinks", linewidth=2)
#ax2.loglog(t30,MSD30,'k-', label=r"Entire chain", linewidth=2)
#ax2.loglog(t50,MSD50,'g-', label=r"Crosslinks", linewidth=2)
#ax2.loglog(t0a,MSD0a,'b--', label=r"$a+bt^{\alpha}$", linewidth=2)
#ax2.loglog(t30a,MSD30a,'k--', label=r"$a+bt^{\alpha}$", linewidth=2)

## Plot for mean-square displacement
ax2.loglog(t0,MSD0,'b-', label=r"$\sigma_0 = 0$", linewidth=2)
ax2.loglog(t20,MSD20,'r-', label=r"$\sigma_0 = 20$", linewidth=2)
ax2.loglog(t30,MSD30,'k-', label=r"$\sigma_0 = 30$", linewidth=2)
ax2.loglog(t50,MSD50,'g-', label=r"$\sigma_0 = 50$", linewidth=2)
#ax2.loglog(t0a,MSD0a,'b--', label=r"$a+bt^{\alpha}$", linewidth=2)
ax2.loglog(t30a,MSD30a,'k--', label=r"$a+bt^{\alpha}$", linewidth=2)

plt.yticks((0.1, 10, 1000), ('$10^{-1}$', r'$10^{1}$', '$10^{3}$'), color='k', size=18)
plt.xticks((0.01, 1, 100), ('$10^{-2}$', r'$10^{1}$', '$10^{2}$'), color='k', size=18)
ax2.locator_params(axis='y', tight=True)
ax2.locator_params(axis='x', tight=True)
#ax2.set_title(r'\TeX\ is Number $\displaystyle\sum_{n=1}^\infty'r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')
#plt.ylabel(r'$\lambda_{max}$',fontsize=24)
#plt.xlabel(r'$r/R_s$',fontsize=24)
plt.ylabel(r'$\langle r^2(t) \rangle$',fontsize=22)
plt.xlabel(r'$t$',fontsize=22)
plt.setp(ax2.get_yticklabels(), fontsize=24)
plt.setp(ax2.get_xticklabels(), fontsize=24)
legend = plt.legend(framealpha=1, frameon=True, loc='upper left',fontsize=14)
legend.get_frame().set_edgecolor('k')
plt.ylim(0.05, 1000)
plt.xlim(0.01, 500)
#plt.xlim(-0.5, 52)
#plt.ylim(0, 4)
plt.savefig('Fig7b.eps')
plt.show()
