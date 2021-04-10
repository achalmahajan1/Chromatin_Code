import numpy as np
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True) #This is included for latex type setting
rc('axes', linewidth=2)#Change the width of axis
import matplotlib.pyplot as plt
plt.gcf().subplots_adjust(bottom=0.15,left=0.15)


Rs = 28
#a = np.loadtxt('Cdr_dt5_E50N23K_003.dat')
#a1 = np.loadtxt('Cdr_dt5_E50N23K_003_Cr6.dat')
#a2 = np.loadtxt('Cdr_dt5_E50N23K_003_Cr6a.dat')
#a3 = np.loadtxt('Cdr_dt5_E50N23K_003_Cr7.dat')
#a4 = np.loadtxt('Cdr_dt5_E50N23K_003_Mo_Cr7a.dat')
#a5 = np.loadtxt('Cdr_dt5_E50N23K_003_Mo_Cr7b.dat')

a1 = np.loadtxt('Cdr_dt2_E0_Phase1b_15.dat')
a2 = np.loadtxt('Cdr_dt2_E0_Phase1b_15.dat')
a3 = np.loadtxt('Cdr_dt05_E0_Phase1b_15.dat')
a4 = np.loadtxt('Cdr_dt1_E0_Phase1b_15.dat')
a5 = np.loadtxt('Cdr_dt2_E0_Phase1b_15.dat')
a6 = np.loadtxt('Cdr_dt5_E0_Phase1b_15.dat')

#a1 = np.loadtxt('Cdr_dt5_E50N23K_003.dat')
#a2 = np.loadtxt('Cdr_dt5_E0N23K_003.dat')
#a3 = np.loadtxt('Cdr_dt5_C50N23K_003_Cr3.dat')
#a4 = np.loadtxt('Cdr_dt5_E50N23K_003_rb.dat')
#a5 = np.loadtxt('Cdr_dt5_E0N23K_003.dat')

#a1 = np.loadtxt('Cdr_dt2_E0N23K_003.dat')
#a2 = np.loadtxt('Cdr_dt2_E50N23K_003_Cr3b.dat')
#a3 = np.loadtxt('Cdr_dt2_E50N23K_003.dat')
#a4 = np.loadtxt('Cdr_dt2_E50N23K_001_Cr3c.dat')
#a5 = np.loadtxt('Cdr_dt2_E50N23K_0008_Cr3d.dat')
#a = np.loadtxt('Cdr_dt001_E50N23K_003_Cr3_Rb001.dat')

#a1 = np.loadtxt('Cdr_dt5_E50N23K_Cr_0.dat')
#a2 = np.loadtxt('Cdr_dt5_E50N23K_003.dat')
#a3 = np.loadtxt('Cdr_dt5_E0N23K_003_Cr3_opprat.dat')
#a4 = np.loadtxt('Cdr_dt5_E0N23K_003_Cr3_opprat1.dat')
#a5 = np.loadtxt('Cdr_dt5_E0N23K_003_Cr3_opprat2.dat')
#Cdr = a[:,1]/a[:,2]
#Cdr = Cdr/Cdr[1]
Cdr001 = a1[:,1]/a1[:,2]
Cdr001 = Cdr001/Cdr001[1]
Cdr01 = a2[:,1]/a2[:,2]
Cdr01 = Cdr01/Cdr01[1]
Cdr05 = a3[:,1]/a3[:,2]
Cdr05 = Cdr05/Cdr05[1]
Cdr1 = a4[:,1]/a4[:,2]
Cdr1 = Cdr1/Cdr1[1]
Cdr2 = a5[:,1]/a5[:,2]
Cdr2 = Cdr2/Cdr2[1]
Cdr5 = a6[:,1]/a6[:,2]
Cdr5 = Cdr5/Cdr5[1]
rad = a1[:,0]/Rs
rad = rad#-rad[2]+0.01

ax2 = plt.subplot(111)
#plt.plot(rad,Cdr001,'b-o', label=r"$\alpha_c  = 0.0, \alpha_t = 0.0$", linewidth=2)
#plt.plot(rad,Cdr02,'r-o', label=r"$\alpha_c  = 0.205, \alpha_t = 0.0$", linewidth=2)
#plt.plot(rad,Cdr05,'g-o', label=r"$\alpha_c  = 0.0, \alpha_t = 0.026$", linewidth=2)
#plt.plot(rad,Cdr1,'c-o', label=r"$\alpha_c  = 0.22, \alpha_t = 0.011$", linewidth=2)
#plt.plot(rad,Cdr5,'k-o', label=r"$\alpha_c  = 0.23, \alpha_t = 0.0285$", linewidth=2)

#plt.plot(rad,Cdr001,'b-o', label=r"$k_{con}/k_{coff} = 0.0/0.0, k_{on}/k_{off} = 200/500$", linewidth=2)
#plt.plot(rad,Cdr02,'r-o', label=r"$k_{con}/k_{coff} = 0.03/0.1, k_{on}/k_{off} = 200/500$", linewidth=2)
#plt.plot(rad,Cdr05,'g-o', label=r"$k_{con}/k_{coff} = 3/10, k_{on}/k_{off} = 2/5$", linewidth=2)
#plt.plot(rad,Cdr1,'c-o', label=r"$k_{con}/k_{coff} = 30/100, k_{on}/k_{off} = 0.2/0.5$", linewidth=2)
#plt.plot(rad,Cdr5,'k-o', label=r"$k_{con}/k_{coff} = 300/1000, k_{on}/k_{off} = 0.02/0.05$", linewidth=2)

#plt.plot(rad,Cdr001,'b-', label=r"$\Delta t = 0.01$", linewidth=2)
#plt.plot(rad,Cdr01,'r-', label=r"$\Delta t = 0.5$", linewidth=2)
plt.plot(rad,Cdr05,'b-o', label=r"$\Delta t = 0.5$", linewidth=2)
plt.plot(rad,Cdr1,'r-o', label=r"$\Delta t = 1.0$", linewidth=2)
plt.plot(rad,Cdr2,'k-o', label=r"$\Delta t = 2.0$", linewidth=2)
plt.plot(rad,Cdr5,'g-o', label=r"$\Delta t = 5.0$", linewidth=2)

#plt.plot(rad,Cdr001,'b-o', label=r"$\sigma_0 = 50$", linewidth=2)
#plt.plot(rad,Cdr02,'r-o', label=r"$\sigma_0 = 0$", linewidth=2)
#plt.plot(rad,Cdr05,'g-o', label=r"$\sigma_0 = -50$", linewidth=2)
#plt.plot(rad,Cdr1,'c-o', label=r"$\Delta t = 1.0$", linewidth=2)
#plt.plot(rad,Cdr5,'k-o', label=r"$\Delta t = 5.0$", linewidth=2)

#plt.plot(rad,Cdr02,'r-o', label=r"$\nu  = 0.0009$", linewidth=2)
#plt.plot(rad,Cdr05,'g-o', label=r"$\nu  = 0.0013$", linewidth=2)
#plt.plot(rad,Cdr1,'c-o', label=r"$\nu  = 0.0021$", linewidth=2)
#plt.plot(rad,Cdr5,'k-o', label=r"$\nu  = 0.003$", linewidth=2)

#plt.plot(rad,Cdr,'r-.', label="$\Delta t  = 0.01(r_h = 0.01)$", linewidth=2)
#plt.plot(rad,Cdr001,'r--', label="$\Delta t  = 0.01 (r_h = 0.05)$", linewidth=2)
#plt.plot(rad,Cdr02,'r-o', label="$\Delta t  = 0.01 (r_h = 0.1)$", linewidth=2)
#plt.plot(rad,Cdr05,'k-.', label="$\Delta t  = 5.0 (r_h = 0.01)$", linewidth=2)
#plt.plot(rad,Cdr1,'k--', label="$\Delta t  = 5.0 (r_h = 0.05)$", linewidth=2)
#plt.plot(rad,Cdr5,'k-o', label="$\Delta t  = 5.0 (r_h = 0.1)$", linewidth=2)

#plt.plot(rad,Cdr001,'g-', label=r"$\alpha_c = 0.203 \; (\sigma_0 = 0)$", linewidth=2)
#plt.plot(rad,Cdr02,'g--', label=r"$\alpha_c = 0.204 \; (\sigma_0 = 50),\textbf{\boldmath$u$} = 0$", linewidth=2)
#plt.plot(rad,Cdr05,'g-o', label=r"$\alpha_c = 0.206 \; (\sigma_0 = 50)$", linewidth=2)
#plt.plot(rad,Cdr1,'r-o', label=r"$\alpha_c = 0.0869 \; (\sigma_0 = 50)$", linewidth=2)
#plt.plot(rad,Cdr5,'k-o', label=r"$\alpha_c = 0.0713 \; (\sigma_0 = 50)$", linewidth=2)

#plt.plot(rad,Cdr,'b-o', label=r"$\leftarrow\kern-8.5pt\hbox{$\circ$} \kern+8.5pt\hbox{$\circ$}\kern-1.5pt\hbox{$\rightarrow$}$", linewidth=2)
#plt.plot(rad,Cdr001,'m-o', label=r"$\leftarrow\kern-8.5pt\hbox{$\bullet$} \kern+8.5pt\hbox{$\bullet$}\kern-1.5pt\hbox{$\rightarrow$}$", linewidth=2)
#plt.plot(rad,Cdr02,'r-o', label=r"$\leftarrow\kern-8.5pt\hbox{$\bullet$} \kern+8.5pt\hbox{$\bullet$}\kern-1.5pt\hbox{$\rightarrow$}$ $(Local)$", linewidth=2)
#plt.plot(rad,Cdr05,'g-o', label=r"$\kern+8.5pt\hbox{$\circ$}\kern-1.5pt\hbox{$\rightarrow$}$", linewidth=2)
#plt.plot(rad,Cdr1,'c-o', label=r"$\kern+8.5pt\hbox{$\bullet$}\kern-1.5pt\hbox{$\rightarrow$}$", linewidth=2)
#plt.plot(rad,Cdr5,'k-o', label=r"$\kern+8.5pt\hbox{$\bullet$}\kern-1.5pt\hbox{$\rightarrow$}$ $(Local)$", linewidth=2)

#plt.setp(ax2.get_xticklabels(), visible=False)
ax2.locator_params(axis='y', tight=True, nbins=4)
ax2.locator_params(axis='x', tight=True, nbins=4)
#ax2.set_title(r'\TeX\ is Number $\displaystyle\sum_{n=1}^\infty'r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')
plt.legend(loc='lower right')
plt.ylabel(r'$C(\Delta r)$',fontsize=26)
plt.xlabel(r'$r/R_s$',fontsize=26)
plt.setp(ax2.get_yticklabels(), fontsize=26)
plt.setp(ax2.get_xticklabels(), fontsize=26)
plt.legend(loc='upper right',fontsize=20)
plt.ylim(-0.01, 1)
plt.xlim(0, 1.15)
plt.savefig('DAC_N1KNc23_Ac0_PhaseSeg1b_15.eps')
#plt.savefig('FigN2KNc23Ac50_003CompareAcN50.eps')
plt.show()
