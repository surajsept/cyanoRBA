import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


def plotLB(y1, y2, y3, ax, model, axes_title):
	ax.set_title(axes_title, fontsize=24, fontweight='bold', loc='left')
	ax.plot(model.p.K_N / y1.Nx, 1 / (86400 * y1.Mu), linewidth=2.0, label='$C_x/K_C = 0.01$')
	ax.plot(model.p.K_N / y2.Nx, 1 / (86400 * y2.Mu), linewidth=2.0, label='$C_x/K_C = 0.05$')
	ax.plot(model.p.K_N / y3.Nx, 1 / (86400 * y3.Mu), linewidth=2.0, label='$C_x/K_C = 1.0$')
	ax.set_xlabel('$K_N/N_x$', fontsize=15)
	ax.set_ylabel('$1/\mu$ (d)', fontsize=15)
	ax.legend(fontsize=12)
	ax.tick_params(labelsize=12)
	plt.tight_layout()


def plotVm_Mu(y, ax, model, axes_title):
	ax.set_title(axes_title, fontsize=24, fontweight='bold', loc='left')
	ax.loglog(y.Nx / model.p.K_N, y.v5, linewidth=2.0, label='uptake flux')
	ax.loglog(y.Nx / model.p.K_N, model.p.k5 * y.Tn, linewidth=2.0, label='$V_{max}$')
	ax.set_xlabel('$N_x/K_N$', fontsize=15)
	ax.set_ylabel('molecules cell$^{-1}$ s$^{-1}$', fontsize=15)
	ax.legend(fontsize=12)
	ax.tick_params(labelsize=12)
	plt.tight_layout()


def twosubplots(y1, y2, y3, model):
	fig, ax = plt.subplots(1, 2, figsize=(10, 4))
	plotVm_Mu(y1, ax[0], model, 'A')
	plotLB(y1, y2, y3, ax[1], model, 'B')
	return plt.show()


def plotTnTy_Mu(X, Y, Y0, Y1, ax, axes_title):
	ax.set_title(axes_title, fontsize=24, fontweight='bold', loc='left')
	ax.semilogx(X, Y, label='$T_N$-strain')
	ax.semilogx(X, Y0, label='$T_Y$-strain')
	ax.semilogx(X, Y1, '--', label='($T_N + T_Y$)-strain')
	ax.set_xlabel('$N_X/K_N$', fontsize=12)
	ax.set_ylabel('Growth rate, $\mu$ ($day^{-1}$)', fontsize=12)
	ax.legend(fontsize=15)
	ax.tick_params(labelsize=12)
	plt.tight_layout()


# return plt.show()

def plotTnTy(X, Y1, Y2, ax, axes_title):
	# fig, ax = plt.subplots(figsize=(5, 4))
	ax.set_title(axes_title, fontsize=24, fontweight='bold', loc='left')
	ax.semilogx(X, Y1, label='$T_N$')
	ax.semilogx(X, Y2, label='$T_Y$')
	ax.set_xlabel('$N_X/K_N$', fontsize=12)
	ax.set_ylabel('relative abundance', fontsize=12)
	ax.legend(fontsize=15)
	ax.tick_params(labelsize=12)
	plt.tight_layout()


# return plt.show()

def subplots_TnTy(y_N, y_Y, y_NY, pF_NY, model):
	fig, ax = plt.subplots(1, 2, figsize=(10, 4))
	plotTnTy_Mu(y_N.Nx / model.p.K_N, 86400*y_N.Mu, 86400*y_Y.Mu, 86400*y_NY.Mu, ax[0], 'A')
	plotTnTy(y_N.Nx / model.p.K_N, pF_NY.T_N2, pF_NY.T_N1, ax[1], 'B')
	return plt.show()


def co_culture(data1, data2, Irr, fsize=15):
	fig, ax = plt.subplots(2, 2, figsize=(10, 5))
	# subplot1
	ax[0, 0].plot(data1.t / 365, data1.nG, 'C0', label='Gleaner', linewidth=2.0)
	ax[0, 0].plot(data1.t / 365, data1.nO, 'C1', label='Opportunist', linewidth=2.0)
	ax[0, 0].set_ylabel('number of cells l$^{-1}$', fontsize=fsize)
	ax[0, 0].tick_params(labelsize=12)
	ax[0, 0].legend(fontsize=12, ncol=2, bbox_to_anchor=(1.02, 1.25))

	# subplot 2
	ax[0, 1].plot((data2.t[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)] / 365)-1.0, data2.nG[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)], 'C0', label='Gleaner', linewidth=2.0)
	ax[0, 1].plot((data2.t[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)] / 365)-1.0, data2.nO[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)], 'C1', label='Opportunist', linewidth=2.0)
	ax[0, 1].tick_params(labelsize=12)
	ax[0, 1].legend(fontsize=12, ncol=2, bbox_to_anchor=(1.02, 1.25))

	# subplot3
	ax[1, 0].plot(data1.t / 365, data1.N, 'C2', label='$N_x$', linewidth=2.0)
	ax[1, 0].set_ylabel('N$_x$ ($\mu$M)', fontsize=fsize)
	ax[1, 0].set_xlabel('years', fontsize=fsize)
	ax[1, 0].tick_params(labelsize=12)
	ax2 = ax[1, 0].twinx()
	ax2.plot(data1.t / 365, 200 * np.ones(len(data1.t)), 'r--', linewidth=2.0)
	ax2.tick_params('y', colors='r', labelsize=12)

	# subplot4
	ax[1, 1].plot((data2.t[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)] / 365)-1.0, data2.N[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)], 'C2', label='$N_x$', linewidth=2.0)
	ax[1, 1].set_xlabel('years', fontsize=fsize)
	ax[1, 1].tick_params(labelsize=12)
	ax2 = ax[1, 1].twinx()
	ax2.plot((data2.t[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)] / 365)-1.0, Irr[findClosestIndex(data2.t, 366):findClosestIndex(data2.t, 4*365)], 'r--', linewidth=2.0)
	ax2.set_ylabel('I ($\mu$ E mu$^{-2}$ s$^{-1}$)', fontsize=15, color='r')
	ax2.tick_params('y', colors='r', labelsize=12)

	plt.tight_layout()
	return plt.show()


def pie(size1, labels, fsize=15, filename='low_Nx', plot='False'):
	fig, ax = plt.subplots(figsize=[4, 4])
	# ax.pie(size1)
	wedges, texts, autotexts = ax.pie(size1, autopct='%1.1f%%')
	ax.axis('equal')
	ax.legend (wedges, labels, ncol=1, fontsize=14,
	            title="Protein fraction",
	            loc="upper left",
	            bbox_to_anchor=(1, 0, 0.5, 1))
	plt.setp(autotexts, weight="bold", color="white", fontsize=15)
	plt.setp (texts, fontsize=12)
	# ax.set_title(title1, fontsize=fsize)
	if not plot == 'False':
		plt.savefig('../Images/cyano_paper/pie_graphs/' + filename + '.png', bbox_inches='tight')
	plt.tight_layout()
	return plt.show()


def piesample(ax, size1):
	wedges, texts, autotexts = ax.pie(size1, autopct='%1.1f%%')
	ax.axis('equal')
	plt.setp(autotexts, weight="bold", color="white", fontsize=15)
	plt.tight_layout()

def findClosestValue(myList, myNumber):
	return min(myList, key=lambda x: abs(x - myNumber))

def findClosestIndex(myList, Value):
	return min(range(len(myList)), key=lambda i: abs(myList[i] - Value))

def pie_subplots(yN, pFN, yC, pFC, yI, pFI):
	fig, ax = plt.subplots(2, 3, figsize=[15, 8])
	piesample(ax[0, 0], pFN[findClosestIndex(yN.Nx, Value=0.22)])
	piesample(ax[1, 0], pFN[findClosestIndex(yN.Nx, Value=10.0)])
	piesample(ax[0, 1], pFC[findClosestIndex(yC.Cx, Value=0.35)])
	piesample(ax[1, 1], pFC[findClosestIndex(yC.Cx, Value=15.0)])
	piesample(ax[0, 2], pFI[findClosestIndex(yI.Irr, Value=20)])
	piesample(ax[1, 2], pFI[findClosestIndex(yI.Irr, Value=400)])
	return plt.show()


def growth_subplots(yN, yC, yI, model, Cx=0.25, Nx=0.5, I=200.0):
	fig, ax = plt.subplots(1, 3, figsize=[12, 4])
	growth_curve(ax[0], yN.Nx / 10.0, 86400 * yN.Mu, xlabel='N$_x$/K$_N$', axes_title='A',
	             ylabel='Growth rate, $\mu$ (d$^{-1}$)',
	             text=['$C_x/K_C$ = %.2f' % round(Cx/model.p.K_C,2), 'I = ' + str(int(I)) + '$\mu E$ m$^{-2}$ s$^{-1}$'])
	growth_curve(ax[1], yC.Cx / 15.0, 86400 * yC.Mu, xlabel='C$_x$/K$_C$', axes_title='B',
	             ylabel='Growth rate, $\mu$ (d$^{-1}$)',
	             text=['$N_x/K_N$ = %.2f' % round(Nx/model.p.K_N,2), 'I = ' + str(int(I)) + '$\mu E$ m$^{-2}$ s$^{-1}$'])
	growth_curve(ax[2], yI.Irr, 86400 * yI.Mu, xlabel='I ($\mu$ E m$^{-2}$ s$^{-1}$)', axes_title='C',
	             ylabel='Growth rate, $\mu$ (d$^{-1}$)',
	             text=['$C_x/K_c$ = %.2f' % round(Cx/model.p.K_C,2), '$N_x/K_n$ = %.2f' % round(Nx/model.p.K_N,2)])
	return plt.show()


def growth_curve(ax, X, Y, xlabel, ylabel, text, axes_title, text_lx=0.3, text_ly=0.1):
	ax.set_title(axes_title, fontsize=20, fontweight='bold')
	ax.plot(X, Y, linewidth=2.0)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.tick_params(labelsize=15)
	plt.tight_layout()
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	textstr = '\n'.join((text))
	# place a text box in upper left in axes coords
	ax.text(text_lx, text_ly, textstr, transform=ax.transAxes, fontsize=14,
	        verticalalignment='bottom', bbox=props)
	
def growth_subplots_new(yN, yC, yI, XN, XC, XI, MuN, MuC, MuI, model, Cx=0.25, Nx=0.5, I=200.0):
	fig, ax = plt.subplots(1, 3, figsize=[12, 4])
	growth_curve_new(ax=ax[0], X1=yN.Nx / 10.0, Y1=86400 * yN.Mu, X2=XN / 10.0, Y2=MuN, xlabel='N$_x$/K$_N$', axes_title='A',
	             ylabel='Growth rate, $\mu$ (d$^{-1}$)',
	             text=['$C_x/K_C$ = %.2f' % round(Cx/model.p.K_C,2), 'I = ' + str(int(I)) + '$\mu E$ m$^{-2}$ s$^{-1}$'])
	growth_curve_new(ax=ax[1], X1=yC.Cx / 15.0, Y1=86400 * yC.Mu, X2=XC / 15.0, Y2=MuC, xlabel='C$_x$/K$_C$', axes_title='B',
	             ylabel='Growth rate, $\mu$ (d$^{-1}$)',
	             text=['$N_x/K_N$ = %.2f' % round(Nx/model.p.K_N,2), 'I = ' + str(int(I)) + '$\mu E$ m$^{-2}$ s$^{-1}$'])
	growth_curve_new(ax=ax[2], X1=yI.Irr, Y1=86400 * yI.Mu,  X2=XI, Y2=MuI, xlabel='I ($\mu$ E m$^{-2}$ s$^{-1}$)', axes_title='C',
	             ylabel='Growth rate, $\mu$ (d$^{-1}$)',
	             text=['$C_x/K_c$ = %.2f' % round(Cx/model.p.K_C,2), '$N_x/K_n$ = %.2f' % round(Nx/model.p.K_N,2)],
	                 label2='Haldane model', text_lx=0.54)
	return plt.show()
	
def growth_curve_new(ax, X1, Y1, X2, Y2, xlabel, ylabel, text, axes_title, label1='BRAM', label2='Monod',
                     text_lx=0.35, text_ly=0.05):
	ax.set_title(axes_title, fontsize=24, fontweight='bold', loc='left')
	ax.plot(X1, Y1, linewidth=2.0, label=label1)
	ax.plot(X2, Y2, 'o', label=label2)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.tick_params(labelsize=15)
	ax.legend(loc=7)
	plt.tight_layout()
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	textstr = '\n'.join((text))
	# place a text box in bottom right in axes coords
	ax.text(text_lx, text_ly, textstr, transform=ax.transAxes, fontsize=14,
	        verticalalignment='bottom', bbox=props)

def basic_plot(X, Y, xlabel, ylabel, text, text_lx=0.3, text_ly=0.1):
	fig, ax = plt.subplots(figsize=[4, 4])
	plt.plot(X, Y, linewidth=2.0)
	plt.xlabel(xlabel, fontsize=15)
	plt.ylabel(ylabel, fontsize=15)
	ax.tick_params(labelsize=15)
	plt.tight_layout()
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	textstr = '\n'.join((text))
	# place a text box in upper left in axes coords
	ax.text(text_lx, text_ly, textstr, transform=ax.transAxes, fontsize=14,
	        verticalalignment='bottom', bbox=props)
	return plt.show()


def plot_gle_opp(Nx1, Nx2, mu1, mu2):
	fig, ax = plt.subplots(figsize=[5, 4])
	plt.plot(Nx1, mu1, label='Gleaner', linewidth=2.0)
	plt.plot(Nx2, mu2, label='Opportunist', linewidth=2.0)
	plt.legend(fontsize=15)
	ax.tick_params(labelsize=14)
	plt.xlabel('External Nitrogen, N$_x$ ($\mu$M)', fontsize=15)
	plt.ylabel('Growth rate, $\mu$ (d$^{-1}$)', fontsize=15)
	plt.tight_layout()
	return plt.show()


def gle_opp_twoaxis(nparray, Iarray):
	legends = ['Gleaner', 'Opportunist', 'N$_x$']
	fig, ax1 = plt.subplots(figsize=(8, 3))
	for i in range(len(legends)):
		ax1.plot(nparray[:, 0], nparray[:, i + 1], label=legends[i])
	ax1.set_xlabel('time (days)', fontsize=15)
	ax1.set_ylabel('concentration', fontsize=15)
	plt.legend(loc=8)

	ax2 = ax1.twinx()
	ax2.plot(nparray[:, 0], Iarray, 'r--')
	ax2.tick_params('y', colors='r')
	ax2.set_ylabel('I ($\mu$ E mu$^{-2}$ s$^{-1}$)', fontsize=15, color='r')
	plt.tight_layout()
	return plt.show()


def nP_ylabel(v):
	if v == 'cyano':
		nP = [7358., 95451., 1681., 2000., 20000., 10000., 1000., 1000.]
		ylabel_P = ['Ribosomes', 'PSU', 'CCM', 'P$_3$', 'M', 'P$_N$', 'Dummy P', 'Q']
		ylabel_f = ['v$_1$', 'v$_2$', 'v$_3$', 'v$_4$', 'v$_5$', 'v$_6$', 'v$_d$']
	else:
		print("check model name")
	return nP, ylabel_P, ylabel_f


def mu_protein(mu1, mu2, mu3, f1, f2, f3):
	nP, ylabel_P, ylabel_f = nP_ylabel('cyano')
	fig = plt.figure(figsize=(12, 6))
	ncol = 3
	for i in range(6):
		fig.add_subplot(2, ncol, i + 1)
		P1 = (nP[i] * f1[:, i]) / np.sum(nP * f1, axis=1)
		P2 = (nP[i] * f2[:, i]) / np.sum(nP * f2, axis=1)
		P3 = (nP[i] * f3[:, i]) / np.sum(nP * f3, axis=1)
		plt.ylabel(ylabel_P[i], fontsize=12)
		plt.xlabel('$\mu$', fontsize=12)
		plt.plot(mu1, P1, label='C$_x$')
		plt.plot(mu2, P2, label='N$_x$')
		plt.plot(mu3, P3, label='I')
		plt.legend()
	plt.tight_layout()
	return plt.show()


def mu_vs_protein(mu, P_conc, xlabel='$\mu$', norm='yes', v='cyano'):
	nP, ylabel_P, ylabel_f = nP_ylabel(v)
	fig = plt.figure(figsize=(12, 6))
	if v == 'cyano' or v == 'glycolysis':
		ncol = 4
	else:
		ncol = 3
	for i in range(len(P_conc[0])):
		fig.add_subplot(2, ncol, i + 1)
		if norm == 'yes':
			P = (nP[i] * P_conc[:, i]) / np.sum(nP * P_conc, axis=1)
			plt.plot(mu, P, '-bo', markersize=4)
		else:
			plt.plot(mu, nP[i] * P_conc[:, i], '-bo', markersize=4)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(-1, 1))
		plt.xlabel(xlabel, fontsize=15)
		plt.ylabel(ylabel_P[i], fontsize=15)
	fig.tight_layout()
	#		fig.savefig('Images/glycolysis/g'+str(v)+'/mu_vs_normProtein'+str(v)+'.png')
	return plt.show()


def nutrients_vs_protein(nutrient_conc, P_conc, xlabel='glc$^x$', v='cyano', norm='yes'):
	nP, ylabel_P, ylabel_f = nP_ylabel(v)
	fig = plt.figure(figsize=(12, 6))
	if v == 'cyano' or v == 'glycolysis':
		ncol = 4
	else:
		ncol = 3
	for i in range(len(P_conc[0])):
		ax = fig.add_subplot(2, ncol, i + 1)
		if norm == 'yes':
			P = (nP[i] * P_conc[:, i]) / np.sum(nP * P_conc, axis=1)
			plt.plot(nutrient_conc, P, '-bo', markersize=3)
		else:
			plt.plot(nutrient_conc, nP[i] * P_conc[:, i], '-bo', markersize=3)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(-1, 1))
		#			ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
		plt.xlabel(xlabel, fontsize=12)
		plt.ylabel(ylabel_P[i], fontsize=12)
	fig.tight_layout()
	#		fig.savefig('Images/glycolysis/g'+str(v)+'/glcx_vs_normProtein'+str(v)+'.png')
	return plt.show()


def plot_surface(X, Y, Z, xlable='[Glc$_x$]', ylable='[O$_2^x$]', zlable='$\mu$'):
	x, y = np.meshgrid(X, Y)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surface = ax.plot_surface(X=x, Y=y, Z=Z, cmap=cm.coolwarm)
	ax.set_xlabel(xlable)
	ax.set_ylabel(ylable)
	ax.set_zlabel(zlable)
	fig.colorbar(surface)
	return plt.show()


def line(x, y, xlabel='external carbon concentration, c$^x$',
         ylabel='growth rate, $\mu$ (day$^{-1}$)', fsize=15, l=5, b=4):
	fig = plt.figure(figsize=(l, b))
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
	plt.plot(x, y, linewidth=3.)
	plt.xlabel(xlabel, fontsize=fsize)
	plt.ylabel(ylabel, fontsize=fsize)
	fig.tight_layout()
	return plt.show()


def plot_N_Mu(Nx, Mu1, Mu2, l=5, b=4, fsize=15):
	fig = plt.figure(figsize=(l, b))
	plt.plot(Nx, Mu1, label='opportunist', linewidth=3.)
	plt.plot(Nx, Mu2, label='gleaner', linewidth=3.)
	plt.legend()
	plt.xlabel('Nitrogen concentration, N$^x$ ($\mu$ M)', fontsize=fsize)
	plt.ylabel('Growth rate, $\mu$ ($day^{-1}$)', fontsize=fsize)
	fig.tight_layout()
	return plt.show()
