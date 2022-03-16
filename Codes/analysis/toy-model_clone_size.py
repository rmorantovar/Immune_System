import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import pickle
import seaborn as sns
from scipy.optimize import curve_fit

colors = ['tab:blue','tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = ['exponential', 'linear']
energy_models = ['MJ', 'MM']
text_pos = np.array([[4e0, 5e-7],[1e2, 6e-2]], dtype = object)
text = [r'$\sim N_{b}^{-1-\frac{\lambda\alpha}{\beta}}$', r'$\sim N_{b}^{-1}$']

path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Python/ensemble/"

alpha = 1
betas = [2.0, 1.0, 0.5]
betas = [2.0, 1.0, 0.5]
ds = [20, 2]
Ns = np.flip([2e2, 2e3, 2e4])
Ns = np.flip([2e3, 2e4])
L = 15
e0 = 4
lamda = 0.62
k = np.linspace(1, L, L*2);

colors_d = [plt.cm.BuPu(np.linspace(0,1,len(Ns)+2)), plt.cm.OrRd(np.linspace(0,1,len(Ns)+2)), plt.cm.GnBu(np.linspace(0,1,len(Ns)+2)),]
markers_beta = ['*', '^', 's', 'o']

fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'bottom':.15})
fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'bottom':.15})
fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'bottom':.15})

for j, d in enumerate(ds):
	geo_exp = np.array([])
	for n, N in enumerate(Ns):
		for i, beta in enumerate(betas):

			exponents = [(-((alpha*lamda)/beta)-1), -1]
			file_clone_sizes = open(path+'clone_size_data_d-%d_beta-%.2f_N-%.1e_'%(d, beta, N)+growth_models[0]+'.pkl','rb')
			#file_clone_sizes = open(path+'clone_size_data_d-%d_beta-%.2f_N-%.1e_'%(d, beta, N)+growth_models[0]+'_'+energy_models[1]+'.pkl','rb')

			clone_sizes = pickle.load(file_clone_sizes)
			#clone_sizes = clone_sizes/np.max(clone_sizes)

			#clone_sizes = clone_sizes/np.max(clone_sizes)
			print(len(clone_sizes))
			file_clone_sizes.close()
			#print(clone_sizes)
			density = True
			clone_size_distribution = np.histogram(clone_sizes, bins = np.logspace(np.log10(np.min(clone_sizes)),np.log10(np.max((clone_sizes))),15), density = density)
			clone_size = clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
			clone_size_counts = clone_size_distribution[0][np.where(clone_size_distribution[0]!=0)]

			delta_clone_size = clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]-clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
			cumsum_clone_size_counts = np.cumsum(clone_size_counts*delta_clone_size)
			#cumsum_clone_size_counts = np.cumsum(clone_size_counts)
			print(cumsum_clone_size_counts[-1])

			popt, pcov = curve_fit(my_linear_func, np.log(clone_size[-3:]), np.log(clone_size_counts[-3:]))
			plaw_fit_csd = clone_size**(popt[1])#*(np.log(clone_size**(1/nu)))**(1-lambd)
			plaw_fit_csd /= (plaw_fit_csd[-2]/(clone_size_counts[-2]))


			#popt, pcov = curve_fit(my_linear_func, np.log(clone_size[-4:-1]), np.log(1-(cumsum_clone_size_counts[:])[-4:-1]))			
			#plaw_fit_csd = clone_size**(popt[1])#*(np.log(clone_size**(1/nu)))**(1-lambd)
			#plaw_fit_csd /= (plaw_fit_csd[-4]/((cumsum_clone_size_counts[-1]-cumsum_clone_size_counts[:])[-4]))
			
			#Distribution from simulations
			ax.plot(clone_size[:], clone_size_counts[:], linestyle = '-', marker = markers_beta[i], ms = 8, linewidth = 2, color = colors_d[j][n+2], alpha = .8)
			#ax.plot(clone_size[:-1], (cumsum_clone_size_counts[:])[:-1], linestyle = '-', marker = markers_beta[i], ms = 6, linewidth = 2, color = colors_d[j][n+2])

			#Fit
			ax.plot(clone_size[:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 8, linewidth = 4, alpha = .8, color = colors_d[j][n+2])
			

			print(popt[1], -(popt[1]+1)*beta)
			geo_exp = np.append(geo_exp, -(popt[1]+1)*beta)
			slope = (np.log(clone_size_counts[1:])-np.log(clone_size_counts[:-1]))/(np.log(clone_size[1:])-np.log(clone_size[:-1]))
			#ax2.scatter(-np.log10(N), -(popt[1]+1)*beta, color = colors_d[j][0+2], marker = markers_beta[i])
			ax2.scatter(d, -(popt[1]+1)*beta, color = colors_d[j][n+2], marker = markers_beta[i])
			#ax2.plot(-np.log(clone_size[1:])/beta, -((slope)+1)*beta, linestyle = '--', marker = markers_beta[i], ms = 8, linewidth = 4, alpha = .8, color = colors_d[j][n+2])


	#ax2.plot(e0*1.2*(k-L), np.ones_like(k)*np.log((L/k)*(d-1)), linewidth = 4 , color = colors_d[j][-1], linestyle = '--');
	#ax2.plot(e0*1.2*(k-L), (4*((d-1)/d) - (4*k)/(L)), linewidth = 4 , color = colors_d[j][-1], linestyle = ':');


d = 20
geo_exp = np.array([])
for n, N in enumerate(Ns):
	for i, beta in enumerate(betas):

		exponents = [(-((alpha*lamda)/beta)-1), -1]
		file_clone_sizes = open(path+'clone_size_data_d-%d_beta-%.2f_N-%.1e_'%(d, beta, N)+growth_models[0]+'_'+energy_models[0]+'.pkl','rb')

		clone_sizes = pickle.load(file_clone_sizes)
		#clone_sizes = clone_sizes/np.max(clone_sizes)

		#clone_sizes = clone_sizes/np.max(clone_sizes)
		print(len(clone_sizes))
		file_clone_sizes.close()
		#print(clone_sizes)
		density = True
		clone_size_distribution = np.histogram(clone_sizes, bins = np.logspace(np.log10(np.min(clone_sizes)),np.log10(np.max((clone_sizes))),15), density = density)
		clone_size = clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
		clone_size_counts = clone_size_distribution[0][np.where(clone_size_distribution[0]!=0)]

		delta_clone_size = clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]-clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
		cumsum_clone_size_counts = np.cumsum(clone_size_counts*delta_clone_size)
		#cumsum_clone_size_counts = np.cumsum(clone_size_counts)
		print(cumsum_clone_size_counts[-1])

		popt, pcov = curve_fit(my_linear_func, np.log(clone_size[-3:]), np.log(clone_size_counts[-3:]))
		plaw_fit_csd = clone_size**(popt[1])#*(np.log(clone_size**(1/nu)))**(1-lambd)
		plaw_fit_csd /= (plaw_fit_csd[-3]/(clone_size_counts[-3]))


		#popt, pcov = curve_fit(my_linear_func, np.log(clone_size[-4:-1]), np.log(1-(cumsum_clone_size_counts[:])[-4:-1]))			
		#plaw_fit_csd = clone_size**(popt[1])#*(np.log(clone_size**(1/nu)))**(1-lambd)
		#plaw_fit_csd /= (plaw_fit_csd[-4]/((cumsum_clone_size_counts[-1]-cumsum_clone_size_counts[:])[-4]))
		
		#Distribution from simulations
		ax.plot(clone_size[:], clone_size_counts[:], linestyle = '-', marker = markers_beta[i], ms = 10, linewidth = 2, color = colors_d[0][n+2])
		#ax.plot(clone_size[:-1], (cumsum_clone_size_counts[:])[:-1], linestyle = '-', marker = markers_beta[i], ms = 6, linewidth = 2, color = colors_d[j][n+2])

		#Fit
		ax.plot(clone_size[:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 8, linewidth = 4, alpha = .8, color = colors_d[0][n+2])
		

		print(popt[1], -(popt[1]+1)*beta)
		geo_exp = np.append(geo_exp, -(popt[1]+1)*beta)
		slope = (np.log(clone_size_counts[1:])-np.log(clone_size_counts[:-1]))/(np.log(clone_size[1:])-np.log(clone_size[:-1]))
		#ax2.scatter(-np.log10(N), -(popt[1]+1)*beta, color = 'darkblue', marker = markers_beta[i], s = 50)
		ax2.scatter(d+2, -(popt[1]+1)*beta, color = 'darkblue', marker = markers_beta[i], s = 50)
		#ax2.plot(-np.log(clone_size[1:])/beta, -((slope)+1)*beta, linestyle = '--', marker = markers_beta[i], ms = 8, linewidth = 4, alpha = .8, color = colors_d[j][n+2])


#ax2.plot(e0*1.2*(k-L), np.ones_like(k)*np.log((L/k)*(d-1)), linewidth = 4 , color = colors_d[j][-1], linestyle = '--');
#ax2.plot(e0*1.2*(k-L), (4*((d-1)/d) - (4*k)/(L)), linewidth = 4 , color = colors_d[j][-1], linestyle = ':');



ds = np.arange(2, 20)

D, K = np.meshgrid(ds, k)

Lamda = np.log((L/K)*(D-1))/e0

CS = ax0.contourf(D, K, Lamda)
cbar = fig0.colorbar(CS)

my_plot_layout(ax = ax0, xlabel=r'$d$', ylabel=r'$k$', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
my_plot_layout(ax = ax, xscale='log', yscale= 'log', xlabel=r'Clone size $N_{b}$', ylabel=r'$p(N_{b})$', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
my_plot_layout(ax = ax2, xlabel=r'$d$', ylabel=r'$\lambda$', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )

#ax.legend(title=r'$\delta$', fontsize = 22, title_fontsize = 22)
fig0.savefig('../../Figures/1_Dynamics/geo_exp_toy0.pdf')
fig.savefig(f'../../Figures/1_Dynamics/CSD_toy_density-{density}.pdf')
fig2.savefig('../../Figures/1_Dynamics/geo_exp_toy.pdf')

