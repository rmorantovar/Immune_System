#---------Import---------
import sys
sys.path.append('../Codes/lib/')
sys.path.append('../Codes/Python/')
from Immuno_models import*
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as mtick
from matplotlib.ticker import PercentFormatter
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d


#------------------------------------

Text_files_path = '../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/'

colors = plt.cm.plasma(np.linspace(0,1,10))
T = 0
dt = 1
Ls = [10, 15, 20, 50]
Types = ["MM"]
M_epitopes = 10

for Type in Types:
    if(Type == "MJ"):
        T = 150
    elif(Type == "MM"):
        T = 22
    for L in Ls:
        fig0, ax0 = plt.subplots(figsize = (12,6), gridspec_kw={'bottom':.15, 'left':.10, 'right':.87})
        fig, ax = plt.subplots(figsize = (12,10), gridspec_kw={'bottom':.12, 'left':.12})
        gain_funct = np.array([])
        for i in np.arange(M_epitopes):

            data = np.loadtxt(Text_files_path + 'Escaping_state_L-%d_L_alphabet-20_N_epitopes-%d_'%(L, i+1)+Type+'.txt')
            time = np.linspace(0, T*dt, T)
            f = interp1d(time, data)
            time_array = np.linspace(0, T*dt, 400)
            gain_funct = np.append(gain_funct, time_array[f(time_array)>.5][-1])
            ax.plot(time, data, color = colors[i], linestyle = '--', label = "%d"%(i+1))
    
        my_plot_layout(ax=ax, xlabel = 'Steps', ylabel = "Retained memory", x_fontsize = 24, y_fontsize = 24)
        ax.set_ylim(-.1, 1.1)
        ax.legend(loc = 5, fontsize = 24)
        fig.savefig('../../Figures/4_Complexity/memory_lost_L-%d_'%(L)+Type+'.png')

        M_epitopes_array = np.linspace(1, M_epitopes, 100)
        ax0.plot(np.arange(1, M_epitopes+1), gain_funct, color = 'darkgreen', lw = 2, ls = 'dashed', alpha = .8)
        ax0.tick_params(axis='y', labelcolor='darkgreen')
        ax02 = ax0.twinx()
        cost_func = my_linear_func(np.arange(M_epitopes), 1, .3)
        ax02.plot(np.arange(1, M_epitopes+1), -cost_func, color = 'darkred', lw = 2, ls = 'dotted', alpha = .8)
        ax02.tick_params(axis='y', labelcolor='darkred')
        fitness = gain_funct - cost_func
        f2 = interp1d(np.arange(1, M_epitopes+1), fitness)
        ax0.plot(np.arange(1, M_epitopes+1), fitness, color = 'darkgrey', label = 'Total fitness', lw = 3)
        ax0.vlines(M_epitopes_array[f2(M_epitopes_array)==np.max(f2(M_epitopes_array))], ax0.get_ylim()[0], np.max(fitness), color = 'black', ls = '-')
        ax0.scatter(M_epitopes_array[f2(M_epitopes_array)==np.max(f2(M_epitopes_array))], np.max(fitness), color = 'black', marker = '*', s = 35, zorder=10)
        my_plot_layout(ax=ax0, xlabel = 'Nº of epitopes', ylabel = 'Time before lost of immunity', x_fontsize = 24, y_fontsize = 24)
        my_plot_layout(ax=ax02, ylabel = 'Cost of Immunity', y_fontsize = 24)
        #ax0.legend(loc = 0, fontsize = 22)
        #fig0.tight_layout()
        fig0.savefig('../../Figures/4_Complexity/Gain_function_L-%d_'%(L)+Type+'.png')
        fig0.savefig('../../Figures/4_Complexity/pdfs/Gain_function_L-%d_'%(L)+Type+'.pdf')



