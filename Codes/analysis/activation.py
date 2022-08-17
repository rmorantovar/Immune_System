import sys
sys.path.append('../library/')
from Immuno_models import*
plt.rcParams['text.usetex'] = True

Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/"


fig_p_a, ax_p_a = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15, 'bottom': .15, 'right':.8})
# Paramters:
N_engagements = 1e2
k_on = 1e6
m_off = 0
k_pr_array = np.array([1e-4, 1])
k_pr_array = np.array([0.1, 1, 3])
markers = ['^', 's', 'o']
colors = ['darkred', 'darkblue', 'green']
#labels = [r'$m_{\mathrm{on}}\ll \k_act$', r'$m_{\mathrm{on}}\approx \k_act$', r'$m_{\mathrm{on}}\gg \k_act$']
#l_off = k_off
l_on = 0
k_act = 0.5
k_act = k_act*24
k_off_array2 = np.logspace(-6, 2, 100)
k_off_array2 = np.logspace(-5, 5, 100)

for j, k_pr in enumerate(k_pr_array):

	file_p_a = open(Text_files_path+'recognition_k_act-%.0e_k_pr-%.0e.pkl'%(k_act/24, k_pr),'rb')
	k_off_array, p_a_array = pickle.load(file_p_a)

	ax_p_a.plot(k_off_array/k_on, p_a_array, marker = markers[j], linestyle = '', ms = 8, color= colors[j], label = r'$%.2f min^{-1}$ '%(k_pr), alpha = .8)
	ax_p_a.plot(k_off_array2/k_on, (k_pr_array[j]*k_act)/((k_pr_array[j]*k_act) + (k_pr_array[j] + k_act)*k_off_array2 +k_off_array2**2), marker = '', linestyle = '-', ms = 8, color = colors[j])
	#ax_p_a.plot(k_off_array2/k_on, np.ones_like(k_off_array2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)
	ax_p_a.plot(k_off_array2/k_on, (k_pr_array[j]*k_act)/((k_pr_array[j]*k_act)+k_off_array2**2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)

	#if(k_pr<k_act):
	#	ax_p_a.plot(k_off_array2/k_on, k_pr_array[j]/(k_pr_array[j]+k_off_array2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)
	#if(k_pr>k_act):
	#	ax_p_a.plot(k_off_array2/k_on, k_act/(k_act+k_off_array2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)
ax_p_a.hlines(.5, np.min(k_off_array2)/k_on, np.max(k_off_array2)/k_on)
ax_p_a.vlines([np.concatenate((k_pr_array/k_on, [0]))], 1e-6, 1, linestyle=':', alpha = .4, color=colors )
my_plot_layout(ax=ax_p_a, yscale = 'log', xscale = 'log', ticks_labelsize = 24, xlabel = r'$k_{\mathrm{off}}$', ylabel = r'$p_a$', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24)
#ax_p_a.set_ylim(bottom = 1e-7)
ax_p_a.legend(loc = 3, fontsize = 24, title = r'$k_{act}=%.2f$'%(k_act), title_fontsize = 24)
fig_p_a.savefig('../../Figures/7_Recognition/K-PR_model.pdf')


