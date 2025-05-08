import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 1
N_r = 1e8
T0 = 0
Tf = 7
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
k_pr = k_pr*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*1000
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
AA = 1
p = 3

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'CMFILVWYAGTSQNEDHRKPFMRTPTRMCW' #L=30
L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
model = 'TCRen'
#model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, model, Text_files_path)
print('min_e_PWM=%.4f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

avg_E = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]) + E_ms
var_E = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.05, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])
beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
print('beta_r = %.1f'%beta_r, 'K_d_r = %.2e'%Kd_r)

S = np.cumsum(betas[:-1]*dE)
Omega = np.sum(np.exp(S)*dE)

for i in np.linspace(3, 10, 20):
    u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*i)/N_A, Es, p, lambda_A, N_c, dE)
    E_peak = Es[:-1][QR == np.max(QR)]

    int_QR = np.sum(QR*dE)
    int_Q01_approx = (Q0[Es[:-1]<E_peak][-1]/(betas[Es<E_peak][-1]))
    #int_Q01_approx = ((Q0[Es[:-1]<E_peak][-1]/(betas[Es<E_peak][-1])) - (Q0[0]/betas[0]))
    int_RQ2_approx = (Q0[Es[:-1]<E_peak][-1]/(p - betas[Es<E_peak][-1]))
    int_QR_approx = int_Q01_approx + int_RQ2_approx
    #int_QR_approx = (Q0[Es[:-1]<E_peak][-1]/(p))
    plt.plot(np.exp(E_peak), int_QR_approx/int_QR, marker = 'o', color = 'red')
    #plt.plot(np.exp(E_peak), int_Q01_approx/int_QR, marker = 'o', color = 'blue')
    #plt.plot(np.exp(E_peak), int_RQ2_approx/int_QR, marker = 'o', color = 'green')
    #plt.plot(np.exp(E_peak), 2*int_Q01_approx/int_QR, marker = 'o', color = 'brown')

plt.vlines([Kd_p, Kd_pr, Kd_r], 1, 2, ls = ['-', '--', ':'])
plt.hlines([.5, 1], 1e-8, 1e-5)
plt.yscale('log')
plt.xscale('log')
plt.show()








# for i in range(-5, 8, 1):

#     E_lim = E_peak + i
#     print(i, E_lim)
#     beta_lim = betas[Es<E_lim][-1]
#     int_Q01 = np.sum(Q0[Es[:-1]<E_lim]*dE[Es[:-1]<E_lim])
#     int_Q01_approx = (Q0[Es[:-1]<E_lim][-1]/(betas[Es<E_lim][-1]))
#     plt.plot(E_lim, int_Q01/int_Q01_approx, marker = 'o', color = 'blue')






