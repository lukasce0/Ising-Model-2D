
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats





"""
Plotting the expectation values as a function of temperature for a 2x2 lattice
"""





expe_val_4b_num = np.loadtxt("expe_val_4b.txt").T

T = expe_val_4b_num[0]




def findExpectationTheoretical(T, L):
    
    N = L**2
    
    Z = 4 * (3 + np.cosh(8 / T))
    
    E_eps_theor = -32 * np.sinh(8 / T) / (Z*N)
    E_eps2_theor = 256 * np.cosh(8 / T) / (Z*N**2)
    E_mag_theor = 8 * (2 + np.exp(8 / T)) / (Z*N)
    E_mag2_theor = 32 * (1 + np.exp(8 / T)) / (Z*N**2)
    
    c_V_theor = N * 1/T**2 * (E_eps2_theor - E_eps_theor**2)
    chi_theor = N * 1/T * (E_mag2_theor - E_mag_theor**2)

    return_vec = [E_eps_theor, E_mag_theor, c_V_theor, chi_theor]

    return return_vec




expe_val_1_theor = findExpectationTheoretical(T, 2)
plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(16,10))

label_list_4b = {0:r"$\langle \epsilon \rangle$ [J]", 1:r"$\langle |m| \rangle$",
                2:r"$\langle c_V \rangle$ $[k_B]$", 3:r"$\langle \chi \rangle$ [1/J]"}

for i in range(4):
    fig = plt.subplot(2,2,i+1)
    
    if i > 1: #x-label only for the bottom plots
        plt.xlabel(r"T [J/$k_B$]")
    plt.ylabel(label_list_4b[i])
    
    plt.plot(T, expe_val_1_theor [i], color="r")
    plt.scatter(T, expe_val_4b_num[i+1])
    plt.grid(linestyle = '--')
    fig.set_facecolor("#EBECF0")


plt.savefig("4b.pdf")










"""
Plottign the error in eps expectation value for two differant temperatures
and as a function of Monte Carlo cycles
"""



E_eps_1_num = np.loadtxt("expe_val_4c_1.txt").T[0]
E_eps_24_num = np.loadtxt("expe_val_4c_24.txt").T[0]


T = 1. #J/k_B
expe_val_1_4c_theor = findExpectationTheoretical(T, 2)
T = 2.4 #J/k_B
expe_val_24_4c_theor = findExpectationTheoretical(T, 2)


abs_error_1 = abs(E_eps_1_num - expe_val_1_4c_theor[0])
abs_error_24 = abs(E_eps_24_num - expe_val_24_4c_theor[0]) 



plt.rcParams.update({'font.size': 17})
fig, ax = plt.subplots(figsize=(8,6))
plt.xlabel(r"Monte Cralo cycles")
plt.ylabel(r"$log_{10}$ of absolute error")
plt.plot(np.log10(abs_error_1), label=r"T = 1.0 [$J/k_B$]")
plt.plot(np.log10(abs_error_24), label=r"T = 2.4 [$J/k_B$]")
ax.set_facecolor("#EBECF0")
plt.grid(linestyle = '--')
plt.legend()

plt.savefig("4c.pdf")













"""
Plotting the expectation values for 20x20 lattice, with two differant temperatures
and at two differant initial states. We want to study burn-in.
"""




expe_val_5a1o = np.loadtxt("expe_val_5a_1_o.txt").T[:2]
expe_val_5a24o = np.loadtxt("expe_val_5a_24_o.txt").T[:2]
expe_val_5a1r = np.loadtxt("expe_val_5a_1_r.txt").T[:2]
expe_val_5a24r = np.loadtxt("expe_val_5a_24_r.txt").T[:2]


fig, ax = plt.subplots(figsize=(8,6))
plt.xlabel(r"Monte Cralo cycles")
plt.ylabel(r"$\langle \epsilon \rangle$ [J]")

plt.plot(expe_val_5a1o[0], label=r"T = 1.0 [$J/k_B$] ordered")
plt.plot(expe_val_5a24o[0], label=r"T = 2.4 [$J/k_B$] ordered")
plt.plot(expe_val_5a1r[0], label=r"T = 1.0 [$J/k_B$] random")
plt.plot(expe_val_5a24r[0], label=r"T = 2.4 [$J/k_B$] random")

plt.legend()
ax.set_facecolor("#EBECF0")
plt.grid(linestyle = '--')
plt.savefig("5e.pdf")




fig, ax = plt.subplots(figsize=(8,6))
plt.xlabel(r"Monte Cralo cycles")
plt.ylabel(r"$\langle |m| \rangle$ [J]")

plt.plot(expe_val_5a1o[1], label=r"T = 1.0 [$J/k_B$] ordered")
plt.plot(expe_val_5a24o[1], label=r"T = 2.4 [$J/k_B$] ordered")
plt.plot(expe_val_5a1r[1], label=r"T = 1.0 [$J/k_B$] random")
plt.plot(expe_val_5a24r[1], label=r"T = 2.4 [$J/k_B$] random")

ax.set_facecolor("#EBECF0")
plt.grid(linestyle = '--')
plt.savefig("5m.pdf")





"""
Here we will wplot the pdf of eps given two differant T for 20x20 lattice.
"""

energy_val_6_1 = np.loadtxt("energy_val_6_1.txt")
energy_val_6_24 = np.loadtxt("energy_val_6_24.txt")

fig, ax = plt.subplots(figsize=(10,8))
plt.xlabel(r"$ \epsilon $ [J]")
plt.ylabel(r"Probability of $\epsilon$ being in bin")

weights_1 = np.ones_like(energy_val_6_1) / len(energy_val_6_1)
weights_24 = np.ones_like(energy_val_6_24) / len(energy_val_6_24)


plt.hist(energy_val_6_24, bins=100, weights=weights_24, histtype = 'bar', log=True, label=r"T = 1.0 [$J/k_B$]")
plt.hist(energy_val_6_1, bins=50, weights=weights_1, log=True,  label=r"T = 2.4 [$J/k_B$]")
#histtype{'bar', 'barstacked', 'step', 'stepfilled'}
plt.legend()

ax.set_facecolor("#EBECF0")
plt.grid(linestyle = '--')
plt.savefig("6.pdf")

var_1 = np.std(energy_val_6_1)**2
var_24 = np.std(energy_val_6_24)**2

print(f"Variance of T = 1.0 is {var_1}")
print(f"Variance of T = 2.4 is {var_24}")








"""
Here we will plot the expectation values as a function of temperature (around critical 
temperature) for differant lattice sizes.
"""





expe_val_40 = np.loadtxt("expe_val_40_8.txt").T
expe_val_60 = np.loadtxt("expe_val_60_8.txt").T
expe_val_80 = np.loadtxt("expe_val_80_8.txt").T
expe_val_100 = np.loadtxt("expe_val_100_8.txt").T
T = expe_val_40[0]

plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(16,10))


for i in range(4):
    fig1 = plt.subplot(2,2,i+1)
    
    if i > 1: #x-label only for the bottom plots
        plt.xlabel(r"T [J/$k_B$]")
    plt.ylabel(label_list_4b[i])
    
    fig1.scatter(T, expe_val_40[i+1])
    fig1.scatter(T, expe_val_60[i+1])
    fig1.scatter(T, expe_val_80[i+1])
    fig1.scatter(T, expe_val_100[i+1])
    fig1.grid(linestyle = '--')
    fig1.set_facecolor("#EBECF0")


plt.savefig("8.pdf")



"""
Plotting the crticial temperatures and finding the L->oo value using
linear regression
"""

TC_list = [[],[],[],[]]


c = 4 #3 for cV 4 for xi

cv_40 = expe_val_40[c]
cv_60 = expe_val_60[c]
cv_80 = expe_val_80[c]
cv_100 = expe_val_100[c]

for i in range(5): #choose over how many points the average will be clculated
    max_40 = max(cv_40)
    max_60 = max(cv_60)
    max_80 = max(cv_80)
    max_100 = max(cv_100)
    
    max_i_40 = np.where(cv_40==max_40)[0][0]
    max_i_60 = np.where(cv_60==max_60)[0][0]
    max_i_80 = np.where(cv_80==max_80)[0][0]
    max_i_100 = np.where(cv_100==max_100)[0][0]
    
    TC_list[0].append(T[max_i_40])
    TC_list[1].append(T[max_i_60])
    TC_list[2].append(T[max_i_80])
    TC_list[3].append(T[max_i_100])

    cv_40[max_i_40] = 0
    cv_60[max_i_60] = 0
    cv_80[max_i_80] = 0
    cv_100[max_i_100] = 0
    

TC_40 = sum(TC_list[0])/len(TC_list[0])
TC_60 = sum(TC_list[1])/len(TC_list[1])
TC_80 = sum(TC_list[2])/len(TC_list[2])
TC_100 = sum(TC_list[3])/len(TC_list[3])

TC_plot = [TC_40, TC_60, TC_80, TC_100]
L_plot = [1/40, 1/60, 1/80, 1/100]

slope, intercept, a, b, c = stats.linregress(L_plot, TC_plot)

x = np.linspace(L_plot[-1], L_plot[0], 100)
lineFind = lambda x: slope*x + intercept

line_y = lineFind(x)

plt.figure(figsize=(10,8))
plt.xlabel(r"1/L")
plt.ylabel(r"$T_C$ $[J/k_B]$")
plt.plot(x,line_y, color = "k", label="Linear regression")
plt.scatter(1/40,TC_40, label="L = 40")
plt.scatter(1/60, TC_60, label="L = 60")
plt.scatter(1/80, TC_80, label="L = 80")
plt.scatter(1/100, TC_100, label="L = 100")
plt.legend()
plt.savefig("9.pdf")



def Delta_l(x,y): #uncertainty in the slope and intercept
    n = len(x)
    D = sum(x**2)-1/n*(sum(x))**2
    E = sum(x*y)-1/n*sum(x)*sum(y)
    F = sum(y**2)-1/n*(sum(y))**2
    da = np.sqrt(1/(n-2)*(D*F-E**2)/D**2)
    db = np.sqrt(1/(n-2)*(D/n+sum(x)/len(x)**2)*(D*F-E**2)/E**2)
    return da, db


diviasion_list = [TC_40-lineFind(1/40), TC_60-lineFind(1/60), TC_80-lineFind(1/80), TC_100-lineFind(1/100)]
satndard_diviasion = np.std(diviasion_list)
d_slope, d_intercept = Delta_l(np.array(L_plot), np.array(TC_plot))

tot_uncertainty = np.sqrt(satndard_diviasion**2 + d_intercept**2) #using variance superpositioning

print(f"The TC is estimated to {intercept}")
print(f"The standard diviasion in the approximation {tot_uncertainty}")


