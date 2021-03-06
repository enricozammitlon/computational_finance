import matplotlib.pyplot as plt
import numpy as np
import csv
import scipy.stats as si

X=47.66
R=2
F=95
T=2.0
C = 1.09
alpha = 0.02
r = 0.0229
T = 2.
sigma = 0.416

def euro_vanilla_call(S, K, T, r, sigma):

    #S: spot price
    #K: strike price
    #T: time to maturity
    #r: interest rate
    #sigma: volatility of underlying asset

    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))

    call = (R*S * si.norm.cdf(d1, 0.0, 1.0) - K * np.exp(-r * T) * si.norm.cdf(d2, 0.0, 1.0))

    return call

variationData=[]
with open('data/varying_imax.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['sMax','iMax','jMax','S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'x':[],'y':[],'jMax':0,'sMax':0}
    for row in reader:
        currentData['x'].append(row['iMax'])
        currentData['y'].append(row['V'])
        currentData['jMax']=row['jMax']
        currentData['sMax']=int(row['sMax']/X)
    variationData.append(currentData)

plt.figure()
plt.grid()
plt.plot(variationData[0]['x'][:40],variationData[0]['y'][:40],label=r'$\beta=0.486,\sigma=3.03,j_{max}=%i,s_{max}=%i$'%(variationData[0]['jMax'],variationData[0]['sMax']),linewidth=2)
plt.xlabel(r'$i_{max}$')
plt.ylabel(r'$V(X,t=0)$')
plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
plt.savefig('images/european_varying_imax.png',bbox_inches='tight', pad_inches=0.2)

with open('data/varying_smax.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['sMax','iMax','jMax','S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'x':[],'y':[],'jMax':[],'iMax':0}
    for row in reader:
        currentData['x'].append(int(row['sMax']/X))
        currentData['y'].append(row['V'])
        currentData['jMax'].append(row['jMax'])
        currentData['iMax']=row['iMax']
    variationData.append(currentData)

fig, ax1 = plt.subplots()
ax1.set_xlabel(r's_{max} (multiples of X)')
ax1.set_ylabel(r'$V(X,t=0)$')
ax1.grid()
ax1.scatter(np.asarray(variationData[1]['x'][:20]),variationData[1]['y'][:20],label=r'$V(X,t=0)$ for $\beta=0.486,\sigma=3.03,i_{max}=%i$'%(variationData[1]['iMax']))
ax2 = ax1.twinx()
ax2.set_ylabel(r'$j_{max}$')
fig.tight_layout()
ax2.plot(np.asarray(variationData[1]['x'][:20]),variationData[1]['jMax'][:20],label=r'$j_{max}$',color="orange")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2,loc='lower right',fancybox=False, framealpha=0.0)
plt.savefig('images/european_varying_smax_zoomed.png',bbox_inches='tight', pad_inches=0.2)

fig, ax1 = plt.subplots()
ax1.set_xlabel(r's_{max} (multiples of X)')
ax1.set_ylabel(r'$V(X,t=0)$')
ax1.grid()
ax1.scatter(np.asarray(variationData[1]['x']),variationData[1]['y'],label=r'$V(X,t=0)$ for $\beta=0.486,\sigma=3.03,i_{max}=%i$'%(variationData[1]['iMax']))
ax2 = ax1.twinx()
ax2.set_ylabel(r'$j_{max}$')
fig.tight_layout()
ax2.plot(np.asarray(variationData[1]['x']),variationData[1]['jMax'],label=r'$j_{max}$',color="orange")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2,loc='lower right',fancybox=False, framealpha=0.0)
plt.savefig('images/european_varying_smax.png',bbox_inches='tight', pad_inches=0.2)

for smax in (10,65):
    currentData={'x':[],'y':[],'iMax':0,'sMax':0}
    with open('data/smax_jmax/'+str(smax)+'_varying_jmax.csv', newline='\n') as csvfile:
        reader = csv.DictReader(csvfile,fieldnames=['sMax','iMax','jMax','S','V'],quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            currentData['x'].append(row['jMax'])
            currentData['y'].append(row['V'])
            currentData['iMax']=row['iMax']
            currentData['sMax']=int(row['sMax']/X)

    plt.figure()
    plt.plot(currentData['x'][:100],currentData['y'][:100],label=r'$\beta=0.486,\sigma=3.03,i_{max}=%i,s_{max}=%i$'%(currentData['iMax'],currentData['sMax']),linewidth=2)
    plt.xlabel(r'$j_{max}$')
    plt.ylabel(r'$V(X,t=0)$')
    plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
    plt.grid()
    plt.savefig('images/smax_jmax/'+str(smax)+'_european_varying_jmax.png',bbox_inches='tight', pad_inches=0.2)
    plt.close()

allData=[]
with open('data/varying_s_beta_1.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'S':[],'V':[]}
    for row in reader:
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
    allData.append(currentData)

with open('data/varying_s_beta_0_4.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'S':[],'V':[]}
    for row in reader:
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
    allData.append(currentData)

plt.figure()
plt.grid()
plt.plot(allData[0]['S'],allData[0]['V'],label=r'$\beta=1,\sigma=0.416$',linewidth=2)
plt.plot(allData[1]['S'],allData[1]['V'],label=r'$\beta=0.486,\sigma=3.03$',linewidth=2)
plt.xlabel(r'$S_0$')
plt.ylabel(r'$V(S,t=0)$')
plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
plt.savefig('images/european_varying_s.png',bbox_inches='tight', pad_inches=0.2)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


currentData={'S':[],'V':[],'beta':[],'sigma':[]}
with open('data/varying_s_sigma_beta.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['beta','sigma','S','V'],quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        if(row['V']==-1):
            continue
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
        currentData['beta'].append(row['beta'])
        currentData['sigma'].append(row['sigma'])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(currentData['sigma'],currentData['beta'],currentData['V'])
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$\beta$')
ax.set_zlabel(r'$V(X,t=0)$')
plt.show()
plt.savefig('images/european_varying_s_varying_sigma_varying_beta.png',bbox_inches='tight', pad_inches=0.2)

plt.figure()
plt.grid()
plt.xlabel(r'maximum iterations')
plt.ylabel(r'$time$ (ms)')

currentData={'x':[],'y':[],'iMax':0,'sMax':0,'time':[]}
with open('data/smax_jmax/10_varying_jmax.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['sMax','iMax','jMax','S','V','time'],quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        currentData['x'].append(row['jMax'])
        currentData['y'].append(row['V'])
        currentData['iMax']=row['iMax']
        currentData['sMax']=int(row['sMax']/X)
        currentData['time'].append(row['time'])


plt.plot(currentData['x'],currentData['time'],label=r'Varying $j_{max}$ const. $i_{max}=40$',linewidth=2)



variationData=[]
currentData={'x':[],'y':[],'jMax':0,'sMax':0,'time':[]}

with open('data/varying_imax.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['sMax','iMax','jMax','S','V','time'],quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        currentData['x'].append(row['iMax'])
        currentData['y'].append(row['V'])
        currentData['jMax']=row['jMax']
        currentData['sMax']=int(row['sMax']/X)
        currentData['time'].append(row['time'])

plt.plot(currentData['x'],currentData['time'],label=r'Varying $i_{max}$ const. $j_{max}=100$',linewidth=2)

plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
plt.savefig('images/european_time.png',bbox_inches='tight', pad_inches=0.2)