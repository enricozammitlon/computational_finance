import matplotlib.pyplot as plt
import numpy as np
import csv
from bisect import bisect_left,bisect_right
F=95.
C = 1.09
R=2.

allData=[]
with open('data/american_varying_s_beta_0_4.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'S':[],'V':[]}
    for row in reader:
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
    allData.append(currentData)

with open('data/no_put_american_varying_s_beta_0_4.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'S':[],'V':[]}
    for row in reader:
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
    allData.append(currentData)

plt.figure()
plt.grid()
parity = [(R*s*C)/F * 100 for s in allData[1]['S'] if (R*s*C)/F <2 ]
x=[ s for s in allData[0]['S'] if s <100]
parity=x
plt.plot(parity,allData[0]['V'][:len(parity)],label=r'American Option Embedded Put',linewidth=2)
plt.plot(parity,allData[1]['V'][:len(parity)],label=r'Embedded Put Removed',linewidth=2,linestyle='--')
equity = [ max(95.,2.*s) for s in allData[1]['S']]
lessC=bisect_right(allData[0]['V'], 100)
plt.annotate(r'$V(S_0=%.2f,t=0)=P_p$'%(parity[lessC]), xy=(parity[lessC],allData[0]['V'][lessC]), xytext=(22, 160),arrowprops=dict(arrowstyle="->"))
lessRS=0
for i in allData[0]['V']:
  if(allData[1]['S'][lessRS]*R>=i):
    break
  lessRS+=1
plt.annotate(r'$V(S_0=%.2f,t=0)=RS_0$'%(parity[lessRS]), xy=(parity[lessRS],allData[0]['V'][lessRS]), xytext=(65, 100),arrowprops=dict(arrowstyle="->"))
plt.plot(allData[1]['S'][:len(parity)],equity[:len(parity)],label=r'Max(F,R$S_0$)',linewidth=1)
plt.xlabel(r'$S_0$')
plt.ylabel(r'$V(S,t=0)$')
plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
plt.savefig('images/american_varying_s.png',bbox_inches='tight', pad_inches=0.2)

allData=[]
with open('data/american_varying_s_kappa_625.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'S':[],'V':[]}
    for row in reader:
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
    allData.append(currentData)
with open('data/american_varying_s_kappa_125.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'S':[],'V':[]}
    for row in reader:
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
    allData.append(currentData)
with open('data/american_varying_s_kappa_187.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'S':[],'V':[]}
    for row in reader:
        currentData['S'].append(row['S'])
        currentData['V'].append(row['V'])
    allData.append(currentData)

plt.figure()
plt.grid()
start=120
end=160
plt.plot(allData[0]['S'][start:end],allData[0]['V'][start:end],label=r'$ \kappa = 0.0625$',linewidth=2)
plt.plot(allData[1]['S'][start:end],allData[1]['V'][start:end],label=r'$ \kappa = 0.125$',linewidth=2)
plt.plot(allData[2]['S'][start:end],allData[2]['V'][start:end],label=r'$ \kappa = 0.1875$',linewidth=2)
plt.xlabel(r'$S_0$')
plt.ylabel(r'$V(S,t=0)$')
plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
plt.savefig('images/american_varying_s_varying_k.png',bbox_inches='tight', pad_inches=0.2)

plt.figure()
plt.grid()
start=100
end=200
plt.plot(allData[0]['S'][start:end],allData[0]['V'][start:end],label=r'$ \kappa = 0.0625$',linewidth=2)
plt.plot(allData[1]['S'][start:end],allData[1]['V'][start:end],label=r'$ \kappa = 0.125$',linewidth=2)
plt.plot(allData[2]['S'][start:end],allData[2]['V'][start:end],label=r'$ \kappa = 0.1875$',linewidth=2)
plt.xlabel(r'$S_0$')
plt.ylabel(r'$V(S,t=0)$')
plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
plt.savefig('images/complete_american_varying_s_varying_k.png',bbox_inches='tight', pad_inches=0.2)
