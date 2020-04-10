import matplotlib.pyplot as plt
import numpy as np
import csv


X=47.66
R=2
F=95

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
plt.plot(variationData[0]['x'],variationData[0]['y'],label=r'$\beta=0.486,\sigma=3.03,jMax=%i,sMax=%i$'%(variationData[0]['jMax'],variationData[0]['sMax']),linewidth=2)
plt.xlabel('iMax')
plt.ylabel(r'$V(S=X,t=0)$')
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
ax1.set_xlabel(r'sMax (multiples of X)')
ax1.set_ylabel(r'$V(S=X,t=0)$')
ax1.grid()
ax1.scatter(np.asarray(variationData[1]['x'][:20]),variationData[1]['y'][:20],label=r'$V(S=X,t=0)$ for $\beta=0.486,\sigma=3.03,iMax=%i$'%(variationData[1]['iMax']))
ax2 = ax1.twinx()
ax2.set_ylabel(r'$jMax$')
fig.tight_layout()
ax2.plot(np.asarray(variationData[1]['x'][:20]),variationData[1]['jMax'][:20],label=r'$jMax$',color="orange")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2,loc='lower right',fancybox=False, framealpha=0.0)
plt.savefig('images/european_varying_smax_zoomed.png',bbox_inches='tight', pad_inches=0.2)

fig, ax1 = plt.subplots()
ax1.set_xlabel(r'sMax (multiples of X)')
ax1.set_ylabel(r'$V(S=X,t=0)$')
ax1.grid()
ax1.scatter(np.asarray(variationData[1]['x']),variationData[1]['y'],label=r'$V(S=X,t=0)$ for $\beta=0.486,\sigma=3.03,iMax=%i$'%(variationData[1]['iMax']))
ax2 = ax1.twinx()
ax2.set_ylabel(r'$jMax$')
fig.tight_layout()
ax2.plot(np.asarray(variationData[1]['x']),variationData[1]['jMax'],label=r'$jMax$',color="orange")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2,loc='lower right',fancybox=False, framealpha=0.0)
plt.savefig('images/european_varying_smax.png',bbox_inches='tight', pad_inches=0.2)

'''
for smax in range(10,101):
    currentData={'x':[],'y':[],'iMax':0,'sMax':0}
    with open('data/smax_jmax/'+str(smax)+'_varying_jmax.csv', newline='\n') as csvfile:
        reader = csv.DictReader(csvfile,fieldnames=['sMax','iMax','jMax','S','V'],quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            currentData['x'].append(row['jMax'])
            currentData['y'].append(row['V'])
            currentData['iMax']=row['iMax']
            currentData['sMax']=int(row['sMax']/X)

    plt.figure()
    plt.plot(currentData['x'],currentData['y'],label=r'$\beta=0.486,\sigma=3.03,iMax=%i,sMax=%i$'%(currentData['iMax'],currentData['sMax']),linewidth=2)
    plt.xlabel('jMax')
    plt.ylabel(r'$V(S=X,t=0)$')
    plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
    plt.grid()
    plt.savefig('images/smax_jmax/'+str(smax)+'_european_varying_jmax.png',bbox_inches='tight', pad_inches=0.2)
    plt.close()
'''
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
plt.plot(allData[1]['S'],np.ones(len(allData[1]['S'])) * F,label=r'Bond Only',linewidth=2)
plt.plot(allData[1]['S'],np.asarray(allData[1]['S'])*R,label=r'Options Only',linewidth=2)

plt.xlabel('S')
plt.ylabel(r'$V(S,t=0)$')
plt.legend(loc='upper center',fancybox=False, framealpha=0.0)
plt.savefig('images/european_varying_s.png',bbox_inches='tight', pad_inches=0.2)