import matplotlib.pyplot as plt 
import csv

with open('test.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['r','p','v'],quoting=csv.QUOTE_NONNUMERIC)
    allData={'r':[],'p':[],'v':[]}
    for row in reader:
        allData['r'].append(row['r'])
        allData['p'].append(row['p'])
        allData['v'].append(row['v'])
    plt.figure()
    plt.plot(allData['r'],allData['p'],label=r'$P(r,t=0,T=2)$',linewidth=2)
    plt.plot(allData['r'],allData['v'],label=r'$V(r,t=0,T=2)$',linewidth=2)
    plt.xlabel('r')
    plt.ylabel('P and V')
    plt.legend(loc='center right',fancybox=False, framealpha=0.0)
    plt.savefig('Solution/plot.png',bbox_inches='tight', transparent="True", pad_inches=0)