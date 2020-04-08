import scipy.stats as si
import numpy as np
import csv
import matplotlib.pyplot as plt

C = 1.09
alpha = 0.02
r = 0.0229
T = 2.
F = 95.
R = 2.
sigma = 0.416
K = 47.66
'''
T=3
C=0.106
alpha=0.01
r=0.0038
R=1
F=56
sigma = 0.369
K=56.47
'''
#Calulate value of coupon
#Through integrating Cexp(-(alpha+r)t)dt from 0 to T
COUPON = C/(alpha+r) * (1- np.exp((-(alpha+r)*T)))

BOND = F*np.exp(-r*T)

variationData=[]
with open('data/analytic.csv', newline='\n') as csvfile:
    reader = csv.DictReader(csvfile,fieldnames=['S','V'],quoting=csv.QUOTE_NONNUMERIC)
    currentData={'x':[],'y':[]}
    for row in reader:
        currentData['x'].append(row['S'])
        currentData['y'].append(row['V'])
    variationData.append(currentData)

def euro_vanilla_call(S, K, T, r, sigma):

    #S: spot price
    #K: strike price
    #T: time to maturity
    #r: interest rate
    #sigma: volatility of underlying asset

    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))

    call = (S * si.norm.cdf(d1, 0.0, 1.0) - K * np.exp(-r * T) * si.norm.cdf(d2, 0.0, 1.0))

    return call

ANALYTIC_PRICE = []
STOCK_PRICE = []
BOND_FLOOR = []
CONV_BOND = []
for s in range(1,140):
    STOCK_PRICE.append(s)
    ANALYTIC_PRICE.append(R*euro_vanilla_call(s, K, T, r, sigma) + BOND +COUPON)

#plt.plot(S,V1, label = " beta = 1")
#plt.plot(STOCK_PRICE,BOND_FLOOR, label = "Bond")
plt.plot(STOCK_PRICE,ANALYTIC_PRICE, label = "Analytic")
plt.plot(variationData[0]['x'],variationData[0]['y'], label = "Crank")
plt.xlabel('Stock price')
plt.ylabel('Eurocall Option')
plt.legend()
plt.savefig('images/analytic.png',bbox_inches='tight', pad_inches=0.2)
