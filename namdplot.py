import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as sc
import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument('-l',metavar="log",type=str,default='cryst.log',help="Name of namd log file")
args=parser.parse_args()

log = args.l
df = pd.DataFrame()
with open(log) as f:
    for line in f:
        tokens=line.split()
        if len(tokens)>0:
            if tokens[0]=='ETITLE:' and df.empty:
                df = pd.DataFrame(columns=tokens[1:])
                i=0
            if tokens[0]=='ENERGY:':
                df.loc[i]=[float(x) for x in tokens[1:]]
                i+=1
            if tokens[0]=='Info:':
                if ' '.join(tokens[1:3]) == 'TOTAL MASS':
                    mass_amu = float(tokens[4])

fac=mass_amu*sc.physical_constants['atomic mass unit-kilogram relationship'][0]*1.e27

df['DENSITY']=fac/df['VOLUME']
ts=2
fig,ax=plt.subplots(2,1,figsize=(6,10),sharex=True)
ax[0].plot(df['TS']*ts/1.e6,df['DENSITY'])
ax[0].set_ylabel('Density (g/cm$^3$)')
ax[1].plot(df['TS']*ts/1.e6,df['POTENTIAL'])
ax[1].set_ylabel('Potential Energy (kcal/mol)')
ax[1].set_xlabel('Time (ns)')
ax[0].set_ylim([1.0,1.6])
plt.savefig('rho-e.png')

