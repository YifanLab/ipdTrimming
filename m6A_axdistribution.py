import re
from turtle import color
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
from lmfit.models import GaussianModel
import scipy.stats
from scipy.signal import savgol_filter
from statsmodels.nonparametric.kernel_regression import KernelReg
from scipy.signal import find_peaks
#from statsmodels.nonparametric.kernel_regression import kernelReg

liney=[]
linex=[]


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def ipddis(ipdcsv):
    linex=[]
    liney=[]
    with open(ipdcsv, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            liney.append(int(row[1]))
        #print(row[0])
            if re.findall(r'\((-?\d+.*),', row[0]):
                seobj = re.findall(r'\((-?.*),', row[0])
                linex.append(float(seobj[0]))
    
    return(linex, liney)



def gaussian(x, am, sigma, mu):
    y = am/(np.sqrt(2*np.pi*(sigma**2)))*np.exp(-0.5*(x-mu)**2/(sigma**2))
    return y



def modelpre(dframe):
    mf_rframe = dframe[dframe['x']>=1.6].reset_index()
    print(mf_rframe)
    modelr = GaussianModel()
    params_r = modelr.guess(mf_rframe['y'], x=mf_rframe['x'])
    result_r = modelr.fit(mf_rframe['y'], params_r, x=mf_rframe['x'])
    sigma_r = result_r.best_values['sigma']
    center_r = result_r.best_values['center']
    amp_r = result_r.best_values['amplitude']
    residen = abs(dframe[(dframe['x']<=center_r) & (dframe['x']>=0)]['y']-gaussian(dframe[(dframe['x']<=center_r) & (dframe['x']>=0)]['x'],amp_r,sigma_r,center_r)-gaussian(dframe[(dframe['x']<=center_r) & (dframe['x']>=0)]['x'],amp_r,sigma_r,center_r))
    resi_cross = dframe['x'][residen.idxmin()]
    #resi_cross = 2.15
    crosspoint = abs(dframe[(dframe['x']<=2) & (dframe['x']>=0)]['y']-gaussian(dframe[(dframe['x']<=2) & (dframe['x']>=0)]['x'],amp_r,sigma_r,center_r))
    cros2x = dframe['x'][crosspoint.idxmin()]
    cros2x=2
    fpr_t = sum(abs((dframe[(dframe['x']>=resi_cross) & (dframe['x']<=cros2x)]['y']-gaussian(dframe[(dframe['x']>=resi_cross) & (dframe['x']<=cros2x)]['x'],amp_r,sigma_r,center_r))))
    nor_t = sum(dframe[dframe['x']>=resi_cross]['y'])
    fnr = scipy.stats.norm(center_r,sigma_r).cdf(resi_cross)
    return(resi_cross, cros2x, fnr, fpr_t, nor_t, amp_r, sigma_r, center_r)



def lmodelpre(dframe):
    mf_rframe = dframe[dframe['x']<=0.5].reset_index()
    #print(mf_rframe)
    modelr = GaussianModel()
    params_r = modelr.guess(mf_rframe['y'], x=mf_rframe['x'])
    result_r = modelr.fit(mf_rframe['y'], params_r, x=mf_rframe['x'])
    sigma_r = result_r.best_values['sigma']
    center_r = result_r.best_values['center']
    amp_r = result_r.best_values['amplitude']
    #print(sigma_r, center_r, amp_r)
    residen = abs(dframe[(dframe['x']>=center_r) & (dframe['x']<=1)]['y']-gaussian(dframe[(dframe['x']>=center_r) & (dframe['x']<=1)]['x'],amp_r,sigma_r,center_r)-gaussian(dframe[(dframe['x']>=center_r) & (dframe['x']<=1)]['x'],amp_r,sigma_r,center_r))
    #print(residen)
    resi_cross = dframe['x'][residen.idxmin()]
    #print(resi_cross)
    #resi_cross = 1.4854
    crosspoint = abs(dframe[(dframe['x']>=-5) & (dframe['x']<=0)]['y']-gaussian(dframe[(dframe['x']>=-2) & (dframe['x']<=0)]['x'],amp_r,sigma_r,center_r))
    cros2x = dframe['x'][crosspoint.idxmin()]
    cros2x=-2
    fpr_t = sum(abs((dframe[(dframe['x']<=resi_cross) & (dframe['x']>=cros2x)]['y']-gaussian(dframe[(dframe['x']<=resi_cross) & (dframe['x']>=cros2x)]['x'],amp_r,sigma_r,center_r))))
    nor_t = sum(dframe[dframe['x']>=resi_cross]['y'])
    fnr = 1 - scipy.stats.norm(center_r,sigma_r).cdf(resi_cross)
    return(resi_cross, cros2x, fnr, fpr_t, nor_t, amp_r, sigma_r, center_r)









negx, negy = ipddis('./sb210_tetrahymena/Rep_2/allmacat_ipd.csv')
negymax = round(np.max(negy)*0.5)
#fig = plt.figure()
negdata = {'x':negx, 'y':negy}
negdata = pd.DataFrame(negdata)








#chrx, chry = ipddis('chrallipd.csv')
#chrymax = round(np.max(chry)*0.3)
#fig = plt.figure()
#chrdata = {'x':chrx, 'y':chry}
#chrdata = pd.DataFrame(chrdata)

posx, posy = ipddis('./allat_brduipd.csv')
posy_smooth = smooth(posy, 10)
ymax = round(np.max(posy)*0.5)
peaks, _ = find_peaks(posy_smooth, height=round(np.max(posy)*0.01))
peaks, _ = find_peaks(-posy_smooth)
print(peaks)
peakv = [posx[x] for x in peaks]
print(peakv)
#fig = plt.figure()
posdata = {'x':posx, 'y':posy}
posdata = pd.DataFrame(posdata)

kr = KernelReg(posy, posx, 'c')
posy_pre, posy_std = kr.fit(posx)


resi_cross, cros2x, fnr, fpr_t, nor_t, amp_r, sigma_r, center_r = modelpre(posdata)
fpt = fpr_t/nor_t
print(fnr, fpt)

#fig, (ax1, ax2) = plt.subplots(2,1)
#print(posdata['x'])
plt.scatter(posdata['x'],posdata['y']-gaussian(posdata['x'],amp_r,sigma_r,center_r),s=0.5, color='blue', alpha=0.5, label='residual')
#plt.scatter(posdata['x'],gaussian(posdata['x'],amp_r,sigma_r,center_r),s=1 ,label='best fit')
plt.plot(posdata['x'],gaussian(posdata['x'],amp_r,sigma_r,center_r),linewidth=0.5, color='red',label='best fit')
plt.vlines([resi_cross], 0, 2000000, linestyles=(0,(5,10)),linewidth=1, colors='black')
    
print(resi_cross)



#plt.set(ylim=(0, ymax))
plt.ylim(0,1500000)
plt.xlim(-2, 4)
#plt.xlim(0,1.8)


##left peak##
#tfpn= np.arange(resi_cross, 2, 1/100)
#plt.fill_between(x=tfpn,y1=gaussian(tfpn,amp_r,sigma_r,center_r), color='gray', alpha=0.2
#)



##left peak##
#tfnn = np.arange(-2, resi_cross, 0.05)
#plt.fill_between(
#    x=posdata['x'], y1=posdata['y']-gaussian(posdata['x'],amp_r,sigma_r,center_r),
#    where=(posdata['x']>0)&(posdata['x']<=resi_cross),
#    color='red',
#    alpha=0.2
#)


##right peak##
tfpn= np.arange(0,resi_cross+0.01,  1/100)
plt.fill_between(x=tfpn,y1=gaussian(tfpn,amp_r,sigma_r,center_r), color='red', alpha=0.2
)

tfnn = np.arange(resi_cross, 2, 0.05)
plt.fill_between(
    x=posdata['x'], y1=posdata['y']-gaussian(posdata['x'],amp_r,sigma_r,center_r),
    where=(posdata['x']>=resi_cross)&(posdata['x']<=2),
    color='gray',
    alpha=0.2
)


plt.style.use('seaborn-whitegrid')
#ax1.set_ylim(0, 400000)
#plt.scatter(negx, negy, s=2, label = 'AN')
#plt.scatter(chrx, chry, s=2, label = 'chr_AN')
plt.scatter(posx, posy, s=2, label = 'Brdu_1.5h Sb210', color='black', alpha=0.8)
#plt.scatter(negx, negy, s=1, label = 'AN')
#plt.plot(chrx, chry,  label = 'chr_AN')
#plt.scatter(posx, posy, s=1, color='red', label = 'AN')
#plt.scatter(posx, posy_pre, s=1, label = 'AN_fit')
#plt.xticks([])

#x,y = ipddis('allat_ipd.csv')
#ax2.vlines([resi_cross], 0, 4000000, linestyles=(0,(5,10)),linewidth=1, colors='black')
#ax2.scatter(x, y, s=1.0, label = 'Mic', color='black', alpha=0.6)
#x,y = ipddis('allac_ipd.csv')
#ax2.scatter(x, y, s=1.0, label = 'Mito', color='blue', alpha=0.6)
#x,y = ipddis('allac_ipd.csv')
#plt.plot(x, y, linewidth=2.0, label = 'AC')
#x,y = ipddis('allag_ipd.csv')
#plt.plot(x, y, linewidth=2.0, label = 'AG')

plt.legend()
plt.show()
