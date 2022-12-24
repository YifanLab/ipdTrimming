from math import inf, ceil
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
from lmfit import models
from lmfit.models import GaussianModel
import scipy.stats
from scipy.signal import savgol_filter
from statsmodels.nonparametric.kernel_regression import KernelReg
from scipy.signal import find_peaks
from matplotlib import gridspec
import random
import kneed
from scipy.signal import peak_widths



def generate_model(spec):
    composite_model = None
    params = None
    x = spec['x']
    y = spec['y']
    x_min = np.min(x)
    x_max = np.max(x)
    x_range = x_max - x_min
    y_max = np.max(y)
    for i, basis_func in enumerate(spec['model']):
        prefix = f'm{i}_'
        model = getattr(models, basis_func['type'])(prefix=prefix)
        if basis_func['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']: # for now VoigtModel has gamma constrained to sigma
            model.set_param_hint('sigma', min=1e-6, max=x_range)
            model.set_param_hint('center', min=x_min, max=x_max)
            model.set_param_hint('height', min=1e-6, max=1.1*y_max)
            model.set_param_hint('amplitude', min=1e-6)
            # default guess is horrible!! do not use guess()
            default_params = {
                prefix+'center': x_min + x_range * random.random(),
                prefix+'height': y_max * random.random(),
                prefix+'sigma': x_range * random.random()
            }
        else:
            raise NotImplemented(f'model {basis_func["type"]} not implemented yet')
        if 'help' in basis_func:  # allow override of settings in parameter
            for param, options in basis_func['help'].items():
                model.set_param_hint(param, **options)
        model_params = model.make_params(**default_params, **basis_func.get('params', {}))
        if params is None:
            params = model_params
        else:
            params.update(model_params)
        if composite_model is None:
            composite_model = model
        else:
            composite_model = composite_model + model
    return composite_model, params






def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))

def _2gaussian(x, am1,mu1,sigma1, am2,mu2,sigma2):
    y = am1/(np.sqrt(2*np.pi*(sigma1**2)))*np.exp(-0.5*(x-mu1)**2/(sigma1**2))+ am2/(np.sqrt(2*np.pi*(sigma2**2)))*np.exp(-0.5*(x-mu2)**2/(sigma2**2))
    #y = am/(np.sqrt(2*np.pi*(sigma**2)))*np.exp(-0.5*(x-mu)**2/(sigma**2))
    return y

def gaussian(x, am, sigma, mu):
    y = am/(np.sqrt(2*np.pi*(sigma**2)))*np.exp(-0.5*(x-mu)**2/(sigma**2))
    return y



def modelpre(dframe):
    mf_rframe = dframe[dframe['x']> -inf].reset_index()
    modelr = GaussianModel()
    params_r = modelr.guess(mf_rframe['y'], x=mf_rframe['x'])
    result_r = modelr.fit(mf_rframe['y'], params_r, x=mf_rframe['x'])
    sigma_r = result_r.best_values['sigma']
    center_r = result_r.best_values['center']
    amp_r = result_r.best_values['amplitude']
    #residen = abs(dframe[(dframe['x']<=center_r) & (dframe['x'])]['y']-gaussian(dframe[(dframe['x']<=center_r) & (dframe['x'])]['x'],amp_r,sigma_r,center_r)-gaussian(dframe[(dframe['x']<=center_r) & (dframe['x']>=0)]['x'],amp_r,sigma_r,center_r))
    #resi_cross = dframe['x'].residen.idxmin()]
    #crosspoint = abs(dframe['y']-gaussian(dframe['x'],amp_r,sigma_r,center_r))
    #cros2x = dframe['x'][crosspoint.idxmin()]
    #fpr_t = sum(dframe[(dframe['x']>=resi_cross) & (dframe['x']<=cros2x)]['y']-gaussian(dframe[(dframe['x']>=resi_cross) & (dframe['x']<=cros2x)]['x'],amp_r,sigma_r,center_r))
    #nor_t = sum(dframe[dframe['x']>=resi_cross]['y'])
    #fnr = scipy.stats.norm(center_r,sigma_r).cdf(resi_cross)
    #return(resi_cross, cros2x, fnr, fpr_t, nor_t, amp_r, sigma_r, center_r)
    return(amp_r, sigma_r, center_r)

posx=[-10.0, -9.8, -9.6, -9.4, -9.2, -9.0, -8.8, -8.6, -8.4, -8.2, -8.0, -7.8, -7.6, -7.4, -7.199999999999999, -7.0, -6.8, -6.6, -6.4, -6.199999999999999, -6.0, -5.8, -5.6, -5.3999999999999995, -5.199999999999999, -5.0, -4.8, -4.6, -4.3999999999999995, -4.199999999999999, -4.0, -3.8, -3.5999999999999996, -3.3999999999999995, -3.1999999999999993, -3.0, -2.8, -2.5999999999999996, -2.3999999999999995, -2.1999999999999993, -2.0, -1.799999999999999, -1.5999999999999996, -1.4000000000000004, -1.1999999999999993, -1.0, -0.7999999999999989, -0.5999999999999996, -0.3999999999999986, -0.1999999999999993, 0.0, 0.20000000000000107, 0.40000000000000036, 0.6000000000000014, 0.8000000000000007, 1.0, 1.200000000000001, 1.4000000000000004, 1.6000000000000014, 1.8000000000000007, 2.0, 2.200000000000001, 2.4000000000000004, 2.6000000000000014, 2.8000000000000007, 3.0, 3.200000000000001, 3.4000000000000004, 3.6000000000000014, 3.8000000000000007, 4.0, 4.200000000000001, 4.4, 4.600000000000001, 4.800000000000001, 5.0, 5.200000000000001, 5.4, 5.600000000000001, 5.800000000000001, 6.0, 6.199999999999999, 6.400000000000002, 6.600000000000001, 6.800000000000001, 7.0, 7.199999999999999, 7.400000000000002, 7.600000000000001, 7.800000000000001, 8.0, 8.2, 8.400000000000002, 8.600000000000001, 8.8, 9.0, 9.200000000000003, 9.400000000000002, 9.600000000000001, 9.8]
posy=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00010048231511254056, 0.00010048231511254056, 0.00010048231511254056, 0.00025120578778135003, 0.0004521704180064311, 0.000803858520900323, 0.001105305466237942, 0.0018589228295819962, 0.004421221864951757, 0.009093649517684894, 0.018890675241157596, 0.0405446141479099, 0.08490755627009655, 0.16001808681671995, 0.25859123794212224, 0.35776728295819976, 0.432475884244373, 0.47352290996784574, 0.48829381028938906, 0.48889670418006437, 0.4803054662379421, 0.45910369774919624, 0.4148914790996784, 0.33968046623794246, 0.24095659163987143, 0.14142885852090004, 0.06641881028938913, 0.024618167202572278, 0.007284967845659184, 0.002009646302250807, 0.0007033762057877793, 0.00025120578778135095, 0.00010048231511254011, 5.024115755627028e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#AGAGCATTGA [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.789272030651358e-05, 4.789272030651358e-05, 4.789272030651358e-05, 0.0002873563218390815, 0.0003831417624521078, 0.0005268199233716486, 0.000670498084291188, 0.0008620689655172423, 0.0009578544061302687, 0.0012452107279693502, 0.0018199233716475133, 0.002969348659003829, 0.004693486590038318, 0.007567049808429105, 0.016139846743295036, 0.031752873563218466, 0.06551724137931024, 0.12437739463601544, 0.21149425287356285, 0.3153256704980843, 0.40086206896551757, 0.45560344827586213, 0.4773946360153258, 0.4836206896551725, 0.47959770114942535, 0.4655651340996168, 0.43242337164750977, 0.3738984674329502, 0.2869252873563222, 0.18323754789272037, 0.09722222222222197, 0.04137931034482764, 0.01786398467432947, 0.008524904214559408, 0.003879310344827591, 0.0021551724137931017, 0.0013888888888888913, 0.0008620689655172419, 0.0006226053639846761, 0.0001915708812260539, 4.7892720306513154e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
posdata = {'x':posx, 'y':posy}
posdata = pd.DataFrame(posdata)
#print(posdata['x'])
#peaks, _ = find_peaks(posy, height=0.1)
#peaks, _ = find_peaks(-y)
#print(peaks)
#peakv = [posx[x] for x in peaks]
#print(peaks, peakv)

spec = {
    'x': posx,
    'y': posy,
    'model':[
        {'type': 'GaussianModel'},
        {'type': 'GaussianModel'}
    ]
}



def plot_to_blog(fig, figure_name):
    filename = os.path.expanduser(f'{image_dir}/{figure_name}')
    fig.savefig(filename)
    return filename
model, params = generate_model(spec)
output = model.fit(spec['y'], params, x=spec['x'])
print(output.best_values)
m0_amp = output.best_values['m0_amplitude']
m0_cen = output.best_values['m0_center']
m0_sigma = output.best_values['m0_sigma']
m1_amp = output.best_values['m1_amplitude']
m1_cen = output.best_values['m1_center']
m1_sigma = output.best_values['m1_sigma']

#output.plot(data_kws={'markersize': 1})
plt.scatter(posdata['x'],gaussian(posdata['x'],m0_amp,m0_sigma,m0_cen),s=2 ,label='best fit m0')
plt.scatter(posdata['x'],gaussian(posdata['x'],m1_amp,m1_sigma,m1_cen),s=2 ,label='best fit m1')
plt.scatter(posdata['x'],_2gaussian(posdata['x'], m0_amp,m0_cen,m0_sigma,m1_amp,m1_cen,m1_sigma),s=2 ,label='best fit m0+m1')
posym = _2gaussian(posdata['x'], m0_amp,m0_cen,m0_sigma,m1_amp,m1_cen,m1_sigma)
peaks, _ =find_peaks(posym, height=max(posym)*0.1)
_, _,_, posyw = peak_widths(posym, peaks, rel_height=1)
print(posyw[0])
kneedle = kneed.KneeLocator(posx[peaks[0]:ceil(posyw[0])], posy[peaks[0]:ceil(posyw[0])], S=1.0, curve="convex", direction="decreasing")
kneepos = kneedle.elbow
print(posyw)
print(peaks)
print(kneepos)

print(np.trapz(posy[0:70], posx[0:70], dx=1))
plt.plot(posx, posy)
plt.vlines(kneepos, 0, 0.1)

#fig.show()
#plot_to_blog(fig, 'xrd-fitting-two-gaussian-noise-lmfit-spec.png')

#popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_2gaussian, posx, posy, p0=[1,0,1,2,0,1], maxfev=5000)
#perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
#pars_1 = popt_2gauss[-inf:0]
#pars_2 = popt_2gauss[0:inf]
#gauss_peak_1 = _1gaussian(posx, *pars_1)
#gauss_peak_2 = _1gaussian(posx, *pars_2)

#fig = plt.figure(figsize=(4,3))
#gs = gridspec.GridSpec(1,1)
#ax1 = fig.add_subplot(gs[0])
#ax1.plot(posx, gauss_peak_1, "g")
#ax1.fill_between(posx, gauss_peak_1.min(), gauss_peak_1, facecolor="green", alpha=0.5)
  
#ax1.plot(posx, gauss_peak_2, "y")
#ax1.fill_between(posx, gauss_peak_2.min(), gauss_peak_2, facecolor="yellow", alpha=0.5)

#resi_cross, cros2x, fnr, fpr_t, nor_t, amp_r, sigma_r, center_r = modelpre(posdata)
#amp_r, sigma_r, center_r = modelpre(posdata)
#print(cros2x)
#plt.scatter(posdata['x'],posdata['y']-gaussian(posdata['x'],amp_r,sigma_r,center_r),s=1 ,label='residual')
#plt.scatter(posdata['x'],gaussian(posdata['x'],amp_r,sigma_r,center_r),s=2 ,label='best fit')
#plt.vlines([resi_cross], 0, max(posy),
#linestyles="dashed",linewidth=0.5, colors='black')
#plt.plot(posx, posy,  label = 'pos_AN')
plt.legend()
plt.show()