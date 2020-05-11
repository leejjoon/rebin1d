
# if True:
#     self.interp[name] = interpolate.interp1d(self.Spectra[name]['lambda'],
#                                              self.Spectra[name]['FLUX'],
#                                              kind='linear',
#                                              bounds_error=False, fill_value = 0.0)

#     wls, trans = self.I.pixel_to_band(row['X'], row['Y'], 
#                                       array=row['DETECTOR'], 
#                                       sparse=True)

#     # integrate the spectrum of the source over the SPHEREx bandpass
#     I = np.sum(trans*Catalog.interp[row['SOURCE_ID']](wls)) / np.sum(trans)

import astropy.io.fits as pyfits

f1 = pyfits.open("bd29d2091_mod_001.fits")
d = f1[1].data
wvl0, flux0 = d['WAVELENGTH'], d['FLUX']
flux_nu0 = flux0*wvl0*wvl0  # unit need to be fixed.
i1, i2 = wvl0.searchsorted([5500, 50000])
sl = slice(i1, i2)
wvl, flux_nu = wvl0[sl], flux_nu0[sl]


from scipy.interpolate import interp1d
intp = interp1d(wvl, flux_nu)
# R = 40
# n= 121
# w0 = 15050.
# dw = 0.5*w0/40.
# ww = np.linspace(w0-dw, w0+dw, n)

# clf()
# plot(wv, flux_nu)
# plot(ww, intp(ww))
# xlim(ww[0], ww[-1])
# fa = intp(ww).mean()
# plot(w0, fa, 'o')
# ii1 = wvl.searchsorted([ww[0], ww[-1]])
# sl1 = slice(*ii1)
# fb = flux_nu[sl1].mean()
# plot(w0, fb, 'x')

from SPHEREx_InstrumentSimulator import smile_lvf
import numpy as np

ia = 1

ix, iy = 1024 + np.zeros((2048, )), np.arange(2048)
lw0, lt_ = smile_lvf(ix, iy, array=ia, central_bandpass_only=True)
lw, lt = smile_lvf(ix, iy, array=ia)  # , central_bandpass_only=True)

simple_interpolation = np.sum(intp(lw*1.e4) * lt, axis=-1)

from intp_integrated import IntP, get_d

intp2 = IntP(wvl, flux_nu)

dd, lw_ = get_d(lw)
vv = intp2.rebin(lw_ * 1.e4) / dd / 1.e4 * lt
with_integration = vv.sum(axis=-1)

# i = 1024
# vv = []
# for lwi, lti in zip(lw, lt):
#     dd, lwi_ = get_d(lwi)
#     vv1 = intp2.rebin(lwi_ * 1.e4) / dd / 1.e4 * lt[i]
#     vv.append(vv1.sum())

clf()

plot(wvl0[sl]/1.e4, flux_nu0[sl])
plot(lw0, simple_interpolation)

plot(lw0, with_integration)

clf()
e_ratio = (simple_interpolation - with_integration)/with_integration

plot(lw0, e_ratio)

clf()
hist(e_ratio * 100., bins=np.linspace(-1.5, 0.3, 40))

# * lt[i] * 1.e4
