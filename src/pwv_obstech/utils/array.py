import numpy as np
from types import FunctionType as function

def tab_interp(tab0, tab, ycol, *, name, description=None, **kwargs):

    x0 = tab0['date'].mjd
    x = tab['date'].mjd
    y = tab[ycol]

    y0 = ma_interp(x0, x, y, **kwargs)
    name = f"{ycol}_{name}"
    tab0.add_column(y0, name=name)
    tab0.unit = y.unit
    tab0.description = y.description if description is None else description

def ma_interp(x0, x, y, xsep_max=0, left=np.nan, right=np.nan, **kwargs):

    x = np.asarray(x)
    x0 = np.asarray(x0)

    y0 = np.interp(x0, x, y, left=left, right=right, **kwargs)

    # mask nan values (out of boundaries)

    mask = np.isnan(y0)

    # mask values that are two distant to interpolating points

    if xsep_max > 0:
        xsep = abs(x0[:,None] - x).min(axis=1)
        mask = np.logical_or(mask, xsep > xsep_max)

    # mask values that are interpolated from a masked value 

    ymask = getattr(y, 'mask', False)

    if ymask is not False:
        ymask0 = np.interp(x0, x, ymask)
        mask = np.logical_or(mask, ymask0 > 0) 
    
    if any(mask):
        y0 = np.ma.masked_array(y0, mask=mask)

    return y0

def ma_bin(t, x, *, step=None, bins=None, max_masked_fraction=0.2):

    n = (bins == None) + (step == None)
    if n == 0 or n == 2:
        raise RuntimeError('Must give either step xor bins')

    if bins is None:
        tmin = int(min(t) / step) * step
        tmax = (int(max(t) / step) + 1) * step
        nsteps = int(0.5 + (tmax - tmin) / step) + 1
        bins = np.linspace(tmin, tmax, nsteps)

    dig = np.digitize(t, bins)

    t_bin = bins[:-1] + step / 2
    x_bin = np.ma.empty_like(t_bin)
    x_bin.mask = True
    for i in np.unique(dig):
        if i == 0 or i == len(t_bin):
            continue
        data = x[dig == i]
        if sum(data.mask) < max_masked_fraction * len(data):
            x_bin[i - 1] = data.mean()

    return t_bin, x_bin

