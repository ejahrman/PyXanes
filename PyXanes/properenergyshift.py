import numpy as np

def labview_calc_angle(xvalue, mmrev, xsteps, loffset, psi):
    rowlanddia = 1000
    netoffset = xvalue * mmrev / xsteps
    offsetmidplane = netoffset + loffset / np.cos(np.deg2rad(psi))
    sin2theta = offsetmidplane * 2 * np.cos(np.deg2rad(psi)) / rowlanddia
    theta = np.arcsin(sin2theta) / 2
    return 90 - np.rad2deg(theta)

def labview_bragg_angle(energy, material, h, k, l):
    hc = 12398.425
    lattice = crystaldata[material]['lattice']
    dspacing = lattice / np.sqrt(h**2 + k**2 + l**2)
    energy0 = hc / (2 * dspacing)
    angle = np.rad2deg(np.arcsin(energy0/energy))
    return angle

def labview_bragg_energy(angle, material, h, k, l):
    hc = 12398.425
    lattice = crystaldata[material]['lattice']
    dspacing = lattice / np.sqrt(h**2 + k**2 + l**2)
    energy0 = hc / (2 * dspacing)
    energy = energy0 / np.sin(np.deg2rad(angle))
    return energy

def labview_calc_steps(angle, mmrev, xsteps, loffset, psi):
    theta = np.deg2rad(90 - angle)
    sin2theta = np.sin(2 * theta)
    rowlanddia = 1000
    offsetmidplane = sin2theta / 2 / np.cos(np.deg2rad(psi)) * rowlanddia
    netoffset = offsetmidplane - loffset / np.cos(np.deg2rad(psi))
    xvalue = netoffset * xsteps / mmrev
    return xvalue

crystaldata = {
    'si':{
        '2d-111':6.2712 #x-ray data booklet
    },
    'ge':{
        '2d-111':6.532 #x-ray data booklet
    }
}
for k in crystaldata.keys():
    crystaldata[k]['lattice'] = crystaldata[k]['2d-111']/2*np.sqrt(3)
    
#overwrite si lattice with value from labview
crystaldata['si']['lattice'] = 5.43095
    
spectrometerconfig = {
    'mmrev':-2.54,
    'xsteps':10000,
    'loffset':86.36,
    'psi':40
}

def proper_energy_shift(spectrum, shiftenergy, amount, material, h, k, l):
    x, y = spectrum
    steps = labview_calc_steps(labview_bragg_angle(x, material=material, h=h, k=k, l=l),**spectrometerconfig)
    cursteps = labview_calc_steps(labview_bragg_angle(shiftenergy, material=material, h=h, k=k, l=l),**spectrometerconfig)
    targetsteps = labview_calc_steps(labview_bragg_angle(shiftenergy+amount, material=material, h=h, k=k, l=l),**spectrometerconfig)
    shift = targetsteps - cursteps
    steps = steps + shift
    return np.array([labview_bragg_energy(labview_calc_angle(steps,**spectrometerconfig), material=material, h=h, k=k, l=l), y])