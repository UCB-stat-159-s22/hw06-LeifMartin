# LIGO-specific readligo.py 
from ligotools import readligo as rl
from ligotools import utils as utils
import unittest as unittest
r"""
____________________________________________________________________________________________

We need a lot of the variables defined in the Ligo notebook in order to test these functions,
hence we define all of them here:
____________________________________________________________________________________________
"""

#-- SET ME   Tutorial should work with most binary black hole events
#-- Default is no event selection; you MUST select one to proceed.
eventname = ''
eventname = 'GW150914' 
#eventname = 'GW151226' 
#eventname = 'LVT151012'
#eventname = 'GW170104'

# want plots?
make_plots = 1
plottype = "png"
#plottype = "pdf"

# Standard python numerical analysis imports:
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import h5py
import json

# the IPython magic below must be commented out in the .py file, since it doesn't work there.
#%matplotlib inline
#%config InlineBackend.figure_format = 'retina'
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# Read the event properties from a local json file
fnjson = "data/BBH_events_v3.json"
try:
    events = json.load(open(fnjson,"r"))
except IOError:
    quit()

# did the user select the eventname ?
try: 
    events[eventname]
except:
    print('You must select an eventname that is in '+fnjson+'! Quitting.')
    quit()


# Extract the parameters for the desired event:
event = events[eventname]
fn_H1 = 'data/'+event['fn_H1']              # File name for H1 data
fn_L1 = 'data/'+event['fn_L1']              # File name for L1 data
fn_template = 'data/'+event['fn_template']  # File name for template waveform
fs = event['fs']                    # Set sampling rate
tevent = event['tevent']            # Set approximate event GPS time
fband = event['fband']              # frequency band for bandpassing signal
#events['fn_H1']='data/H-H1_LOSC_4_V2-1126259446-32.hdf5'
r"""
________________________________

And here begins the first actual tests:
________________________________
"""
def test_loaddata_1():
	assert (isinstance(rl.loaddata(fn_H1, 'H1')[0], (np.ndarray, np.generic))), 'The first output of loaddata should be an array.'

def test_loaddata_2():
	assert (isinstance(rl.loaddata(fn_H1, 'H1')[1], (np.ndarray, np.generic))), 'The second output of loaddata should be an array.'

def test_loaddata_3():
	assert (isinstance(rl.loaddata(fn_H1, 'H1')[2], dict)), 'The third output of loaddata should be a dictionary.'

def test_loaddata_4():
	assert (len(rl.loaddata(fn_H1, 'H1'))==3), 'Loaddata should hade three outputs.'

r"""
___________________________________________________________________________________________________________________________

Now as loaddata is definitly 100% good we can use that to make more variables that we will need to test our other functions.
____________________________________________________________________________________________________________________________
"""

strain_H1, time_H1, chan_dict_H1 = rl.loaddata(fn_H1, 'H1')
strain_L1, time_L1, chan_dict_L1 = rl.loaddata(fn_L1, 'L1')


time = time_H1
dt = time[1] - time[0]

deltat = 5
indxt = np.where((time >= tevent-deltat) & (time < tevent+deltat))

# number of sample for the fast fourier transform:
NFFT = 4*fs
Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT = NFFT)
Pxx_L1, freqs = mlab.psd(strain_L1, Fs = fs, NFFT = NFFT)
# We will use interpolations of the ASDs computed above for whitening:
psd_H1 = interp1d(freqs, Pxx_H1)
psd_L1 = interp1d(freqs, Pxx_L1)
# Here is an approximate, smoothed PSD for H1 during O1, with no lines. We'll use it later.    
Pxx = (1.e-22*(18./(0.1+freqs))**2)**2+0.7e-23**2+((freqs/2000.)*4.e-23)**2
psd_smooth = interp1d(freqs, Pxx)

BNS_range = 1
if BNS_range:
    #-- compute the binary neutron star (BNS) detectability range

    #-- choose a detector noise power spectrum:
    f = freqs.copy()
    # get frequency step size
    df = f[2]-f[1]

    #-- constants
    # speed of light:
    clight = 2.99792458e8                # m/s
    # Newton's gravitational constant
    G = 6.67259e-11                      # m^3/kg/s^2 
    # one parsec, popular unit of astronomical distance (around 3.26 light years)
    parsec = 3.08568025e16               # m
    # solar mass
    MSol = 1.989e30                      # kg
    # solar mass in seconds (isn't relativity fun?):
    tSol = MSol*G/np.power(clight,3)     # s
    # Single-detector SNR for detection above noise background: 
    SNRdet = 8.
    # conversion from maximum range (horizon) to average range:
    Favg = 2.2648
    # mass of a typical neutron star, in solar masses:
    mNS = 1.4

    # Masses in solar masses
    m1 = m2 = mNS    
    mtot = m1+m2  # the total mass
    eta = (m1*m2)/mtot**2  # the symmetric mass ratio
    mchirp = mtot*eta**(3./5.)  # the chirp mass (FINDCHIRP, following Eqn 3.1b)

    # distance to a fiducial BNS source:
    dist = 1.0                           # in Mpc
    Dist =  dist * 1.0e6 * parsec /clight # from Mpc to seconds

    # We integrate the signal up to the frequency of the "Innermost stable circular orbit (ISCO)" 
    R_isco = 6.      # Orbital separation at ISCO, in geometric units. 6M for PN ISCO; 2.8M for EOB 
    # frequency at ISCO (end the chirp here; the merger and ringdown follow) 
    f_isco = 1./(np.power(R_isco,1.5)*np.pi*tSol*mtot)
    # minimum frequency (below which, detector noise is too high to register any signal):
    f_min = 20. # Hz
    # select the range of frequencies between f_min and fisco
    fr = np.nonzero(np.logical_and(f > f_min , f < f_isco))
    # get the frequency and spectrum in that range:
    ffr = f[fr]

    # In stationary phase approx, this is htilde(f):  
    # See FINDCHIRP Eqns 3.4, or 8.4-8.5 
    htilde = (2.*tSol/Dist)*np.power(mchirp,5./6.)*np.sqrt(5./96./np.pi)*(np.pi*tSol)
    htilde *= np.power(np.pi*tSol*ffr,-7./6.)
    htilda2 = htilde**2

    # loop over the detectors
    dets = ['H1', 'L1']
    for det in dets:
        if det is 'L1': sspec = Pxx_L1.copy()
        else:           sspec = Pxx_H1.copy()
        sspecfr = sspec[fr]
        # compute "inspiral horizon distance" for optimally oriented binary; FINDCHIRP Eqn D2:
        D_BNS = np.sqrt(4.*np.sum(htilda2/sspecfr)*df)/SNRdet
        # and the "inspiral range", averaged over source direction and orientation:
        R_BNS = D_BNS/Favg
        print(det+' BNS inspiral horizon = {0:.1f} Mpc, BNS inspiral range   = {1:.1f} Mpc'.format(D_BNS,R_BNS))

r"""
________________________________

And here begins the remaining tests:
________________________________
"""

#Utils tests:

def test_whiten():
	assert (isinstance(utils.whiten(strain_H1, psd_H1, dt),(np.ndarray, np.generic))), 'The output of whiten should be an array.'

def test_reqshift():
	"""
	This test assumes whiten is in order, hence why whiten itself is tested
	earlier in an attempt to make sure this test does not return a false pass.
	"""
	strain_H1_whiten = utils.whiten(strain_H1,psd_H1,dt)
	strain_L1_whiten = utils.whiten(strain_L1,psd_L1,dt)
	# We need to suppress the high frequency noise (no signal!) with some bandpassing:
	bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
	normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
	strain_H1_whitenbp = filtfilt(bb, ab, strain_H1_whiten) / normalization
	strain_L1_whitenbp = filtfilt(bb, ab, strain_L1_whiten) / normalization
	assert (isinstance(utils.reqshift(strain_H1_whitenbp),(np.ndarray, np.generic))), 'Reqshift should return an array.'

def test_write_wavfile():
	#unittest.assertRaises(TypeError,utils.write_wavfile(1,1,1)),'goog'
	assert(callable(utils.write_wavfile)), 'write_wavfile is a function, and should hence be callable.'

def test_plotmasta():
	#unittest.assertRaises(TypeError,utils.write_wavfile(1,1,1)),'goog'
	assert(callable(utils.plotmasta)), 'plotmasta is a function, and should hence be callable.'

def test_naive():
	"""
	Just to make sure..
	"""
	assert (1==1), 'I think you need to delete your computer and leave.'

