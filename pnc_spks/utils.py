import numpy as np

def unpackbits(x,num_bits = 16):
    '''
    unpacks numbers in bits.
    '''
    xshape = list(x.shape)
    x = x.reshape([-1,1])
    to_and = 2**np.arange(num_bits).reshape([1,num_bits])
    return (x & to_and).astype(bool).astype(int).reshape(xshape + [num_bits])

from numpy.lib.stride_tricks import as_strided

def _check_arg(x, xname):
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('%s must be one-dimensional.' % xname)
    return x

def xcorr(x, y, maxlag):
    """
    Cross correlation with a maximum number of lags.

    `x` and `y` must be one-dimensional numpy arrays with the same length.

    This computes the same result as
        numpy.correlate(x, y, mode='full')[len(a)-maxlag-1:len(a)+maxlag]

    The return vaue has length 2*maxlag + 1.

    Adapted from code from Warren Weckesser (stackoverflow).
    
    Joao Couto - January 2016
    """
    x = _check_arg(x, 'x')
    y = _check_arg(y, 'y')
    py = np.pad(y.conj(), 2*maxlag, mode='constant')
    T = as_strided(py[2*maxlag:], shape=(2*maxlag+1, len(y) + 2*maxlag),
                   strides=(-py.strides[0], py.strides[0]))
    px = np.pad(x, maxlag, mode='constant')
    return T.dot(px)

def convolve_with_kernel(fr,sigma = 5,dt = 1,kernel = 'exp'):
    ''' 
    Convolves with a kernel.
        - kernel can be exp or gauss if it is not a string you should provive a 1d vector smaller than fr (the kernel to convolve with).
    ''' 
    dt *= 1e-3
    sigma *= 1e-3
    xx = np.arange(-3*sigma,3*sigma,dt)
    if not type(kernel) is str:
        s = kernel
    elif kernel == 'gauss':
        s =  np.exp(-(xx)**2/(2*sigma)**2)/(np.sqrt(2.*np.pi)*sigma)
        s = s/np.sum(s)
    elif kernel == 'exp':
        s =  np.exp(-1. * np.sqrt(2)*np.abs(xx/sigma))/(np.sqrt(2.)*sigma)
        s = s/np.sum(s)
    return np.convolve(fr,s,'same')

def find_spikes(dat,thresh = None,wpre=16,wpos=20,threshstd=6):
    ''' Find spikes and extract sample times and waveforms '''
    tmp = np.zeros(shape=dat.shape)
    if thresh is None:
        thresh = compute_spike_threshold(dat,threshstd)
    tmp[dat<-thresh] = 1
    tstamps = np.where(np.diff(tmp)>0)[0]
    # align...
    for i,t in enumerate(tstamps):
        tmp = dat[t-wpre:t+wpos]
        tmpmax = np.argmin(tmp)
        tstamps[i] = t-wpre+tmpmax
        #extract waveforms
        waves = np.zeros(shape=(len(tstamps),wpre+wpos))
        for i,t in enumerate(tstamps):
            try:
                waves[i,:] = dat[t-wpre:t+wpos]
            except:
                print('Failed for spike {0}'.format(t))
    return tstamps,waves

def compute_spike_threshold(x,stdmin=4):
    ''' Compute spike threshold from filtered raw trace.
    Uses the formula from R. Quiroga, Z. Nadasdy, and Y. Ben-Shaul:
       thr = stdmin*sigma_n ,with
       sigma_n = median(|x|/0.6745)
 NOTE: 
   Default stdmin is 4.
   In the WaveClus batch scripts stdmin is set to 5.
    Joao Couto - January 2016    
    '''
    return stdmin * np.median(np.abs(x))/0.6745;

def whitening_matrix(x, fudge=1e-18):
    """
    wmat = whitening_matrix(dat, fudge=1e-18)
    Compute the whitening matrix.
        - dat is a matrix nsamples x nchannels
    Apply using np.dot(dat,wmat)
    Adapted from phy
    """
    assert x.ndim == 2
    ns, nc = x.shape
    x_cov = np.cov(x, rowvar=0)
    assert x_cov.shape == (nc, nc)
    d, v = np.linalg.eigh(x_cov)
    d = np.diag(1. / np.sqrt(d + fudge))
    w = np.dot(np.dot(v, d), v.T)
    return w


def alpha_function(N, amplitude = 1, t_rise = 2, t_decay = 250, srate = 1000.,norm = True):

    t_rise = t_rise/srate;
    t_decay = t_decay/srate;
    
    fun_max  = (t_rise*t_decay/(t_decay-t_rise)) * np.log(t_decay-t_rise);
    normalization_factor = 1; #%(exp(-fun_max/t_rise) - exp(-fun_max/t_decay))/(t_rise-t_decay);
    ii = np.arange(0,N)
    kernel = np.hstack([np.zeros(N),
                        amplitude*(1.0/(normalization_factor*(t_decay-t_rise))) * (np.exp(-((ii/srate)/t_decay))
                                                                                   - np.exp(-(ii/srate)/t_rise))])
    if norm:
        kernel /= np.sum(kernel)
    return kernel

def binary_spikes(spks,edges,sigma=None,kernel = None):
    ''' Create a vector of binary spikes.

    binsize = 0.001
    edges = np.arange(0,5,binsize)
    bspks = binary_spikes(spks,edges,5)/binsize

    Joao Couto - March 2016
    '''
    bins = [np.histogram(sp,edges)[0] for sp in spks]
    if not sigma is None:
        # Convolve spk trains
        x = np.arange(np.floor(-3*sigma),np.ceil(3*sigma))
        kernel = np.exp(-(x/sigma)**2/2)/(sigma*np.sqrt(2*np.pi))
    if not kernel is None:
        bins = [np.convolve(a,kernel,'same') for a in bins]
    return np.vstack(bins)

from scipy.interpolate import interp2d
from scipy.signal import ellip, filtfilt,butter

def bandpass_filter(X,srate,band=[3,300]):
    b, a = ellip(4, 0.1, 40, np.array(band)/(srate/2.),btype='bandpass')
    return filtfilt(b, a, X,axis = 0)#, method="gust"

def current_source_density(lfp,chmap, chanspacing=60, interpolate=False):
    # Interpolate so that we get even sampling
    selchannels = np.array(chmap.ichan)
    ux = np.unique(chmap.x)
    ix = np.argmax([np.sum(chmap.x==u) for u in ux])
    chidx = chmap.x==ux[ix]
    y = np.array(chmap[chidx].y)
    duration = lfp.shape[1]
    x = np.arange(duration)
    z = lfp[chmap[chidx].ichan,:]
    f = interp2d(x,y,z)
    ny = np.arange(np.min(y)-chanspacing,np.max(y)+chanspacing,chanspacing)
    nlfp = f(x,ny)
    # duplicate outmost channels
    csd = np.empty((nlfp.shape[0]-2,nlfp.shape[1]))
    smoothed_lfp = np.empty_like(nlfp)
    for i in range(csd.shape[0]):
        smoothed_lfp[i+1,:] = (1./4.) *(nlfp[i,:] + 2.*nlfp[i+1,:] + nlfp[i+2,:])
    smoothed_lfp[0,:] = (1./4.) *(3.*nlfp[0,:] + nlfp[1,:])
    smoothed_lfp[-1,:] = (1./4.) *(3.*nlfp[-1,:] + nlfp[-2,:])
    smoothed_lfp = smoothed_lfp
    for i in range(csd.shape[0]):
        csd[i,:] = -(1./(chanspacing*1.e-3)**2.)*(smoothed_lfp[i,:]-2.*smoothed_lfp[i+1,:]+smoothed_lfp[i+2,:])
    f = interp2d(x,np.linspace(np.min(y)-chanspacing,np.max(y)+chanspacing,csd.shape[0]),csd)
    ny = np.arange(np.min(y)-chanspacing,np.max(y)+chanspacing,5.)
    return f(x,ny),smoothed_lfp[:,1:-1]
