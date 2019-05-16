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

#######################################################################################

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
