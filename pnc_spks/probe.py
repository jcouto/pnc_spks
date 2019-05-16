#!/usr/bin/env python
# Functions and scripts to generate probe files.
# Also functions to read channelmaps.
import numpy as np
import pandas as pd
import argparse

def read_channelmap_file(filename):
    ''' Reads a channelmap file
    Channel map files are csv files with a header on the first row.
'''
    chmap = pd.DataFrame.from_csv(filename, index_col=None)
    return chmap.sort_values(by=['Probe'],ascending=[True])

def chmap_get_active_eletrodes(chmap):
    return np.array(chmap.Data[~np.array(chmap.Disconnect,dtype=bool)], dtype=int)

def write_prb_file(probe, filename, nchannels = None, radius = None):
    prb = {}
    for ii,pp in enumerate(probe):
        prb[ii] = pp
    with open(filename,'w') as f:
        if not nchannels is None:
            f.write('total_nb_channels = {0}\n'.format(int(nchannels)))
        if not radius is None:
            f.write('radius = {0}\n'.format(int(radius)))
        f.write( 'channel_groups = ' + repr(prb) + '\n' )
    print('Probe written to file: '+ filename)

def make_probe(chmap, maxdist = 35, append_unconnected=True):
    '''
    Uses the probe geometry to create a probe (kwikteam)
    Some of this is inspired in Nick Steinmetz's matlab code.
    '''
    use = ~np.array(chmap.Disconnect,dtype=bool)
    chans = np.array(chmap.Data)
    shank = np.array(chmap.Shank)[use]
    chan = np.array(chmap.Data)[use]
    xcoord = np.array(chmap.X)[use]
    ycoord = np.array(chmap.Y)[use]

    [xx,yy] = np.meshgrid(xcoord,ycoord)
    alldist = ((xx-xx.T)**2+(yy-yy.T)**2)**.5
    [s1, s2] = np.meshgrid(shank, shank);

    probe = []
    for s in np.unique(shank):
        tmp = {}
        # "Nums" variables reflect channels in the data file
        # "Inds" variables reflect entries in the order you passed them
        chanNums = chan[shank==s]
        chanInds = np.where(shank==s)[0];
        nChans = np.sum(shank==s)
        xc = xcoord[shank==s]
        yc = ycoord[shank==s]
        tmp['channels'] = [c for c in chanNums]
        conn = (alldist<maxdist) & (s1==s) & (s2==s)
        [ind1, ind2] = np.where(conn);
        ind1 = ind1 - np.min(ind1)
        ind2 = ind2 - np.min(ind2)

        tmp['graph'] = [[chanNums[x],chanNums[y]] for x,y in zip(ind1[ind2>ind1],ind2[ind2>ind1])]
        tmp['geometry'] = {chanNums[i]:[xc[i],yc[i]] for i in range(len(chanNums))}
        # This is NOT working very well!!!!!
        # Known bug... does not connect the 2 ends....
        # Since it is not used anywhere other than in phy so I am not worrying about it...
        if append_unconnected:
            unconn = np.where(~use)[0]
            for ch in range(len(unconn)):
                surr = []
                mm = np.mod(unconn[ch],2)==0
                if not unconn[ch] == 0:
                    try:
                        lowerchans = np.where(use[unconn[ch]:-1])[0]
                        tmp0 = lowerchans[np.where(np.mod(lowerchans,2) == mm)[0][-1]]
                        tmp1 = lowerchans[np.where(np.mod(lowerchans,2) != mm)[0][-1]]
                        if tmp0 in tmp['channels']:
                            surr.append(tmp0)
                        if tmp1 in tmp['channels']:
                            surr.append(tmp1)
                    except:
                        # In the bottom of probe
                        pass

                if not unconn[ch]==len(use):
                    try:
                        upperchans = np.where(use[:unconn[ch]])[0]
                        tmp0 = upperchans[np.where(
                            mod(upperchans,2) == mm)[0][0]] + unconn[ch]
                        tmp1 = upperchans[np.where(
                            mod(upperchans,2) != mm)[0][0]] + unconn[ch]
                        if tmp0 in tmp['channels']:
                            surr.append(tmp0)
                        if tmp1 in tmp['channels']:
                            surr.append(tmp1)
                    except:
                        # On the top of the probe
                        pass
                surr = np.array(surr)[(surr > np.min(chan)) & (surr < np.max(chan))]
                surr = [ i if (not i in unconn) else [] for i in surr]
                for c1 in range(len(surr)):
                    for c2 in range(c1+1,len(surr)):
                        toappend = [chans[surr[c1]], chans[surr[c2]]]
                        toappend_flipped = [chans[surr[c2]], chans[surr[c1]]]
                        if (not toappend in tmp['graph'] and
                            not toappend_flipped in tmp['graph']):
                            tmp['graph'].append(toappend)
        probe.append(tmp.copy())
    return probe

def plot_probe_sites(prb,site_distance=0):
    ''' Plot probe sites and connections.
    Joao Couto - February 2016
    '''
    import pylab as plt
    fig = plt.figure(figsize = (len(prb)*3.5,10),facecolor='w')
    axes = []
    for s,probe in enumerate(prb):
        if len(axes) == 0:
            axes.append(plt.subplot(1,len(prb),s+1))
        else:
            axes.append(plt.subplot(1,len(prb),s+1,sharex=axes[-1],sharey = axes[-1]))
        ax = axes[-1]
        for ii,chan in enumerate(probe['channels']):
            xy = probe['geometry'][chan]
            plt.plot(xy[0],xy[1],'ko',markersize=20,markerfacecolor=[.7,.7,.7],
                     markeredgecolor=[.3,.3,.3])
        for conn in probe['graph']:
            orig = probe['geometry'][conn[0]]
            dest = probe['geometry'][conn[1]]
            plt.plot([orig[0],dest[0]],[orig[1],dest[1]],color=[.8,.3,.3])

        for ii,chan in enumerate(probe['channels']):
            xy = probe['geometry'][chan]
            plt.text(xy[0]+site_distance*.05,
                     xy[1]+site_distance*.05,
                     '{0}'.format(chan),
                     color='k')
        ax.axis('tight')
        tmp = ax.get_xlim()
        ax.set_xlim([tmp[0]-np.abs(np.diff(tmp))*0.2,tmp[1]+np.abs(np.diff(tmp))*0.2])
        tmp = ax.get_ylim()
        ax.set_ylim([tmp[0]-np.abs(np.diff(tmp))*0.2,tmp[1]+np.abs(np.diff(tmp))*0.2])
        ax.grid(True)
    ax = axes[0]
    ax.set_xlabel('X micrometer')
    ax.set_ylabel('Y micrometer')


def script_plot_channelmap():
    # TO DO: Add more parameters to control site size and so on...
    parser = argparse.ArgumentParser(description='''Helper to visualize channelmap files.''')
    parser.add_argument('filename',type=str,help='''channelmap file (read docs)''')
    parser.add_argument('-d','--max-distance',type=float,
                        help='max distance to connect sites (irrelevant)',
                        default=50)
    args = parser.parse_args()

    chmap = read_channelmap_file(args.filename)
    prb = make_probe(chmap, maxdist = args.max_distance, append_unconnected=True)
    plot_probe_sites(prb)
    import pylab as plt
    plt.show()


def script_create_prb():
    parser = argparse.ArgumentParser(description='''Creates a prb file from a channelmap.''')
    parser.add_argument('filename',type=str,help='''channelmap file (read docs)''')
    parser.add_argument('-d','--max-distance',type=float,
                        help='max distance to connect sites (irrelevant)',
                        default=50)
    parser.add_argument('-r','--radius',type=float,
                        help='radius for Circus',
                        default=100)
    args = parser.parse_args()

    chmap = read_channelmap_file(args.filename)
    probe = make_probe(chmap, maxdist = args.max_distance)
    nchannels = np.max(chmap.Data) + 1

    write_prb_file(probe,args.filename.replace('.channelmap','.prb'),
                   nchannels,args.radius)

