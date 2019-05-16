# Adapted from the vispy examples...
# Joao Couto - February 2016

from vispy import gloo
from vispy import app
from vispy.util.event import EmitterGroup

from glob import glob
import numpy as np
import os
import sys
import math
from ..io import *
from ..probe import *

# OPENGL stuff...
VERT_SHADER = """
#version 120
attribute vec2 a_position;
attribute float a_index;
varying float v_index;

attribute vec3 a_color;
varying vec3 v_color;

uniform vec2 u_pan;
uniform vec2 u_scale;

void main() {
    vec2 position_tr = u_scale * (a_position + u_pan);
    gl_Position = vec4(position_tr, 0.0, 1.0);
    v_color = a_color;
    v_index = a_index;
}
"""

FRAG_SHADER = """
#version 120
varying vec3 v_color;
varying float v_index;
void main() {
    gl_FragColor = vec4(v_color, 1.0);
    if ((fract(v_index) > .00001) && (fract(v_index) < .99999))
        gl_FragColor.a = 0.;
}
"""
keys_description = '''

Description of keyboard shorcuts:
    [F] - move forward on file
    [B] - move backwards on file
    [R] - Reset view
    [left arrow] - move view to the left
    [right arrow] - move view to the right
    [D] - spatial derivative
    [H] - highpass filter
    [M] - medial subtraction
    [K] - show keys
    [Q] - exit

'''

class BinaryViewer(app.Canvas):
    def __init__(self,fd, chtoplot = None,
                 srate = None, N = 30000, trace_separation = 0.3,
                 timeoffset = 0,highpass = False):
        app.Canvas.__init__(self, keys='interactive')
        self.size = (1200,800)
        self.title = 'binary file viewer'
        self.program = gloo.Program(VERT_SHADER, FRAG_SHADER)
        self.offset = trace_separation
        self.highpass = highpass
        self.mediansub = False
        self.fd = fd
        self.srate = srate
        self.chtoplot = chtoplot
        self.spatial_filter = False
        self.n = N
        self.m = len(self.chtoplot)
        self.ymean = None
        self.noffset = float(timeoffset*srate)
        print(self.n,self.m)
        self.data = np.zeros(int(self.n*self.m), dtype=[
            ('a_position', np.float32, 2),
            ('a_color', np.float32, 3),
            ('a_index', np.float32, 1),
        ])
        self.dark = False
        if self.dark:
            self.data['a_color'] = np.repeat(np.random.uniform(
                size=(self.m, 3),
                low=.0, high=.3),self.n, axis=0)
        else:
            self.data['a_color'] = np.repeat(np.random.uniform(
                size=(self.m, 3),
                low=.3, high=.9),self.n, axis=0)

        self.data['a_index'] = np.repeat(np.arange(self.m), self.n)

        self.read_from_binary()
        self.dbuf = gloo.VertexBuffer(self.data)
        self.program.bind(self.dbuf)

        self.program['u_pan'] = (0, 0.)
        self.program['u_scale'] = (1., 1.)

        gloo.set_viewport(1, 1, *self.physical_size)

        gloo.set_state(clear_color='white', blend=True,
                       blend_func=('src_alpha', 'one_minus_src_alpha'))
        self.show()

    def read_from_binary(self):
        n = int(self.n)
        m = int(self.m)
        sampleoffset = int(self.noffset)
        print('T = {0}'.format(sampleoffset/float(self.srate)))
        x = np.tile(np.linspace(-1., 1., n), m)
        y = np.array(self.fd[int(0)+sampleoffset:n + sampleoffset,
                             self.chtoplot.astype(int)].astype(np.float32)).T
        if self.ymean is None:
            self.ymean = np.repeat(np.mean(
                y,axis=1), y.shape[1]).reshape(y.shape)
            self.ymin = np.min(y[:])
            self.ymax = np.max(y[:])
        y -= self.ymean

        if self.highpass:
            import scipy.signal as signal
            b, a = signal.butter(3,[i/(self.srate/2.)
                                    for i in [300,0.95*self.srate/2.]], 'bandpass')
            for ch in range(y.shape[0]):
                y[ch,:] = signal.filtfilt(b, a, y[ch,:])

        if self.mediansub:
            ymedian = np.tile(np.median(
                y,axis=0), (y.shape[0],1))
            y -= ymedian
        if self.spatial_filter:
            y[0:-2,:] = np.diff(y,2,0)
            y[-2:,:] = np.zeros((2,y.shape[1]))
#        y = (y-self.ymin)/self.ymax
        y +=  np.repeat(np.arange(y.shape[0]),
                        y.shape[1]).reshape(y.shape)*self.offset


        # Update data
        self.data['a_position'] = np.zeros((n*m, 2), dtype=np.float32)
        self.data['a_position'][:, 0] = x.ravel()
        self.data['a_position'][:, 1] = .9*(y.ravel()/y.max()*2-1)
        srate = self.srate
    def on_resize(self, event):
        gloo.set_viewport(0, 0, *event.physical_size)

    def on_draw(self, event):
        if self.dark:
            gloo.clear(color=(0.8, 0.8, 0.8, 1.0))
        else:
            gloo.clear(color=(0.1, 0.1, 0.1, 1.0))
        self.program.draw('line_strip')

    def _normalize(self, x_y):
        x, y = x_y
        w, h = float(self.size[0]), float(self.size[1])
        return x/(w/2.)-1., y/(h/2.)-1.

    def on_mouse_move(self, event):
        if event.is_dragging:
            x0, y0 = self._normalize(event.press_event.pos)
            x1, y1 = self._normalize(event.last_event.pos)
            x, y = self._normalize(event.pos)
            dx, dy = x - x1, -(y - y1)
            button = event.press_event.button

            pan_x, pan_y = self.program['u_pan']
            scale_x, scale_y = self.program['u_scale']

            if button == 1:
                self.program['u_pan'] = (pan_x+dx/scale_x, pan_y+dy/scale_y)
            elif button == 2:
                scale_x_new, scale_y_new = (scale_x * math.exp(2.5*dx),
                                            scale_y * math.exp(2.5*dy))
                self.program['u_scale'] = (scale_x_new, scale_y_new)
                self.program['u_pan'] = (pan_x -
                                         x0 * (1./scale_x - 1./scale_x_new),
                                         pan_y +
                                         y0 * (1./scale_y - 1./scale_y_new))
            self.update()
        x,y = self._normalize(event.pos)
        self.currenttimepos = float(((x + 1.0)/2.0)*self.n+self.noffset)/float(self.srate)

    def on_mouse_wheel(self, event):
        dx = np.sign(event.delta[1])*.05
        scale_x, scale_y = self.program['u_scale']
        scale_x_new, scale_y_new = (scale_x * math.exp(2.5*dx),
                                    scale_y * math.exp(2.5*dx))
        self.program['u_scale'] = (scale_x_new, scale_y_new)
        self.update()

    def on_key_press(self, event):
        scale_x, scale_y = self.program['u_scale']
        pan_x, pan_y = self.program['u_pan']
        n = self.n
        if event.key=='F':
            self.noffset += n
            self.read_from_binary()
            self.dbuf.set_data(self.data)
        elif event.key=='K':
            print(keys_description)
        elif event.key=='B':
            self.noffset -= n
            self.read_from_binary()
            self.dbuf.set_data(self.data)
        elif event.key=='H':
            self.highpass = not self.highpass
            self.read_from_binary()
            self.dbuf.set_data(self.data)
        elif event.key == 'R':
            self.program['u_pan'] = (0, 0.)
            self.program['u_scale'] = (1., 1.)
        elif event.key=='Left':
            if 'Shift' in event.modifiers:
                dx = 0.4
            else:
                dx = 0.02
            dy = 0.
            self.program['u_pan'] = (pan_x+dx/scale_x, pan_y+dy/scale_y)
        elif event.key=='H':
            self.highpass = not self.highpass
            self.read_from_binary()
            self.dbuf.set_data(self.data)
        elif event.key=='M':
            self.mediansub = not self.mediansub
            self.read_from_binary()
            self.dbuf.set_data(self.data)
        elif event.key=='Right':
            if 'Shift' in event.modifiers:
                dx = -0.4
            else:
                dx = -0.02
            dy = 0.
            self.program['u_pan'] = (pan_x+dx/scale_x, pan_y+dy/scale_y)
        elif event.key=='Q':
            sys.exit()
        elif event.key=='D':
            self.spatial_filter = not self.spatial_filter
            self.read_from_binary()
            self.dbuf.set_data(self.data)

        else:
           pass
        self.update()
usagemsg = '''

plot-bin [-h] [--chmap CHANNEL_MAP]
         [--channels-to-discard CHANNELS_TO_DISCARD]
         [--srate SRATE] [--time-offset TIME_OFFSET]
         [--duration DURATION] [--offset OFFSET]
         [--nchannels NCHANNELS] [--dtype DTYPE]
         [--no-spikeglx] [--high-pass]
         filename
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
{0}
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
'''.format(keys_description)

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Plot binary file with electrophysiology data.', usage = usagemsg)
    parser.add_argument('fname',metavar='filename',
                        type=str,help='data filename',
                        nargs=1)
    parser.add_argument('--chmap',
                        action='store',
                        help='path to channel map file',
                        default=None,type=str,required=False)
    parser.add_argument('--channels-to-discard',
                        action='store',
                        help='channels to discard...',
                        default='',type=str,required=False)
    parser.add_argument('--srate',
                        action='store',
                        help='sampling rate (Hz)',
                        default=None,type=float,required=False)
    parser.add_argument('--time-offset',
                        action='store',
                        help='time offset (s)',
                        default=0,type=float,required=False)
    parser.add_argument('--duration',
                        action='store',
                        help='Duration to display (in seconds)',
                        default=0.5,type=float,required=False)
    parser.add_argument('--offset',
                        action='store',
                        help='Offset in the Y axis (Y units)',
                        default=80, type=float, required=False)
    parser.add_argument('--nchannels',
                        action='store',
                        help='number of channels',
                        default=None, type=int, required=False)
    parser.add_argument('--dtype',
                        action='store',
                        help='data type',
                        default='int16',type=str,required=False)
    parser.add_argument('--no-spikeglx',
                        action='store_false',
                        help='file is not recorded with spikeglx',
                        default=False)
    parser.add_argument('--high-pass',
                        action='store_true',
                        help='highpass data',
                        default=False)

    opts = parser.parse_args()
    # Handle filename
    if len(opts.fname):
        # For not take only one filename
        filename = os.path.realpath(opts.fname[0])
    else:
        print('Check arguments.')
        sys.exit(-1)
    if not os.path.isfile(filename):
        print('File not found: ' + filename)
        sys.exit(-1)

    # Handle channelmap
    channels_to_discard = []
    if len(opts.channels_to_discard):
        channels_to_discard = [int(i) for i in
                               opts.channels_to_discard.split(',')]
    if not opts.chmap is None:
        channelmap = os.path.realpath(opts.chmap)
        if not os.path.isfile(channelmap):
            print('File not found: ' + channelmap)
            sys.exit(-1)
        channelmap = read_channelmap_file(channelmap)
        chtoplot = chmap_get_active_eletrodes(channelmap)
    else:
        print('Plotting channels with no order...')
        chtoplot = None

    if not opts.no_spikeglx:
        fd,meta = load_spikeglx_binary(filename)
        srate = meta['sRateHz']
    else:
        if not opts.srate is None and not opts.nchannels is None:
            fd = load_dat(filename, opts.nchannels,dtype=opts.dtype)
            srate = opts.srate
        else:
            print('Missing some paramenters there (srate,nchannels)...')
            sys.exit(-1)
    if chtoplot is None:
        chtoplot = np.arange(fd.shape[1])
        channs = []
        for ch in chtoplot:
            if not ch in channels_to_discard:
                channs.append(ch)
        chtoplot = np.array(channs)
    # Load additional data from behaviour file.
    c = BinaryViewer(fd = fd,
                     srate=srate,
                     N=np.floor(srate*opts.duration),
                     chtoplot=chtoplot,
                     trace_separation = opts.offset,
                     timeoffset = opts.time_offset,
                     highpass = opts.high_pass)
    app.run()



if __name__ == '__main__':
    main()
