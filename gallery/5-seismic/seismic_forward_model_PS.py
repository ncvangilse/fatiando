import numpy as np
import matplotlib.pyplot as plt
# import ipywidgets as ipw
from fatiando.seismic import RickerWavelet, GaussianWavelet, FDElasticPSV, FDAcoustic2D
from matplotlib import animation
from fatiando.seismic.wavefd.utils import apply_damping, lame_lamb, lame_mu, xz2ps
import matplotlib.colors as colors


def animate_s(self, every=1, plottype=['wavefield'], cutoff=None,
            cmap=plt.cm.seismic, scale=1, every_particle=5,
            ax=None, interval=100, embed=False, blit=False,
            fps=10, dpi=70, writer='ffmpeg', **kwargs):
    nz, nx = self.shape
    mx, mz = nx * self.dx, nz * self.dz
    if ax is None:
        plt.figure(facecolor='white')
        ax = plt.subplot(111)
        # fig, axes = plt.subplots(2)
        # for ax in axes:
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.set_xlim(0, mx)
        ax.set_ylim(0, mz)
        ax.invert_yaxis()
    fig = ax.get_figure()
    wavefield = None
    particles = None
    vectors = None
    if 'wavefield' in plottype:
        extent = [0, mx, mz, 0]
        p = np.empty(self.shape, dtype=np.float32)
        s = np.empty(self.shape, dtype=np.float32)
        # imshow_args = dict(cmap=cmap, extent=extent)
        imshow_args = dict(extent=extent)
        # if cutoff is not None:
        #     imshow_args['vmin'] = -cutoff
        #     imshow_args['vmax'] = cutoff
        # wavefield_p = ax.imshow(np.zeros(self.shape), **imshow_args)
        wavefield = ax.imshow(np.zeros(self.shape), **imshow_args)
        fig.colorbar(wavefield, pad=0, aspect=30).set_label(
            'Divergence + Curl')
    if 'particles' in plottype or 'vectors' in plottype:
        xs = np.linspace(0, mx, nx)[::every_particle]
        zs = np.linspace(0, mz, nz)[::every_particle]
        x, z = np.meshgrid(xs, zs)
    if 'particles' in plottype:
        markersize = kwargs.get('markersize', 1)
        style = kwargs.get('style', '.k')
        particles, = plt.plot(x.ravel(), z.ravel(), style,
                              markersize=markersize)
    if 'vectors' in plottype:
        linewidth = kwargs.get('linewidth', 0.1)
        vectors = plt.quiver(x, z, np.zeros_like(x),
                             np.zeros_like(z),
                             scale=1 / scale, linewidth=linewidth,
                             pivot='tail', angles='xy',
                             scale_units='xy')
    # Check the aspect ratio of the plot and adjust figure size to match
    aspect = min(self.shape) / max(self.shape)
    try:
        aspect /= ax.get_aspect()
    except TypeError:
        pass
    if nx > nz:
        width = 10
        height = width * aspect * 0.8
    else:
        height = 8
        width = height * aspect * 1.5
    fig.set_size_inches(width, height)

    def plot(i):
        ax.set_title('iteration: {:d}'.format(i * every))
        ux, uz = self[i * every]
        if wavefield is not None:
            p = np.empty(self.shape, dtype=np.float32)
            s = np.empty(self.shape, dtype=np.float32)
            xz2ps(ux, uz, p, s, nx, nz, self.dx, self.dz)
            # wavefield.set_array(p + s)
            norm = colors.DivergingNorm(vmin=-cutoff, vmax=cutoff, vcenter=0)
            # p = np.clip(p, -cutoff, cutoff)
            # s = np.clip(s, -cutoff, cutoff)
            # p[p>cutoff] = cutoff
            # p[p<-cutoff] = -cutoff
            # s[s>cutoff] = cutoff
            # s[s<-cutoff] = -cutoff
            _p = plt.cm.bwr(norm(p))
            _s = plt.cm.BrBG(norm(s))
            wavefield.set_array(_p*_s)
            # print(_s)
        if particles is not None or vectors is not None:
            ux = ux[::every_particle, ::every_particle]
            uz = uz[::every_particle, ::every_particle]
        if particles is not None:
            particles.set_data(x.ravel() + scale * ux.ravel(),
                               z.ravel() + scale * uz.ravel())
        if vectors is not None:
            vectors.set_UVC(ux, uz)
        return wavefield, particles, vectors

    frames = self.simsize // every
    anim = animation.FuncAnimation(fig, plot, frames=frames, blit=blit,
                                   interval=interval)
    plt.show()
    return anim

FDElasticPSV.animate_s = animate_s

shape = (200, 600)
spacing = 5

vp1_solido = 3000
vp2_solido = 4000
vs1_solido = 2000
vs2_solido = 3000


shape2 = (300, 600)
l = 150
twosolid_density = np.zeros(shape2, dtype='float32') + 1800
twosolid_density[l:, :] = 2100
twosolid_vs = np.zeros(shape2, dtype='float32') + vs1_solido
twosolid_vs[l:, :] = vs2_solido
twosolid_vp = np.zeros(shape2, dtype='float32') + vp1_solido
twosolid_vp[l:, :] = vp2_solido


fonte_explosiva = (70, shape[1]//4)




ps_twosolid = FDElasticPSV(twosolid_vp, twosolid_vs, twosolid_density, spacing=spacing)
ps_twosolid.add_point_source(position=fonte_explosiva, dip=45, wavelet=GaussianWavelet(1, 100))

ps_twosolid.run(1200)
# ps_twosolid = FDElasticPSV.from_cache('C:\\Users\\nchvg\\repos\\geophysics\\SCPT\\FDElasticPSV-ttman7zy.h5')

ps_twosolid.animate_s(every=5, plottype=['wavefield',],
                    cutoff=1e-6, scale=1e7, every_particle=10,
                    dpi=60, fps=6, embed=False)


times = np.linspace(0, ps_twosolid.dt*ps_twosolid.simsize, ps_twosolid.simsize)
