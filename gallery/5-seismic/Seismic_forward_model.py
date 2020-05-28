import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as ipw
from fatiando.seismic import RickerWavelet, GaussianWavelet, FDElasticPSV, FDAcoustic2D

shape = (200, 600)
spacing = 5
extent = [0, shape[1]*spacing, shape[0]*spacing, 0]

vp1 = 1500
vp2 = 2000

increase_density = np.zeros(shape, dtype='float32') + 1000
increase_velocity = np.zeros(shape, dtype='float32') + vp1
increase_density[100:,:] = 1500
increase_velocity[100:,:] = vp2

p_increase = FDAcoustic2D(increase_velocity, increase_density, spacing=spacing)
fonte = (0, shape[1]//4)
p_increase.add_point_source(fonte, RickerWavelet(1, 60))

p_increase.run(1500)

p_increase.animate(every=20, embed=False, dpi=60, cutoff=0.4, fps=7)


def plot_with_p_rays_increasing(tempo, incidencia):
    fig = plt.figure()
    ax = plt.subplot(111)
    p_increase.snapshot(frame=tempo, ax=ax, cutoff=0.2, cmap='Greys')
    fig.set_size_inches(14, 5.5)
    y_bottom = shape[0] * spacing
    y_interface = 100 * spacing
    y_source = fonte[0] * spacing
    x_source = fonte[1] * spacing
    x_incidence = (np.tan(np.radians(incidencia)) * (y_interface - y_source)
                   + x_source)
    x_reflect = 2 * x_incidence - x_source
    arg = (vp2 / vp1) * np.sin(np.radians(incidencia))
    if arg <= 1:
        refract = np.arcsin(arg)
        x_refract = (np.tan(refract) * (y_bottom - y_interface) + x_incidence)
        ax.plot([x_incidence, x_refract], [y_interface, y_bottom], '-r', linewidth=2)
    ax.plot([x_source, x_incidence], [y_source, y_interface], '-k', linewidth=2)
    ax.plot([x_incidence, x_reflect], [y_interface, 0], '-b', linewidth=2)
    ax.hlines(y_interface, 0, shape[1] * spacing, colors='grey')


ipw.interactive(plot_with_p_rays_increasing,
                tempo=ipw.IntSlider(value=0, min=0, max=p_increase.it, step=50),
                incidencia=ipw.FloatSlider(value=45, min=0, max=90, step=0.5))

