"""
Seismic: 2D finite difference simulation of scalar wave propagation.

Difraction example in cylindrical wedge model. Based on:
R. M. Alford, K. R. Kelly and D. M. Boore -
Accuracy of finite-difference modeling of the acoustic wave equation.
Geophysics  1974
"""
import numpy as np
from matplotlib import animation
from fatiando.seismic import wavefd
from fatiando.vis import mpl
from fatiando.seismic.wavefd import FDAcoustic2D
from fatiando.seismic import RickerWavelet, GaussianWavelet

# Set the parameters of the finite difference grid
from matplotlib.colors import hsv_to_rgb

x_len = 10.
z_len = 10.
ds = 0.1
shape = (int(z_len / ds), int(x_len / ds))
# ds = 100.  # spacing
area = [0, x_len, 0, z_len]
# Set the parameters of the finite difference grid
velocity = np.zeros(shape) + 100.
# velocity[:50,:] = 100
velocity[25:40,:] = 500
# velocity[100:, 100:] = 0.
fc = 150. #Frequency of source
# sources = [wavefd.GaussSource(125 * ds, 75 * ds, area, shape,  1., fc)]
max_peak_amplitude = 1
amplitude = max_peak_amplitude * 1e4 / 8.1 #By experimentation, if a receiver is placed on top of a source, this factor will yield an amplitude equal to the max_peak_amplitude
sources = [wavefd.GaussianWavelet(x_len, 0.1, area, shape, amplitude, fc)]
dt = wavefd.scalar_maxdt(area, shape, np.max(velocity))
# dt = 0.05
duration = 0.5
maxit = int(duration / dt)
# stations = [[0, i*5] for i in range(6)]  # x, z coordinate of the seismometer
stations = [[0, .1]]  # x, z coordinate of the seismometer
snapshots = 3  # every 3 iterations plots one
simulation = wavefd.scalar(
    velocity, area, dt, maxit, sources, stations, snapshots)

# This part makes an animation using matplotlibs animation API
velocity_span = velocity.max() - velocity.min()
background = (velocity - velocity.min()) / velocity_span - 0.5
# background = velocity*0
fig = mpl.figure(figsize=(8, 6))
mpl.subplots_adjust(right=0.98, left=0.11, hspace=0.5, top=0.93)
mpl.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=3)

H = 0.34 * (velocity - velocity.min()) / velocity_span
if np.isnan(H).any():
    H =np.ones_like(velocity) * 0.34
S = np.ones_like(velocity) * 0.5
V = np.ones_like(velocity) * 0.5

HSV = np.dstack((H, S, V))
RGB = hsv_to_rgb(HSV)

wavefield = mpl.imshow(RGB[::-1], extent=area)


# wavefield = mpl.imshow(np.zeros_like(velocity), extent=area,
#                        cmap=mpl.cm.gray_r, vmin=-1, vmax=1)
mpl.points(stations, '^b', size=8)
mpl.ylim(area[2:][::-1])
mpl.xlabel('x (km)')
mpl.ylabel('z (km)')
mpl.m2km()
mpl.subplot2grid((4, 3), (3, 0), colspan=3)
seismogram1, = mpl.plot([], [], '-k')
mpl.xlim(0, duration)
rx_amplitude = max_peak_amplitude * 0.54 / np.sqrt(5**2+5**2)
mpl.ylim(-1*rx_amplitude, rx_amplitude)
mpl.ylabel('Amplitude')
times = np.linspace(0, dt * maxit, maxit)
# This function updates the plot every few timesteps


def animate(i):
    t, u, seismogram = next(simulation)
    seismogram1.set_data(times[:t + 1], seismogram[0][:t + 1])

    # H = np.zeros_like(velocity)
    u_span = u.max() - u.min()
    V = 0.5 + u/u_span

    # u_span = u.max() - u.min()
    # _u = u/u_span
    #
    #
    # log_u = np.log10(np.abs(_u))
    # # f = -2
    # # log_u[log_u<f] = f
    # log_u_span = log_u.max() - log_u.min()
    # log_u = (log_u - log_u.min()) / log_u_span
    # neg_mask = u <= 0
    # u = log_u.copy()
    # u[neg_mask] *= -1
    # V = 0.5 + 0.5 * u
    #
    # if ~np.isnan(V).any():
    #     pass
    # # V = 0.5 + 0.5 * u / np.abs(u).max()

    # V = (u-u.min()) / u_span
    print(V.min(), V.mean(), V.max())
    # S = np.ones_like(velocity) * 0.5
    # V = np.ones_like(velocity) * 0.5

    HSV = np.dstack((H, S, V))
    # HSV = np.dstack((V, S, H/0.34))
    RGB = hsv_to_rgb(HSV)

    wavefield.set_data(RGB[::-1])
    # print(u.max())
    return wavefield, seismogram1


anim = animation.FuncAnimation(
    fig, animate, frames=int(maxit / snapshots), interval=1)
mpl.show()

#
# """
# Seismic: 2D finite difference simulation of scalar wave propagation.
#
# Difraction example in cylindrical wedge model. Based on:
# R. M. Alford, K. R. Kelly and D. M. Boore -
# Accuracy of finite-difference modeling of the acoustic wave equation.
# Geophysics  1974
# """
# import numpy as np
# from matplotlib import animation
# from fatiando.seismic import wavefd
# from fatiando.vis import mpl
#
# # Set the parameters of the finite difference grid
# shape = (200, 200)
# ds = 100.  # spacing
# area = [0, shape[0] * ds, 0, shape[1] * ds]
# # Set the parameters of the finite difference grid
# velocity = np.zeros(shape) + 6000.
# velocity[100:, 100:] = 0.
# fc = 15.
# sources = [wavefd.GaussSource(125 * ds, 75 * ds, area, shape,  1., fc)]
# dt = wavefd.scalar_maxdt(area, shape, np.max(velocity))
# duration = 2.5
# maxit = int(duration / dt)
# stations = [[75 * ds, 125 * ds]]  # x, z coordinate of the seismometer
# snapshots = 3  # every 3 iterations plots one
# simulation = wavefd.scalar(
#     velocity, area, dt, maxit, sources, stations, snapshots)
#
# # This part makes an animation using matplotlibs animation API
# background = (velocity - 4000) * 10 ** -1
# fig = mpl.figure(figsize=(8, 6))
# mpl.subplots_adjust(right=0.98, left=0.11, hspace=0.5, top=0.93)
# mpl.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=3)
# wavefield = mpl.imshow(np.zeros_like(velocity), extent=area,
#                        cmap=mpl.cm.gray_r, vmin=-1000, vmax=1000)
# mpl.points(stations, '^b', size=8)
# mpl.ylim(area[2:][::-1])
# mpl.xlabel('x (km)')
# mpl.ylabel('z (km)')
# mpl.m2km()
# mpl.subplot2grid((4, 3), (3, 0), colspan=3)
# seismogram1, = mpl.plot([], [], '-k')
# mpl.xlim(0, duration)
# mpl.ylim(-200, 200)
# mpl.ylabel('Amplitude')
# times = np.linspace(0, dt * maxit, maxit)
# # This function updates the plot every few timesteps
#
#
# def animate(i):
#     t, u, seismogram = simulation.next()
#     seismogram1.set_data(times[:t + 1], seismogram[0][:t + 1])
#     wavefield.set_array(background[::-1] + u[::-1])
#     return wavefield, seismogram1
#
#
# anim = animation.FuncAnimation(
#     fig, animate, frames=maxit / snapshots, interval=1)
# mpl.show()