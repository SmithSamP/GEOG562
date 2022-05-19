import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


# The parametrized function to be plotted

def fx(u, A, epsilon, h):
    return h*(np.cos(u)+epsilon)

def fy(u, A, epsilon, h):
    return h*A*np.sin(u)

def fA(eps):
    return np.sqrt((1-eps)/(1+eps))

u = np.linspace(0, 2*np.pi, 100)
epsilon = 0.2
A = np.sqrt((1-epsilon)/(1+epsilon))

h = 9

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()

fig.set_figwidth(7)
fig.set_figheight(10)

line, = plt.plot(fx(u,epsilon,A,h), fy(u,epsilon,A,h), lw=2)
ax.set_xlabel('Time [s]')
plt.axhline(y = 0, color="black", linestyle="--")
plt.axvline(x = 0, color="black", linestyle="--")
plt.xlim([-10, 10])
plt.ylim([-10, 10])

# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control the frequency.
axh = plt.axes([0.25, 0.1, 0.65, 0.03])
h_slider = Slider(
    ax=axh,
    label='Frequency [Hz]',
    valmin=-10,
    valmax=10,
    valinit=h,
)

# Make a vertically oriented slider to control the amplitude
axepsilon = plt.axes([0.1, 0.25, 0.0225, 0.63])
epsilon_slider = Slider(
    ax=axepsilon,
    label="Amplitude",
    valmin=0,
    valmax=1,
    valinit=epsilon,
    orientation="vertical"
)


# The function to be called anytime a slider's value changes
# def update(val):
#     line.set_ydata(f(t, epsilon_slider.val, h_slider.val))
#     fig.canvas.draw_idle()

def update(val):
    line.set_xdata(fx(u, epsilon_slider.val, fA(epsilon_slider.val), h_slider.val))
    line.set_ydata(fy(u, epsilon_slider.val, fA(epsilon_slider.val), h_slider.val))
    fig.canvas.draw_idle()


# register the update function with each slider
h_slider.on_changed(update)
epsilon_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    h_slider.reset()
    epsilon_slider.reset()
button.on_clicked(reset)

plt.show()