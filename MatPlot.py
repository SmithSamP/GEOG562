import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button

# Prepare arrays x, y, z
u = np.linspace(0, 2*np.pi, 100)
epsilon = 0.0
A = np.sqrt((1-epsilon)/(1+epsilon))
h = 9
x =  h*(np.cos(u)+epsilon)
y = h*A*np.sin(u)

f = plt.figure()
f.set_figwidth(7)
f.set_figheight(10)



plt.plot(x, y, label='parametric elipse')
plt.subplots_adjust(left=0.25, bottom=0.25)
plt.axhline(y = 0, color="black", linestyle="--")
plt.axvline(x = 0, color="black", linestyle="--")
plt.xlim([-10, 10])
plt.ylim([-10, 10])


# axes = plt.axes([0.25, 0.1, 0.65, 0.03])
# h_slider = Slider(
#     ax=axes,
#     label='h',
#     valmin=0,
#     valmax=30,
#     valinit=3,
# )
#
# def update(val):
#     line.set_ydata(f(t, epsilon_slider.val, h_slider.val))
#     plt.canvas.draw_idle()
#
#
# # register the update function with each slider
# slider.on_changed(update)
plt.show()