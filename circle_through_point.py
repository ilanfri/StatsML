# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

# Given any 2D Cartesian input coordinate (xinput, yinput) compute and plot
#       a semi-circle which goes both through the point and the origin
#       in the (+,+) quadrant


def yfunc(x, x0):
    return np.sqrt(x0 ** 2 - (x - x0) ** 2)


def radius(xinput, yinput):
    xinput = float(xinput)
    yinput = float(yinput)
    return 1 / (2 * xinput) * (xinput ** 2 + yinput ** 2)


xinput = 2
yinput = 3
r = radius(xinput, yinput)
print "Radius:", r

xvalues = np.linspace(0, 2 * r, 100)

fig = plt.figure()
plt.plot(xvalues, yfunc(xvalues, r))
plt.plot(xinput, yinput, 'o')
plt.show()

# DERIVATION OF RELEVANT EQUATIONS:
# The equation of a circle with centre at (x0,y0) is
# (x-x0)**2 + (y-y0)**2 = r**2
# But our circle is always centred on the x-axis, so y0 = 0
# (x-x0)**2 + y**2 = r**2                    (***)
# Applying the constraint that the circle go through the origin (0,0) gives
# r = x0
# Using the substitution r=x0 on equation (***) and applying the constraint that the circle must
#      also go through our specified point X = (x1,y1) gives
# (x1 - x0)**2 + y1**2 = x0**2
#    and solving this for x0 gives
# x0 = 1/(2*x1) * (x1**2 + y1**2)
# This is both our x-coordinate centre for the circle, as well as its radius (since r = x0)
# Inserting this result for x0 and r into equation (***) and solving for y gives us
#      the desired equation for the circle
