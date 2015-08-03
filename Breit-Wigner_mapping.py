import numpy as np
#import pylab as pyl
import matplotlib.pyplot as plt
import math
import random

# DEMONSTRATES THE BREIT-WIGNER TRANSFORM TO EFFICIENTLY SAMPLE FROM A BREIT-
# WIGNER DISTRIBTUION USING A SIMPLE VARIABLE TRANSFORMATION

mass = 500
width = 50

generated_s = []

histwidth = 8 * width
binwidth = 1
bins = []

# Number of generated points per bin
nperbin = 10000

N = nperbin * (histwidth / binwidth)

print "Generating " + str(N) + " invariant mass points.."

for i in range(0, histwidth / binwidth):
    bins.append((mass - 2 * width) + i * binwidth)

# random.random() generates [0,1] therefore all invariant masses are larger than the mass! If both tails of the width are to be generated random numbers [-pi/2,pi/2] are required using random.uniform(-pi/2,pi/2)

# Used to avoid the divergences in tan
epsilon = 0.1

for i in range(0, N):
    #s = mass**2 + mass*width*math.tan(random.random())
    s = mass ** 2 + mass * width * math.tan(
        random.uniform(-math.pi / 2 + epsilon, math.pi / 2 - epsilon))
    generated_s.append(math.sqrt(s))
    #if math.sqrt(s)>mass+width or math.sqrt(s)<mass-width:
    #    print "Point generated beyond the width!: ",math.sqrt(s)

plt.axis([mass - 2 * width, mass + 2 * width, 0, 10 * nperbin])
plt.hist(generated_s, bins=histwidth / binwidth)
plt.title('Breit-Wigner resonance with mass ' + str(mass) + ' GeV and width ' +
          str(width) + ' GeV')
plt.xlabel('Invariant mass (GeV)')
plt.ylabel('Event count')
plt.savefig('breit_wigner_peak.pdf')
plt.show()
