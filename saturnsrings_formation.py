"""
This project is a last year bachelor (L3) project for the computer course. 
The goal is to simulate the "formation" of Saturn's ring when adding particles of 
density similar to those surrounding Saturn -- mostly ice -- with other corresponding 
parameters (gravity, Saturn's mass, etc.). Another part of the project, mostly conducted 
by student Matis Melquiond is focused on Saturn's moon effects, using a lot of the "Rebound" 
module. If you use this code, think about changing file paths. The structure of this code 
has been inspired by several other Github codes involving moon formation around a jovian planet.

For more information, contact me on Github.

2022
"""

# IMPORTS
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib.patches import Circle
import rebound
from IPython.display import Image

# INITIALIZATION OF SIMULATION
sim = rebound.Simulation()        # simu
sim.units = ('km', 's', 'kg')     # units
sim.dt = 10                       # time step in corresponding units
sim.softening = 0.2               # softens small scales
sim.gravity    = "basic"          # chosen gravity
sim.add("699")       # add Saturn using its code 699
sim.move_to_com()    # move to new centre of mass
sim.status()         # print status of the simu

# SIMULATING RINGS

def powerlaw(slope, min_v, max_v): 
    """
    This powerlaw is meant to help picking a random particle radius
    """
    y = np.random.uniform()
    pow_max = pow(max_v, slope+1.)
    pow_min = pow(min_v, slope+1.)
    return pow((pow_max-pow_min)*y + pow_min, 1./(slope+1.))

particle_density = .1 # kg/m^3; particles are mostly ice
print('Adding ring particles')

def ring(ds):
    """
    Generating a "cloud" of particles at a chosen orbit. 
    The structure is similar to many codes generating moons. 
    We assume particles are spheres and we use use np.random.uniform()
    to add a number between -5000 and +5000 to the x and y coords to 
    spread particles out so they're not all on top of each other.
    """
    count = 0
    while count < 1000:
        radius = 50*powerlaw(slope=-4, min_v=1, max_v=4)/1000  # m get the radius
        mass = particle_density*4./3.*np.pi*(radius**3)  # kg
        rs = ds/2 # our particles will be 273,496/2 km away from the centre of Saturn
        theta = np.random.uniform(0,2*np.pi)
        # Converting coordinates from polar to cartesian
        x = rs*np.cos(theta)
        y = rs*np.sin(theta)
        x += np.random.uniform(low=-1, high=1)*5000 # spreading particles
        y += np.random.uniform(low=-1, high=1)*5000 
        G = 6.67428e-11 # gravitional constant
        # Converting the final velocity into km/s and times by .75 so particles fall towards Saturn
        v = np.sqrt(G * sim.particles[0].m/ (rs*1000)) / 1000 *.75# km/s;  v = sqrt(GM/r) particles speed
        # Add particles
        sim.add(
            m=mass, # mass of Saturn (stored in sim.particles[0].m)
            r=radius, # distance rs in m (not km !!)
            x=x,
            y=y,
            z=np.random.normal(),
            vx = -v*np.sin(theta),
            vy = v*np.cos(theta),
            vz = 0.)
        count += 1
    print('Finished adding ring particles')


def plotParticles(sim, k):
    """
    Takes an input simulation and integration number k
    Outputs a png of the current system with Saturn scaled to real size but particles enlarged
    """ 
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111, aspect='equal')
    ax.set_ylabel("y / km")
    ax.set_xlabel("x / km")
    ax.set_xlim(-200000, 200000)
    ax.set_ylim(-200000, 200000)
    ax.set_aspect('equal')
    ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
    ax.set_facecolor('black')
    for i, p in enumerate(sim.particles):
        if i == 0:
            fc, ec, a, r = "navajowhite", "None", 1, 58232
        else:
            fc, ec, a, r = "white", "None", 1, p.r*5000
        circ = Circle((p.x, p.y), r, facecolor=fc, edgecolor=ec, alpha=a)
        ax.add_patch(circ)
    plt.savefig('./myimagestest3/dynamics_'+str(k)+'.png', dpi=100)
    fig.clf()
    plt.close()

# EXECUTION (piling up images, computing time may increase as you add new rings)

# Ring A
ring(273496)#configurating new particles depending of ring diameter (so distance to the planet)
plotParticles(sim, 0)  # Add it on the figure : plot initial setup
for i in range(50): # jump 50 times
    sim.integrate(sim.t + 100) # gets the current sime time (sim.t) + jumps an extra 100 s
    plotParticles(sim, i+1) # plot again
sim.save("./myimagestest3/myfirstsim.bin") # save particles last state in simu

# Ring B
ring(253166)
plotParticles(sim, 0)
for i in range(50):
    sim.integrate(sim.t + 100)
    plotParticles(sim, i+1)
sim.save("./myimagestest3/myfirstsim.bin")

# Ring C
ring(184059) 
plotParticles(sim, 0)
for i in range(50):
    sim.integrate(sim.t + 100)
    plotParticles(sim, i+1)
sim.save("./myimagestest3/myfirstsim.bin")

# Ring D
ring(148983)
plotParticles(sim, 0)
for i in range(50):
    sim.integrate(sim.t + 100)
    plotParticles(sim, i+1)
sim.save("./myimagestest3/myfirstsim.bin")

# Ring E
ring(964288)
plotParticles(sim, 0)
for i in range(50):
    sim.integrate(sim.t + 100)
    plotParticles(sim, i+1)
sim.save("./myimagestest3/myfirstsim.bin")

# Ring F
ring(280367)
plotParticles(sim, 0)
for i in range(50):
    sim.integrate(sim.t + 100)
    plotParticles(sim, i+1)
sim.save("./myimagestest3/myfirstsim.bin")

# Ring G
ring(349554)
plotParticles(sim, 0)
for i in range(50):
    sim.integrate(sim.t + 100)
    plotParticles(sim, i+1)
sim.save("./myimagestest3/myfirstsim.bin")

"""
Then, you can just add those images to create a gif/video that shows the evolution/formation !!!
"""
