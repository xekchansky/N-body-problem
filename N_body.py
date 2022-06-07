import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

#G = 6.6743 * 10**(-15) #mass in kg, dist in meters, time in seconds
G = 1 #???

class dot:
    def __init__(self, mass=1, position=[0,0,0], velocity=[0,0,0]):
        self.mass = mass
        self.radius = ((3*mass)/(4*np.pi))**0.3
        self.position = np.array(position).astype(float)
        self.velocity = np.array(velocity).astype(float)
    
    def update_position(self):
        self.position += self.velocity
        
    def update_velocity(self, force, dt):
        #a = force / self.mass
        #v += at
        self.velocity += np.array(force) * dt / self.mass
        #print(self.velocity)
        
class group:
    def __init__(self, N=None, generator=None, dots=None):
        
        if N is not None:
            self.N = N
        elif dots is not None:
            self.N = len(dots)
        else:
            raise RuntimeError('group init error')
        
        self.dots = []
        if generator == None:
            self.dots=dots
        else:
            self.dots=generator.generate(self.N)
            
        self.history = []
        for i in range(self.N): self.history.append([])
        
    def update(self, max_t, max_iter):
        self.ts = np.linspace(0, max_t, max_iter)
        self.dt = self.ts[1] - self.ts[0]
        
        for t in self.ts:
            for i in range(len(self.dots)):
                for j in range(i+1, len(self.dots)):
                    d1 = self.dots[i]
                    d2 = self.dots[j]
                    update_dots(G, d1, d2, self.dt)
                    
            for i, dot in enumerate(self.dots):
                dot.update_position()
                self.history[i].append(dot.position.copy())
    
    def visualize_trajectories(self, colors=None):
        #sns.set(style = "darkgrid")
        fig = plt.figure()
        ax = fig.gca(projection='3d', adjustable='box')

        for i, dot_history in enumerate(self.history):
            ax.scatter(np.array(dot_history).T[0], np.array(dot_history).T[1], np.array(dot_history).T[2], c=colors[i])

        #ax.set_box_aspect([1,1,1])
        # ax.set_proj_type('ortho')
        set_axes_equal(ax) 
        plt.show()
        
    def visualize_dots(self, i, ax, colors=None):
        ax.clear()
        for c, dot_history in enumerate(self.history):
            if colors is not None and len(colors) == len(self.history):
                ax.scatter(np.array(dot_history).T[0][i], np.array(dot_history).T[1][i], np.array(dot_history).T[2][i], c=colors[c])
            else:
                ax.scatter(np.array(dot_history).T[0][i], np.array(dot_history).T[1][i], np.array(dot_history).T[2][i], s=self.dots[c].radius)
            _set_axes_radius(ax, [0,0,0], 200)
        
    def animate(self, colors=[]):
        #fig = plt.figure()
        #ax = fig.gca(projection='3d', adjustable='box')
        #ani = FuncAnimation(fig, self.visualize_dots, frames=1000, interval=10, repeat=True)
        #set_axes_equal(ax) 
        #plt.show()
        
        fig = plt.figure()
        ax = fig.gca(projection='3d', adjustable='box')
        ani = FuncAnimation(fig, self.visualize_dots, frames=1000, interval=5, repeat=True)
        plt.show(block=False)
        
class generator():
    def __init__(self, m_loc, m_scale, p_loc, p_scale, v_loc, v_scale):
        self.m_loc = m_loc
        self.m_scale = m_scale
        self.p_loc = p_loc
        self.p_scale = p_scale
        self.v_loc = v_loc
        self.v_scale = v_scale
        
    def generate(self, size):
        
        masses = np.abs(np.random.normal(loc=self.m_loc, scale=self.m_scale, size=size))
        
        xs = np.random.normal(loc=self.p_loc, scale=self.p_scale, size=size)
        ys = np.random.normal(loc=self.p_loc, scale=self.p_scale, size=size)
        zs = np.random.normal(loc=self.p_loc, scale=self.p_scale, size=size)
        coords = np.array([xs,ys,zs]).T
        
        v_xs = np.random.normal(loc=self.v_loc, scale=self.v_scale, size=size)
        v_ys = np.random.normal(loc=self.v_loc, scale=self.v_scale, size=size)
        v_zs = np.random.normal(loc=self.v_loc, scale=self.v_scale, size=size)
        velocities = np.array([v_xs,v_ys,v_zs]).T
        
        dots = []
        for i in range(size):
            d = dot(mass=masses[i], position=coords[i], velocity=velocities[i])
            dots.append(d)
        return dots
        
        
def update_dots(G, dot1, dot2, dt):
    vec = dot1.position - dot2.position
    r_pow_2 = sum((vec)**2)
    
    #check colision
    if r_pow_2 < (dot1.radius + dot2.radius)**2:
        F = 0 #(ignore colision)
    else:
        F = G*dot1.mass*dot2.mass/r_pow_2
        force_vec = (vec * F)/(r_pow_2**0.5)
        dot1.update_velocity(-force_vec, dt)
        dot2.update_velocity(force_vec, dt)
    
def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])
    
def animate_dots(g):
    fig = plt.figure()
    ax = fig.gca(projection='3d', adjustable='box')
    #g.visualize_dots(5)
    ani = FuncAnimation(fig, g.visualize_dots, frames=1000, interval=10, repeat=True)
    #set_axes_equal(ax)
    plt.show(block=False)