import numpy as np
from tqdm import tqdm

from Law_Of_Motion import LOM_Force
from  BH_Tree import OctTree
    
#standart algorithm
def update_dots(dot1, dot2, dt):
    vec = dot1.position - dot2.position
    r = sum((vec)**2) ** 0.5
    
    #check colision
    if r < dot1.radius + dot2.radius:
        F = 0 #(ignore colision)
    else:
        F = LOM_Force(dot1.mass, dot2.mass, r)
        force_vec = (vec * F)/(r)
        dot1.update_velocity(-force_vec, dt)
        dot2.update_velocity(force_vec, dt)

class dot:
    def __init__(self, mass=1, position=[0,0,0], velocity=[0,0,0]):
        self.mass = mass
        self.radius = ((3*mass)/(4*np.pi))**0.3
        self.position = np.array(position).astype(float)
        self.velocity = np.array(velocity).astype(float)
    
    def update_position(self, dt):
        self.position += self.velocity * dt
        
    def update_velocity(self, force, dt):
        #a = force / self.mass
        #v += at
        self.velocity += np.array(force) * dt / self.mass
        
class group:
    def __init__(self, N=None, generator=None, dots=None):
        
        #init group size
        if N is not None:
            self.N = N
        elif dots is not None:
            self.N = len(dots)
        else:
            raise RuntimeError('group init error')
        
        #init dots
        self.dots = []
        if generator == None:
            self.dots=dots
        else:
            self.dots=generator.generate(self.N)
            
        #init history
        self.history = []
        for i in range(self.N): self.history.append([])
            
        #init center of mass
        self.COM_position, self.COM_velocity = self.calculate_COM()
        self.COM_history = []
        self.COM_history.append((self.COM_position, self.COM_velocity))
        self.true_COM_history = []
        self.true_COM_history.append((self.COM_position, self.COM_velocity))
        
    def update_standart(self, run_time, num_iterations):
        #updates every dot in group for given time duration and iterations with standart algorithm
        #input: 
        #    run_time - time duration for simulation
        #    num_iterations - number of iterarions to which split given time duration
        #return:
        #    self.dots.position
        #    self.history
        #    self.COM_history
        
        #create timestamps
        self.ts = np.linspace(0, run_time, num_iterations+1).astype(float)
        self.dt = self.ts[1] - self.ts[0]
        
        #iterate over timestamps
        for t in tqdm(self.ts):
            for i in range(len(self.dots)):
                for j in range(i+1, len(self.dots)):
                    d1 = self.dots[i]
                    d2 = self.dots[j]
                    update_dots(d1, d2, self.dt)
                    
            #update history of every dot
            for i, dot in enumerate(self.dots):
                dot.update_position(self.dt)
                self.history[i].append(dot.position.copy())
                
            #update history of COM
            self.COM_history.append(self.calculate_COM())
            
            #update true history of COM
            COM_init_position, COM_init_velocity = self.true_COM_history[0]
            true_COM_position = COM_init_position + (COM_init_velocity * (t + self.dt))
            self.true_COM_history.append((true_COM_position, COM_init_velocity))
            
    def update_barnes_hut(self, run_time, num_iterations, theta=0.0001):
        #updates every dot in group for given time duration and iterations with Barnes-Hut algorithm
        #input: 
        #    run_time - time duration for simulation
        #    num_iterations - number of iterarions to which split given time duration
        #return:
        #    self.dots.position
        #    self.history
        #    self.COM_history
        #create timestamps
        self.ts = np.linspace(0, run_time, num_iterations+1).astype(float)
        self.dt = self.ts[1] - self.ts[0]
        
        #iterate over timestamps
        for t in tqdm(self.ts):
            
            #get bounds for tree:
            m = max(self.get_max_coords()) + 1
            
            #build octo tree
            self.Q=OctTree(self.dots,m,-m,m,-m,m,-m)
            
            #update dots velocities
            for dot in self.dots:
                acceleration=self.Q.forces(dot,theta)
                dot.velocity += acceleration * self.dt
            
            #update history of every dot
            for i, dot in enumerate(self.dots):
                dot.update_position(self.dt)
                self.history[i].append(dot.position.copy())
                
            #update history of COM
            self.COM_history.append(self.calculate_COM())
            
            #update true history of COM
            COM_init_position, COM_init_velocity = self.true_COM_history[0]
            true_COM_position = COM_init_position + (COM_init_velocity * (t + self.dt))
            self.true_COM_history.append((true_COM_position, COM_init_velocity))
                
    def calculate_COM(self):
        #calculates current center of mass
        #input:
        #    self.dots
        #return: 
        #    COM_position - position of center of masses
        #    COM_velocity - velocity of center of masses
        
        COM_position = np.zeros(self.dots[0].position.shape)
        COM_velocity = np.zeros(self.dots[0].velocity.shape)
        total_mass = 0
        for dot in self.dots:
            COM_position += dot.mass * dot.position
            COM_velocity += dot.mass * dot.velocity
            total_mass += dot.mass
        COM_position /= total_mass
        COM_velocity /= total_mass
        return COM_position, COM_velocity
    
    def calculate_COM_residuals(self):
        #calculates differences between actual and true COM for every timestamp
        #input:
        #    self.COM_history
        #    self.true_COM_history
        #output:
        #    self.COM_residuals
        
        self.COM_position_residuals = []
        self.COM_velocity_residuals = []
        for i in range(len(self.ts)):
            self.COM_position_residuals.append(np.linalg.norm(self.COM_history[i][0] - self.true_COM_history[i][0]))
            self.COM_velocity_residuals.append(np.linalg.norm(self.COM_history[i][1] - self.true_COM_history[i][1]))
            
    def get_max_coords(self):
        max_x = 0
        max_y = 0
        max_z = 0
        for dot in self.dots:
            x, y, z = dot.position
            max_x = max(abs(x), max_x)
            max_y = max(abs(y), max_y)
            max_z = max(abs(z), max_z)
        return max_x, max_y, max_z
        
class generator:
    def __init__(self, m_loc, m_scale, p_loc, p_scale, v_loc, v_scale):
        #mass
        self.m_loc = m_loc
        self.m_scale = m_scale
        #position
        self.p_loc = p_loc
        self.p_scale = p_scale
        #velocity
        self.v_loc = v_loc
        self.v_scale = v_scale
        
    def generate(self, size):
        #generates array of dots with given size and gaussian distribution
        
        #generate masses
        masses = np.abs(np.random.normal(loc=self.m_loc, scale=self.m_scale, size=size))
        
        #generate positions
        xs = np.random.normal(loc=self.p_loc, scale=self.p_scale, size=size)
        ys = np.random.normal(loc=self.p_loc, scale=self.p_scale, size=size)
        zs = np.random.normal(loc=self.p_loc, scale=self.p_scale, size=size)
        coords = np.array([xs,ys,zs]).T
        
        #generate velocities
        v_xs = np.random.normal(loc=self.v_loc, scale=self.v_scale, size=size)
        v_ys = np.random.normal(loc=self.v_loc, scale=self.v_scale, size=size)
        v_zs = np.random.normal(loc=self.v_loc, scale=self.v_scale, size=size)
        velocities = np.array([v_xs,v_ys,v_zs]).T
        
        #construct dots
        dots = []
        for i in range(size):
            d = dot(mass=masses[i], position=coords[i], velocity=velocities[i])
            dots.append(d)
        return dots
        