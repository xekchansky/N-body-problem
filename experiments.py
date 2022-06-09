import numpy as np
from time import time

from N_body import dot, group, generator

class Experiment:
    def __init__(self):
        self.generator = generator(m_loc=0, 
                                   m_scale=100,
                                   p_loc=0, 
                                   p_scale=100,
                                   v_loc=0, 
                                   v_scale=0.1)
        
    def run(self, N, run_time, num_iterations):
        #runs experiment
        #input: 
        #    N - number of bodies
        #    run_time - time duration for simulation
        #    num_iterations - number of iterarions to which split given time duration
        #return:
        #    position_residuals
        #    velocity_residuals
        #    timestamps
        #    elapsed_time
        
        #init group
        self.group = group(N=N, generator=self.generator)
        
        #run experiment
        start_time = time()
        #self.group.update_standart(run_time=run_time, num_iterations=num_iterations)
        self.group.update_barnes_hut(run_time=run_time, num_iterations=num_iterations, theta=0.0001)
        elapsed_time = time() - start_time
        
        #agregate metrics
        self.group.calculate_COM_residuals()
        
        return (np.array(self.group.COM_position_residuals),
                np.array(self.group.COM_velocity_residuals), 
                np.array(self.group.ts), 
                elapsed_time)