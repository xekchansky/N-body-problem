import pickle
import numpy as np
from time import time

from N_body import dot, group, generator

class Experiment:
    def __init__(self,seed=None):
        self.generator = generator(m_loc=0, 
                                   m_scale=100,
                                   p_loc=0, 
                                   p_scale=100,
                                   v_loc=0, 
                                   v_scale=0.1,
                                   seed=seed)
        
    def run(self, method, N, run_time, num_iterations, theta=0.0001, g=None, workers=1, time_test=False):
        #runs experiment
        #input: 
        #    method - method of computation
        #    N - number of bodies
        #    run_time - time duration for simulation
        #    num_iterations - number of iterarions to which split given time duration
        #    
        #return:
        #    position_residuals
        #    velocity_residuals
        #    timestamps
        #    elapsed_time
        
        #init group
        self.group = group(N=N, generator=self.generator)
        
        #run experiment
        start_time = time()
        if method == 'standart':
            self.group.update_standart(run_time=run_time, num_iterations=num_iterations)
        else:
            self.group.update_barnes_hut(run_time=run_time, num_iterations=num_iterations, theta=theta, workers=workers, time_test=time_test)
        elapsed_time = time() - start_time
        
        if time_test:
            return elapsed_time
        else:
            #agregate metrics
            self.group.calculate_COM_residuals()
            return (np.array(self.group.COM_position_residuals),
                    np.array(self.group.COM_velocity_residuals), 
                    np.array(self.group.ts), 
                    elapsed_time)
    
class Experiments:
    def __init__(self):
        self.exp1 = Experiment(seed=0)
        self.exp2 = Experiment(seed=0)
    
    def run(self, max_N, max_iter, max_t=10, max_t_iter=100):
        res = dict()
        res['Ns'] = np.linspace(1, max_N, max_iter).astype(int)
        res['standart_times'] = []
        res['barneshut_times'] = []
        for N in res['Ns']:
            print(N)
            a1, a2, a3, elap_time = self.exp1.run(method='standart', N=N, run_time=max_t, num_iterations=max_t_iter)
            res['standart_times'].append(elap_time/max_t_iter)
            a1, a2, a3, elap_time = self.exp2.run(method='barnes_hut', N=N, run_time=max_t, num_iterations=max_t_iter)
            res['barneshut_times'].append(elap_time/max_t_iter)
            
            #save res
            with open('res.pickle', 'wb') as handle:
                pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)
                
    def run_parallel(self, max_N, max_iter, max_t=10, max_t_iter=100, workers=8):
        res = dict()
        res['Ns'] = []
        #res['standart_times'] = []
        res['barneshut_times'] = []
        for N in np.linspace(1, max_N, max_iter).astype(int):
            print(N)
            
            #a1, a2, a3, elap_time = self.exp1.run(method='standart', 
            #                                      N=N, run_time=max_t, 
            #                                      num_iterations=max_t_iter)
            #
            #res['standart_times'].append(elap_time/max_t_iter)
            
            elap_time = self.exp2.run(method='barnes_hut', 
                                                  N=N, run_time=max_t, 
                                                  num_iterations=max_t_iter, 
                                                  workers=workers, time_test=True)
            
            res['barneshut_times'].append(elap_time/max_t_iter)
            res['Ns'].append(N)
            
            #save res
            with open('res_parallel.pickle', 'wb') as handle:
                pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)