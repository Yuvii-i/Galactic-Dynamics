import numpy as np

def leapfrog_step(p, v, a, dt):
    
    hv = v + 0.5*dt*a        
    p = p + dt*hv
             
    return p, hv