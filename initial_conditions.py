import numpy as np
import matplotlib.pyplot as plt
import random

G = 1

def make_disc(N,Rd,M, sigma):
    
    p = lambda r: (1/Rd**2) * r * np.exp(-r/Rd)
    
    r_max = 5 * Rd
    p_max = 1/(np.exp(1) * Rd)
    
    r = []
    
    while len(r) < N:
        
        r_random = random.uniform(0,r_max)
        p_random = random.uniform(0,p_max)
        if p_random < p(r_random):
            
            r.append(r_random)
            
    theta = np.random.uniform(0,2*np.pi,N)
    
    r = np.array(r)
    
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    
    v = 1.2*np.sqrt(G*M/r)
    vx = - v * np.sin(theta)
    vy = v * np.cos(theta)
    
    
    vx += np.random.normal(0, sigma, N)
    vy += np.random.normal(0, sigma, N)
    
    center_pos = np.array([0.0,0.0])
    center_vel = np.array([0.0,0.0])
    
    pos = np.vstack([center_pos, np.transpose(np.array([x,y]))])
    vel = np.vstack([center_vel, np.transpose(np.array([vx,vy]))])
    
    mass = np.concatenate([[M], np.full(N, M/N)])
    
    return pos , vel, mass
