import numpy as np
from barnes_hut import *
from initial_conditions import make_disc
from leapfrog import leapfrog_step

def run_simulation(p, v, M, dt, steps, theta, save_every=10):
    
    N = np.shape(p)[0]
    
    saved_pos = []
    saved_vel = []
    
    saved_pos.append(p.copy())
    saved_vel.append(v.copy())
    
    x_min,x_max,y_min,y_max = np.min(p[:,0]), np.max(p[:,0]), np.min(p[:,1]), np.max(p[:,1])
    
    a = []
        
    Q = QuadNode(x_min,x_max,y_min,y_max)
        
    for j in range(N):
            
        particle = Particle(p[j,0],p[j,1], M/N)
            
        Q.insert(particle)
            
    for j in range(N):
        
        particle = Particle(p[j,0],p[j,1], M/N)
            
        ax, ay = Q.compute_acceleration(particle, theta)
        a.append((ax,ay))
    
    for i in range(1,steps):
        
        p, vhalf = leapfrog_step(p, v, np.array(a), dt)
        
        x_min,x_max,y_min,y_max = np.min(p[:,0]), np.max(p[:,0]), np.min(p[:,1]), np.max(p[:,1])
        
        a = []
        
        Q = QuadNode(x_min,x_max,y_min,y_max)
        
        for j in range(N):
            
            particle = Particle(p[j,0],p[j,1], M/N)
            
            Q.insert(particle)
            
        for j in range(N):
            
            particle = Particle(p[j,0],p[j,1], M/N)
            
            ax, ay = Q.compute_acceleration(particle, theta)
            a.append((ax,ay))
            
        v = vhalf + 0.5*dt*np.array(a)
        
        if i % save_every == 0:
            saved_pos.append(p)
            saved_vel.append(v)
            
    return np.array(saved_pos), np.array(saved_vel)
        
        
positions, velocities = make_disc(100, 1.0, 1.0)
saved_pos, saved_vel = run_simulation(positions, velocities, M=1.0, dt=0.01, steps=100, theta=0.5)

import cProfile
cProfile.run('run_simulation(positions, velocities, M=1.0, dt=0.01, steps=10, theta=0.5)')        