import numpy as np
from barnes_hut import *
from initial_conditions import make_disc
from leapfrog import leapfrog_step
from animation import *
from energy import *

def run_simulation(p, v, mass, dt, steps, theta, save_every=1):
    
    N = np.shape(p)[0]
    spiral_strength = []
    saved_pos = []
    saved_vel = []
    
    saved_pos.append(p.copy())
    saved_vel.append(v.copy())
    spiral_strength.append(measure_spiral_strength(p))
    x_min,x_max,y_min,y_max = np.min(p[:,0]), np.max(p[:,0]), np.min(p[:,1]), np.max(p[:,1])
    
    a = []
        
    Q = QuadNode(x_min,x_max,y_min,y_max)
        
    for j in range(N):
            
        particle = Particle(p[j,0],p[j,1],mass[j])
            
        Q.insert(particle)
            
    for j in range(N):
        
        particle = Particle(p[j,0],p[j,1], mass[j])
            
        ax, ay = Q.compute_acceleration(particle, theta)
        a.append((ax,ay))
    
    for i in range(1,steps):
        
        p, vhalf = leapfrog_step(p, v, np.array(a), dt)
        
        x_min,x_max,y_min,y_max = np.min(p[:,0]), np.max(p[:,0]), np.min(p[:,1]), np.max(p[:,1])
        
        a = []
        
        Q = QuadNode(x_min,x_max,y_min,y_max)
        
        for j in range(N):
            
            particle = Particle(p[j,0],p[j,1], mass[j])
            
            Q.insert(particle)
            
        for j in range(N):
            
            particle = Particle(p[j,0],p[j,1], mass[j])
            
            ax, ay = Q.compute_acceleration(particle, theta)
            a.append((ax,ay))
            
        v = vhalf + 0.5*dt*np.array(a)
        
        if i % save_every == 0:
            saved_pos.append(p)
            saved_vel.append(v)
            spiral_strength.append(measure_spiral_strength(p))
            
    return np.array(saved_pos), np.array(saved_vel), np.array(spiral_strength)

def measure_surface_density(positions, masses, rbins):
    
    mass_per_particle = masses[1] 
    r = np.linalg.norm(positions, axis=1)
    
    particle_count, bin_edges = np.histogram(r, rbins)
     
    dr = (bin_edges[1:] - bin_edges[:-1])
    r_mid = (bin_edges[1:] + bin_edges[:-1])/2
    
    mass_in_bin = mass_per_particle*particle_count
    
    area = 2*np.pi*r_mid*dr
    
    sigma = mass_in_bin/area
    
    return r_mid, sigma

def measure_spiral_strength(position):
    
    phi = np.arctan2(position[:,1],position[:,0])
    
    A2 = np.sum(np.exp(2j * phi)) / N
    
    return np.abs(A2)
    
def measure_Q(position, velocity, mass, r_bins):
    
    r_mid, sigma = measure_surface_density(position, mass, r_bins)
    
    r = np.linalg.norm(position, axis=1)
    v = np.linalg.norm(velocity, axis=1)
    
    particle_count, bin_edges = np.histogram(r, r_bins)
    
    cs = np.empty(len(bin_edges)-1) 
    
    for i in range(len(bin_edges)-1):
        r_left = bin_edges[i]
        r_right = bin_edges[i+1]
        mask = (r >= r_left) & (r < r_right)
        
        if np.sum(mask) < 2:
            cs[i] = np.nan
        else:
            cs[i] = np.std(v[mask])
        
    k = np.sqrt(G*mass[0]/(r_mid**(3/2)))
    
    return r_mid, k*cs/ (np.pi*G*sigma)
    
N = 500
Rd = 1.0
M = 1.0

pos1, vel1, mass1 = make_disc(500, 1.0, 10.0, 0.1)
pos2, vel2, mass2 = make_disc(500, 1.0, 10.0, 0.5)
pos3, vel3, mass3 = make_disc(500, 1.0, 10.0, 2.0)


p1,v1,s1 = run_simulation(pos1, vel1, mass1, dt=0.01, steps=500, theta=0.5)
p2,v2,s2 = run_simulation(pos2, vel2, mass2, dt=0.01, steps=500, theta=0.5)
p3,v3,s3 = run_simulation(pos3, vel3, mass3, dt=0.01, steps=500, theta=0.5)

t = np.arange(500)

plt.plot(t,s1, label="Sigma=0.1")
plt.plot(t,s2, label="Sigma=0.5")
plt.plot(t,s3, label="Sigma=2.0")
plt.legend()
plt.show()