import numpy as np
from barnes_hut import *
from initial_conditions import make_disc
from leapfrog import leapfrog_step
from animation import *
from energy import *
from analysis import *

def run_simulation(p, v, mass, dt, steps, theta, save_every=1):
    
    N = np.shape(p)[0]
    spiral_strength = []
    saved_pos = []
    saved_vel = []
    
    saved_pos.append(p.copy())
    saved_vel.append(v.copy())
    spiral_strength.append(measure_spiral_strength(p, N))
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
            spiral_strength.append(measure_spiral_strength(p, N))
            
    return np.array(saved_pos), np.array(saved_vel), np.array(spiral_strength)

def main():
    N = 500
    Rd = 1.0
    M = 10.0
    sigma_v1 = 0.1
    sigma_v2 = 0.5
    sigma_v3 = 2.0
    dt = 0.01
    steps = 500
    theta = 0.5

    pos1, vel1, mass1 = make_disc(N, Rd, M, sigma_v1)
    pos2, vel2, mass2 = make_disc(N, Rd, M, sigma_v2)
    pos3, vel3, mass3 = make_disc(N, Rd, M, sigma_v3)
    
    r_bins = np.linspace(0, 5, 25)
    r_mid, Q = measure_Q(pos2, vel2, mass2, r_bins)
    plt.figure()
    plt.plot(r_mid, Q)
    plt.axhline(y=1, color='r', linestyle='--', label='Q=1 threshold')
    plt.xlabel('r')
    plt.ylabel('Toomre Q')
    plt.title('Toomre Q Profile')
    plt.ylim(0, 10)
    plt.legend()
    plt.show()

    p1,v1,s1 = run_simulation(pos1, vel1, mass1, dt, steps, theta)
    p2,v2,s2 = run_simulation(pos2, vel2, mass2, dt, steps, theta)
    p3,v3,s3 = run_simulation(pos3, vel3, mass3, dt, steps, theta)

    t = np.arange(steps)

    plt.plot(t,s1, label="Sigma=0.1")
    plt.plot(t,s2, label="Sigma=0.5")
    plt.plot(t,s3, label="Sigma=2.0")
    plt.xlabel("Time")
    plt.ylabel("Spiral strength")
    plt.legend()
    plt.show()
    
if __name__ == '__main__':
    main()