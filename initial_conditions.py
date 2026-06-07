import numpy as np
import matplotlib.pyplot as plt
import random

G = 1

def make_disc(N,Rd,M):
    
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
    
    v = np.sqrt(G*M/r)
    vx = - v * np.sin(theta)
    vy = v * np.cos(theta)
    
    return np.transpose(np.array([x,y])) , np.transpose(np.array([vx,vy]))

x, v = make_disc(10000, 1.0, 1.0)

plt.scatter(x[:, 0], x[: , 1])
plt.show()