import numpy as np

G = 1

def calculate_energy(p, v, m, epsilon=0.5):
    
    ke = 0.5 * np.sum(m*(np.linalg.norm(v, axis=2))**2, axis=1)
    
    r = p[:,np.newaxis,:,:] - p[:,:,np.newaxis,:] 
    rnorm = np.linalg.norm(r, axis=3)  + epsilon
    
    mass = np.triu(m[:,np.newaxis]*m[np.newaxis,:])
    
    pe = np.triu(-G*mass[np.newaxis,:,:]/rnorm).sum(axis=(1,2))
    
    return ke+pe, ke, pe