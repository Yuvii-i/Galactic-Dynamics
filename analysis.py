import numpy as np

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

def measure_spiral_strength(position, N):
    
    phi = np.arctan2(position[:,1],position[:,0])
    
    A2 = np.sum(np.exp(2j * phi)) / N
    
    return np.abs(A2)
    
def measure_Q(position, velocity, mass, r_bins, G=1):
    
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