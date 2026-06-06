import numpy as np
from collections import namedtuple

G = 1
epsilon = 0.1

Particle = namedtuple('Particle', ['x', 'y', 'mass'])

class QuadNode:
    
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.mass = 0
        self.cx = 0
        self.cy = 0
        self.children = None
        self.particle = None
        
    def is_empty(self):
        return self.particle is None and self.children is None
    
    def is_leaf(self):
        return self.particle is not None and self.children is None
    
    def subdivide(self):
        
        x_mid = (self.x_max + self.x_min)/2
        y_mid = (self.y_max + self.y_min)/2
        
        self.children = [ QuadNode(x_mid, self.x_max, y_mid, self.y_max),   #Top Right
                         QuadNode(self.x_min, x_mid, y_mid, self.y_max),    #Top Left
                         QuadNode(x_mid, self.x_max, self.y_min, y_mid),    #Bottom Right
                         QuadNode(self.x_min, x_mid, self.y_min, y_mid)     #Bottom Left 
                         ]
        
    def get_child(self, particle):
        
        x_mid = (self.x_max + self.x_min)/2
        y_mid = (self.y_max + self.y_min)/2
        
        if particle.x >= x_mid and particle.y >= y_mid:
            return self.children[0]
        
        elif particle.x < x_mid and particle.y > y_mid:
            return self.children[1]
        
        elif particle.x > x_mid and particle.y < y_mid:
            return self.children[2]
        
        else:
            return self.children[3]
        
    def insert(self, particle):
        
        if self.is_empty():
            self.particle = particle
            
        
        elif self.is_leaf():
            if self.particle is particle:
                return
            
            else:
                old = self.particle
                self.particle = None
                
                self.subdivide()
                
                self.get_child(old).insert(old)
                self.get_child(particle).insert(particle)
                
        else:
            self.get_child(particle).insert(particle)
            
        self.cx = (self.mass*self.cx + particle.mass*particle.x)/(self.mass + particle.mass)
        self.cy = (self.mass*self.cy + particle.mass*particle.y)/(self.mass + particle.mass)   
        self.mass += particle.mass

    def compute_force(self, particle, theta):
        
        M = self.mass
        m = particle.mass
        r = np.array([self.cx - particle.x, self.cy - particle.y])
        d = self.x_max - self.x_min
        
        if self.is_empty():
            return np.array([0.0,0.0])
        
        elif self.is_leaf():
            if self.particle is particle:
                return np.array([0.0,0.0])
            
            else:
                return G*M*m*r/(np.linalg.norm(r)**2 + epsilon**2)**(3/2)
        
        elif d/np.linalg.norm(r) < theta:
            return G*M*m*r/(np.linalg.norm(r)**2 + epsilon**2)**(3/2)
        
        else:
            F = np.zeros((4,2))
            
            for i in range(4):
                F[i] =self.children[i].compute_force(particle, theta)
                
            return np.sum(F, axis=0)
        
Q = QuadNode(-2,2,-2,2)

p1 = Particle(-1.0, 0.0, 1.0)
p2 = Particle(1.0, 0.0, 1.0)

[Q.insert(i) for i in [p1,p2]]

print(Q.compute_force(p1,0.5))

 