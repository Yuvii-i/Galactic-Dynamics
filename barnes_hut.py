import numpy as np
from collections import namedtuple

G = 1
epsilon = 0.1
zero_array = np.array([0.0,0.0])

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

    def compute_acceleration(self, particle, theta=0.5):
        
        M = self.mass
        rx, ry = self.cx - particle.x, self.cy - particle.y
        rnorm = rx**2 + ry**2 + epsilon**2
        d = self.x_max - self.x_min
        
        if self.is_empty():
            return zero_array
        
        elif self.is_leaf():
            if self.particle is particle:
                return zero_array
            
            else:
                return G*M*rx/(rnorm)**(3/2), G*M*ry/(rnorm)**(3/2)
        
        elif d/np.sqrt(rx**2 + ry**2) < theta:
            return G*M*rx/(rnorm)**(3/2), G*M*ry/(rnorm)**(3/2)
        
        else:
            Fx , Fy = 0.0, 0.0
            
            for child in self.children:
                fx, fy =child.compute_acceleration(particle, theta)
                Fx += fx
                Fy += fy
                
            return Fx, Fy
        
        
        