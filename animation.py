import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def show_animation(p, N, steps, save_every=1, trail=10):
    
    fig, ax = plt.subplots()
    
    scatter = ax.scatter(p[0,:,0], p[0,:,1], s=6)
    
    
    def update(frame):
        
        scatter.set_offsets(p[frame])      
    
    ani = FuncAnimation(fig, update, frames=int(steps/save_every), interval=20)
    
    plt.show()