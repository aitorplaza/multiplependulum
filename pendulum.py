import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate

N = 4
g = 9.8

lengths = np.arange(1,1.3-(0.3/N),0.3/(N))
# ~ lengths = np.arange(1,1.4,0.04)
print(lengths)

theta_0 = 0.9*np.ones((1,N))
dtheta_0 = np.zeros((1,N))

time = 0
dt = 0.01
nSteps = 10000;

t_series = np.zeros((nSteps,1))
theta_series = np.zeros((nSteps,N))
dtheta_series = np.zeros((nSteps,N))
Etotal_series = np.zeros((nSteps,1))

theta = theta_0
dtheta = dtheta_0

m=1
# ~ for i in range(0,nSteps):

    # ~ # calculate accelerations
    # ~ ddtheta = -np.multiply((g/lengths),np.sin(theta))

    # ~ #integrate    
    # ~ dtheta = dtheta + ddtheta*dt
    # ~ theta = theta + dtheta*dt + 0.5*ddtheta*dt**2
    
    # ~ # save the values
    # ~ t_series[i] = dt*i
    # ~ theta_series[i,:] = theta
    # ~ dtheta_series[i,:] = theta
    
    # ~ Etotal = 0;
    # ~ for j,l in enumerate(lengths):
        # ~ Etotal=Etotal+0.5*m*l*dtheta[0,j]**2 + m*g*l*(1-np.cos(theta[0,j]))
        
    # ~ Etotal_series[i,:] = Etotal


def derivs(state, t):
    dydx = np.zeros((2*N,1))
    state=state[0]
    print(state)
    for j,l in enumerate(lengths):
        dydx[j] = state[N+j]
        dydx[N+j] = -(g/lengths[j])*np.sin(state[j])
    return dydx.reshape(1,2*N)

# ~ print(derivs([1,1,1,1,1,2,3,4],1))

# create a time array from 0..100 sampled at 0.05 second steps
dt = 0.05
t = np.arange(0.0, 0.1, dt)

# initial state
state = np.concatenate((theta_0,dtheta_0)).reshape(1,2*N)
print("state = ", state[0])

print(derivs(state,1))

y = integrate.odeint(derivs, state, t)


print(y)
'''
    
plt.figure(0)
plt.plot(t_series,Etotal_series)
plt.xlabel('Time (s)')
plt.ylabel('Total Energy of the system(J)')
plt.savefig('plot1.png', bbox_inches='tight')

plt.figure(2)
plt.plot(t_series,theta_series[:,0])
plt.xlabel('Time (s)')
plt.ylabel('theta 1(J)')
plt.savefig('theta1.png', bbox_inches='tight')

        
# Create gif:
x_series = np.zeros((nSteps,N))
y_series = np.zeros((nSteps,N))
for i,l in enumerate(lengths):
    x_series[:,i] = l*np.sin(theta_series[:,i])
    y_series[:,i] = -l*np.cos(theta_series[:,i])



fig = plt.figure(1)
ax = plt.axes(xlim=(-1.3, 1.3), ylim=(-1.4, 0.2))
plt.plot( 0.0, 0.0,'o',markersize=3,color='k') 
plt.axis('off')
lines = [plt.plot([], [], 'k-', linewidth = 1 )[0] for _ in range(N)] #lines to animate
balls = [plt.plot([], [], 'o' ,markersize=12)[0] for _ in range(N)] #balls to animate

patches = lines + list(balls)#things to animate

def init():
    #init lines
    for line in lines:
        line.set_data([], [])

    #init balls
    for ball in balls:
        ball.set_data([], [])
        
    return patches 
    
def animate(i):
    #animate lines
    for j,line in enumerate(lines):
        x = [0.0, x_series[i,j]]
        y = [0.0, y_series[i,j]]
        line.set_data(x, y)
        
    #animate balls
    for j,ball in enumerate(balls):
        x = [x_series[i,j]]
        y = [y_series[i,j]]
        ball.set_data(x, y)
    
    return patches

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=5000, interval=10, blit=True)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=100, metadata=dict(artist='Me'), bitrate=1800)
anim.save('im.mp4', writer=writer)
#~ plt.show()


'''




