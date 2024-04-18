"""
Created on Thu Apr 18 12:17:17 2024

@author: amissharma
"""

import numpy as np 
import matplotlib.pyplot as plt
import math
from matplotlib.animation import FuncAnimation
#parameters for analytical solution plot of the orbit 
x0 = 10 
v0 = 1/10  
h = x0*v0 
p = h**2 
E0 = (1/2)*(v0)**2 - 1/x0 
a = -1/(2*E0)
e = math.sqrt(1 - p/a) 
P = 2*math.pi*a**(3/2) 
t = np.arange(0, 100, 0.1)
n=100


#initial conditions

itr=1500
dt = P/500
#Cromer
rc= np.zeros((2,itr+1))
vc=np.zeros((2,itr+1))
rc[0,0]=10
rc[1,0]=0
vc[0,0]=0
vc[1,0]=0.1
#Newton
rn= np.zeros((2,itr+1))
vn=np.zeros((2,itr+1))
rn[0,0]=10
rn[1,0]=0
vn[0,0]=0
vn[1,0]=0.1
#2nd order r-kutta 
rk= np.zeros((2,itr+1))
vk=np.zeros((2,itr+1))
rk[0,0]=10
rk[1,0]=0
vk[0,0]=0
vk[1,0]=0.1
ax = - rk[0,0] / (rk[0,0]**2 + rk[1,0]**2)**(3/2)
ay =  - rk[1,0] / (rk[0,0]**2 + rk[1,0]**2)**(3/2) 
aax= (rk[0,0]+dt*vk[0,0])/ ((rk[0,0] + dt*vk[0,0])**2 + (rk[1,0]+dt*vk[1,0])**2)**(3/2) 
aay= (rk[1,0]+dt*vk[0,0])/ ((rk[0,0] + dt*vk[0,0])**2 + (rk[1,0]+dt*vk[1,0])**2)**(3/2) 
#position verlet 
rp= np.zeros((2,itr+1))
vp=np.zeros((2,itr+1))
rp[0,0]=10
rp[1,0]=0
vp[0,0]=0
vp[1,0]=0.1



for i in range(0,itr):
    #Cromer update
    vc[0,i+1]= vc[0,i]- rc[0,i]*dt/(rc[0,i]**2 + rc[1,i]**2)**(3/2)
    vc[1,i+1]=vc[1,i] - rc[1,i]*dt/(rc[0,i]**2 + rc[1,i]**2)**(3/2) 
    rc[0,i+1]= rc[0,i] + vc[0,i+1]*dt
    rc[1,i+1]=rc[1,i] + vc[1,i+1]*dt
    #Newton update
    rn[0,i+1]=rn[0,i]+vn[0,i]*dt 
    rn[1,i+1]=rn[1,i]+vn[1,i]*dt 
    vn[0,i+1]= vn[0,i]- rn[0,i+1]*dt/(rn[0,i+1]**2 + rn[1,i+1]**2)**(3/2)
    vn[1,i+1]= vn[1,i]- rn[1,i+1]*dt/(rn[0,i+1]**2 + rn[1,i+1]**2)**(3/2)
    #2nd Order RK update
    rk[0,i+1]= rk[0,i]+ vk[0,i]*dt +0.5*ax*dt**2
    rk[1,i+1]= rk[1,i]+ vk[1,i]*dt +0.5*ay*dt**2
    vk[0,i+1]= vk[0,i]+0.5*dt*(ax+aax)
    vk[1,i+1]= vk[1,i]+0.5*dt*(ay+aay)
    ax= - rk[0,i+1] / (rk[0,i+1]**2 + rk[1,i+1]**2)**(3/2)
    ay= - rk[1,i+1] / (rk[0,i+1]**2 + rk[1,i+1]**2)**(3/2) 
    aax= -(rk[0,i+1]+dt*vk[0,i+1])/ ((rk[0,i+1] + dt*vk[0,i+1])**2 + (rk[1,i+1]+dt*vk[1,i+1])**2)**(3/2) 
    aay= -(rk[1,i+1]+dt*vk[0,i+1])/ ((rk[0,i+1] + dt*vk[0,i+1])**2 + (rk[1,i+1]+dt*vk[1,i+1])**2)**(3/2) 
    #Position Verlet update
    rp[0,i+1]=rp[0,i]+0.5*dt*vp[0,i] 
    rp[1,i+1]=rp[1,i]+0.5*dt*vp[1,i] 
    apx= rp[0,i+1]/(rp[0,i+1]**2 + rp[1,i+1]**2)**(3/2)
    apy= rp[1,i+1]/(rp[0,i+1]**2 + rp[1,i+1]**2)**(3/2)
    vp[0,i+1] = vp[0,i] - dt*apx
    vp[1,i+1] = vp[1,i] - dt*apy
    rp[0,i+1]=rp[0,i+1]+0.5*dt*vp[0,i+1] 
    rp[1,i+1]=rp[1,i+1]+0.5*dt*vp[1,i+1] 
    

fig, ax = plt.subplots()
ax.set_xlim(-2, 11)
ax.set_ylim(-6,6)
line1, = ax.plot([], [], 'ro-', lw=1,markevery=5, label='Cromer')
line2, = ax.plot([], [], 'gd-', lw=1, markevery=5,label='Newton')
line3, = ax.plot([], [], 'ko-', lw=1, markevery=5,ms=4,label='Runge-Kutta O(2)')
line4, = ax.plot([], [], 'b*-', lw=0.5,markevery=5, label='Position Verlet')
plt.plot(p*np.cos(t*2*math.pi/n)/(1-e*np.cos(t*2*math.pi/n)),p*np.sin(t*2*math.pi/n)/(1-e*np.cos(t*2*math.pi/n)), label='Analytical plot')

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    return line1,line2, line3,line4

def update(frame):
    xc=rc[0, :frame+1]
    yc=rc[1, :frame+1]
    xn=rn[0, :frame+1]
    yn=rn[1, :frame+1]
    xk=rk[0, :frame+1]
    yk=rk[1, :frame+1]
    xp=rp[0, :frame+1]
    yp=rp[1, :frame+1]
    line1.set_data(xc, yc)
    line2.set_data(xn, yn)
    line3.set_data(xk, yk)
    line4.set_data(xp, yp)
    
    return line1,line2,line3,line4

ani = FuncAnimation(fig, update, frames=np.arange(itr), init_func=init, blit=True, repeat=False,interval=5)
plt.title('Kepler solution with various methods' )
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()
plt.show() 


#modify plotmarkers to get plots to your satisfaction. One can also change I.C to 
#get different orbits  

