"""
A stationary damped pendulum
For fun.

"""
import numpy as np
import matplotlib.pyplot as plt
import pendulum

### Pendulum swings with amplitude smaller than $4 v_c$

#We assume a pendulum length of 20 m and a trolley speed of vc = 1 m/s. The load needs to be transported by 6 m, hence deltat=  3 s. 

#The initial speed of the pendulum is v0 = 1 m/s

# `tau` is calculated according to eq. 25 in the paper. After that value is otained, t1 according to eq. 22

l  = 20
g  = 9.807 
w0 =np.sqrt(g/l)
T0 = 2*np.pi/w0
vt = 1
deltat = 3
v0 = 1
tau=T0/np.pi*np.arccos(v0/(-4*vt*np.sin(w0*deltat*0.5))) 
t1 = T0/4-tau/2-deltat/2 
while t1<0:
    t1=t1+T0/2
t1= t1+T0   # You can delete this line, I just like to observe the pendulum for a period
print('tau = {0:6.1f}'.format(tau))
print('t1= {0:6.1f}'.format(t1))


tl1 = [0,t1,t1+deltat,t1+tau,t1+tau+deltat]  #times when a change in velocity  happens
vl1 = [0,vt,0        ,vt    ,0]       #the velocity it changes to at that time
move1 = pendulum.const_vel(tl1,vl1)
pend1 = pendulum.simplependulum(x0=0,v0=v0,w0=w0,xi=0.0,trolley_pos=move1.xc)
pend1.go(24)
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend1.out_t,pend1.out_vm,'b-')
ax[0].plot(pend1.out_t,pend1.out_vc,'r:')
ax[1].plot(pend1.out_t,pend1.out_xm,'b-')
ax[1].plot(pend1.out_t,pend1.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 1: Damping a pendulum with v0<4vc")


## Pendulum swings with amplitude larger than 4 vc
#We assume a pendulum length of 20 m and a trolley speed of vc = 1 m/s. The load needs to be transported by 6 m, hence deltat=  3 s. 

#The initial speed of the pendulum is v0 = 5 m/s

#tau is 2*T0. With this value, t1  is calculated according to eq. 22
l  = 20
g  = 9.807 
w0 =np.sqrt(g/l)
T0 = 2*np.pi/w0
vt = 1
deltat = 3
v0 = 5
tau=T0*2

t1 = T0/4-tau/2-deltat/2 
while t1<0:
    t1=t1+T0/2
t1= t1+T0   # You can delete this line, I just like to observe the pendulum for a period
print('tau = {0:6.1f}'.format(tau))
print('t1= {0:6.1f}'.format(t1))
tl2 = [0,t1,t1+deltat,t1+tau,t1+tau+deltat]  #times when a change in velocity  happens
vl2 = [0,vt,0        ,vt    ,0]       #the velocity it changes to at that time
move2 = pendulum.const_vel(tl2,vl2)
pend2 = pendulum.simplependulum(x0=0,v0=v0,w0=w0,xi=0.0,trolley_pos=move2.xc)
pend2.go(40)
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend2.out_t,pend2.out_vm,'b-')
ax[0].plot(pend2.out_t,pend2.out_vc,'r:')
ax[1].plot(pend2.out_t,pend2.out_xm,'b-')
ax[1].plot(pend2.out_t,pend2.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 1: Damping a pendulum with v0>4vc")

print("Press Enter to quit")
input()