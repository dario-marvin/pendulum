"""
A stationary damped pendulum
For fun.

"""
import numpy as np
import matplotlib.pyplot as plt
import pendulum

T0 = 8.97
pend = pendulum.simplependulum(x0=1,v0=0,w0=2*np.pi/T0,xi=0.1)
pend.go(30)
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend.out_t,pend.out_vm,'b-')
ax[0].plot(pend.out_t,pend.out_vc,'r:')
ax[1].plot(pend.out_t,pend.out_xm,'b-')
ax[1].plot(pend.out_t,pend.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 1: A damped crane pendulum")


#Trick 1 with an undamped pendulum
# start moving and stop a period later

tl1 = [0,10,10+T0]  #times when a change in velocity  happens
vl1 = [0,1,0]       #the velocity it changes to at that time
move1 = pendulum.const_vel(tl1,vl1)
pend1 = pendulum.simplependulum(x0=0,v0=0,w0=2*np.pi/T0,xi=0.0,trolley_pos=move1.xc)
pend1.go(30)
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend1.out_t,pend1.out_vm,'b-')
ax[0].plot(pend1.out_t,pend1.out_vc,'r:')
ax[1].plot(pend1.out_t,pend1.out_xm,'b-')
ax[1].plot(pend1.out_t,pend1.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()

print("Press Enter to quit")
input()