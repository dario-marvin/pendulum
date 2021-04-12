"""
A stationary damped pendulum
For fun.

"""
import numpy as np
import matplotlib.pyplot as plt
import pendulum

T0 = 8.97
pend = pendulum.simplependulum(x0=1,v0=0,w0=2*np.pi/T0,xi=0.1)
pend.go(30) # This simulates for 30 seconds, adjust the number to get more
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend.out_t,pend.out_vm,'b-')
ax[0].plot(pend.out_t,pend.out_vc,'r:')
ax[1].plot(pend.out_t,pend.out_xm,'b-')
ax[1].plot(pend.out_t,pend.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 1: A stationary damped pendulum")


#Trick 1 with an undamped pendulum
# start moving and stop a period later

tl1 = [0,10,10+T0]  #times when a change in velocity  happens
vl1 = [0,1,0]       #the velocity it changes to at that time
move1 = pendulum.const_vel(tl1,vl1)
pend1 = pendulum.simplependulum(x0=0,v0=0,w0=2*np.pi/T0,xi=0.0,trolley_pos=move1.xc)
pend1.go(30) # This simulates for 30 seconds, adjust the number to get more
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend1.out_t,pend1.out_vm,'b-')
ax[0].plot(pend1.out_t,pend1.out_vc,'r:')
ax[1].plot(pend1.out_t,pend1.out_xm,'b-')
ax[1].plot(pend1.out_t,pend1.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 2: Trick 1 start moving and stop a full period later ")


### Trick 2 with an undamped pendulum. 
#Go to v0/2 and then half a period later to v0
n=1 # you can choose 
tl2 = [0,10,10+T0/2+n*T0]
vl2 = [0,0.5,1]
move2 = pendulum.const_vel(tl2,vl2)
pend2 = pendulum.simplependulum(x0=0,v0=0,w0=2*np.pi/T0,xi=0.0,trolley_pos=move2.xc)
pend2.go(15+T0/2+n*T0) # This simulates for 15+T0/2+n*T0 seconds, adjust the argument to get more or less
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend2.out_t,pend2.out_vm,'b-')
ax[0].plot(pend2.out_t,pend2.out_vc,'r:')
ax[1].plot(pend2.out_t,pend2.out_xm,'b-')
ax[1].plot(pend2.out_t,pend2.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 3: Trick 2 start moving half speed and  go full speed an odd mutiple of half a period later")

### Trick 2 with stopping. 
#Go to v0/2 and then half a period later to v0
n=0 # Here we choose n=0
# for other n the 20 belwo must be adjusted
t1 = 5 
t2 = t1+T0/2+n*T0
t3 = t2+ 23 #glide for 23 seconds
t4 = t3+T0/2+n*T0
tl3 = [0,t1,t2,t3,t4]
vl3 = [0,0.5,1,0.5,0]
move3 = pendulum.const_vel(tl3,vl3)
pend3 = pendulum.simplependulum(x0=0,v0=0,w0=2*np.pi/T0,xi=0.0,trolley_pos=move3.xc)
pend3.go(t4 +10)   # This simulates for t4+10 seconds, adjust the argument to get more or less
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend3.out_t,pend3.out_vm,'b-')
ax[0].plot(pend3.out_t,pend3.out_vc,'r:')
ax[1].plot(pend3.out_t,pend3.out_xm,'b-')
ax[1].plot(pend3.out_t,pend3.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 4: Trick 2 used twice")

### Trick 1 with acceleration. 
delt = 4
t1 = 10
t2 = t1+delt
t3 = t1+T0
t4 = t3+delt
tl4 = [0,t1,t2,t3,t4]  #times when a change in velocity  happens
al4 = [0,1 ,0 ,-1,0]       #the velocity it changes to at that time
move4= pendulum.const_acc(tl4,al4)
pend4 = pendulum.simplependulum(x0=0,v0=0,w0=2*np.pi/T0,xi=0.0,trolley_pos=move4.xc)
pend4.go(t4+10)  # This simulates for t4+10 s, adjust the argument to get more or less
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend4.out_t,pend4.out_vm,'b-')
ax[0].plot(pend4.out_t,pend4.out_vc,'r:')
ax[1].plot(pend4.out_t,pend4.out_xm,'b-')
ax[1].plot(pend4.out_t,pend4.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 5: Trick 1 with finite acceleration")

### Trick 2 with acceleration. 
delt = 4
t2 = t1+delt
t3 = t1+T0/2
t4 = t3+delt
tl5 = [0,t1,t2,t3,t4]  #times when a change in velocity  happens
al5 = [0,1 ,0 ,1,0]       #the velocity it changes to at that time
move5= pendulum.const_acc(tl5,al5)
pend5 = pendulum.simplependulum(x0=0,v0=0,w0=2*np.pi/T0,xi=0.0,trolley_pos=move5.xc)
pend5.go(t4+10) # This simulates for t4+10 s, adjust the argument to get more or less
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(pend5.out_t,pend5.out_vm,'b-')
ax[0].plot(pend5.out_t,pend5.out_vc,'r:')
ax[1].plot(pend5.out_t,pend5.out_xm,'b-')
ax[1].plot(pend5.out_t,pend5.out_xc,'r:')
ax[1].set_xlabel('t/s')
ax[0].set_ylabel('v/(m/s)')
ax[1].set_ylabel('x/m')
fig.show()
print("Figure 6: Trick 2 with finite acceleration")



print("Press Enter to quit")
input()