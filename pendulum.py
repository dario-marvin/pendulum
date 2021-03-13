# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 16:59:34 2021

@author: schlammi
"""


class simplependulum:
    def __init__(self,x0=0,v0=0,l=20,t1=0,deltat=0,tau=0,vt=1):
        """
        x0,v0 : Pendulum styarting positions
        l = length of the pendulum
        t1 start of the first move
        deltat duration of the move
        t1+tau start of the second move
        trolley speed.
        everythin in si units m, s
        
        """
        self. g    = 9.807
        self.l = l
        
        self.omega = np.sqrt(self.g/self.l)
        self.t1 = t1
        self.deltat =  deltat
        self.tau =  tau
        self.vt = vt
        self.M=np.matrix(((0,1),(-self.omega**2,0)))
        self.dt = 0.005  # time step in the simulation
        self.y = np.matrix((x0,v0)).T
        self.t =0
        
        # output arrays 
        self.out_t  = []
        self.out_x  = []
        self.out_v  = []
        self.out_xt  =[]
        self.out_vt  =[]
        
        self.out_t.append(self.t)
        self.out_x.append(self.y[0,0])
        self.out_v.append(self.y[1,0])
        self.out_xt.append(self.trolley(self.t))
        self.out_vt.append(0)
        

    def go(self,tstop):
        while self.t<tstop:
            self.iterate()
    
    
    def trolley(self,t):
        if t<self.t1:
            return 0
        if t>=self.t1 and t<self.t1+self.deltat:
            return self.vt*(t-self.t1)
        if t>=self.t1+self.deltat  and t<self.t1+self.tau:
            return self.vt*self.deltat
        if t>=self.t1+self.tau and t<self.t1+self.tau+self.deltat:
            return self.vt*self.deltat+self.vt*(t-(self.t1+self.tau))
        if t>=self.t1+self.tau+self.deltat:
            return 2*self.vt*self.deltat
    
    def N(self,t):
        n = self.trolley(t)
        return np.matrix((0,n )).T

    def f(self,t,y):
        return self.omega**2*self.N(t)+self.M*y
    
    def iterate(self):
        F1 = self.dt * self.f(self.t         ,self.y)
        F2 = self.dt * self.f(self.t+self.dt/2,self.y +F1/2)
        F3 = self.dt * self.f(self.t+self.dt/2,self.y +F2/2)
        F4 = self.dt * self.f(self.t+self.dt,self.y +F3)
        self.y = self.y + 1.0/6.0*(F1+2*F2+2*F3+F4)
        self.t = self.t +self.dt

        self.out_t.append(self.t)
        self.out_x.append(self.y[0,0])
        self.out_v.append(self.y[1,0])
        self.out_xt.append(self.trolley(self.t))
        self.out_vt.append((self.trolley(self.t)-self.trolley(self.t-self.dt))/self.dt)
        