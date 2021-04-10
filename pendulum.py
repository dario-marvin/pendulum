"""
pendulum.py 

implements simplependulum
"""
import numpy as np


class const_vel:
    def __init__(self,ts,vs):
        """
        ts and vs are both lists, the class sorts the lists.
        while maintainint the pairs for the sorted pairs tp,vp
        if t>tp[i] then vp[i] will be returned
        if t<tp[0] 0 will be returned
        """
        self.tp,self.vp = zip(* sorted(zip(ts,vs),key = lambda x: x[0]))
        self.mint = min(self.tp)
    
    def xc(self,t):
        if t<self.mint:
            return 0
        x = 0  # we always start at 0
        ot =0
        ov=0
        for t_,v_ in zip(self.tp,self.vp):
            if t>=t_ :
                    x=x+(t_-ot)*ov
            else:
                return x+(t-ot)*ov
            ot=t_  
            ov=v_
        return x+(t-t_)*v_

class const_acc:
    def __init__(self,ts,a):
        """
        ts and a are both lists, the class sorts the lists.
        while maintainint the pairs for the sorted pairs tp,ap
        if t>ap[i] then ap[i] will be returned
        if t<tp[0] 0 will be returned
        """
        self.tp,self.ap = zip(* sorted(zip(ts,a),key = lambda x: x[0]))
        self.mint = min(self.tp)
    
    def xc(self,t):
        if t<self.mint:
            return 0
        x = 0  # we always start at 0
        ot =0
        oa =0
        v =0
        for t_,a_ in zip(self.tp,self.ap):
            if t>=t_ :
                    x=x+0.5*oa*(t_-ot)**2+(t_-ot)*v
                    v=v+(t_-ot)*oa
            else:
                return x+0.5*oa*(t-ot)**2+(t-ot)*v 
            ot=t_  
            oa=a_
        return x+0.5*a_*(t-ot)**2+(t-ot)*v 


    
    
class simplependulum:
    def __init__(self,x0=0,v0=0,w0=0.700249,xi=0,trolley_pos=None,dt =0.005):
        """
        x0,v0       : pendulum initial conditions x0 in m, v0 in m/s
        w0          : angular frequncy of the pendulum sqrt(g/l)
        xi          : damping ratio xi = c/(2m w0), xi =1 critically damped
        trolley_pos : a function f(t) that returns the positon of the trolley
                      as a function of t
        """
        self.omega = w0
        self.xi    = xi
        self.M=np.matrix(((0,1),(-self.omega**2,-2*self.omega*self.xi)))
        self.dt =dt  # time step in the simulation
        self.xm = np.matrix((x0,v0)).T
        self.t =0
        self.xc = trolley_pos
        
        # output arrays 
        self.out_t    = []    # time
        self.out_xm   = []    # x positon of mass
        self.out_vm   = []    # vel of mass
        self.out_xc   = []    # x position of trolley
        self.out_vc   = []    # velocity of trolley
        
        xc0 = 0
        if self.xc!=None:
            xc0 = xc0 = self.xc(self.t)
        self.out_t.append(self.t)
        self.out_xm.append(self.xm[0,0])
        self.out_vm.append(self.xm[1,0])
        self.out_xc.append(xc0)
        self.out_vc.append(0)
        

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
    
    def get_xc(self,t):
        if self.xc!=None:
            return self.xc(t)
        else:
            return 0
   

    def f(self,t,y):
        F=np.matrix((0,self.get_xc(t) )).T # Matrix for the ext. force.
        return self.omega**2*F+self.M*y
    
    def iterate(self):
        F1 = self.dt * self.f(self.t         ,self.xm)
        F2 = self.dt * self.f(self.t+self.dt/2,self.xm +F1/2)
        F3 = self.dt * self.f(self.t+self.dt/2,self.xm +F2/2)
        F4 = self.dt * self.f(self.t+self.dt,self.xm +F3)
        self.xm = self.xm + 1.0/6.0*(F1+2*F2+2*F3+F4)
        self.t = self.t +self.dt

        xcn = self.get_xc(self.t)
        xcb = self.get_xc(self.t-self.dt)
        self.out_t.append(self.t)
        self.out_xm.append(self.xm[0,0])
        self.out_vm.append(self.xm[1,0])
        self.out_xc.append(xcn)
        self.out_vc.append((xcn-xcb)/self.dt)
        
        
def two_cos(t,times,amp):
    """
    Simulates two cosine shaped moves. Each one is amp/2
    times[0]: start of first move
    times[1]: end of first move
    times[2]: start of second move
    times[3]: end of second move
    """
    
    def f1(t,tau=28):
        return 1.0*(t>0)

    def f2(t,tau=28):
        return (t>0)*np.cos(t/tau*np.pi/2)
    
    def f3(t,tau=28):
        return (t>0)*np.sin(t/tau*np.pi/2)


    tau1 = times[1]-times[0]
    tau2 = times[3]-times[2]
    ret= amp/2*(f1(t-times[0],tau=tau1)+f1(t-times[3],tau=tau2)
               -f2(t-times[0],tau=tau1)-f2(t-times[3],tau=tau2)
               -f3(t-times[1],tau=tau1)+f3(t-times[2],tau=tau2))
            
    return ret
    
    
        
        
class torsionpendulum:
    def __init__(self,times,funct,x0=0,v0=0,amp=1.4e-8,oramp=True,t3=0,t4=0):
        """
        Note: Input and Output (ay amd fr_a) in arcsecs, internal calculation in rads
        times is a vector of times that desribes certain things
        funct is a fcuntion that describes the external torque, it has the following arguyments
        funct(t,times,amp)
        
        """
        self.oramp = oramp
        self.rad2arcsecs = 180.0/np.pi*3600
        self.I =0.076233 # in kgm^2
        self.omega = 2*np.pi/120.0
        self.times = times
        self.funct = funct
        self.amp0 =amp # in Nm
        self.M=np.matrix(((0,1),(-self.omega**2,0)))
        self.dt = 1
        self.y = np.matrix((x0/self.rad2arcsecs,v0/self.rad2arcsecs)).T
        self.t =0
        self.ay  =[]
        self.avelo =[]
        self.at =[]
        self.at.append(self.t)
        self.ay.append(self.y[0,0]*self.rad2arcsecs)
        self.avelo.append(self.y[1,0]*self.rad2arcsecs) 
        
        
        self.aN  =[ self.torque(self.t)*1e9]
        self.sLI = [0]
        self.cLI = [0]
        self.alc = 0
        self.als = 0
        self.fit_t = []
        self.fit_y = []
        self.fr_t =[]
        self.fr_a =[]
        
    def fita(self,t,y):
        if len(self.fit_t)>0:
            if t-self.fit_t[0]>1.0/self.omega/2/np.pi:
                self.fit_t.pop(0)
                self.fit_y.pop(0)
        self.fit_t.append(t)
        self.fit_y.append(y)
        if len(self.fit_t)>3:
            wt = self.omega*np.array(self.fit_t)
            O = np.ones(len(wt))
            X = np.matrix(np.vstack ((O,np.cos(wt),np.sin(wt) )) ).T
            d = np.matrix(self.fit_y).T
            result = (X.T*X).I*(X.T*d)
            self.fr_t.append(self.fit_t[-1])
            self.fr_a.append(np.sqrt(result[1,0]**2+result[2,0]**2))
            
            

    def go(self,tstop):
        while self.t<tstop:
            self.iterate()
        self.ay = np.array(self.ay)
        self.av = (self.ay[1:]-self.ay[0:-1])/self.dt
        self.av = np.hstack((self.av[0], self.av))
    
    def f1(self,t,tau=28):
        return 1.0*(t>0)

    def f2(self,t,tau=28):
        return (t>0)*np.cos(t/tau*np.pi/2)
    
    def f3(self,t,tau=28):
        return (t>0)*np.sin(t/tau*np.pi/2)

    
    def torque(self,t):
        return self.funct(t,self.times,self.amp0)

    
    def N(self,t):
        n = self.torque(t)
        return np.matrix((0,n/self.I)).T

    def f(self,t,y):
        return self.N(t)+self.M*y
    
    def iterate(self):
        F1 = self.dt * self.f(self.t         ,self.y)
        F2 = self.dt * self.f(self.t+self.dt/2,self.y +F1/2)
        F3 = self.dt * self.f(self.t+self.dt/2,self.y +F2/2)
        F4 = self.dt * self.f(self.t+self.dt,self.y +F3)
        self.y = self.y + 1.0/6.0*(F1+2*F2+2*F3+F4)
        self.t = self.t +self.dt
        self.ay.append(self.y[0,0]*self.rad2arcsecs)
        self.avelo.append(self.y[1,0]*self.rad2arcsecs) 
        self.at.append(self.t)
        self.fita(self.t,self.y[0,0]*self.rad2arcsecs)
        self.aN.append(self.torque(self.t)*1e9)
                
        
        
        
        