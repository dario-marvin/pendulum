"""
pendulum.py 

implements simplependulum

simulation classes: 
class  simplependulum(x0=0,v0=0,w0=0.700249,xi=0,trolley_pos=None,dt =0.005)
class  torisonpendulum(theta0=0,theta_dot0=0,w0=0.0523598,xi=0,\
                       I=0.076233,ext_torque=None,dt =1))

stimulus classes:
class const_vel(ts,vs)
class const_acc(ts,as)
cass  cos_tor(ts,amp) #ts: four long, amp: scalar


"""

import numpy as np
import matplotlib.pyplot as plt
import os

def segment(p1, p2):
    """
    Helpful for drawing system
    """
    return [p1[0], p2[0]], [p1[1], p2[1]]

class const_vel:
    def __init__(self, ts, vs):
        """
        ts and vs are both lists, the class sorts the lists.
        while maintainint the pairs for the sorted pairs tp,vp
        if t>tp[i] then vp[i] will be returned
        if t<tp[0] 0 will be returned
        """
        self.tp, self.vp = zip(*sorted(zip(ts, vs), key=lambda x: x[0]))
        self.mint = min(self.tp)

    def xc(self, t):
        if t < self.mint:
            return 0
        x = 0  # we always start at 0
        ot = 0
        ov = 0
        for t_, v_ in zip(self.tp, self.vp):
            if t >= t_:
                x = x + (t_ - ot) * ov
            else:
                return x + (t - ot) * ov
            ot = t_
            ov = v_
        return x + (t - t_) * v_


class const_acc:
    def __init__(self, ts, a):
        """
        ts and a are both lists, the class sorts the lists.
        while maintainint the pairs for the sorted pairs tp,ap
        if t>ap[i] then ap[i] will be returned
        if t<tp[0] 0 will be returned
        """
        self.tp, self.ap = zip(*sorted(zip(ts, a), key=lambda x: x[0]))
        self.mint = min(self.tp)

    def xc(self, t):
        if t < self.mint:
            return 0
        x = 0  # we always start at 0
        ot = 0
        oa = 0
        v = 0
        for t_, a_ in zip(self.tp, self.ap):
            if t >= t_:
                x = x + 0.5 * oa * (t_ - ot) ** 2 + (t_ - ot) * v
                v = v + (t_ - ot) * oa
            else:
                return x + 0.5 * oa * (t - ot) ** 2 + (t - ot) * v
            ot = t_
            oa = a_
        return x + 0.5 * a_ * (t - ot) ** 2 + (t - ot) * v


class cos_tor:
    def __init__(self, ts, amp):
        """
        Simulates two cosine shaped moves. Each one is amp/2
        ts[0]: start of first move
        ts[1]: end of first move
        ts[2]: start of second move
        ts[3]: end of second move
        """
        self.ts = ts
        self.amp = amp

    def f1(self, t):
        return 1.0 * (t > 0)

    def f2(self, t, tau):
        return (t > 0) * np.cos(t / tau * np.pi / 2)

    def f3(self, t, tau):
        return (t > 0) * np.sin(t / tau * np.pi / 2)

    def n(self, t):
        tau1 = self.ts[1] - self.ts[0]
        tau2 = self.ts[3] - self.ts[2]
        return self.amp / 2 * (self.f1(t - self.ts[0]) + self.f1(t - self.ts[3])
                               - self.f2(t - self.ts[0], tau=tau1) - self.f2(t - self.ts[3], tau=tau2)
                               - self.f(t - self.ts[1], tau=tau1) + self.f3(t - self.ts[2], tau=tau2))


class simplependulum:
    def __init__(self, x0=0, v0=0, w0=0.700249, xi=0, trolley_pos=None, dt=0.005):
        """
        x0,v0       : pendulum initial conditions x0 in m, v0 in m/s
        w0          : angular frequncy of the pendulum sqrt(g/l)
        xi          : damping ratio xi = c/(2m w0), xi =1 critically damped
        trolley_pos : a function f(t) that returns the positon of the trolley
                      as a function of t
        """
        self.omega = w0
        self.xi = xi
        self.M = np.matrix(((0, 1), (-self.omega ** 2, -2 * self.omega * self.xi)))
        self.dt = dt  # time step in the simulation
        self.xm = np.matrix((x0, v0)).T
        self.t = 0
        self.xc = trolley_pos
        self.L = 9.80665 / (self.omega * self.omega)

        # output arrays
        self.out_t = []  # time
        self.out_xm = []  # x positon of mass
        self.out_vm = []  # vel of mass
        self.out_xc = []  # x position of trolley
        self.out_vc = []  # velocity of trolley

        xc0 = 0
        if self.xc != None:
            xc0 = self.xc(self.t)
        self.out_t.append(self.t)
        self.out_xm.append(self.xm[0, 0])
        self.out_vm.append(self.xm[1, 0])
        self.out_xc.append(xc0)
        self.out_vc.append(0)

    def go(self, tstop):
        while self.t < tstop:
            self.iterate()

    def get_xc(self, t):
        if self.xc != None:
            return self.xc(t)
        else:
            return 0

    def f(self, t, y):
        F = np.matrix((0, self.get_xc(t))).T  # Matrix for the ext. force.
        return self.omega ** 2 * F + self.M * y

    def iterate(self):
        F1 = self.dt * self.f(self.t, self.xm)
        F2 = self.dt * self.f(self.t + self.dt / 2, self.xm + F1 / 2)
        F3 = self.dt * self.f(self.t + self.dt / 2, self.xm + F2 / 2)
        F4 = self.dt * self.f(self.t + self.dt, self.xm + F3)
        self.xm = self.xm + 1.0 / 6.0 * (F1 + 2 * F2 + 2 * F3 + F4)
        self.t = self.t + self.dt

        xcn = self.get_xc(self.t)
        xcb = self.get_xc(self.t - self.dt)
        self.out_t.append(self.t)
        self.out_xm.append(self.xm[0, 0])
        self.out_vm.append(self.xm[1, 0])
        self.out_xc.append(xcn)
        self.out_vc.append((xcn - xcb) / self.dt)

    def draw(self, prefix: str = '', title: str = None, save=True, ax=None):

        min_x = min(np.min(self.out_xc), np.min(self.out_xm)) - 1
        max_x = max(np.max(self.out_xc), np.max(self.out_xm)) + 1
        max_y = 1
        min_y = -self.L - 2
        skip = 40

        for i in range(0, len(self.out_xc), skip):
            print(f'{int(i/skip)}/{int(len(self.out_xc)/skip)}')
            fig, ax = plt.subplots()
            ax.set_title(title)
            ax.axis('off')
            ax.set_aspect('equal')

            ax.set_xlim((min_x, max_x))
            ax.set_ylim((min_y, max_y))

            ax.plot(*segment((min_x, 0), (max_x, 0)), color='black')
            lcl_t = self.out_t[i]
            lcl_xm = self.out_xm[i]  # x position of pendulum mass
            lcl_xc = self.out_xc[i]  # x position of trolley
            lcl_vc = self.out_vc[i]  # x speed of trolley
            lcl_theta = np.arcsin((lcl_xm - lcl_xc) / self.L)

            lcl_ym = -self.L * np.cos(lcl_theta)

            pend = plt.Circle((lcl_xm, lcl_ym), 0.3, color='red', fill=True)
            trolley = plt.Rectangle((lcl_xc-0.3, -0.2), 0.6, 0.4, color='blue', fill=True)
            ax.add_patch(pend)
            ax.add_patch(trolley)

            ax.plot(*segment((lcl_xc, 0), (lcl_xm, lcl_ym)), color='black')
            ax.text(min_x, min_y + 0.5, f't: {lcl_t:.1f} s', color='black', size=10)
            ax.text(lcl_xc-1, 0.6, f'{lcl_vc:.1f} m/s', color='black', size=10)

            fig.tight_layout()
            fig.savefig(f'../output/{prefix + str(int(i/skip)).zfill(5)}.png')
            plt.clf()
            plt.cla()
            plt.close()

    def make_mp4(self, fn: str = 'system', remove_images=False, framerate=40):
        os.chdir('../output')
        os.system(
            f"ffmpeg -pattern_type glob -framerate {framerate} -i '*.png' -c:v libx264 -crf 18 -pix_fmt yuv420p {fn}.mp4")
        if remove_images:
            os.system('rm *.png')


class torsionpendulum:
    def __init__(self, theta0=0, theta_dot0=0, w0=0.0523598, xi=0, \
                 I=0.076233, ext_torque=None, dt=1):
        """
        simulates a torsion pendulum
        theta,theta_dot0  : pendulum initial conditions x0 in rad, v0 in rad/s
        w0                : angular frequncy of the pendulum sqrt(kappa/I)
        xi                : damping ratio xi = c/(2I w0), xi =1 critically damped
        I                 : moment of inertia in kgm^2
        ext_troque        : a function f(t) that returns the external torque
                      as a function of t in Nm
        """

        self.omega = w0
        self.xi = xi
        self.M = np.matrix(((0, 1), (-self.omega ** 2, -2 * self.omega * self.xi)))
        self.dt = dt  # time step in the simulation
        self.I = I  # moment of inertia
        self.theta = np.matrix((theta0, theta_dot0)).T
        self.t = 0
        self.n = ext_torque

        # output arrays
        self.out_t = []  # time
        self.out_theta = []  # positon of the pendulum, theta
        self.out_theta_dot = []  # velocity of the pendulum
        self.out_n = []  # ext. torque
        self.out_n_dot = []  # derivative of the ext. torque

        n0 = 0
        if self.n != None:
            n0 = self.n(self.t)
        self.out_t.append(self.t)
        self.out_theta.append(self.theta[0, 0])
        self.out_theta_dot.append(self.theta[1, 0])
        self.out_n.append(n0)
        self.out_n_dot.append(0)

    def go(self, tstop):
        while self.t < tstop:
            self.iterate()

    def torque(self, t):
        return self.funct(t, self.times, self.amp0)

    def get_n(self, t):
        if self.n != None:
            return self.n(t)
        else:
            return 0

    def f(self, t, y):
        N = np.matrix((0, self.get_n(t) / self.I)).T
        return N + self.M * y

    def iterate(self):
        F1 = self.dt * self.f(self.t, self.theta)
        F2 = self.dt * self.f(self.t + self.dt / 2, self.theta + F1 / 2)
        F3 = self.dt * self.f(self.t + self.dt / 2, self.theta + F2 / 2)
        F4 = self.dt * self.f(self.t + self.dt, self.theta + F3)
        self.theta = self.theta + 1.0 / 6.0 * (F1 + 2 * F2 + 2 * F3 + F4)
        self.t = self.t + self.dt

        nn = self.get_n(self.t)
        nb = self.get_n(self.t - self.dt)
        self.out_t.append(self.t)
        self.out_theta.append(self.theta[0, 0])
        self.out_theta_dot.append(self.theta[1, 0])
        self.out_n.append(nn)
        self.out_n_dot.append((nn - nb) / self.dt)


class find_move:
    def __init__(self, deltat=19.2, w0=0.0523598, na=3.11e-08, I=0.076233):

        """
        class to find the moves to damp a torsion pendulum to any desired
        amplitude
        deltat : time it takes for a move
        w0     : resonance frequency
        na     : amplitude of the external torque
        I      : moment of inertia
        """
        self.I = I
        self.deltat = deltat
        self.w0 = w0
        self.na = na

    def vel(self, t1, t3, v0=1.0):
        """
        calculates the velocity amplitude of the pendulum after both moves
        first   move from t1 to t1+deltat
        second move from t3 to t3+deltat
        v0 is the initial velocity amplitude (when the pendulum swings
            through the equilibrium pos)
        """
        pi = np.pi  # make the eq. easier to read
        w0 = self.w0
        I = self.I
        na = self.na

        ka = w0 * w0 * I  # kappa
        c1 = np.cos(t1 * w0)
        s1 = np.sin(t1 * w0)
        c3 = np.cos(t3 * w0)
        s3 = np.sin(t3 * w0)

        w0t = w0 * self.deltat  # pi*w0*2*dt/(2*pi)
        c1t = np.cos(t1 * w0 + w0t)
        s1t = np.sin(t1 * w0 + w0t)
        c3t = np.cos(t3 * w0 + w0t)
        s3t = np.sin(t3 * w0 + w0t)

        cc = v0 - na * w0 * pi / (2 * ka * (pi * pi - 4 * w0t * w0t)) * (2 * w0t * (-c3 + c1t) + pi * (s1 + s3t))
        sc = na * w0 * pi / (2 * ka * (pi * pi - 4 * w0t * w0t)) * (pi * (c1 + c3t) + 2 * w0t * (s3 - s1t))
        return np.sqrt(cc * cc + sc * sc)

    def find_t1_t3(self, desired_vamp=1e-9, v0=1):
        """
        returns the best t1 and t3  when the moves start to
        obtain the desired velocity amp
        """
        best_damp = 9e99
        best_amp = 9e99
        tol_amp = 1e-9
        T0 = np.pi * 2 / self.w0
        for t1f in np.arange(0, T0 - self.deltat, 0.4):
            for t3f in np.arange(t1f + self.deltat, T0, 0.4):
                amp = self.vel(t1f, t3f, v0)
                damp = np.sqrt((amp - desired_vamp) ** 2)
                if damp < best_damp:
                    best_t1 = t1f
                    best_t3 = t3f
                    best_damp = damp
                    best_amp = amp
                if damp < tol_amp:
                    break
        return best_t1, best_t3, best_amp
