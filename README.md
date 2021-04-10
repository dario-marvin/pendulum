# pendulum

These files go with the article *The crane operator's trick an other shenanigans with a pendulum*

## pendulum.py
The file pendulum.py is a python module with several classes.

### simplependulum
The class simplependulum can be used to simulate a simple pendulum, the crane in the article.
The following differential equation is solve numerically.


$m\ddot{x}_m + c\dot{x}_m + m \frac{g}{l}(x_m-x_c) = 0$ or  $\ddot{x}_m +\frac{c}{m} \dot{x}_m+ \omega_o^2 (x_m-x_c)=0$.

With  $\omega_o^2=\frac{g}{l}$ and $\xi = \frac{c}{2 m\omega_o}$, it reads,

$$
\ddot{x}_m =- 2\xi\omega_o \dot{x}_m -\omega_o^2 x_m + \omega_o^2 x_c.
$$

This second order differential equation can be reduced to a first order differential equation by using the variable

$x=\left( \begin{array}{c} x_m\\\dot{x}_m\end{array} \right)$

Hence,
$$\left( \begin{array}{c} \dot{x}_m\\\ddot{x}_m\end{array} \right) = 
\left( \begin{array}{cc} 0 &1\\ -\omega_o^2 & -2\xi\omega_o\end{array} \right) 
\left( \begin{array}{c} x_m\\\dot{x}_m\end{array} \right) + 
\omega^2 \left( \begin{array}{c} 0\\ x_c \end{array} \right) $$

In other words,
$$
\dot{x} =M x + \omega^2 F.
$$

This equation is solved with the Runge-Kutta method, see [Wiki](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)


The simulation is setup by

`pend = pendulum.simplependulum(x0=1,v0=0,w0=2*np.pi/T0,xi=0.1,trolley_pos=func)`

where x0 and v0 give the initial conditionof the pendulum at t=0 and w0 and xi are the paramters discussed above. The argument trolley_pos points at a fucntion that returns the trolley position as a function of time. If it's left out the trolley is assumed to be at 0 at all times.


One the simulation is setup, it can be run, by calling

`pend.go(30)`

where the argument is the number of simulated seconds the simulation runs. In this case the trajectory of the pendulum is calculated for 30 s.The function does not return anything. The ouptut of the simulation is stored in five lists owned by the instance. They are


- `pend.out_t`  : time
- `pend.out_xm`   : position of the mass
- `pend.out_vm`   : velocity of the mass
- `pend.out_xc`   : position of the trolley
- `pend.out_vc`   : velocity of the trolley



### Jupyter notebooks

Several jupyter notebooks demonstrate the pendulum class.


#### section2.ipynb
Shows simulations of the topics that are discussed in section 2 of the article.


#### section3.ipynb
Shows simulations of the topics that are discussed in section 3 of the article.