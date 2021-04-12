# README.md
The following are supplemental materials for the article entitled *The Crane Operator's Tricks and other Shenanigans with a Pendulum*.


**pendulum.py**  Python module with several classes.

**section2.ipynb**  Jupyter notebook that shows the physics in section 2 of the paper.

**section3.ipynb**  Jupyter notebook that shows the physcis in section 3 of the paper.

**section4.ipynb**  Jupyter notebook that shows the physics in section 4 of the paper.

**section2.py**  Python script that shows the physics in section 2 of the paper.

**section3.py**  Python script that shows the physics in section 3 of the paper.

**section4.py**  Python script that shows the physics in section 4 of the paper.

**supplemental_material.pdf**  More math on the torsion pendulum



## Description of the classes in pendulum.py

### simplependulum
This class is used to simulate the dynamics of a simple pendulum, specifically the crane in the article. The math used in this simulation is summarized in Section II of supplemental_material.pdf.

The setup of the simulation is as follows:

`pend = pendulum.simplependulum(x0=1,v0=0,w0=2*np.pi/T0,xi=0.1,trolley_pos=func)`

where x0 and v0 are the initial position (m) and velocity (m/s) of the pendulum, respectively, at t=0, w0 is the natural frequency (rad/s), and xi is the unitless damping ratio. The argument trolley_pos points at a fucntion that returns the trolley position as a function of time. If omitted the trolley is assumed to be at 0 at all times.

To execute the simulation, run the following:

`pend.go(30)`

where the argument is the amount of time the pendulum is simulated for in seconds. In this case the trajectory of the pendulum is calculated for 30 s. The function does not return anything and the ouptut of the simulation is stored in five lists owned by the instance:

- `pend.out_t`    : time
- `pend.out_xm`   : position of the mass
- `pend.out_vm`   : velocity of the mass
- `pend.out_xc`   : position of the trolley
- `pend.out_vc`   : velocity of the trolley

### torsionpendulum

The setup of this class is similar to simplependulum and is as follows:

`pend = pendulum.torsionpendulum(theta0=0,theta_dot0=0,w0=0.0523598,xi=0,I=0.076233,ext_torque=None,dt=1)`

where theta0 and theta_dot0 are the initial angle (rad) and angular velocity (rad/s), respectively, at t=0, w0 is the natural frequency (rad/s), xi is the unitless damping ratio, and I is the moment of inertia (kg m^2).

Similar to simplependulum, to execute the simulation, run the following:

`pend.go(30)`

The five output lists are:

- `pend.out_t`         : time
- `pend.out_theta`     : theta (angular excursion of the torsion pendulum)
- `pend.out_theta_dot` : time derivative of theta
- `pend.out_n`         : external torque
- `pend.out_n_dot`     : time derivative of the external torque


