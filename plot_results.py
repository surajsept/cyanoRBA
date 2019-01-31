import simulation
import parameters
import plot as pl
import numpy as np
s = simulation.simulate_cyano()
p = parameters.cyano()

### Growth laws ####################################################################################
# Growth rate versus External nitrogen
yN = s.vary_Nx(model=s.m, Cx=0.2, ub_Nx=20.0) # returns a dataframe with variable values
pl.basic_plot (yN.Nx / p.K_N, yN.Mu, xlabel='N$_x$/K$_N$', ylabel='Growth rate, $\mu$ (d$^{-1}$)',
               text=['C$_x$ = 0.2 $\mu$ M', 'I = 200 $\mu$ E m$^{-2}$ s$^{-1}$'])

# Growth rate versus External carbon
yC = s.vary_Cx(model=s.m1, Cx=0.2, ub_Nx=20.0)  # returns a dataframe with variable values
pl.basic_plot (yC.Cx / p.Km_C, yC.Mu, xlabel='C$_x$/K$_N$', ylabel='Growth rate, $\mu$ (d$^{-1}$)',
               text=['N$_x$ = 0.2 $\mu$ M', 'I = 200 $\mu$ E m$^{-2}$ s$^{-1}$'])

# Growth rate versus irradiance
yI = s.vary_irradiance(model=s.m, Cx=0.2, Nx=0.2, ub_Irr=1000.0) # returns a dataframe with variable values
pl.basic_plot (yI.Irr, yI.Mu, xlabel='I ($\mu$ E m$^{-2}$ s$^{-1}$)', text_lx=0.5,
               ylabel='Growth rate, $\mu$ (d$^{-1}$)',
               text=['C$_x$ = 0.2 $\mu$ M', 'N$_x$ = 0.2 $\mu$ M'])

### Resource allocation ############################################################################
# protein allocation
pF_N = s.p_fracs(yN, s.m)
pF_C = s.p_fracs(yC, s.m)
pF_I = s.p_fracs(yI, s.m)
pl.pie_subplots(pF_N, pF_C, pF_I)

# Vmax vs uptake flux
y1 = s.vary_Nx(s.m, steps=20, lb_Nx=0.1, Cx=0.25, ub_Nx=1000)
y2 = s.vary_Nx(s.m, steps=20, lb_Nx=0.1, Cx=0.75, ub_Nx=1000)
y3 = s.vary_Nx(s.m, steps=20, lb_Nx=0.1, Cx=15.0, ub_Nx=1000)
pl.twosubplots(y1, y2, y3, s.m)


### Gleaner-Opportunist ############################################################################
# Growth rate versus external nitrogen
y1, y2 = s.twostrains(model=s.m, I=200.0, Cx=0.2)
pl.plot_gle_opp (Nx1=y1.Nx, Nx2=y2.Nx, mu1=y1.Mu, mu2=y2.Mu)

# Metabolic diversity and cost of regulation
y_Y = s.vary_Nx(s.m0, lb_Nx=0.1, stepsize=0.001, steps=40)
y_N = s.vary_Nx(s.m, lb_Nx=0.1, stepsize=0.001, steps=40)
y_NY = s.vary_Nx(s.m1, lb_Nx=0.1, stepsize=0.001, steps=40)
pF_NY = s.protein_fracs(y_NY, s.m)
pl.subplots_TnTy(y_N, y_Y, y_NY, pF_NY, s.m)

# Time-dependent oscillations
# uncomment the line below to run time-dependent simulation.
# t, y = s.culture(y0=[1, 1, 1], tfinal=1000.0, t0=0.0)   # uses a CVODE (Sundials) integrator
t = np.linspace(0, 10, 100)
x = s.rk2a(s.dummyF, [1, 1, 1], t)  # Second-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0
