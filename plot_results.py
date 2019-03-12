import simulation
import parameters
import plot as pl
import numpy as np
import pandas as pd
s = simulation.simulate_cyano()
p = parameters.cyano_pars()

### Growth laws ####################################################################################
### Growth rate versus External nitrogen
#yN = s.vary_Nx(model=s.m, Cx=0.25, ub_Nx=20.0, stepsize=0.001) # returns a dataframe with variable values
yN = pd.read_csv('results/yN.csv')	#load simulation results from csv file

### Growth rate versus External carbon
#yC = s.vary_Cx(model=s.m, Nx=0.5, ub_cix=20.0, stepsize=0.001)  # returns a dataframe with variable values
yC = pd.read_csv('results/yC.csv')	#load simulation results from csv file


### Growth rate versus irradiance
#yI = s.vary_irradiance(model=s.m, Cx=0.25, Nx=0.5, ub_Irr=1000.0) # returns a dataframe with variable values
yI = pd.read_csv('results/yI.csv')	#load simulation results from csv file

pl.growth_subplots(yN, yC, yI)


### Resource allocation ############################################################################
### protein allocation
#pF_N = s.p_fracs(yN, s.m)
#pF_C = s.p_fracs(yC, s.m)
#pF_I = s.p_fracs(yI, s.m)
#pl.pie_subplots(pF_N, pF_C, pF_I)

### Vmax vs uptake flux
#y1 = s.vary_Nx(s.m, steps=20, lb_Nx=0.1, Cx=0.25, ub_Nx=1000)
y1 = pd.read_csv('results/y1.csv')	#load simulation results from csv file
#y2 = s.vary_Nx(s.m, steps=20, lb_Nx=0.1, Cx=0.75, ub_Nx=1000)
y2 = pd.read_csv('results/y2.csv')	#load simulation results from csv file
#y3 = s.vary_Nx(s.m, steps=20, lb_Nx=0.1, Cx=15.0, ub_Nx=1000)
y3 = pd.read_csv('results/y3.csv')	#load simulation results from csv file
pl.twosubplots(y1, y2, y3, s.m)


### Gleaner-Opportunist ############################################################################
### Growth rate versus external nitrogen
#y1, y2 = s.twostrains(model=s.m, I=200.0, Cx=0.2, ub_Nx=1.0)
yG = pd.read_csv('results/yG.csv')	#load simulation results from csv file
yO = pd.read_csv('results/yO.csv')	#load simulation results from csv file
pl.plot_gle_opp (Nx1=yG.Nx, Nx2=yO.Nx, mu1=86400*yG.Mu, mu2=86400*yO.Mu)

### Metabolic diversity and cost of regulation
#y_Y = s.vary_Nx(s.m0, lb_Nx=0.1, stepsize=0.001, steps=40)
#y_Y = pd.read_csv('results/y_Y.csv')	#load simulation results from csv file
#y_N = s.vary_Nx(s.m, lb_Nx=0.1, stepsize=0.001, steps=40)
#y_N = pd.read_csv('results/y_N.csv')	#load simulation results from csv file
#y_NY = s.vary_Nx(s.m1, lb_Nx=0.1, stepsize=0.001, steps=40)
#y_NY = pd.read_csv('results/y_NY.csv')	#load simulation results from csv file
#pF_NY = s.protein_fracs(y_NY, s.m)
#pl.subplots_TnTy(y_N, y_Y, y_NY, pF_NY, s.m)

### Time-dependent oscillations
# y = s.rk4(f=s.chemostat, x0=[1., 1., 10., 200.], t=np.linspace(0, 1000, 10000)    # uses 4th order runge-kutta method to solve the ODEs
y1 = pd.read_csv('results/chemostat_fixedI.csv')	#load simulation results from csv file
y2 = pd.read_csv('results/chemostat_variableI.csv')	#load simulation results from csv file
pl.co_culture(y1, y2, s.I(y2.t))
