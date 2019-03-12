import numpy as np
import pandas as pd
import time

import cyano_v1
import cyano_v1_0
import cyano_v1_1


### Aggregates basic functions to simulate the cyano model ###################################
class simulate_cyano():
    def __init__(self):
        self.m = cyano_v1  # RBA model of (T_N)-strain
        self.m0 = cyano_v1_0  # RBA model of (T_Y)-strain
        self.m1 = cyano_v1_1  # RBA model of (T_N + T_Y)-strain
        return print('Cyano model is ready to simulate')

    ####################################### Other simulations ######################################

    def vary_Nx(self, model, I=200.0, Cx=0.25, lb_Nx=0.01, ub_Nx=20.0, steps=20, stepsize=0.1,
                newpars={}):
        '''
		To simulate specific growth rate at diffrent conecntration of external Nitrogen
		:param model: model instance of the class (e.g. instatiating the current class returns
		cy.m, cy.m0 and cy.m1)
		:param I: light intensity. A float value
		:param Cx: external inorganic carbon concentration (microM). A float value
		:param lb_Nx: lowest concentration of external nitrogen. A float value
		:param ub_Nx: highest concentration of exteral nitrogen. A float value
		:param steps: number of points for which the solution is expected. Integer value
		:param stepsize: incremental concentration value. A float value
		:param newpars: new parameters if necessary. A python dictionary
		:return: df: A dataframe. DataFrame is a 2-dimensional labeled data structure with
		columns of potentially different types.
		'''
        if ub_Nx <= 20.0:
            Nx = np.hstack((np.arange(lb_Nx, 0.2 * ub_Nx, stepsize * 0.2 * ub_Nx),
                            np.linspace(0.2 * ub_Nx, ub_Nx, steps)))
        else:
            Nx = np.hstack((np.linspace(lb_Nx, 1.0, steps),
                            np.linspace(1.0, 10.0, steps),
                            np.linspace(10.0, 100.0, steps),
                            np.linspace(100.0, ub_Nx, steps)))
        Mu, Flux = [], []
        for i in range(len(Nx)):
            mu, flux = model.binary_search(mu=0.01, Cx=Cx, Nx=Nx[i], I=I, newpars=newpars)
            Mu.append(mu)
            Flux.append(flux)
        y = np.hstack((np.transpose(np.vstack((Nx, Mu))), np.array(Flux)))
        df = pd.DataFrame(y, columns=['Nx', 'Mu'] + model.x_name)
        return df

    def vary_Cx(self, model, Nx=0.5, I=200.0, lb_cix=0.1, ub_cix=15.0, steps=20, newpars={},
                stepsize=0.1):
        '''
		To simulate the phototrophic growth at different concentration of external inorganic carbon
		:param model: the model instance
		:param Nx: concentration of external nitrogen (microM). A float value
		:param I: light intensity. A float value
		:param lb_cix: lowest concentration of external carbon. A float value
		:param ub_cix: highest concentration of external carbon. A float value
		:param steps: number of points for which the solution is expected. Integer value
		:param stepsize: incremental concentration value. A float value
		:param newpars: new parameters if necessary. A python dictionary
		:return: A data frame
		'''
        cix = np.hstack((np.arange(lb_cix, 0.2 * ub_cix, stepsize * 0.2 * ub_cix),
                         np.linspace(0.2 * ub_cix, ub_cix, steps)))
        Mu, Flux = [], []
        for i in range(len(cix)):
            mu, flux = model.binary_search(mu=0.01, Cx=cix[i], Nx=Nx, I=I, newpars=newpars)
            Mu.append(mu)
            Flux.append(flux)
        y = np.hstack((np.transpose(np.vstack((cix, Mu))), np.array(Flux)))
        df = pd.DataFrame(y, columns=['Cx', 'Mu'] + model.x_name)
        return df

    def vary_irradiance(self, model, Cx=0.3, Nx=0.2, lb_Irr=0.1, ub_Irr=1000.0, steps=20,
                        newpars={}):
        '''
		To simulate phototrophic growth at different light intensities
		:param model: model instance
		:param Cx: concentration of external inorganic carbon
		:param Nx: concentration of external inorganic carbon
		:param lb_Irr: lowest light intensity. A float
		:param ub_Irr: highest light intensity. A float
		:param steps: number of points for which the solution is expected. Integer value
		:param newpars: new parameters if necessary. A python dictionary
		:return: A dataframe
		'''
        Irr = np.hstack((np.linspace(lb_Irr, 0.2 * ub_Irr, steps),
                         np.linspace(0.2 * ub_Irr, ub_Irr, steps)))
        Mu, Flux = [], []
        for i in range(len(Irr)):
            mu, flux = model.binary_search(mu=0.01, Cx=Cx, Nx=Nx, I=Irr[i], newpars=newpars)
            Mu.append(mu)
            Flux.append(flux)
        y = np.hstack((np.transpose(np.vstack((Irr, Mu))), np.array(Flux)))
        df = pd.DataFrame(y, columns=['Irr', 'Mu'] + model.x_name)
        return df

    def p_fracs(self, y, model):
        '''
		To get relative concentrations of different metabolic proteins
		:param y: a dataframe with solutions of phototrotrophic growth
		:param model: model instance
		:return: numpy n-dimensional array of size similar to the number of columns in the dataframe
		'''
        a = np.transpose(np.array([y.Tc, y.Tn, y.PSU, y.R, y.Pq, y.CB, y.Mc, y.Mq]))
        alpha_P = np.array([model.p.nTc, model.p.nTn, model.p.nPSU, model.p.nR, model.p.nPq,
                            model.p.nCB, model.p.nMc, model.p.nMq])
        aa = a * alpha_P
        pFrac = np.transpose(np.vstack((np.transpose(aa[:, :-2]), np.sum(aa[:, -2:], axis=1))))
        return pFrac

    def protein_fracs(self, y, model):
        '''
		To get relative concentrations of different metabolic proteins
		:param y: a dataframe with solutions of phototrotrophic growth
		:param model: model instance
		:return: a dataframe with relative concentrations of metabolic proteins
		'''
        labels = ['T_C', 'T_N1', 'T_N2', 'PSU', 'R', 'P_Q', 'CB', 'M_c', 'M_q']
        a = np.transpose(np.array([y.Tc, y.Tn1, y.Tn2, y.PSU, y.R, y.Pq, y.CB, y.Mc, y.Mq]))
        alpha_P = np.array([model.p.nTc, model.p.nTn1, model.p.nTn2, model.p.nPSU, model.p.nR,
                            model.p.nPq, model.p.nCB, model.p.nMc, model.p.nMq])
        aa = a * alpha_P
        pFrac = (aa.T / aa.sum(axis=1)).T
        df = pd.DataFrame(pFrac, columns=labels)
        return df

    ####################################### Opportunist vs Gleaner #################################

    def twostrains(self, model, I=200.0, ub_Nx=1.0, steps=20, Cx=0.3, lb_Nx=1e-4):
        '''
		Simulates the specific growth rates of gleaner and opportunist
		:param model: model instance
		:param I: light intensity
		:param ub_Nx: highest concentration of external nitrogen
		:param steps:
		:param Cx: concentration of external inorganic carbon
		:param lb_Nx: lowest concentration of external nitrogen
		:return: y1, y2: dataframes with solutions in regard to phototrophic growth of gleaner
		and opportunist
		'''
        K_N2 = 50.0  # micro mol per litre
        k5_2 = 20.0  # per second
        k2_2 = 200.0  # per second
        y1 = self.vary_Nx(model, I=I, Cx=Cx, lb_Nx=lb_Nx, ub_Nx=ub_Nx,
                          newpars={}, steps=steps)
        y2 = self.vary_Nx(model, I=I, Cx=Cx, lb_Nx=lb_Nx, ub_Nx=ub_Nx,
                          newpars={'Km_N': K_N2, 'k2': k2_2, 'k5': k5_2}, steps=steps)
        return y1, y2

    ####################################### Chemostat simulation ###################################

    def I(self, t, I0=100):
        return abs(I0 + 0.01 + I0 * np.sin(2 * np.pi * t / 365))

    ### simulating co-existence using RBA model of cyanobacteria
    def strains(self, Nx, I):
        '''
		Two strains of cyanobacteria referred to as gleaner and opportunist in the Manuscript
		:param Nx: External nitrogen concentration (in microM). Float value
		:param I: Light intensity (in microE per m^2 per second). Float value.
		:return: mu1 (per second), mu2 (per second) = specific growth rate of gleaner, specific
		growth rate of opportunist
		'''
        Cx = 0.3  # micro mol per litre
        K_N2 = 50.0  # micro mol per litre
        k5_2 = 20.0  # per second
        k2_2 = 200.0  # per second
        # strain 1: gleaner
        mu1, flux1 = self.m.binary_search(mu=0.01, Cx=Cx, Nx=Nx, I=I, newpars={})

        # strain 2: opportunist
        mu2, flux2 = self.m.binary_search(mu=0.01, Cx=Cx, Nx=Nx, I=I,
                                          newpars={'Km_N': K_N2, 'k2': k2_2, 'k5': k5_2})
        return mu1, mu2, flux1[4], flux2[4]

    def chemostat(self, y, t):
        """ODEs used to calculate the time-dependent changes
		INPUT:
			t	- time point
			y	- initial conditions. y can be an array of floats or integers.
		OUTPUT:
			dydt	- y(t) at time point t. dydt is an array of floats or integers.

		"""
        Y = 2.786e10
        d = 0.25  	# per day
        Nx = 10.0  	# micro mol per litre
        # Irr = abs (100 + 0.01 + 100 * np.sin (2 * np.pi * t / 24))    # daily cycle
        # Irr = abs(100 + 0.01 + 100 * np.sin(2 * np.pi * int(t) / 365))  # yearly cycle
        Irr = 200.0
        mu1, mu2, vn1, vn2 = self.strains(Nx=y[2], I=Irr)
        mu1 = 86400 * mu1			# per day
        mu2 = 86400 * mu2			# per day
        vn1 = 86400 * vn1 / 6e17		# micro mol per day
        vn2 = 86400 * vn2 / 6e17		# micro mol per day
        dRho1_dt = mu1 * y[0] - d * y[0]	# number of cells
        dRho2_dt = mu2 * y[1] - d * y[1]	# number of cells
        dNdt = d * (Nx - y[2]) - vn1 * y[0] - vn2 * y[1]  # micro mol per litre
        dydt = [dRho1_dt, dRho2_dt, dNdt]
        return dydt

    def rk2a(self, f, x0, t):
        """Second-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

		USAGE:
			x = rk2a(f, x0, t)

		INPUT:
			f     - function of x and t equal to dx/dt.  x may be multivalued,
					in which case it should a list or a NumPy array.  In this
					case f must return a NumPy array with the same dimension
					as x.
			x0    - the initial condition(s).  Specifies the value of x when
					t = t[0].  Can be either a scalar or a list or NumPy array
					if a system of equations is being solved.
			t     - list or NumPy array of t values to compute solution at.
					t[0] is the the initial condition point, and the difference
					h=t[i+1]-t[i] determines the step size h.

		OUTPUT:
			x     - NumPy array containing solution values corresponding to each
					entry in t array.  If a system is being solved, x will be
					an array of arrays.

		NOTES:
			This version is based on the algorithm presented in "Numerical
			Analysis", 6th Edition, by Burden and Faires, Brooks-Cole, 1997.
		"""
        start_time = time.time()
        n = len(t)
        x = np.array([x0] * n)
        for i in range(n - 1):
            h = t[i + 1] - t[i]
            k1 = h * np.array(f(x[i], t[i])) / 2.0
            x[i + 1] = x[i] + h * np.array(f(x[i] + k1, t[i] + h / 2.0))

        SimTime = (time.time() - start_time)
        summary = 'Total time of simulation: %s seconds' % (SimTime)
        print(summary)
        df = pd.DataFrame(x, columns=['nG', 'nO', 'N'])
        df['t'] = t
        df.to_csv('results/chemostat.csv', index=False)
        return df

    def rk4(self, f, x0, t):
        """Fourth-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

		USAGE:
			x = rk4(f, x0, t)

		INPUT:
			f     - function of x and t equal to dx/dt.  x may be multivalued,
					in which case it should a list or a NumPy array.  In this
					case f must return a NumPy array with the same dimension
					as x.
			x0    - the initial condition(s).  Specifies the value of x when
					t = t[0].  Can be either a scalar or a list or NumPy array
					if a system of equations is being solved.
			t     - list or NumPy array of t values to compute solution at.
					t[0] is the the initial condition point, and the difference
					h=t[i+1]-t[i] determines the step size h.

		OUTPUT:
			x     - NumPy array containing solution values corresponding to each
					entry in t array.  If a system is being solved, x will be
					an array of arrays.
		"""
        start_time = time.time()
        n = len(t)
        x = np.array([x0] * n)
        for i in range(n - 1):
            h = t[i + 1] - t[i]
            k1 = h * np.array(f(x[i], t[i]))
            k2 = h * np.array(f(x[i] + 0.5 * k1, t[i] + 0.5 * h))
            k3 = h * np.array(f(x[i] + 0.5 * k2, t[i] + 0.5 * h))
            k4 = h * np.array(f(x[i] + k3, t[i + 1]))
            x[i + 1] = x[i] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0
        SimTime = (time.time() - start_time)
        summary = 'Total time of simulation: %s seconds' % (SimTime)
        print(summary)
        df = pd.DataFrame(x, columns=['nG', 'nO', 'N'])
        df['t'] = t
        df.to_csv('results/chemostat.csv', index=False)
        return df


####################################################################################################
if __name__ == '__main__':
    cy = simulate_cyano()
