from gurobipy import *
import numpy as np
import parameters

p = parameters.cyano_pars()

### Create variables:
# flux vector and species concentration
rxns = ['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'vd']
gammas = ['gR', 'gP1', 'gP2', 'gP3', 'gP4', 'gP5', 'gP6', 'gP']
RP = ['R', 'PSU', 'Tc', 'CB', 'Mc', 'Tn', 'Mq', 'Pq']
Misc = ['Q',]

x_name = rxns + gammas + RP + Misc
idx = np.arange (0, len (x_name))
obj = np.zeros (len (x_name))

species = ['Ci', 'e', 'C3', 'Q', 'aa', 'Ni', 'R', 'PSU', 'Tc', 'CB', 'Mc', 'Tn', 'Mq', 'Pq']

stoic_mat = np.array([[0.0, 1.0, -3.0, 0.0, 0.0, 0.0, 0.0],
					  [8.0, -1.0, -10.0, -35.0, -1.0, 0.0, 0.0],
					  [0.0, 0.0, 1.0, -2.0, 0.0, -1.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
					  [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, p.nPSU],
					  [0.0, 0.0, 0.0, -2.0, 1.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
					  ])

alpha_P = np.array([p.nR, p.nPSU, p.nTc, p.nCB, p.nMc, p.nTn, p.nMq, p.nPq])

########################################## MODEL ###################################################
#### Create model for a given growth rate and external carbon
def cyano_model(mu, Cx, Nx, I, newpars):
	p = parameters.cyano_pars(newpars)

	### Create the minimal model in gurbobi
	model = Model('cyano_model')
	# lower bounds, all reactions are irreversible
	lb = np.zeros(len(x_name))
	# upper bounds
	ub = GRB.INFINITY * np.ones(len(x_name))
	# add all variables to the model
	x = model.addVars(idx, lb=lb, ub=ub, obj=obj, vtype=GRB.CONTINUOUS, name='x')
	
	### Mass conservation
	mass_cons = np.zeros ([len (species), len (x_name)])
	
	# add stoichiometic coefficients
	for i in range(len(species)):
		for j in range(len(rxns)):
			mass_cons[species.index (species[i]), x_name.index (rxns[j])] = stoic_mat[i, j]
		
	# add translation coefficients
	mass_cons[species.index ('e'), len (rxns):len (rxns) + len (gammas)] = -3.0 * alpha_P
	mass_cons[species.index ('aa'), len (rxns):len (rxns) + len (gammas)] = -alpha_P
	np.fill_diagonal (mass_cons[-len (RP):, len (rxns):], 1.0)
	
	# dilution of R, P, Q by growth
	np.fill_diagonal (mass_cons[-len (RP):, -(len (RP) + len (Misc)):], -mu)
	mass_cons[species.index ('Q'), x_name.index ('Q')] = -mu

	# add constraint for mass conservation
	B_eq = np.zeros(len(species))
	model.addConstrs((quicksum(mass_cons[i, j] * x[j] for j in range(len(x_name))) ==
	                  B_eq[i] for i in range(len(B_eq))), 'mass_conservation')

	### Constant density
	A_eq2 = np.hstack((np.zeros(len(rxns)+len(gammas)), alpha_P, np.zeros(len(Misc))))
	model.addConstr(quicksum(A_eq2[i] * x[i] for i in range(len(A_eq2))) == 1.4e10, 'fixed_density')

	### Constant concentration of Q protein
	A_eq3 = np.zeros (len (x_name))
	A_eq3[x_name.index ('Q')] = p.nQ
	B_eq3 = np.hstack ((np.zeros (len (rxns) + len (gammas)), alpha_P, [p.nQ,]))
	model.addConstr (quicksum (A_eq3[i] * x[i] for i in range (len (A_eq3))) == 0.5 * quicksum (
		B_eq3[i] * x[i] for i in range (len (A_eq3))), 'fixed_Q_conc_ratio')

	### Constant concentration of dummy P protein
	A_eq4 = np.zeros(len(x_name))
	A_eq4[x_name.index('Pq')] = p.nPq
	B_eq4 = np.hstack((np.zeros(len(rxns)+len(gammas)), alpha_P, np.zeros(len(Misc))))
	model.addConstr(quicksum(A_eq4[i] * x[i] for i in range(len(A_eq4))) == 0.2 *
	                quicksum(B_eq4[i] * x[i] for i in range(len(A_eq4))), 'fixed_P_conc_ratio')
	
	### Constant concentration of Ribosomes
	A_eq4 = np.zeros (len (x_name))
	A_eq4[x_name.index ('R')] = p.nR
	B_eq4 = np.hstack ((np.zeros (len (rxns) + len (gammas)), alpha_P, np.zeros (len (Misc))))
	model.addConstr (quicksum (A_eq4[i] * x[i] for i in range (len (A_eq4))) >= 0.05 *
	                 quicksum (B_eq4[i] * x[i] for i in range (len (A_eq4))), 'fixed_R_conc')

	### Constant concentration of Rubisco
	A_eq4 = np.zeros (len (x_name))
	A_eq4[x_name.index ('CB')] = p.nCB
	A_eq4[x_name.index ('Mc')] = p.nMc
	A_eq4[x_name.index ('Mq')] = p.nMq
	B_eq4 = np.hstack ((np.zeros (len (rxns) + len (gammas)), alpha_P, np.zeros (len (Misc))))
	model.addConstr (quicksum (A_eq4[i] * x[i] for i in range (len (A_eq4))) >= 0.1 *
	                 quicksum (B_eq4[i] * x[i] for i in range (len (A_eq4))), 'fixed_CB_conc')

	### Enzyme capacity R, P1, P2, P3, P4, P5, P6, P7, P
	enzyme_capacity = np.zeros ([len (RP), len (x_name)])
	
	enzyme_capacity[RP.index ('R'), len (rxns):len (rxns) + len (gammas)] = alpha_P / p.g_max
	np.fill_diagonal (enzyme_capacity[RP.index ('PSU'):, :], 1.0)
	
	v_PSU = (p.k1 * p.sigma * I) / (p.sigma * I + p.k1 + p.kd * p.sigma * I)
	v_Tc = (p.k2 * Cx) / (p.K_C + Cx)
	v_Tn = (p.k5 * Nx) / (p.K_N + Nx)
	v_d = (p.kd * np.square (p.sigma * I)) / (p.sigma * I + p.k1 + p.kd * p.sigma * I)
	
	rates = [1.0, v_PSU - v_d, v_Tc, p.k3, p.k4, v_Tn, p.k6, p.kP]
	for i in range (len (enzyme_capacity)):
		enzyme_capacity[i, len (rxns) + len (gammas) + i] = - rates[i]
	
	B_eq5 = np.zeros (len (enzyme_capacity))
	model.addConstrs ((quicksum (enzyme_capacity[i, j] * x[j] for j in range (len (x_name))) <=
	                   B_eq5[i] for i in range (len (B_eq5))), 'enzyme_capacity')

	### optimize
	model.setParam('OutputFlag', False)
	model.optimize()
	return model

####################################### BINARY SEARCH ##############################################
def binary_search(mu, Cx=10.0, Nx=10.0, I=200 * 8 * 60, newpars={}, mu_lb=0.0, mu_ub=0.0,
                  mu_prev=0.0):
	# if prev_model is set to 1, calculate mu for the previous model for which a solution exist
	prev_model = 0.0
	m = cyano_model(mu, Cx, Nx, I, newpars)
	# if optimization terminated successfully
	if m.status == GRB.status.OPTIMAL:
		mu_lb = mu
		if mu_ub == 0.0:
			mu_new = 2.0*mu_lb
		else:
			mu_new = (mu_lb+mu_ub)/2.0
			mu_prev = mu
	# and if not
	else:
		prev_model = 1.0
		mu_ub = mu
		if mu_lb == 0.0:
			mu_new = mu_ub/2.0
		else:
			mu_new = (mu_ub+mu_lb)/2.0

	if 100.0*(np.abs(mu_new-mu)/mu) > 0.01:
		return binary_search(mu_new, Cx, Nx, I, newpars, mu_lb, mu_ub, mu_prev)
	else:
		if prev_model == 0.0:
			return [mu, m.getAttr('x')]
		else:
			m = cyano_model(mu_prev, Cx, Nx, I, newpars)
			return [mu_prev, m.getAttr('x')]