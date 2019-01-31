
##################################### MODEL PARAMETERS #############################################
class cyano_pars():
	defaultpara={
		'kd': 0.56, 'sigma':0.2,
		'K_C':15.0, 'K_N':10.0, 'K_N1':10.0, 'K_N2':10.0, # microM
		'k1':250.0, 'k2':20.0, 'k2_2':200.0, 'k3': 1.0, 'k4': 10.0, 'k5': 50.0, 'k6': 100.0,
		'k7': 100.0, 'kQ': 100.0, 'kP': 100.0, 'k5_1':30.0, 'k5_2':50.0,
		'nR': 7358.0, 'nPSU':95451.0, 'nTc':1681.0, 'nCB':2000.0, 'nMc':20000.0, 'nTn':10000.0,
		'nTn1':20000.0, 'nTn2':10000.0,
		'nMq': 1000.0, 'nP7': 1000.0, 'nPq':1000.0, 'nQ':100.0, 'nC2':1.0, 'g_max': 22.0,
    }
	
	def __init__(self, pars={}):
		mypars = pars.copy()
		for k in cyano_pars.defaultpara.keys():
			mypars.setdefault(k,cyano_pars.defaultpara[k])
		for k,v in mypars.items():
			setattr(self, k, mypars[k])