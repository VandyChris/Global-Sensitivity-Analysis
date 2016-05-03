import numpy as np
##------First Order index, double loop MCS-------------------------------------------------------
def FirstOrder(Input, func, n, *var):
	'''
	this function conducts the sensitivity analysis using the first-order index
	here Input is a N*D array, where N is the number of sampling points, 
	and D is the dimension of the Input
	func is a function to generate Output, which can transfer N*D Input
	to N*1 Outputs, such as a GP model
	n is the number of Input samples you'd like to use. For example, N can be
	10000 but only 5000 are used in the analysis if you set n=5000
	*var: the argument of func besides Input
	'''
	(N, D) = np.shape(Input)

	if N<n:
		print "Error: not enough sampling pints"
		return

	Input = Input[:n, :]
	var_Output = np.var(func(Input, *var), ddof = 1)

	index = np.zeros(D)
	for i in range(D):
		Output_xi = np.zeros((n, n))
		for j in range(n):
			Input_temp = Input.copy()
			Input_temp[:, i] = np.ones(n) * Input[j, i]
			Output_xi[:, j] = func(Input_temp, *var)
		index[i] = np.var(np.mean(Output_xi, axis = 0), ddof = 1)/var_Output
		print i+1

	return index

#-------------------Total Effects index, double-loop MCS-----------------------------------------
def TotalEffect(Input, func, n, *var):
	'''
	this function conducts the sensitivity analysis using the total effects index
	here Input is a N*D array, where N is the number of sampling points, 
	and D is the dimension of the Input
	func is a function to generate Output, which can transfer N*D Input
	to N*1 Outputs, such as a GP model
	n is the number of Input samples you'd like to use. For example, N can be
	10000 but only 5000 are used in the analysis if you set n=5000
	*var: the argument of func besides Input
	'''
	(N, D) = np.shape(Input)

	if N<n:
		print "Error: not enough sampling pints"
		return

	Input = Input[:n, :]
	var_Output = np.var(func(Input, *var), ddof = 1)

	index = np.zeros(D)
	for i in range(D):
		Output_xi = np.zeros((n, n))
		for j in range(n):
			Input_temp = np.tile(Input[j, :], (n, 1))
			Input_temp[:, i] = Input[:, i].copy()
			Output_xi[:, j] = func(Input_temp, *var)

		index[i] = np.mean(np.var(Output_xi, axis = 0), ddof = 1)/var_Output
		print i+1

	return index
 
#-------------------------First-order, MGSA -----------------------------------
def MGSA_FirstOrder(Input, Output, ndomain):
    '''
    input: nsample * nd matrix, where nsample is the number of sample, and nd
    is the input dimension
    output: nsample * 1 array
    ndomain: number of sub-domain the to divide a single input
    This algorithm is proposed by me. Please cite the following paper if you use this code.
    Li, Chenzhao, and Sankaran Mahadevan. "An efficient modularized sample-based method to estimate the first-order Sobol’index." Reliability Engineering & System Safety (2016).
    '''
    (nsample, nd) = np.shape(Input);
    
    # convert the input samples into cdf domains
    U = np.linspace(0.0, 1.0, num = ndomain + 1)
    cdf_input = np.zeros((nsample, nd))
    cdf_values = np.linspace(1.0/nsample, 1.0, nsample)
    
    j = 1
    for i in range(nd):
        IX = np.argsort(Input[:, i])
        IX2 = np.argsort(IX)
        cdf_input[:, i] = cdf_values[IX2]
        
    # compute the first-order indices
    VY = np.var(Output, ddof = 1)
    VarY_local = np.zeros((ndomain, nd))
    for i in range(nd):
        cdf_input_i = cdf_input[:, i]
        output_i = Output
        U_i = U
        for j in range(ndomain):
            sub = cdf_input_i < U_i[j+1]
            VarY_local[j, i] = np.var(output_i[sub], ddof = 1)
            inverse_sub = ~sub
            cdf_input_i = cdf_input_i[inverse_sub]
            output_i = output_i[inverse_sub]
            
    index = 1.0 - np.mean(VarY_local, axis = 0)/VY
    
    return index
    
    
    
def Sobol07_FirstOrder(x, z, f, *var):
    '''this function computes the first-order index by Sobol's single-loop method
    the function is V_i = mean((f(x)-c)*(f(xi, z_-i)-f(z))), see my MSA paper'''
    
    '''
    n: number of sample
    k: number of dimension
    x: n*k matrix
    z: n*k matrix
    f: function handle, so y_x = f(x), y_z = f(z)
    '''
    
    y_x = f(x, *var)
    y_z = f(z, *var)
    k = np.shape(x)[1]
    y_all = np.concatenate((y_x, y_z))
    Vy = np.var(y_all)
    c = np.mean(y_all)
    sen_vector = np.zeros(k)
    
    for i in range(k):
        xz = z.copy()
        xz[:, i] = x[:, i]
        y_xz = f(xz, *var)
        V_i = np.mean((y_x - c)*(y_xz - y_z))
        sen_vector[i] = V_i/Vy
        print i+1
        
    return sen_vector
    
def Sobo07_TotalEffects(x, z, f, *var):
    '''this function computes the first-order index by Sobol's single-loop method
    the function is V_i = mean((f(x)-f(x_-i, zi)^2)/2, see my MSA paper'''
    
    y_x = f(x, *var)
    y_z = f(z, *var)
    k = np.shape(x)[1]
    y_all = np.concatenate((y_x,y_z))
    Vy = np.var(y_all)
    sen_vector = np.zeros(k)
    
    for i in range(k):
        zx = x.copy()
        zx[:, i] = z[:, i]
        y_zx = f(zx, *var)
        VT_i = np.mean((y_x - y_zx)**2)/2.0
        sen_vector[i] = VT_i/Vy
        print i+1
        
    return sen_vector