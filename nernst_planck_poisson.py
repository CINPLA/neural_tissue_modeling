import numpy as np
from dolfin import *
import types
from ufl.algebra import Sum
from ufl.tensors import ListTensor
from time import sleep



PETScOptions.set('pc_type', 'asm') 
PETScOptions.set('sub_pc_type', 'lu')
PETScOptions.set('pc_asm_overlap', '10')

class System:
	"""
	Electro diffusion solver
	"""
	def __init__(self, geometry, temperature=1.):
		self.geometry = geometry
		self.concentration_list = []
		R = 1. # Rayleighs constant
		F = 1. # Faradays constant 
		T = temperature;
		self.eps = 1.
		self.factor = F/(R*T)
		self.a_list = []
		self.L_list = []
		
		self.Phi = TrialFunction(self.geometry.V)
		self.Phi_p = Function(self.geometry.V)
		self.E = nabla_grad(self.Phi_p)
		self.w_Phi = TestFunction(self.geometry.V)
		self.rho = Function(self.geometry.V)
		self.plotE = Function(self.geometry.VV)


	def add_concentration(self, concentration):
		self.concentration_list.append(concentration)


	def electric_field(self):
		
		rho = Function(self.geometry.V)
		for i in range(len(self.concentration_list)):
			rho += self.concentration_list[i].charge*self.concentration_list[i].c_n
		self.rho.assign(rho)
		solve(self.a_Phi == self.L_Phi, self.Phi_p, solver_parameters={"linear_solver": "gmres", "symmetric": True}, \
      			form_compiler_parameters={"optimize": True})



		# E = Function(geometry.VV)

		# for i in range(len(self.concentration_list)):
		# 	E += self.concentration_list[i].electric_field()
		# self.E.assign(E)

		



class Geometry: 
	"""
	The geometry the solver is used on
	"""
	def __init__(self, mesh, space='CG', order=1): 
		domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
		self.meshtype = mesh
		if isinstance(mesh, list): 
			self.dim = len(mesh)
			if len(mesh) == self.dim:
				self.mesh = domain_type[self.dim-1](*mesh)
			else:
				print 'dimension mismatch in set_geometry! mesh does not match dimension'
				print str(self.dim)
				print str(len(mesh))
				sys.exit()

		elif isinstance(mesh, str):
			print 'interpreting mesh input as filename...'
			try:
				self.mesh = Mesh(mesh)
				self.dim = self.mesh.geometry().dim()
			except IOError: 
				print "Could not find the file spesified, exiting...."
				sys.exit(1)

		elif isinstance(mesh,Mesh):
			self.mesh = Mesh(mesh)
			self.dim = self.mesh.geometry().dim()
		else:
			print "input not understood! Exiting..."
			sys.exit(1)

		self.V = FunctionSpace(self.mesh, space, order)
		self.VV = VectorFunctionSpace(self.mesh, space, order)
		#self.V_t = TensorFunctionSpace(self.mesh, space, order)
		#self.V_v = VectorFunctionSpace(self.mesh, space, order)
		#self.vertex_to_dof_map = self.V.dofmap().vertex_to_dof_map(self.V.mesh())


class Time_solver:
	"""
	a class for keeping track of different finite difference 
	time solvers, such as forward/backward euler, Jack Nicholson etc. 
	Uses theta rule for the above! 
	"""
	def __init__(self, system, method, dt=0.01):
		self.dt = dt
		self.system = system
		accepted_methods = ['FE', 'BE', 'CN']
		if method == 'FE': 
			print 'choosing FD time scheme: Forward Euler...'
			self.theta = 0.; 

		elif method == 'BE':
			print 'choosing FD time scheme: Backward Euler...'
			self.theta = 1.; 

		elif method == 'CN':
			print 'choosing FD time scheme: Crank-Nicolson...'
			self.theta = 0.5
		elif isinstance(method, float):
			self.theta = method
		else:
			print 'Unknown method!'
			sys.exit()

	def set_form(self): 
		for i in range(len(self.system.concentration_list)):
			con = self.system.concentration_list[i]
			theta = self.theta
			Dt_v_k_n = con.c-con.c_p
			c_mid = theta*con.c + (1-theta)*con.c_p
			dt = Constant(self.dt)
			# M_grad_v_p = self.M*nabla_grad(self.v_p)
			# inner_prod = project(inner(M_grad_v_p, nabla_grad(self.v_p)),self.V)
			# plot(inner_prod)
			# interactive()
			# self.system.electric_field()
			cc = nabla_grad(c_mid)
			form = (Dt_v_k_n*con.w + dt*inner(con.D*(cc- con.charge*self.system.factor*c_mid*self.system.E), nabla_grad(con.w)))*dx
			(a, L) = system(form)
			self.system.a_list.append(a)
			self.system.L_list.append(L)
		form = (inner(nabla_grad(self.system.Phi), nabla_grad(self.system.w_Phi)) + self.system.rho/self.system.eps*self.system.w_Phi)*dx
 		(self.system.a_Phi, self.system.L_Phi) = system(form)

	def solve_for_time_step(self):
		self.system.electric_field()
		for i in range(len(self.system.concentration_list)):
			con = self.system.concentration_list[i]
			a = self.system.a_list[i]
			L = self.system.L_list[i]
			self.system.electric_field()
			solve(a == L, con.c_n, solver_parameters={"linear_solver": "gmres", "symmetric": True}, \
      			form_compiler_parameters={"optimize": True, "cpp_optimize": True})
			con.c_p.assign(con.c_n)
			con.c_p.vector().apply("insert")

	def solve(self, T, t0=0, plot_interactive=True, save_result=False):
		t = t0
		#p = plot(self.system.E, mode='glyphs', interactive=False)
		while t < T:
			self.solve_for_time_step()
			self.system.plotE.assign(project(self.system.E,self.system.geometry.VV))
			
			plot(self.system.plotE, mode="glyphs")
			#plot(self.system.concentration_list[0].c_p, range_min=0., range_max=1.0)
			#plot(self.system.concentration_list[1].c_p, range_min=0., range_max=1.0)


			t += self.dt

class Concentration:
	"""
	Class keeping track of the concentration object
	"""

	def __init__(self, con_system, charge, initial_condition,diffusion_tensor=Constant(1.0)): 
		self.system = con_system
		self.charge = charge
		self.V = self.system.geometry.V
		self.VV = self.system.geometry.VV
		self.vertex_to_dof_map = vertex_to_dof_map(self.V)
		self.dof_to_vertex_map = dof_to_vertex_map(self.V)
		self.set_initial_condition(initial_condition)
		self.concentration_p = initial_condition
		self.D = diffusion_tensor
	
	def set_initial_condition(self, c0):
		print 'setting initial conditions... ',
		if isinstance(c0, np.ndarray):
			print "initial condition as array, assuming sorted properly... ",

			self.c_p = project(Expression('exp(x[0])'), self.V)
			print self.c_p.vector().array().shape, c0.shape
			
			self.c_p.vector().set_local(c0)
			self.c_p.vector().apply("insert")
		else:
			# self.u = []
			# self.u.append(u0)
			self.c_p = project(c0, self.V)
		self.initial_condition_set = True
		self.E = Function(self.VV)
		self.c = TrialFunction(self.V)
		self.w = TestFunction(self.V)
		self.c_n = Function(self.V)
		self.c_n.assign(self.c_p)

		print 'inital conditions set!'

	# def electric_field(self):
		
	# 	coor = self.V.mesh().coordinates()
		
	# 	E = Function(self.VV)
	# 	v2d_map_VV = vertex_to_dof_map(self.VV)
	# 	d2v_map_VV = dof_to_vertex_map(self.VV)
	# 	v2d_map_V = vertex_to_dof_map(self.V)
	# 	d2v_map_V = dof_to_vertex_map(self.V)
	# 	E_vertices = E.vector().array();
	# 	density = self.c_p.vector().array()
	# 	density = density[v2d_map_V]
	# 	eps = 1.
	# 	num_vertices, dim = np.shape(coor)
	# 	R = VectorFunctionSpace(self.VV.mesh(), 'R', 0)
	# 	c = TestFunction(R)
	# 	for i in range(num_vertices):
	# 		integrand = Function(self.VV)
	# 		integrand_vertices = integrand.vector().array()
	# 		for j in range(num_vertices):
	# 			#print i,j
	# 			if i != j: 
	# 				r = coor[j] - coor[i]
	# 				rn = np.linalg.norm(r)
	# 				integrand_vertices[2*j:2*j+2] = -1./(4*np.pi*eps)*self.charge*density[j]*r/rn**3
	# 		integrand_dof = integrand_vertices[d2v_map_VV]
	# 		integrand.vector()[:] = integrand_dof
	# 		E_local = assemble(dot(integrand,c)*dx).array()
	# 		E_vertices[2*i: 2*i+2] = E_local
	# 	E_dof = E_vertices[d2v_map_VV]
	# 	E.vector()[:] = E_dof
	# 	# E = interpolate(Expression(('10*x[0]', '0')),self.VV)
	# 	return E


def default_f(c, mesh, space, time): 
	"""
	a default proposition for the function f used by monodomain solver
	"""
	return -c

def gaussian_u0_2d():
	"""
	creates a gaussian bell curve centered at the origin
	"""
	c0 = Expression('exp(-(x[0]*x[0] + x[1]*x[1]))')
	return c0

def gaussian_u0_2d_2():
	"""
	creates a gaussian bell curve centered at the origin
	"""
	c0 = Expression('exp(-((x[0]-1)*(x[0]-1) + x[1]*x[1]))')
	return c0


if __name__ == '__main__':
	geometry = Geometry([100,100])
	con_system = System(geometry)
	time_solver = Time_solver(con_system,'BE',dt=0.001)

	g1 = gaussian_u0_2d()
	g2 = gaussian_u0_2d_2()

	D = Expression('((x[0]<0.47) || (x[0]>0.53)) || ((x[1] > 0.47) && (x[1] < 0.53))')
	D3 = Expression('x[0] < 0.5')
	D2 = Constant(1.)
	anions = Concentration(con_system, -5., Expression('x[0]<0.4'), D2)
	kations = Concentration(con_system, 5., Expression('x[0]<0.4'), D3)


	#kations = Concentration(con_system, 1., Expression('0'), D2)

	con_system.add_concentration(anions)
	con_system.add_concentration(kations)

	time_solver.set_form()
	time_solver.solve(1)

	# plot(con_system.E, mode='glyphs', interactive=True)


