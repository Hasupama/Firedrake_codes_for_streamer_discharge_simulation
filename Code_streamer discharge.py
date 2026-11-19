"""streamer simulation"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["PYOP2_CACHE_DIR"] = "/.../.firedrakecache"
os.environ["FIREDRAKE_TSFC_KERNEL_CACHE_DIR"] = "/.../.firedrakecache"

# import all packages from firedrake into cache
from firedrake import *
import time as comp_time
from firedrake.petsc import PETSc
import warnings
warnings.filterwarnings("ignore")

solver_N_parameters ={"pc_factor_mat_solver_type": "mumps","pc_type": "hypre",'pc_hypre_type': 'boomeramg',
                          "snes_rtol": 1.e-4,
                          "snes_atol": 1.e-5,
                     	 "maximum_iterations": 100,
                    		"report": True,
                     	 	"error_on_nonconvergence": True,
                     	 									 }# use vinewtonrsls for direct approach


solver_phi_parameters ={"pc_factor_mat_solver_type": "mumps","pc_type": "hypre",'pc_hypre_type': 'boomeramg',
                          "snes_rtol": 1.e-3,
                          "snes_atol": 1.e-3,
                     	 "maximum_iterations": 100,
                    		"report": True,
                     	 	"error_on_nonconvergence": True,
                     	 									 }	
															 

mesh = Mesh('domain_4mm.msh')

# x[0] , x[1]
x = SpatialCoordinate(mesh)

t = 0.0  #STARTING TIME
tMax =18e-9# FINISHING TIME
Dt = Constant(1e-12)  #TIME STEP
Dt_old=Constant(1e-12)

scalar_lin_space = FunctionSpace(mesh, 'CG', 1)
vector_lin_space = VectorFunctionSpace(mesh, 'CG', 1)#for electric field

Space_N = MixedFunctionSpace([scalar_lin_space,scalar_lin_space])#for charged particles

Space_phi = FunctionSpace(mesh, 'CG', 2)#for phi

dunk_N = TrialFunction(Space_N)
del_unk_N = TestFunction(Space_N)
del_Ne, del_Np = split(del_unk_N)

dphi = TrialFunction(Space_phi)
del_phi = TestFunction(Space_phi)

unk_N = Function(Space_N, name=r'$N$')
unk0_N = Function(Space_N)
unk00_N = Function(Space_N)
unk000_N = Function(Space_N)

Ne, Np = split(unk_N)
Ne0, Np0 = split(unk0_N)
Ne00, Np00 = split(unk00_N)

phi = Function(Space_phi, name=r'$\phi$ in V')
phi0 = Function(Space_phi)

Ne_initial = Function(scalar_lin_space)
Np_initial = Function(scalar_lin_space)
Ne_initial.interpolate( ln(1E13* exp( -(1e-20)* ( (x[0])**2 + (x[1])**2 ) ) ))# omit ln () for direct approach
Np_initial.interpolate(ln (1E13  + 5E18 * exp( -((1/4e-4)**2)* ( (x[0]-0./1000.)**2 + (x[1]-10/1000.)**2 ) ) ))

#ASSIGNE Ne ans Np INITIALLLY
unk_N.sub(0).assign(Ne_initial)
unk0_N.sub(0).assign(Ne_initial)
unk00_N.sub(0).assign(Ne_initial)
unk000_N.sub(0).assign(Ne_initial)

unk_N.sub(1).assign(Np_initial)
unk0_N.sub(1).assign(Np_initial)
unk00_N.sub(1).assign(Np_initial)
unk000_N.sub(1).assign(Np_initial)

e = 1.602176634E-19 #charge of electron in  C
eps_r = 1#dielectric constant air
eps_0 = 8.854187E-12 #in F/m
E = -grad(phi)
sign_e=-1
sign_p=+1
E_mag = ( E[i]*E[i] )**0.5 #in V/m
D = (4.3628E-3*E_mag**0.22) #in m2/s  
eta = 340.75 #in 1/m
alpha =  (1.1944E6 + 4.3666E26 * E_mag**(-3))*exp(-2.73E7/E_mag)
mu_e=2.3987*E_mag**-0.26
We1=mu_e*E_mag #in m/s
We=sign_e*(E/E_mag)*We1
We_abs = ( We[i]*We[i] )**0.5 # We[0]*We[0] + We[1]*We[1]
Si = exp(Ne)*alpha*We_abs   
tr = Dt / Dt_old
trp1 = 1.0 + tr
tr2p1 = 1.0 + 2.0 * tr

F_Ne =  x[0]*( exp(Ne)*((Ne * tr2p1 - trp1**2.0 * Ne0 + tr**2.0 * Ne00) / trp1/Dt)*del_Ne +  (-exp(Ne)*We[i] + (D*exp(Ne)).dx(i))*del_Ne.dx(i) - Si*del_Ne + exp(Ne)*eta*We_abs*del_Ne  )*dx
F_Np =  x[0]*( exp(Np)*((Np * tr2p1 - trp1**2.0 * Np0 + tr**2.0 * Np00) / trp1/Dt)*del_Np  - Si*del_Np  + exp(Ne)*eta*We_abs*del_Np )*dx
F_phi =x[0]*( phi.dx(i)*del_phi.dx(i) - (e/(eps_r*eps_0))*(exp(Np)-exp(Ne))*del_phi )*dx

bc_list_N, bc_list_phi = [],[]

'''
"AllRestBoundaries", 10
"TopAnode", 11
"BottomCathode", 12
'''

bc_list_N.append( DirichletBC(Space_N.sub(0), ln(1e13), 11) )
bc_list_phi.append( DirichletBC(Space_phi, 0.0, 12) )
bc_list_phi.append( DirichletBC(Space_phi, 18.75E3, 11) )

Form_N = F_Ne  +F_Np
Form_phi = F_phi 

problem_N = NonlinearVariationalProblem(Form_N, unk_N, bcs=bc_list_N)
solver_N = NonlinearVariationalSolver(problem_N, solver_parameters=solver_N_parameters)

problem_phi = NonlinearVariationalProblem(Form_phi, phi, bcs=bc_list_ph)
solver_phi = NonlinearVariationalSolver(problem_phi, solver_parameters=solver_N_parameters)

sol_file = File("...pvd")

E_sol = Function(vector_lin_space,name=r'E field in V/m')
Ne_sol = Function(scalar_lin_space,name=r'Electrons exp')
Np_sol = Function(scalar_lin_space,name=r'positive ions exp')
Ne_sol1 = Function(scalar_lin_space,name=r'Electrons ')
Np_sol1 = Function(scalar_lin_space,name=r'positive ions ')

start_time = comp_time.time()

t+=Dt(0)
elapsed = int(comp_time.time() - start_time)
e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
PETSc.Sys.Print('time: %.3f ns solved in %.0f h %.0f min %.0f s ' % (t*1E9, e_h, e_m, e_s)  )

Dt_max=1e-12
Dt_min=1e-15
ii=0
ij=0
total_charge=[]
while t <= tMax:
	Dt = Constant(1e-12)

	solver_phi.solve()
	E_sol.interpolate( E )
	if ii==12:
		PETSc.Sys.Print('sad sad')
		break

	while ii<13:
		try:

			if Dt(0)>Dt_max:
				Dt.assign(Dt_max)
			if Dt(0)<Dt_min:
				Dt.assign(Dt_min)
	
			solver_N.solve()#use bounds for direct approach
			
			Dt_old.assign(Dt)
						
			Ne_sol, Np_sol = unk_N.split()
			Ne_sol1.interpolate(exp(Ne_sol))
			Np_sol1.interpolate(exp(Np_sol))
				
			t+=Dt(0)
			
			if t>ij:
				sol_file.write(E_sol,Ne_sol1,Np_sol1, time=t)
				ij=ij+1e-9
				rho = e * (exp(Np) - exp(Ne))  
				Q = assemble(2*pi * x[0] * rho * dx)
				total_charge.append(Q)
				PETSc.Sys.Print(total_charge)

			elapsed = int(comp_time.time() - start_time)
			e_h, e_m, e_s = int(elapsed/3600), int(elapsed % 3600 / 60), int( ( elapsed % 3600 ) % 60 )
			PETSc.Sys.Print('time: %.3f ns solved in %.0f h %.0f min %.0f s ' % (t*1E9, e_h, e_m, e_s)  )
	
			unk000_N.assign(unk00_N)
			unk00_N.assign(unk0_N)
			unk0_N.assign(unk_N)
			
			break
		except:
			PETSc.Sys.Print('not converged at all ',float(Dt), ii,t)
			unk_N.assign(unk0_N)
			unk0_N.assign(unk00_N)
			unk00_N.assign(unk000_N)	
			Dt.assign(Dt/2)	
			ii+=1
			PETSc.Sys.Print(Dt(0))
			pass
		PETSc.Sys.Print(total_charge)
		if ii==12:
			PETSc.Sys.Print('not converged')
			break
