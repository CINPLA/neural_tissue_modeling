import numpy as np

##### Astrocyte model

### PHYSICAL CONSTANTS ################## Units: mol, m, V, s, C, ohm, A
R = 8.315 # Gas constant  J*mol^-1*K^-1
F = 96500 # Faraday's constant (C/mol)
T = 298 # Temperature, K
psifac = R*T/F #Standard potential (V)

#### GEOMETRICAL PARAMETERS
# Notation here: o means outside (ECS), i means inside (intracellular)
l = 3.00e-4 #Length of cell (m)
A_i = 0.4 # Chen&Nicholson (Tissue volume fraction being astrocytes)
A_o = 0.2 #Chen&Nicholson (Tissue volume fraction being ECS)
O_m = A_i/5.00e-8 # (Astrocytic membrane area per tissue volume)

#### INPUT PARAMETERS
tstart = 100 # Input starts
tstop = 400 # and ends (s)
#tstart = 10 tstop = 45 # Used with compare_models.m
xstart = 0 
xstop = l*0.1 # Region for input (m)

Imax = 5.5e-7 # Input amplitude K/Na-(exchange) (mol/(s m^2))
Kdec = 2.9e-8 # Output factor (m/s)

### DISCRETISATION
simt = 600 # Simulate for 600 s
xres = 1*100
x = np.linspace(0,l,xres) # 100 x-points
tres = 1*100
t = np.linspace(0,simt,tres) # 100 t-points


### MEMBRANE PARAMETERS
g_Na = 1 # Na+ conductance (S*m^-2)
g_Cl=0.5 # Cl- concuctance (S*m^-2)
g_K0=16.96 # K+ conductance (S*m^-2)
C_m = 1.00e-2 # Membrane capacitance (Farad/m^2)

#### ELECTRODIFFUSION PARAMETERS (from Grodzinski-book)
D_K = 1.96e-9 #Diffusion coefficients (m^2/s)
D_Cl=2.03e-9
D_Na=1.33e-9
lambda_o = 1.6 # ECS tortuousity (1), C&N 
lambda_i=3.2 # Intracellular tortuousity (1), C&N


### STATIC CHARGES GIVING THE CORRECT CONCENTRATION/V_M-relation
c_Ko0=3.0+0.082
c_Ki0=100.0 - 0.041
c_Nao0=145.0 - 0.338
c_Nai0=15.0 + 0.189
c_Clo0=134 - 0.29
c_Cli0=5 + 0.145
v_m0 = -83.6e-3 # Membrane pot. (V)
u0 = [c_Ko0, c_Ki0, c_Nao0, c_Nai0, c_Clo0, c_Cli0, v_m0, v_m0]
X_iZ_i = -(F*A_i/O_m)*(c_Ki0 + c_Nai0 - c_Cli0) + C_m*v_m0 #Intracellular static charge (C/m^2)
X_oZ_o = -(F*A_o/O_m)*(c_Ko0 + c_Nao0 - c_Clo0) - C_m*v_m0 #ECS static charge (C/m^2)

