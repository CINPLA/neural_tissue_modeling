import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

size_x = 10
size_y = 10

grid_map = np.zeros([size_x, size_y]) + 1


grid_map[:,0] = 0
grid_map[:,-1] = 0
grid_map[0,:] = 0
grid_map[-1,:] = 0

grid_map[:,3:4] = 3
grid_map[:,6:7] = 3

grid_map[0,3:4] = 6
grid_map[-1,3:4] = 6
grid_map[0,6:7] = 6
grid_map[-1,6:7] = 6

fig = plt.figure()


im = plt.imshow(grid_map, extent=[0,1,0,1], vmin=0, vmax=1)
im.set_interpolation('nearest')
plt.colorbar()
#fig.set_title('Default')


z = np.zeros([size_x, size_y])
z[3:4, 3:4] = 5

dt = 0.001
dx = 0.1
T = 2.
t0 = 0.
N = int((T-t0)/dt)
t = np.linspace(t0, T, N)

D = 1

for i in range(N):
	z[0,1:-1] = z[2,1:-1]
	z[1:-1,0] = z[1:-1,2]
	z[-1,1:-1] = z[-3,1:-1]
	z[1:-1,-1] = z[1:-1,-3]
	z_new = np.zeros([size_x, size_y])
			
	for j in range(1,size_x-1):
		for k in range(1,size_y-1):
			z_new[j,k] = z[j,k] + D*(dt/dx**2)*((z[j,k+1] - 2*z[j,k] + z[j,k-1]) \
			 + (z[j+1,k] - 2*z[j,k] + z[j-1,k]))

	im.set_data(z)
	z = np.copy(z_new)
	

	print np.sum(z[1:-1, 1:-1])

plt.show()