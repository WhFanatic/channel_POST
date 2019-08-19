#!/work1/cuigx2_work/whn/anaconda_install/anaconda3/bin/python
import numpy as np
from numpy.fft import ifft as fft # in fortran program, ifft and fft are inversed
from numpy.fft import fft as ifft # in fortran program, ifft and fft are inversed
from numpy.polynomial.chebyshev import chebfit, chebder, chebval
from struct import pack,unpack


# class Data:
# 	def __init__(self, workpath="workpath.txt"):
# 		with open(workpath) as fp:
# 			for line in fp.readlines():
# 				line = line.strip()
# 				if "SOURCE_PATH" == line[:11]: self.data_path = line.split()[-1].strip()
# 				if "TARGET_PATH" == line[:11]: self.postdata_path = line.split()[-1].strip()

# 		self.getPara()

# 	def getPara():
# 		pass


# a = Data()

# exit()


## user input
data_path = '../data/data_F180/'#_ofw_b5_arti/'
postdata_path = data_path + 'postdata/'



## computation parameters
H = 1.0 # half channel height is 1.0 by default

file_name = 'XINDAT'
with open(data_path + file_name) as fp:
	lines = fp.readlines()

info = lines[1].split()
nu = 1.0/float(info[0])
dt = float(info[3])
Nx, Ny, Nz = int(info[4]), int(info[5]), int(info[6])
Lx, Ly, Lz = 2*np.pi/float(info[1]), 2.0*H, 2*np.pi/float(info[2])

info = lines[3].split()
nsp = int(info[1])
nsprb = int(info[4])

info = lines[17].split()
j_probes = [int(n) for n in info]
# j_probes = [int(n) for n in info[: 1+info[1:].index('0')]]

## get data set
file_name = 'datainfo'
with open(data_path + 'fielddata/' + file_name) as fp:
	lines = fp.readlines()

tstep_start, tstep_end = [int(n) for n in lines[0].split()]
tsteps = range(tstep_start, tstep_end+1, nsp)
print 'steps from %i to %i every %i, total %i samples.\n'%(tstep_start, tstep_end, nsp, len(tsteps))

## non-dimensional parameters
file_name = 'XUTP.DAT'
with open(data_path + 'runtimedata/' + file_name) as fp:
	Re_tau = []
	for line in fp:
		info = line.split()
		if float(info[0]) >= dt * min(tsteps):
			Re_tau.append( float(info[1])**0.5 * H / nu )
		elif float(info[0]) > dt * max(tsteps):
			break
Re_tau = np.mean(Re_tau)

try:
	file_name = 'XFR.DAT'
	with open(data_path + 'runtimedata/' + file_name) as fp:
		Re = []
		for line in fp:
			info = line.split()
			if float(info[0]) >= dt * min(tsteps):
				Re.append( 0.5 * float(info[1]) / nu )
			elif float(info[0]) > dt * max(tsteps):
				break
	Re = np.mean(Re) # Re is not valid in off-wall cases
except IOError:
	Re = 1.0 / nu

## get coordinates
try:
	file_name = 'INITFIELD.DAT'
	with open(data_path + 'runtimedata/' + file_name) as fp:
		y_mesh = []
		for line in fp:
			y_mesh.append(float(line.split()[0]))
	y_mesh = np.array(y_mesh)
except IOError:
	y_mesh = - np.cos( np.arange(Ny+1) * np.pi / Ny )
x_mesh = float(Lx) * np.arange(Nx) / Nx
z_mesh = float(Lz) * np.arange(Nz) / Nz
k_x = np.array( range(Nx/2) + range(-Nx/2,0) ) * (2*np.pi/Lx)
k_z = np.array( range(Nz/2) + range(-Nz/2,0) ) * (2*np.pi/Lz)
lambda_x = 2*np.pi / (k_x + 1e-20)
lambda_z = 2*np.pi / (k_z + 1e-20)

## non-dimensionalize
u_tau = Re_tau * nu / H
tau_w = u_tau**2
delta_nu = H / Re_tau
t_nu = delta_nu / u_tau
y_plus = (y_mesh + 1.0) / delta_nu
lambda_x_plus = lambda_x / delta_nu
lambda_z_plus = lambda_z / delta_nu


print 'Re_tau = %.2f'%Re_tau
print 'Re = %.2f'%Re
print 'dt =', dt
print 'Nx Ny Nz =', Nx, Ny, Nz
print 'Lx Ly Lz = %.1f pi , %.1f , %.1f pi\n'%(Lx/np.pi, Ly, Lz/np.pi)
print 'u_tau = %.3f, delta_nu = %.3f, t_nu = %.3f, tau_w = %.3f'%(u_tau, delta_nu, t_nu, tau_w)
print 'probe y+ =', [int(y_plus[max(min(n,Ny),0)]*10)/10.0 for n in j_probes]
print ''




## centeral difference coefficients for 1st & 2nd order derivative in y direction (upwind difference at wall)
diff1 = np.zeros([Ny+1, 3])
diff2 = np.zeros([Ny+1, 3])
for j in range(Ny+1):
	if 0 < j < Ny:
		dy1 = y_mesh[j-1] - y_mesh[j]
		dy2 = y_mesh[j+1] - y_mesh[j]
		diff1[j] = np.array([dy2/dy1/(dy2-dy1), -(1/dy1+1/dy2), -dy1/dy2/(dy2-dy1)])
		diff2[j] = 2/(dy1*dy2) * np.array([-dy2/(dy2-dy1), 1, dy1/(dy2-dy1)])
	elif j == 0:
		dy1 = y_mesh[1] - y_mesh[0]
		dy2 = y_mesh[2] - y_mesh[0]
		diff1[j] = np.array([-(1/dy1+1/dy2), dy2/dy1/(dy2-dy1), -dy1/dy2/(dy2-dy1)])
		diff2[j] = 2/(dy1*dy2) * np.array([1, -dy2/(dy2-dy1), dy1/(dy2-dy1)])
	elif j == Ny:
		dy1 = y_mesh[-2] - y_mesh[-1]
		dy2 = y_mesh[-3] - y_mesh[-1]
		diff1[j] = np.array([-dy1/dy2/(dy2-dy1), dy2/dy1/(dy2-dy1), -(1/dy1+1/dy2)])
		diff2[j] = 2/(dy1*dy2) * np.array([dy1/(dy2-dy1), -dy2/(dy2-dy1), 1])




try:
	U_mean, V_mean, W_mean, P_mean, dUdy = list( np.loadtxt(open(postdata_path + 'means.txt')).T [1:] )
except IOError:
	U_mean, V_mean, W_mean, P_mean = [ np.zeros(Ny+1) for n in range(4) ]
	for tstep in tsteps:
		print 'Reading mean values of step %i ...'%tstep
		with open( data_path+'fielddata/U/U'+str(tstep).zfill(8)+'.BIN', 'rb' ) as fp:
			for j in range(Ny+1):
				fp.seek((j+1) * Nz*Nx/2 * 8)
				U_mean[j] += unpack('f', fp.read(4))
		with open( data_path+'fielddata/V/V'+str(tstep).zfill(8)+'.BIN', 'rb' ) as fp:
			for j in range(Ny+1):
				fp.seek((j+1) * Nz*Nx/2 * 8)
				V_mean[j] += unpack('f', fp.read(4))
		with open( data_path+'fielddata/W/W'+str(tstep).zfill(8)+'.BIN', 'rb' ) as fp:
			for j in range(Ny+1):
				fp.seek((j+1) * Nz*Nx/2 * 8)
				W_mean[j] += unpack('f', fp.read(4))
		with open( data_path+'fielddata/P/P'+str(tstep).zfill(8)+'.BIN', 'rb' ) as fp:
			for j in range(Ny+1):
				fp.seek((j+1) * Nz*Nx/2 * 8)
				P_mean[j] += unpack('f', fp.read(4))
	U_mean = (U_mean + U_mean[::-1]) / 2 / len(tsteps) # average among time steps and between upper & lower half of the channel
	V_mean = (V_mean - V_mean[::-1]) / 2 / len(tsteps)
	W_mean = (W_mean + W_mean[::-1]) / 2 / len(tsteps)
	P_mean = (P_mean + P_mean[::-1]) / 2 / len(tsteps)
	dUdy = diff1[:,0] * U_mean[[0]+range(0,Ny-1)+[Ny-2]] + diff1[:,1] * U_mean[[1]+range(1,Ny)+[Ny-1]] + diff1[:,2] * U_mean[[2]+range(2,Ny+1)+[Ny]]
	np.savetxt(postdata_path + 'means.txt', np.array([y_mesh,U_mean,V_mean,W_mean,P_mean,dUdy]).T)

dPdy = diff1[:,0] * P_mean[[0]+range(0,Ny-1)+[Ny-2]] + diff1[:,1] * P_mean[[1]+range(1,Ny)+[Ny-1]] + diff1[:,2] * P_mean[[2]+range(2,Ny+1)+[Ny]]
dUdyy = diff2[:,0] * U_mean[[0]+range(0,Ny-1)+[Ny-2]] + diff2[:,1] * U_mean[[1]+range(1,Ny)+[Ny-1]] + diff2[:,2] * U_mean[[2]+range(2,Ny+1)+[Ny]]
dPdyy = diff2[:,0] * P_mean[[0]+range(0,Ny-1)+[Ny-2]] + diff2[:,1] * P_mean[[1]+range(1,Ny)+[Ny-1]] + diff2[:,2] * P_mean[[2]+range(2,Ny+1)+[Ny]]
means = {'U':U_mean, 'V':V_mean, 'W':W_mean, 'P':P_mean, 'dUdy':dUdy}
print 'Check V, W mean: %f, %f \n' %( np.sum(V_mean**2), np.sum(W_mean**2) )



# space time parameters
C_ave_plus = 15 # average convection velocity ~= 15 * u_tau
Nt = (Lx / (C_ave_plus*u_tau)) / (nsprb * dt) # number of time steps in a sample
for n in np.sort(np.ravel([[(2**i * 3**j) for i in range(1,20)] for j in range(20)])):
	if n >= Nt: # find the next integet that is greater than or equal to Nt and is power of 2 and 3
		Nt = n
		break
Lt = Nt * (nsprb * dt)
t_mesh = float(Lt) * np.arange(Nt) / Nt
k_t = np.array( range(Nt/2) + range(-Nt/2,0) ) * (2*np.pi/Lt)

delta_x = np.sort(k_x) / (2*np.pi/Lx) * (Lx/Nx)
delta_z = np.sort(k_z) / (2*np.pi/Lz) * (Lz/Nz)
delta_t = np.sort(k_t) / (2*np.pi/Lt) * (Lt/Nt)
delta_x_plus = delta_x / delta_nu
delta_z_plus = delta_z / delta_nu
delta_t_plus = delta_t / t_nu
k_x_plus = k_x * delta_nu
k_z_plus = k_z * delta_nu
k_t_plus = k_t * t_nu

print 'Lt = %.2f, Lt^+ = %.2f, Nt = %i\n'%(Lt, Lt/t_nu, Nt)




## math functions
def spec(q):
	return fft(fft(q, axis=-2), axis=-1)
def phys(q):
	return ifft(ifft(q, axis=-2), axis=-1).real
def normalize(q):
	return q / ( np.max(abs(q)) + 1e-20 )
def convol(u, v=0): # pseudo spectral convolution using 3/2 method
	u2 = np.zeros([Nz*3/2, Nx*3/2], dtype=complex)
	uv = np.zeros([Nz, Nx], dtype=complex)
	u2[ :Nz/2,  :Nx/2 ] = u[ :Nz/2,  :Nx/2 ]
	u2[ :Nz/2, 1-Nx/2:] = u[ :Nz/2, 1-Nx/2:]
	u2[1-Nz/2:, :Nx/2 ] = u[1-Nz/2:, :Nx/2 ]
	u2[1-Nz/2:,1-Nx/2:] = u[1-Nz/2:,1-Nx/2:]

	if type(v) == int and v == 0:
		u2v2 = spec(phys(u2)**2)
	else:
		v2 = np.zeros([Nz*3/2, Nx*3/2], dtype=complex)
		v2[ :Nz/2,  :Nx/2 ] = v[ :Nz/2,  :Nx/2 ]
		v2[ :Nz/2, 1-Nx/2:] = v[ :Nz/2, 1-Nx/2:]
		v2[1-Nz/2:, :Nx/2 ] = v[1-Nz/2:, :Nx/2 ]
		v2[1-Nz/2:,1-Nx/2:] = v[1-Nz/2:,1-Nx/2:]
		u2v2 = spec(phys(u2) * phys(v2))

	uv[ :Nz/2 ,  :Nx/2 ] = u2v2[ :Nz/2 ,  :Nx/2 ]
	uv[ :Nz/2 , 1-Nx/2:] = u2v2[ :Nz/2 , 1-Nx/2:]
	uv[1-Nz/2:,  :Nx/2 ] = u2v2[1-Nz/2:,  :Nx/2 ]
	uv[1-Nz/2:, 1-Nx/2:] = u2v2[1-Nz/2:, 1-Nx/2:]

	return uv



## operation functions
def get_path_name(file_id, file_type='U', path='fielddata/'):
	file_path = data_path + path + file_type + '/'
	file_name = file_type + str(file_id).zfill(8) + '.BIN'
	return  file_path + file_name

def parse_path_name(file_path_name):
	file_path = file_path_name.split('/')[:-1]
	file_name = file_path_name.split('/')[-1]
	file_id = int(file_name[-12:-4])
	file_type = file_name[:-12]
	path = file_path[-2] + '/'
	return file_id, file_type, path

def read_channel_infosec(file_path_name):
	recl = Nx/2 * Nz * 8
	with open(file_path_name, 'rb') as fp:
		return np.array( unpack( recl/4*'f', fp.read(recl) ) )

def read_channel_layer(file_path_name, layer_id):
	q = np.zeros([Nz,Nx], dtype=complex)
	recl = Nx/2 * Nz * 8
	with open(file_path_name, 'rb') as fp:
		fp.seek(recl * (layer_id + 1)) # +1 for skipping info section
		data = np.array( unpack( recl/4*'f', fp.read(recl) ) ).reshape([Nz,Nx/2,2])

	asemble = np.array([ complex(1,0), complex(0,1) ])
	q[:Nz/2,:Nx/2] = np.sum( data[:Nz/2] * asemble, axis=-1 )
	q[Nz/2+1:,:Nx/2] = np.sum( data[Nz/2+1:] * asemble, axis=-1 )
	asemble = np.array([ complex(1,0), complex(0,-1) ])
	q[:Nz/2,Nx/2+1:] = np.sum( data[range(0,-Nz/2,-1), Nx/2-1:0:-1] * asemble, axis=-1 )
	q[Nz/2+1:,Nx/2+1:] = np.sum( data[Nz/2-1:0:-1, Nx/2-1:0:-1] * asemble, axis=-1 )

	# q[:,:Nx/2] = data[:,:,0] + complex(0,1) * data[:,:,1]
	# q[:,Nx/2+1:] = np.conj( q[[0]+range(Nz-1,0,-1)][:,range(Nx/2-1,0,-1)] )

	return q

def read_channel_cut(file_path_name, layer_ids, z): # read an xt or xy plane
	q = np.zeros([len(layer_ids),Nx/2], dtype=complex)
	recl = Nx/2 * 8
	with open(file_path_name, 'rb') as fp:
		for j in range(len(layer_ids)):
			fp.seek( recl * ((layer_ids[j]+1)*Nz + z) )
			data = np.array( unpack( recl/4*'f', fp.read(recl) ) ).reshape([Nx/2,2])
			q[j] = data[:,0] + complex(0,1) * data[:,1]
	return q

def read_channel_column(file_path_name, layer_ids, x, z):
	q = np.zeros(len(layer_ids), dtype=complex)
	with open(file_path_name, 'rb') as fp:
		for j in range(len(layer_ids)):
			fp.seek( ((layer_ids[j]+1) * Nz*Nx/2 + z * Nx/2 + x) * 8 )
			data = unpack('2f', fp.read(8))
			q[j] = data[0] + complex(0,1) * data[1]
	return q

def read_channel_layer_fluc(file_path_name, layer_id, Qm='default1'):
	q = read_channel_layer(file_path_name, layer_id)

	if type(Qm) != str:
		q[0,0] -= Qm
	elif Qm == 'default0':
		q[0,0] = 0
	elif Qm == 'default1':

		file_id, file_type, path = parse_path_name(file_path_name)
		if file_type in ('U', 'V','W', 'P'):
			q[0,0] -= means[file_type][layer_id]
		elif file_type in ('UU','VV','WW','UV','VW','WU'):
			ft1, ft2 = file_type[0], file_type[1]
			if ft1 == ft2:
				Qm2 = Qm1 = means[ft1][layer_id]
				q2 = q1 = read_channel_layer_fluc(get_path_name(file_id,ft1), layer_id)
			else:
				Qm1 = means[ft1][layer_id]
				Qm2 = means[ft2][layer_id]
				q1 = read_channel_layer_fluc(get_path_name(file_id,ft1), layer_id)
				q2 = read_channel_layer_fluc(get_path_name(file_id,ft2), layer_id)
			q -= Qm1 * q2 + Qm2 * q1
			q[0,0] -= Qm1 * Qm2
		else:
			print 'fluc read of file type %s not supported !' %file_type
			exit()

	return q

def write_channel_infosec(file_path_name, info):
	data = np.ravel(info)
	with open(file_path_name, 'wb') as fp: # !!! note: with 'wb', file cleared if already exsit
		fp.seek(0,0)
		fp.write( pack( len(data)*'f', *data ) )

def write_channel_layer(file_path_name, layer_id, q):
	data = q[:,:Nx/2]
	data = np.ravel( np.transpose( [data.real, data.imag], [1,2,0] ) )
	recl = Nx/2 * Nz * 8
	with open(file_path_name, 'ab') as fp: # note: with 'ab', all writes append to the end of the file regardless of the current seek position
		write_pos = recl * (layer_id + 1)
		fp.seek(0,2) # put the file pointer at the end of file, this must be explicitly done on windows
		blank_size = write_pos - fp.tell()
		if blank_size > 0: # if first time write this file, add 0 for blanks in the info section
			fp.write( pack( blank_size/4*'f', *[0 for n in range(blank_size/4)] ) )
		fp.write( pack( recl/4*'f', *data ) )

def read_channel_fluc_easy(file_path_name, x=range(Nx), y=range(Ny+1), z=range(Nz)):
	q = []
	for j in y:
		qq = read_channel_layer( file_path_name, j )
		qq[0,0] = 0 # mean fluctuation in time not considered
		q.append(phys(qq)[z][:,x])
	return np.squeeze(q)

def read_channel_layer_dy(file_path_name, layer_id):
	if 0 < layer_id < Ny:
		layer_ids = [layer_id-1, layer_id, layer_id+1]
	elif layer_id == 0:
		layer_ids = [layer_id, layer_id+1, layer_id+2]
	elif layer_id == Ny:
		layer_ids = [layer_id-2, layer_id-1, layer_id]

	q = np.array([ read_channel_layer(file_path_name, j) for j in layer_ids ])
	dqdy = np.sum(diff1[layer_id].reshape([3,1,1]) * q, axis=0)
	dqdyy= np.sum(diff2[layer_id].reshape([3,1,1]) * q, axis=0)

	return dqdy, dqdyy

def compute_channel_layer_dt(file_id, layer_id):
	if layer_id in [0, Ny]:
		return [ np.zeros([Nz,Nx], dtype=complex) for n in range(3) ]

	kx, kz = k_x, k_z.reshape([Nz,1])
	im = complex(0,1)

	u, v, w, p = [read_channel_layer(get_path_name(file_id,ft), layer_id) for ft in ['U','V','W','P']]
	dudy, dudyy = read_channel_layer_dy(get_path_name(file_id,'U'), layer_id)
	dvdy, dvdyy = read_channel_layer_dy(get_path_name(file_id,'V'), layer_id)
	dwdy, dwdyy = read_channel_layer_dy(get_path_name(file_id,'W'), layer_id)
	dpdy, dpdyy = read_channel_layer_dy(get_path_name(file_id,'P'), layer_id)

	uu,vv,ww,uv,vw,wu = [read_channel_layer(get_path_name(file_id,ft), layer_id) for ft in ['UU','VV','WW','UV','VW','WU']]
	duvdy, duvdyy = read_channel_layer_dy(get_path_name(file_id,'UV'), layer_id)
	dvvdy, dvvdyy = read_channel_layer_dy(get_path_name(file_id,'VV'), layer_id)
	dvwdy, dvwdyy = read_channel_layer_dy(get_path_name(file_id,'VW'), layer_id)

	# NS equation for whole values
	dudt = im*kx*uu - duvdy + im*kz*wu + im*kx*p - 1.0/Re * (kx**2*u - dudyy + kz**2*u)
	dvdt = im*kx*uv - dvvdy + im*kz*vw - dpdy	 - 1.0/Re * (kx**2*v - dvdyy + kz**2*v)
	dwdt = im*kx*wu - dvwdy + im*kz*ww + im*kz*p - 1.0/Re * (kx**2*w - dwdyy + kz**2*w)
	dudt[0,0] += tau_w

	# # NS equation for fluctuations
	# dudt = im*Um*kx*u + (im*kx*uu - duvdy + im*kz*wu) + im*kx*p - 1.0/Re * (kx**2*u - dudyy + kz**2*u) - Uy*v
	# dvdt = im*Um*kx*v + (im*kx*uv - dvvdy + im*kz*vw) - dpdy	- 1.0/Re * (kx**2*v - dvdyy + kz**2*v)
	# dwdt = im*Um*kx*w + (im*kx*wu - dvwdy + im*kz*ww) + im*kz*p - 1.0/Re * (kx**2*w - dwdyy + kz**2*w)
	# dudt[0,0] += duvdy[0,0]
	# dvdt[0,0] += dvvdy[0,0]
	# dwdt[0,0] += dvwdy[0,0]

	return dudt, dvdt, dwdt


def read_channel_flow_animatioin(file_id, file_type='U'):
	import matplotlib.pyplot as plt
	plt.ion()

	file_path_name = get_path_name(file_id, file_type, 'probedata/')
	time_start, time_end = list( read_channel_infosec(file_path_name)[[0,1]] )
	prblen = int(round( (time_end - time_start) / (dt*nsprb) )) + 1 # total number of time steps probed

	for t in range(0, prblen, int(np.ceil(0.1/(dt*nsprb)))):
		plt.contourf(
			x_mesh, z_mesh,
			phys( read_channel_layer(file_path_name, t) ) - means[file_type][file_id],
			levels=np.linspace(-0.5,0.5,15)	)
		plt.axis('scaled')
		plt.xlabel('x')
		plt.ylabel('z')
		plt.title('t = %.2f'%(t*nsprb*dt))
		plt.pause(0.01)
		plt.close()


## for debug
def debug_div():
	import matplotlib.pyplot as plt

	tstep = tsteps[0]
	j = 21

	u, v, w = [ read_channel_layer(get_path_name(tstep,ft), j) for ft in ('U','V','W') ]
	dudx = - complex(0,1) * k_x * u
	dvdy, dvdyy = read_channel_layer_dy(get_path_name(tstep,'V'), j)
	dwdz = - complex(0,1) * k_z.reshape([Nz,1]) * w
	div = phys( dudx + dvdy + dwdz )

	# dudt, dvdt0, dwdt = compute_channel_layer_dt(tstep, j-1)
	# dudt, dvdt1, dwdt = compute_channel_layer_dt(tstep, j+1)
	# dudt, dvdt, dwdt = compute_channel_layer_dt(tstep, j)
	# dudtdx = - complex(0,1) * k_x * dudt
	# dvdtdy = np.sum(diff1[j].reshape([3,1,1]) * np.array([dvdt0, dvdt, dvdt1]), axis=0)
	# dwdtdz = - complex(0,1) * k_z.reshape([Nz,1]) * dwdt
	# div = phys( dudtdx + dvdtdy + dwdtdz ) * dt


	plt.subplot(211)
	plt.contourf(x_mesh, z_mesh, div)
	plt.subplot(212)
	plt.plot([np.min(div), np.min(div)], [0,Nx*Nz])
	plt.plot([np.max(div), np.max(div)], [0,Nx*Nz])
	plt.plot([np.mean(div), np.mean(div)], [0,Nx*Nz])
	plt.plot([dUdy[j]/10, dUdy[j]/10], [0,Nx*Nz])
	plt.legend(['min div', 'max div', 'mean div', '0.1 mean shear'], loc='best')
	plt.hist(np.ravel(div), bins=100, density=False)
	plt.yscale('log')
	plt.show()


def debug_dudt():
	import matplotlib.pyplot as plt

	tstep = tsteps[0]
	j = 21

	dudx, dvdx, dwdx, dpdx = [ phys(-complex(0,1) * k_x * read_channel_layer(get_path_name(tstep,ft), j)) for ft in ('U','V','W','P') ]
	dudt, dvdt, dwdt, dpdt = [ phys( read_channel_layer(get_path_name(tstep,ft), j) ) for ft in ('UT','VT','WT','PT') ]

	fig, axs = plt.subplots(3,2,squeeze=True)
	axs[0,0].contourf(x_mesh, z_mesh,-dudt, levels=np.linspace(-2.5,2.5,10)*np.var(dudt))
	axs[0,1].contourf(x_mesh, z_mesh, dudx, levels=np.linspace(-1.5,1.5,10)*np.var(dudx))
	axs[1,0].contourf(x_mesh, z_mesh,-dvdt, levels=np.linspace(-2.5,2.5,10)*np.var(dvdt))
	axs[1,1].contourf(x_mesh, z_mesh, dvdx, levels=np.linspace(-1.5,1.5,10)*np.var(dvdx))
	axs[2,0].contourf(x_mesh, z_mesh,-dwdt, levels=np.linspace(-2.5,2.5,10)*np.var(dwdt))
	axs[2,1].contourf(x_mesh, z_mesh, dwdx, levels=np.linspace(-1.5,1.5,10)*np.var(dwdx))
	axs[0,0].set_title(r"$\partial_tu$")
	axs[0,1].set_title(r"$\partial_xu$")
	axs[1,0].set_title(r"$\partial_tv$")
	axs[1,1].set_title(r"$\partial_xv$")
	axs[2,0].set_title(r"$\partial_tw$")
	axs[2,1].set_title(r"$\partial_xw$")
	plt.tight_layout()
	plt.show()


# read_channel_flow_animatioin(21, 'U')

# debug_dudt()
# debug_div()
