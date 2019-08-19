#!/work/cuigx2_work/whn/anaconda_install/anaconda2/bin/python
from post_channel import *

pressure_path = postdata_path + 'pressure/'

# chasing method for tri-diag equation
def chasing(a,b,c,d):
	N = len(a)
	l = np.zeros(N)
	r = np.zeros(N)
	m = np.zeros(N, dtype=complex)
	q = np.zeros(N, dtype=complex)

	r[0] = b[0]
	m[0] = d[0]
	for j in range(1,N):
		l[j] = a[j] / r[j-1]
		r[j] = b[j] - l[j]*c[j-1]
		m[j] = d[j] - l[j]*m[j-1]
	q[-1] = m[-1] / r[-1]
	for j in range(N-1)[::-1]:
		q[j] = (m[j] - c[j]*q[j+1]) / r[j]

	return q

def bp_generate(file_id, ifwall=False):
	bp = np.zeros([2,Nz,Nx/2], dtype=complex)

	if ifwall:	# it is tested that, for boundaries at the wall, the 'ifwall' flag does not make discernible difference (relative error 1e-9). The only source of difference is from dvvdy
		dvdyy0 = np.sum( [diff2[0,j] * read_channel_layer(get_path_name(file_id, 'V'), j) for j in (0,1,2)], axis=0 )
		dvdyy1 = np.sum( [diff2[-1,j]* read_channel_layer(get_path_name(file_id, 'V'), j+Ny-2) for j in (0,1,2)], axis=0 )
		bp[0] = nu * dvdyy0[:,:Nx/2]
		bp[1] = nu * dvdyy1[:,:Nx/2]
		bp[:,0,0] = 0 # Dirichlet BC at wall for (0,0) mode
		return bp

	im = complex(0,1)
	kx, kz = k_x, k_z.reshape([Nz,1])

	for j in (0, Ny):
		js = np.arange(3, dtype=int) + (Ny-2) * (j/Ny)

		uv = read_channel_layer(get_path_name(file_id,'UV'), j)
		vw = read_channel_layer(get_path_name(file_id,'VW'), j)
		vt = read_channel_layer(get_path_name(file_id,'VT'), j)

		vv0= read_channel_layer(get_path_name(file_id,'VV'), js[0])
		vv1= read_channel_layer(get_path_name(file_id,'VV'), js[1])
		vv2= read_channel_layer(get_path_name(file_id,'VV'), js[2])

		v0 = read_channel_layer(get_path_name(file_id, 'V'), js[0])
		v1 = read_channel_layer(get_path_name(file_id, 'V'), js[1])
		v2 = read_channel_layer(get_path_name(file_id, 'V'), js[2])

		dvvdy = diff1[j,0] *vv0 + diff1[j,1] *vv1 + diff1[j,2] *vv2
		dvdyy = diff2[j,0] * v0 + diff2[j,1] * v1 + diff2[j,2] * v2

		bp[j/Ny] = ( nu * dvdyy + im * (kx*uv + kz*vw) - dvvdy - vt )[:,:Nx/2]

	bp[:,0,0] = - read_channel_column(get_path_name(file_id,'VV'), [0,Ny], 0, 0) # Dirichlet BC at wall for (0,0) mode

	return bp

def rhs_generate(file_id):
	rhs = np.zeros([Ny+1,Nz,Nx/2], dtype=complex)

	im = complex(0,1)
	kx, kz = k_x, k_z.reshape([Nz,1])

	vv0,uv0,vw0 = [read_channel_layer(get_path_name(file_id,ft), 0) for ft in ('VV','UV','VW')]
	vv1,uv1,vw1 = [read_channel_layer(get_path_name(file_id,ft), 1) for ft in ('VV','UV','VW')]

	for j in range(1,Ny):

		uu, ww, wu = [read_channel_layer(get_path_name(file_id,ft), j) for ft in ('UU','WW','WU')]
		vv2,uv2,vw2= [read_channel_layer(get_path_name(file_id,ft), j+1) for ft in ('VV','UV','VW')]

		duvdy = diff1[j,0] * uv0 + diff1[j,1] * uv1 + diff1[j,2] * uv2
		dvwdy = diff1[j,0] * vw0 + diff1[j,1] * vw1 + diff1[j,2] * vw2
		dvvdyy= diff2[j,0] * vv0 + diff2[j,1] * vv1 + diff2[j,2] * vv2

		rhs[j] = ( kx**2 * uu + kz**2 * ww - dvvdyy + 2 * ( im*kx*duvdy + im*kz*dvwdy + kx*wu*kz ) )[:,:Nx/2]

		if j == 1:
			uu, ww, wu = [read_channel_layer(get_path_name(file_id,ft), 0) for ft in ('UU','WW','WU')]
			duvdy = diff1[0,0] * uv0 + diff1[0,1] * uv1 + diff1[0,2] * uv2
			dvwdy = diff1[0,0] * vw0 + diff1[0,1] * vw1 + diff1[0,2] * vw2
			dvvdyy= diff2[0,0] * vv0 + diff2[0,1] * vv1 + diff2[0,2] * vv2
			rhs[0] = ( kx**2 * uu + kz**2 * ww - dvvdyy + 2 * ( im*kx*duvdy + im*kz*dvwdy + kx*wu*kz ) )[:,:Nx/2]

		if j == Ny-1:
			uu, ww, wu = [read_channel_layer(get_path_name(file_id,ft), Ny) for ft in ('UU','WW','WU')]
			duvdy = diff1[-1,0] * uv0 + diff1[-1,1] * uv1 + diff1[-1,2] * uv2
			dvwdy = diff1[-1,0] * vw0 + diff1[-1,1] * vw1 + diff1[-1,2] * vw2
			dvvdyy= diff2[-1,0] * vv0 + diff2[-1,1] * vv1 + diff2[-1,2] * vv2
			rhs[-1] = ( kx**2 * uu + kz**2 * ww - dvvdyy + 2 * ( im*kx*duvdy + im*kz*dvwdy + kx*wu*kz ) )[:,:Nx/2]

		vv0,uv0,vw0 = vv1,uv1,vw1
		vv1,uv1,vw1 = vv2,uv2,vw2

	return rhs

# Poisson equation for pressure
def poisson(file_id, bp, rhs, file_type='P'): # bp: [2,Nz,Nx/2], rhs: [Ny+1,Nz,Nx/2]
	p = np.zeros([Ny+1, Nz, Nx/2], dtype=complex)
	pj= np.zeros([Nz, Nx], dtype=complex)

	BP = bp.copy()[:2,:Nz,:Nx/2] * complex(1,0)
	RHS = rhs.copy()[:Ny+1,:Nz,:Nx/2] * complex(1,0)
	a, b, c = list( diff2.copy().T )
	d, e, f = list( diff1.copy()[[0,-1]].T )

	a[0] = a[0]*f[0] - c[0]*d[0]
	b[0] = b[0]*f[0] - c[0]*e[0]
	temp = c[0]
	c[0] = b[0]
	b[0] = a[0]
	a[0] = temp

	c[-1] = c[-1]*d[-1] - a[-1]*f[-1]
	b[-1] = b[-1]*d[-1] - a[-1]*e[-1]
	temp = a[-1]
	a[-1] = b[-1]
	b[-1] = c[-1]
	c[-1] = temp

	RHS[0] = RHS[0] * f[0] - BP[0] * a[0]
	RHS[-1]= RHS[-1]* d[-1]- BP[-1]* c[-1]

	for k in range(Nz):
		print file_type + str(file_id) + ' solving progress %.2f%%'%(100.0*k/Nz)
		for i in range(Nx/2):

			if not ( i == 0 and k == 0 ):	# Neumann BC at wall
				B = - ( k_x[i]**2 + k_z[k]**2 ) * np.ones(Ny+1)
				B[0] *= f[0]
				B[-1]*= d[-1]

				p[:,k,i] = chasing( a, b+B, c, RHS[:,k,i] )

			elif i == 0 and k == 0:	# Dirichlet BC at wall for (0,0) mode
				A, B, C = [ np.zeros(Ny+1) for n in range(3) ]
				B[0] = a[0] * ( d[0] - f[0]*a[1]/c[1] ) - 1
				C[0] = a[0] * ( e[0] - f[0]*b[1]/c[1] )
				A[-1]= c[-1]* ( e[-1]- d[-1]*b[-2]/a[-2] )
				B[-1]= c[-1]* ( f[-1]- d[-1]*c[-2]/a[-2] ) - 1

				RHS[0,0,0] += BP[0,0,0] * ( a[0] - 1 ) - RHS[1,0,0] * a[0]*f[0]/c[1]
				RHS[-1,0,0]+= BP[-1,0,0]* ( c[-1]- 1 ) - RHS[-2,0,0]* c[-1]*d[-1]/a[-2]

				p[:,0,0] = chasing( a+A, b+B, c+C, RHS[:,0,0] )

	file_name = file_type + str(file_id).zfill(8) + '.BIN'
	fpn = pressure_path + file_name
	write_channel_infosec(fpn, [0])
	for j in range(Ny+1):
		pj[:,:Nx/2] = p[j]
		pj[:,Nx/2+1:] = np.conj( p[j,[0]+range(Nz-1,0,-1)][:,range(Nx/2-1,0,-1)] )
		write_channel_layer(fpn, j, pj)

	return fpn


def green_function(k):
	G = np.zeros([Ny+1,Ny+1])
	dGdy = np.zeros([Ny+1,Ny+1])
	for j in range(Ny+1):
		y = y_mesh[j]

		eta = y_mesh[:j]
		if k == 0:
			G[j,:j] = (y-eta) / 2
			dGdy[j,:j] = 1
		else:
			G[j,:j] = - np.cosh(k*(y-1)) * np.cosh(k*(eta+1)) / ( 2*k * np.cosh(k) * np.sinh(k) )
			dGdy[j,:j] = - np.sinh(k*(y-1)) * np.cosh(k*(eta+1)) / ( 2 * np.cosh(k) * np.sinh(k) )

		eta = y_mesh[j:]
		if k == 0:
			G[j,j:] = (eta-y) / 2
			dGdy[j,j:] = -1
		else:
			G[j,j:] = - np.cosh(k*(y+1)) * np.cosh(k*(eta-1)) / ( 2*k * np.cosh(k) * np.sinh(k) )
			dGdy[j,j:] = - np.sinh(k*(y+1)) * np.cosh(k*(eta-1)) / ( 2 * np.cosh(k) * np.sinh(k) )

	return G, dGdy

def green_function_int(G,R):
	p = []
	a, b, c = list( diff1.copy().T )
	e1, e2 = c[0]/c[1], a[-1]/a[-2]

	a[0] -= a[1] * e1
	b[0] -= b[1] * e1
	c[0] = b[0]
	b[0] = a[0]
	a[0] = 0

	b[-1] -= b[-2] * e2
	c[-1] -= c[-2] * e2
	a[-1] = b[-1]
	b[-1] = c[-1]
	c[-1] = 0

	b[0] = 1e20 # Dirichlet BC: q(eta=-1) = 0

	for d in G*R:
		d[0] -= d[1] * e1
		d[-1] -= d[-2] * e2
		q = chasing(a,b,c,d)
		p.append(q[-1] - q[0])

	return np.array(p, dtype=complex)

def pressure_poisson(file_id):
	return poisson( file_id, bp_generate(file_id,True), rhs_generate(file_id), 'P' )

def pressure_green(file_id, file_type='GP'):
	p = np.zeros([Ny+1, Nz, Nx/2], dtype=complex)
	pj= np.zeros([Nz, Nx], dtype=complex)
	rhs = rhs_generate(file_id)
	for k in range(Nz):
		print file_type + str(file_id) + ' solving progress %.2f%%'%(100.0*k/Nz)
		for i in range(Nx/2):
			G = green_function((k_x[i]**2+k_z[k]**2)**0.5)[0]
			R = rhs[:,k,i]
			p[:,k,i] = green_function_int(G,R)

	file_name = file_type + str(file_id).zfill(8) + '.BIN'
	fpn = pressure_path + file_name
	write_channel_infosec(fpn, [0])
	for j in range(Ny+1):
		pj[:,:Nx/2] = p[j]
		pj[:,Nx/2+1:] = np.conj( p[j,[0]+range(Nz-1,0,-1)][:,range(Nx/2-1,0,-1)] )
		write_channel_layer(fpn, j, pj)

	return fpn

# decompose pressure into rapid term and slow term
def pressure_decompose(file_id):
	im = complex(0,1)
	kx, kz = k_x, k_z.reshape([Nz,1])
	plane = np.ones([Nz,Nx])

	# rapid term
	v = np.array([ read_channel_layer_fluc(get_path_name(file_id, 'V'), j) for j in range(Ny+1) ])
	rhs = ( 2 * im * kx * dUdy.reshape([Ny+1,1,1]) * v )[:,:,:Nx/2]
	bp = np.zeros([2,Nz,Nx/2])
	fpn1 = poisson(file_id, bp=bp, rhs=rhs, file_type='PR')

	# slow term
	rhs = np.zeros([Ny+1,Nz,Nx/2], dtype=complex)

	vv0,uv0,vw0 = [read_channel_layer_fluc(get_path_name(file_id,ft), 0) for ft in ('VV','UV','VW')]
	vv1,uv1,vw1 = [read_channel_layer_fluc(get_path_name(file_id,ft), 1) for ft in ('VV','UV','VW')]

	for j in range(1,Ny):
		uu, ww, wu = [read_channel_layer_fluc(get_path_name(file_id,ft), j) for ft in ('UU','WW','WU')]
		vv2,uv2,vw2= [read_channel_layer_fluc(get_path_name(file_id,ft), j+1) for ft in ('VV','UV','VW')]

		duvdy = diff1[j,0] * uv0 + diff1[j,1] * uv1 + diff1[j,2] * uv2
		dvwdy = diff1[j,0] * vw0 + diff1[j,1] * vw1 + diff1[j,2] * vw2
		dvvdyy= diff2[j,0] * vv0 + diff2[j,1] * vv1 + diff2[j,2] * vv2
		rhs[j] = ( kx**2 * uu + kz**2 * ww - dvvdyy + 2 * ( im*kx*duvdy + im*kz*dvwdy + kx*wu*kz ) )[:,:Nx/2]

		if j == 1:
			uu, ww, wu = [read_channel_layer_fluc(get_path_name(file_id,ft), 0) for ft in ('UU','WW','WU')]
			duvdy = diff1[0,0] * uv0 + diff1[0,1] * uv1 + diff1[0,2] * uv2
			dvwdy = diff1[0,0] * vw0 + diff1[0,1] * vw1 + diff1[0,2] * vw2
			dvvdyy= diff2[0,0] * vv0 + diff2[0,1] * vv1 + diff2[0,2] * vv2
			rhs[0] = ( kx**2 * uu + kz**2 * ww - dvvdyy + 2 * ( im*kx*duvdy + im*kz*dvwdy + kx*wu*kz ) )[:,:Nx/2]	

		if j == Ny-1:
			uu, ww, wu = [read_channel_layer_fluc(get_path_name(file_id,ft), Ny) for ft in ('UU','WW','WU')]
			duvdy = diff1[-1,0] * uv0 + diff1[-1,1] * uv1 + diff1[-1,2] * uv2
			dvwdy = diff1[-1,0] * vw0 + diff1[-1,1] * vw1 + diff1[-1,2] * vw2
			dvvdyy= diff2[-1,0] * vv0 + diff2[-1,1] * vv1 + diff2[-1,2] * vv2
			rhs[-1] = ( kx**2 * uu + kz**2 * ww - dvvdyy + 2 * ( im*kx*duvdy + im*kz*dvwdy + kx*wu*kz ) )[:,:Nx/2]

		vv0,uv0,vw0 = vv1,uv1,vw1
		vv1,uv1,vw1 = vv2,uv2,vw2

	bp = np.zeros([2,Nz,Nx/2], dtype=complex)
	for j in (0, Ny):
		vv = read_channel_layer_fluc(get_path_name(file_id,'VV'), j)
		bp[j/Ny,0,0] = - vv[0,0]

	fpn2 = poisson(file_id, bp=bp, rhs=rhs, file_type='PS')

	# slow term + stokes term
	temp = bp[:,0,0]
	bp = np.zeros([2,Nz,Nx/2], dtype=complex)

	for j in (0, Ny):
		js = np.arange(3, dtype=int) + (Ny-2) * (j/Ny)

		Um = U_mean[j]

		uv = read_channel_layer_fluc(get_path_name(file_id,'UV'), j)
		vw = read_channel_layer_fluc(get_path_name(file_id,'VW'), j)
		vt = read_channel_layer(get_path_name(file_id,'VT'), j)

		vv0= read_channel_layer_fluc(get_path_name(file_id,'VV'), js[0])
		vv1= read_channel_layer_fluc(get_path_name(file_id,'VV'), js[1])
		vv2= read_channel_layer_fluc(get_path_name(file_id,'VV'), js[2])

		v0 = read_channel_layer_fluc(get_path_name(file_id, 'V'), js[0])
		v1 = read_channel_layer_fluc(get_path_name(file_id, 'V'), js[1])
		v2 = read_channel_layer_fluc(get_path_name(file_id, 'V'), js[2])

		dvvdy = diff1[j,0] *vv0 + diff1[j,1] *vv1 + diff1[j,2] *vv2
		dvdyy = diff2[j,0] * v0 + diff2[j,1] * v1 + diff2[j,2] * v2
		v = (v0, v2)[j/Ny]

		bp[j/Ny] = ( nu * dvdyy + im * ( kx*(uv+Um*v) + kz*vw ) - dvvdy - vt )[:,:Nx/2]

	bp[:,0,0] = temp

	fpn3 = poisson(file_id, bp=bp, rhs=rhs, file_type='PSST')

	# stokes term
	rhs = np.zeros([Ny+1,Nz,Nx/2])
	bp[:,0,0] = 0
	fpn4 = poisson(file_id, bp=bp, rhs=rhs, file_type='PST')


	return fpn1, fpn2, fpn3, fpn4


def rhs_generate_decomp(file_id, layer_id, scale=1):
	rhs = np.zeros([6,Nz,Nx], dtype=complex)

	im = complex(0,1)
	kx, kz = k_x, k_z.reshape([Nz,1])

	u,v,w = [ read_channel_layer(get_path_name(file_id,ft), layer_id) * scale for ft in ('U','V','W') ]
	dudy,dvdy,dwdy = [ read_channel_layer_dy(get_path_name(file_id,ft), layer_id) * scale for ft in ('U','V','W') ]

	rhs[0] = convol( -im*kx*u )
	rhs[1] = convol( dvdy )
	rhs[2] = convol( -im*kz*w )
	rhs[3] = convol( dudy, -im*kx*v )
	rhs[4] = convol( -im*kz*v, dwdy )
	rhs[5] = convol( -im*kx*w, -im*kz*u )

	return rhs



# for test, output only the divergence related terms in RHS of pressure Poisson eq 
def rhs2_generate(file_id):
	rhs = np.zeros([Ny+1,Nz,Nx/2], dtype=complex)

	im = complex(0,1)
	kx, kz = k_x, k_z.reshape([Nz,1])

	v0 = read_channel_layer(get_path_name(file_id,'V'), 0)
	v1 = read_channel_layer(get_path_name(file_id,'V'), 1)

	for j in range(1,Ny):
		u, w = [read_channel_layer(get_path_name(file_id,ft), j) for ft in ('U','W')]
		v2 = read_channel_layer(get_path_name(file_id,'V'), j+1)
		dvdy = diff1[j,0] * v0 + diff1[j,1] * v1 + diff1[j,2] * v2
		rhs[j] = ( convol(-kx*im*u) + convol(dvdy) + convol(-kz*im*w) )[:,:Nx/2]

		if j == 1:
			u, w = [read_channel_layer(get_path_name(file_id,ft), 0) for ft in ('U','W')]
			dvdy = diff1[0,0] * v0 + diff1[0,1] * v1 + diff1[0,2] * v2
			rhs[0] = ( convol(-kx*im*u) + convol(dvdy) + convol(-kz*im*w) )[:,:Nx/2]

		if j == Ny-1:
			u, w = [read_channel_layer(get_path_name(file_id,ft), Ny) for ft in ('U','W')]
			dvdy = diff1[-1,0] * v0 + diff1[-1,1] * v1 + diff1[-1,2] * v2
			rhs[-1] = ( convol(-kx*im*u) + convol(dvdy) + convol(-kz*im*w) )[:,:Nx/2]

		v0 = v1
		v1 = v2

	return rhs
