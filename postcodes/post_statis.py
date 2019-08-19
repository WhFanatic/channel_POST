#!/work1/cuigx2_work/whn/anaconda_install/anaconda2/bin/python
from post_channel import *


def compute_statis():
	RS_uu, RS_vv, RS_ww, RS_uv = [ np.zeros(Ny+1) for n in range(4) ]
	CP_pp, CP_up, CP_vp, CP_wp = [ np.zeros(Ny+1) for n in range(4) ]
	VX_xx, VX_yy, VX_zz, VX_xy = [ np.zeros(Ny+1) for n in range(4) ]
	kx, kz = k_x, k_z.reshape([Nz,1])

	for tstep in tsteps:
		print 'reading statis of time step %i ...'%tstep

		UU_mean = read_channel_column( get_path_name(tstep,'UU'), range(Ny+1), 0, 0 ).real
		VV_mean = read_channel_column( get_path_name(tstep,'VV'), range(Ny+1), 0, 0 ).real
		WW_mean = read_channel_column( get_path_name(tstep,'WW'), range(Ny+1), 0, 0 ).real
		UV_mean = read_channel_column( get_path_name(tstep,'UV'), range(Ny+1), 0, 0 ).real

		for j in range(Ny+1):
			u, v, w, p = [read_channel_layer(get_path_name(tstep,ft), j) for ft in ['U','V','W','P']]
			dudy,dudyy = read_channel_layer_dy(get_path_name(tstep,'U'), j)
			dwdy,dwdyy = read_channel_layer_dy(get_path_name(tstep,'W'), j)

			u[0,0] -= U_mean[j]
			p[0,0] -= P_mean[j]
			dudy[0,0] -= dUdy[j]
			dudyy[0,0] -= dUdyy[j]

			vtx_x = phys( complex(0,1) * kz * v + dwdy		)
			vtx_y = phys(-complex(0,1) *(kz * u - kx * w)	)
			vtx_z = phys(-complex(0,1) * kx * v - dudy		)
			u, v, w, p = phys(u), phys(v), phys(w), phys(p)

			RS_uu[j] += UU_mean[j] - U_mean[j]**2
			RS_vv[j] += VV_mean[j]
			RS_ww[j] += WW_mean[j]
			RS_uv[j] += UV_mean[j]

			CP_pp[j] += np.mean(p**2)
			CP_up[j] += np.mean(u*p)
			CP_vp[j] += np.mean(v*p)
			CP_wp[j] += np.mean(w*p)

			VX_xx[j] += np.mean(vtx_x**2)
			VX_yy[j] += np.mean(vtx_y**2)
			VX_zz[j] += np.mean(vtx_z**2)
			VX_xy[j] += np.mean(vtx_x*vtx_y)

	# average between upper & lower half and among time steps
	RS_uu = (RS_uu + RS_uu[::-1]) / 2 / len(tsteps)
	RS_vv = (RS_vv + RS_vv[::-1]) / 2 / len(tsteps)
	RS_ww = (RS_ww + RS_ww[::-1]) / 2 / len(tsteps)
	RS_uv = (RS_uv - RS_uv[::-1]) / 2 / len(tsteps) # v reverses sign across center of channel

	CP_pp = (CP_pp + CP_pp[::-1]) / 2 / len(tsteps)
	CP_up = (CP_up + CP_up[::-1]) / 2 / len(tsteps)
	CP_vp = (CP_vp - CP_vp[::-1]) / 2 / len(tsteps) # v reverses sign across center of channel
	CP_wp = (CP_wp + CP_wp[::-1]) / 2 / len(tsteps)

	VX_xx = (VX_xx + VX_xx[::-1]) / 2 / len(tsteps)
	VX_yy = (VX_yy + VX_yy[::-1]) / 2 / len(tsteps)
	VX_zz = (VX_zz + VX_zz[::-1]) / 2 / len(tsteps)
	VX_xy = (VX_xy - VX_xy[::-1]) / 2 / len(tsteps) # v reverses sign across center of channel

	statis = np.array([
		RS_uu, RS_vv, RS_ww, RS_uv,
		CP_pp, CP_up, CP_vp, CP_wp,
		VX_xx, VX_yy, VX_zz, VX_xy	])

	return statis

def compute_ES2D(layer_id):
	es2D = np.zeros([4, Nz, Nx]) # Euu, Evv, Eww, Epp
	for tstep in tsteps:
		es2D += np.array([	abs( read_channel_layer(get_path_name(tstep,ft),	layer_id) )**2	for ft in ['U','V','W','P']])
		es2D += np.array([	abs( read_channel_layer(get_path_name(tstep,ft), Ny-layer_id) )**2	for ft in ['U','V','W','P']])
	es2D /= len(tsteps) * 2
	es2D[:,0,0] = 0
	return es2D # [4, Nz, Nx]

def compute_BGT2D(layer_id):
	bgt = np.zeros([3, 7, Nz, Nx]) # budgets of Reynolds stress transportation equations
	flx = np.zeros([3, 3, Nz, Nx]) # energy flux contributing to diffusion terms
	kx, kz = k_x, k_z.reshape([Nz,1])

	# get means of this layer
	Um = U_mean[layer_id]
	Pm = P_mean[layer_id]
	Py = dPdy[layer_id]
	Uy = dUdy[layer_id]
	Uyy = dUdyy[layer_id]
	Pyy = dPdyy[layer_id]

	for tstep in tsteps:
		# read an instantaneous field
		u, v, w, p = [read_channel_layer(get_path_name(tstep,ft), layer_id) for ft in ['U','V','W','P']]
		dudy, dudyy = read_channel_layer_dy(get_path_name(tstep,'U'), layer_id)
		dvdy, dvdyy = read_channel_layer_dy(get_path_name(tstep,'V'), layer_id)
		dwdy, dwdyy = read_channel_layer_dy(get_path_name(tstep,'W'), layer_id)
		dpdy, dpdyy = read_channel_layer_dy(get_path_name(tstep,'P'), layer_id)

		uu,vv,ww,uv,vw,wu = [read_channel_layer(get_path_name(tstep,ft), layer_id) for ft in ['UU','VV','WW','UV','VW','WU']]
		duvdy, duvdyy = read_channel_layer_dy(get_path_name(tstep,'UV'), layer_id)
		dvvdy, dvvdyy = read_channel_layer_dy(get_path_name(tstep,'VV'), layer_id)
		dvwdy, dvwdyy = read_channel_layer_dy(get_path_name(tstep,'VW'), layer_id)

		# get fluctuations
		u[0,0] -= Um
		p[0,0] -= Pm
		dudy[0,0] -= Uy
		dpdy[0,0] -= Py
		dudyy[0,0] -= Uyy
		dpdyy[0,0] -= Pyy

		uu -= 2 * Um * u
		uu[0,0] -= Um**2
		uv -= Um * v
		wu -= Um * w
		duvdy -= Uy * v + Um * dvdy
		duvdyy -= Uyy * v + 2 * Uy * dvdy + Um * dvdyy

		# production
		bgt[0,0] += -2 * dUdy[layer_id] * (np.conj(u) * v).real

		# turbulent flux & diffusion
		flx[0,0] += - (np.conj(u) * uv).real
		flx[1,0] += - (np.conj(v) * vv).real
		flx[2,0] += - (np.conj(w) * vw).real
		bgt[0,1] += - (np.conj(u) * duvdy + np.conj(dudy) * uv).real
		bgt[1,1] += - (np.conj(v) * dvvdy + np.conj(dvdy) * vv).real
		bgt[2,1] += - (np.conj(w) * dvwdy + np.conj(dwdy) * vw).real

		# turbulent transportation among scales
		bgt[0,2] += 2*kx * (np.conj(uu) * u).imag + 2*kz * (np.conj(wu) * u).imag - (np.conj(u) * duvdy + np.conj(dudy) * uv).real + 2 * (np.conj(uv) * dudy).real
		bgt[1,2] += 2*kx * (np.conj(uv) * v).imag + 2*kz * (np.conj(vw) * v).imag - (np.conj(v) * dvvdy + np.conj(dvdy) * vv).real + 2 * (np.conj(vv) * dvdy).real
		bgt[2,2] += 2*kx * (np.conj(wu) * w).imag + 2*kz * (np.conj(ww) * w).imag - (np.conj(w) * dvwdy + np.conj(dwdy) * vw).real + 2 * (np.conj(vw) * dwdy).real

		# pressure flux & diffusion
		flx[1,1] += -2 * (np.conj(p) * v).real
		bgt[1,3] += -2 * (v * np.conj(dpdy) + dvdy * np.conj(p)).real

		# pressure redistribution
		bgt[0,4] += 2 * kx * (np.conj(p) * u).imag
		bgt[1,4] += 2 * (np.conj(p) * dvdy).real
		bgt[2,4] += 2 * kz * (np.conj(p) * w).imag

		# viscous flux & diffusion
		flx[0,2] += 2.0*nu * (np.conj(u) * dudy).real
		flx[1,2] += 2.0*nu * (np.conj(v) * dvdy).real
		flx[2,2] += 2.0*nu * (np.conj(w) * dwdy).real
		bgt[0,5] += 2.0*nu * ( abs(dudy)**2 + (np.conj(u) * dudyy).real )
		bgt[1,5] += 2.0*nu * ( abs(dvdy)**2 + (np.conj(v) * dvdyy).real )
		bgt[2,5] += 2.0*nu * ( abs(dwdy)**2 + (np.conj(w) * dwdyy).real )

		# viscous dissipation
		bgt[0,6] += -2.0*nu * ( (kx*abs(u))**2 + (kz*abs(u))**2 + abs(dudy)**2 )
		bgt[1,6] += -2.0*nu * ( (kx*abs(v))**2 + (kz*abs(v))**2 + abs(dvdy)**2 )
		bgt[2,6] += -2.0*nu * ( (kx*abs(w))**2 + (kz*abs(w))**2 + abs(dwdy)**2 )

	bgt /= len(tsteps)
	flx /= len(tsteps)
	print 'budget balance', np.sum(bgt), 'at layer j =', layer_id

	return bgt, flx




def get_statis():
	try:
		return np.loadtxt(open(postdata_path + 'statis.txt')).T [1:] # drop y_mesh record
	except IOError:
		statis = compute_statis()
		np.savetxt(postdata_path + 'statis.txt', np.array([y_mesh]+list(statis)).T)
		return statis # [12, Ny+1]

def get_ES2D(layer_id):
	file_path_names = [ postdata_path+'ES2D/ES2D_'+ft+'.BIN' for ft in ('UU','VV','WW','PP') ]
	try:
		j = (layer_id, Ny-layer_id) [int(layer_id > Ny/2)]
		return np.array([ abs( read_channel_layer(fpn, j) ) for fpn in file_path_names ]) # [4, Nz, Nx]
	except IOError:
		for fpn in file_path_names:
			write_channel_infosec(fpn, [0])
		for j in range(Ny/2+1):
			print 'Computing energy spectra of slice %i ...'%j
			es2Ds = list( compute_ES2D(j) )
			for fpn in file_path_names:
				write_channel_layer( fpn, j, es2Ds.pop(0) * complex(1,0) )
		return get_ES2D(layer_id) # [4, Nz, Nx]

def get_ES1D():
	es1D_x = np.zeros([4, Ny+1, Nx])
	es1D_z = np.zeros([4, Ny+1, Nz])
	for j in range(Ny/2+1):
		es2D = get_ES2D(j)
		es1D_x[:,j,:] = es1D_x[:,Ny-j,:] = np.sum(es2D, axis=-2)
		es1D_z[:,j,:] = es1D_z[:,Ny-j,:] = np.sum(es2D, axis=-1)
	return es1D_x, es1D_z # [4, Ny+1, Nx], [4, Ny+1, Nz]

def get_BGT2D(layer_id):
	return compute_BGT2D(layer_id) # bgt:[3, 7, Nz, Nx], flx:[3, 3, Nz, Nx]

def get_BGT():
	try:
		bgtflx = np.loadtxt(open(postdata_path + 'budgets.txt')).T [1:] # drop y_mesh record
		bgtflx = bgtflx.reshape([3, 10, Ny+1])
	except IOError:
		bgtflx = np.zeros([3, 10, Ny+1])
		for j in range(Ny+1):
			bgt2D, flx2D = get_BGT2D(j)
			bgtflx[:,:7,j] = np.sum(bgt2D, axis=(-1,-2))
			bgtflx[:,7:,j] = np.sum(flx2D, axis=(-1,-2))
		np.savetxt(postdata_path + 'budgets.txt', np.array([y_mesh]+list(bgtflx.reshape([30,Ny+1]))).T)

	bgt = bgtflx[:,:7,:]
	flx = bgtflx[:,7:,:]
	return bgt, flx # bgt:[3, 7, Ny+1], flx:[3, 3, Ny+1]







