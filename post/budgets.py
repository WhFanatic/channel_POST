from basic import *
from operators import Operators


class Budgets:
	def __init__(self, para):
		self.para = para

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

