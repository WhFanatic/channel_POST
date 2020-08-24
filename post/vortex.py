from basic import *
from operators import Operators


class Vortex:
	def __init__(self, para):
		self.para = para


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

				vtx_x = phys( 1j * kz * v + dwdy		)
				vtx_y = phys(-1j *(kz * u - kx * w)	)
				vtx_z = phys(-1j * kx * v - dudy		)
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

		