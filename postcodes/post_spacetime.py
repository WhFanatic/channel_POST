#!/work1/cuigx2_work/whn/anaconda_install/anaconda2/bin/python
from post_channel import *
from post_statis import *


def compute_ES2D_DT(layer_id):
	ududt, dudtdudt = np.zeros([Nz,Nx], dtype=complex), np.zeros([Nz,Nx])
	vdvdt, dvdtdvdt = np.zeros([Nz,Nx], dtype=complex), np.zeros([Nz,Nx])
	wdwdt, dwdtdwdt = np.zeros([Nz,Nx], dtype=complex), np.zeros([Nz,Nx])
	pdpdt, dpdtdpdt = np.zeros([Nz,Nx], dtype=complex), np.zeros([Nz,Nx])

	for tstep in tsteps:
		for j in (layer_id, Ny-layer_id):
			u, v, w, p = [ read_channel_layer(get_path_name(tstep,ft), j) for ft in ('U','V','W','P') ]
			dudt, dvdt, dwdt, dpdt = [ read_channel_layer(get_path_name(tstep,ft), j) for ft in ('UT','VT','WT','PT') ]

			u[0,0] -= U_mean[j]
			p[0,0] -= P_mean[j]

			ududt += np.conj(u) * dudt
			vdvdt += np.conj(v) * dvdt
			wdwdt += np.conj(w) * dwdt
			pdpdt += np.conj(p) * dpdt

			dudtdudt += abs(dudt)**2
			dvdtdvdt += abs(dvdt)**2
			dwdtdwdt += abs(dwdt)**2
			dpdtdpdt += abs(dpdt)**2

	ududt /= len(tsteps) * 2
	vdvdt /= len(tsteps) * 2
	wdwdt /= len(tsteps) * 2
	pdpdt /= len(tsteps) * 2
	dudtdudt /= len(tsteps) * 2
	dvdtdvdt /= len(tsteps) * 2
	dwdtdwdt /= len(tsteps) * 2
	dpdtdpdt /= len(tsteps) * 2

	return 	np.array([ududt, vdvdt, wdwdt, pdpdt, dudtdudt, dvdtdvdt, dwdtdwdt, dpdtdpdt])

def get_ES2D_DT(layer_id):
	file_types = ( 'UDUDT','VDVDT','WDWDT','PDPDT', 'DUDTDUDT','DVDTDVDT','DWDTDWDT','DPDTDPDT'	)
	file_path_names = [ postdata_path+'ES2D/ES2D_'+ft+'.BIN' for ft in file_types ]
	try:
		j = (layer_id, Ny-layer_id) [int(layer_id > Ny/2)]
		es2Ds = np.array([ read_channel_layer(fpn, j) for fpn in file_path_names ])
		return es2Ds[:4], abs(es2Ds[4:]) # es2Ds_DT: [4, Nz, Nx] complex;  es2Ds_DTDT: [4, Nz, Nx] float
	except IOError:
		for fpn in file_path_names:
			write_channel_infosec(fpn, [0])
		for j in range(Ny/2+1):
			print 'Computing dt spectra of slice %i ...'%j
			es2Ds = list( compute_ES2D_DT(j) )
			for fpn in file_path_names:
				write_channel_layer(fpn, j, es2Ds.pop(0) * complex(1,0))
		return get_ES2D_DT(layer_id) # es2Ds_DT: [4, Nz, Nx] complex;  es2Ds_DTDT: [4, Nz, Nx] float


def get_CONV2D(layer_id):
	return get_ES2D_DT(layer_id)[0].imag / (k_x * get_ES2D(layer_id) + 1e-20) # [4, Nz, Nx]

def get_CONV1D(layer_id):
	es2D = get_ES2D(layer_id) # [4, Nz, Nx)]
	conv2D = get_CONV2D(layer_id) # [4, Nz, Nx)]
	conv1D_x = np.sum(conv2D * es2D, axis=-2) / ( np.sum(es2D, axis=-2) + 1e-20 )
	conv1D_z = np.sum(conv2D * es2D * k_x**2, axis=-1) / ( np.sum(es2D * k_x**2, axis=-1) + 1e-20 )
	return conv1D_x, conv1D_z # [4, Nx], [4, Nz]


def get_ESTS(layer_id):
	file_path_names_1 = [ get_path_name(layer_id, ft, 'probedata/') for ft in ('U','V','W','P') ]
	file_path_names_2 = [ postdata_path+'ESTS/ESTS_'+ft+str(layer_id).zfill(8)+'.BIN' for ft in ('UU','VV','WW','PP') ]
	try:
		esTS_xt = np.zeros([4, Nt, Nx])
		esTS_zt = np.zeros([4, Nt, Nz])
		for t in range(Nt):
			esTS_xz = np.array([ abs( read_channel_layer(fpn, t) ) for fpn in file_path_names_2 ])
			t1, t2 = t, (Nt-t)%Nt # the complemented half of Nx should be positioned at the symetric t, except for t==0
			esTS_xt[:, t1,:Nx/2 ] = np.sum(esTS_xz[:,:,:Nx/2], axis=-2)
			esTS_xt[:, t2, Nx/2:] = np.sum(esTS_xz[:,:,Nx/2:], axis=-2)
			esTS_zt[:, t1, :] += np.sum(esTS_xz[:,:,:Nx/2], axis=-1)
			esTS_zt[:, t2, :] += np.sum(esTS_xz[:,:,Nx/2:], axis=-1)
		return esTS_xt, esTS_zt # [4, Nt, Nx], [4, Nt, Nz]
	except IOError:
		for fpn_1, fpn_2 in zip(file_path_names_1, file_path_names_2):
			time_start, time_end = list( read_channel_infosec(fpn_1)[[0,1]] )
			prblen = int(round( (time_end - time_start) / (dt*nsprb) )) + 1 # total number of time steps probed
			Ns = (prblen - Nt/2) / (Nt/2) # number of samples using 50% overlap window
			prblen = (Ns + 1) * Nt/2 # total number of time steps used in computing spectra
			print 'Computing TS spectra of %s, total length %i, sample number %i, sample length %i ...'%(fpn_1[-13:], prblen, Ns, Nt)
			with open(fpn_2, 'wb') as fp:
				# compute space time spectra from DNS data
				for k in range(Nz):
					print '\tprogress %.2f%% ...'%(100.0*k/Nz)

					q = read_channel_cut(fpn_1, range(prblen), k) # [prblen, Nx/2]
					if k == 0:
						q[:,0] -= np.mean(q[:,0]) # get fluctuation values

					qq = np.zeros([Nt, Nx/2])
					for s in range(Ns):
						winrange = range(s * Nt/2, s * Nt/2 + Nt)
						window = np.hanning(Nt).reshape([Nt,1]) / (np.sum(np.hanning(Nt)**2)/Nt)**0.5
						qq += abs( fft( q[winrange] * window, axis=0 ) )**2
					qq /= Ns
					qq = qq * np.mean(abs(q)**2, axis=0) / ( np.sum(qq, axis=0) + 1e-20 ) # normalization to regain total energy after applying window function

					for t in range(len(qq)):
						recl = Nx/2 * 8
						data = np.ravel( np.transpose( [qq[t], np.zeros(Nx/2)], [1,0] ) )
						fp.seek( recl * ((t+1)*Nz + k) )
						fp.write( pack( recl/4*'f', *data ) )	

		return get_ESTS(layer_id) # [4, Nt, Nx], [4, Nt, Nz]

def get_CORTS(layer_id):
	file_path_names_1 = [ get_path_name(layer_id, ft, 'probedata/') for ft in ('U','V','W','P') ]
	file_path_names_2 = [ postdata_path+'ESTS/CORTS_'+ft+str(layer_id).zfill(8)+'.BIN' for ft in ('UU','VV','WW','PP') ]

	try:
		corTS_xt = np.zeros([4, Nt, Nx])
		corTS_zt = np.zeros([4, Nt, Nz])
		for t in range(Nt/2):
			corTS_xz = np.array([ read_channel_layer(fpn, t) for fpn in file_path_names_2 ])

			corTS_xz[:,:,Nx/2:] = corTS_xz.imag[:,:,:Nx/2][:,:,::-1]
			corTS_xz = corTS_xz.real[:,np.argsort(k_z)][:,:,np.argsort(k_x)]

			t1, t2 = Nt/2+t, Nt/2-t
			corTS_xt[:, t2, :] = corTS_xz[:, Nz/2, -np.arange(Nx)]
			corTS_zt[:, t2, :] = corTS_xz[:, -np.arange(Nz), Nx/2]
			corTS_xt[:, t1, :] = corTS_xz[:, Nz/2, :]
			corTS_zt[:, t1, :] = corTS_xz[:, :, Nx/2]

		corTS_xt = np.array([ normalize(cor) for cor in corTS_xt ])
		corTS_zt = np.array([ normalize(cor) for cor in corTS_zt ])

		return corTS_xt, corTS_zt # [4, Nt, Nx], [4, Nt, Nz]
	except IOError:
		# the following algorithm is verified by comparing the results with (Kim 1993 Pof) and (Kim 1989 JFM)
		for fpn_1, fpn_2 in zip(file_path_names_1, file_path_names_2):
			time_start, time_end = list( read_channel_infosec(fpn_1)[[0,1]] )
			prblen = int(round( (time_end - time_start) / (dt*nsprb) )) + 1 # total number of time steps probed
			Ns = (prblen - Nt/2) / (Nt/2) # number of samples using 50% overlap window
			prblen = (Ns + 1) * Nt/2 # total number of time steps used in computing spectra

			qm = np.mean( read_channel_column(fpn_1, range(prblen), 0, 0) )
			
			Nt_divide = 6
			if ( (Nt/2) % Nt_divide ):
				print '\nNt_divide not valid !\n'
				exit()
			Ntt = (Nt/2) / Nt_divide
			Nss = Ns * Nt_divide

			if Nt_divide >= 1:
				print 'Computing TS correlation of %s, total length %i, sample number %i, sample length %i ...'%(fpn_1[-13:], prblen, Nss, Ntt)
				
				write_channel_infosec(fpn_2, [0])
			
				for tt in range(Nt_divide):
					qq = np.zeros([Ntt,Nz,Nx])

					for ss in range(tt+1):
						q2 = [ read_channel_layer_fluc(fpn_1, Ntt*ss+t, Qm=qm) for t in range(Ntt) ]

						for s in range(Nss):
							print '\tprogress %.2f%% ...'%(100.0*(((1+tt)*tt/2+ss)*Nss+s)/((Nt_divide+1)*Nt_divide/2*Nss))

							q1_range = range( Ntt * (ss+s*(1+tt)),		Ntt * (ss+s*(tt+1)+1) )
							q2_range = range( Ntt * (ss+s*(1+tt)+tt),	Ntt * (ss+(s+1)*(tt+1)) )
							if q1_range[0]/Ntt >= Nss:
								break

							q1 = np.array(q2)

							if tt != 0:
								q2 = [ read_channel_layer_fluc(fpn_1, t, Qm=qm) for t in q2_range ]

							for t in range(Ntt):
								qq[t] += np.sum( phys( np.conj(q1) * q2 ), axis=0 ) / ( Nss * Ntt )
								q2.pop(0)
								q2.append( read_channel_layer_fluc(fpn_1, q2_range[-1]+t, Qm=qm) )

					for t in range(Ntt):
						write_channel_layer(fpn_2, tt*Ntt+t, qq[t]+complex(0,1)*qq[t,:,::-1])

			elif Nt_divide == 1:
				# for Nt_divide = 1, the above codes are equivalent to the following
				print 'Computing TS correlation of %s, total length %i, sample number %i, sample length %i ...'%(fpn_1[-13:], prblen, Ns, Nt)
				
				qq = np.zeros([Nt/2,Nz,Nx])
				q2 = [ read_channel_layer_fluc(fpn_1, t, Qm=qm) for t in range(Nt/2) ]
				for s in range(Ns):
					print '\tprogress %.2f%% ...'%(100.0*s/Ns)
					q1 = np.array(q2)
					for t in range(Nt/2):
						qq[t] += np.sum( phys( np.conj(q1) * q2 ), axis=0 ) / ( Ns * Nt/2 )
						q2.pop(0)
						q2.append( read_channel_layer_fluc(fpn_1, (s+1)*Nt/2+t, Qm=qm) )

				write_channel_infosec(fpn_2, [0])
				for t in range(Nt/2):
					write_channel_layer(fpn_2, t, qq[t]+complex(0,1)*qq[t,:,::-1])


		return get_CORTS(layer_id) # [4, Nt, Nx], [4, Nt, Nz]


def get_ESTS_PARA(layer_id):
	kx, kt = k_x, k_t.reshape([Nt,1])

	esTSs = get_ESTS(layer_id)[0] # [4, Nt, Nx]
	es1Ds = np.sum(esTSs, axis=-2) + 1e-20
	dRdxx = np.sum( -kx**2 * esTSs, axis=(-1,-2) ) + 1e-20
	dRdtt = np.sum( -kt**2 * esTSs, axis=(-1,-2) )
	dRdxt = np.sum( -kx*kt * esTSs, axis=(-1,-2) )

	elip_U = - dRdxt / dRdxx
	elip_V = ( dRdtt / dRdxx - (dRdxt / dRdxx)**2 )**0.5
 
	ridge = -C_ave_plus*u_tau * k_x	# calculate 1st and 2nd order moments around the ridge of the energy spectrum, to avoid the effect of sidelobe due to limited time resolution
	kt = np.array([ k_t-omega for omega in ridge ]).T	# [Nt, Nx]
	kt[kt<min(k_t)] += max(k_t) - min(k_t)
	kt[kt>max(k_t)] -= max(k_t) - min(k_t)

	lamw_omega_c = np.sum( kt * esTSs, axis=-2 ) / es1Ds	# [4, Nx]
	lamw_B = np.sum( [ (kt-omegas)**2 for omegas in lamw_omega_c ] * esTSs, axis=-2 ) / es1Ds
	lamw_omega_c += ridge

	return elip_U, elip_V, lamw_omega_c, lamw_B # [4], [4], [4,Nx], [4,Nx]



## reconstruction: Elliptic model

def compute_ELIP():
	elip_U = np.zeros([4, Ny+1]) # uu, vv, ww, pp
	elip_V = np.zeros([4, Ny+1]) # uu, vv, ww, pp

	for j in range(1, Ny/2+1):
		qq = get_ES2D(j) # [4, Nz, Nx]
		qdqdt, dqdtdqdt = get_ES2D_DT(j)

		dRdxx = np.sum( -k_x**2 * qq, axis=(-1,-2) ) # [4]
		dRdxt = np.sum( -k_x*complex(0,1) * qdqdt, axis=(-1,-2) ).real
		dRdtt = np.sum( -dqdtdqdt, axis=(-1,-2) )

		elip_U[:,j] = elip_U[:,Ny-j] = - dRdxt / dRdxx
		elip_V[:,j] = elip_V[:,Ny-j] = ( dRdtt / dRdxx - (dRdxt / dRdxx)**2 )**0.5
		print 'Elliptic parameters of slice %i are: U+ = %.2f, V+ = %.2f'%(j, np.mean(elip_U[:,j])/u_tau, np.mean(elip_V[:,j])/u_tau)

	return np.array( list(elip_U) + list(elip_V) ) # [8, Ny+1]


def get_ELIP():
	try:
		elips = np.loadtxt(open(postdata_path + 'elip.txt')).T [1:]
		return elips[:4], elips[4:] # elip_U: [4,Ny+1];  elip_V: [4,Ny+1]
	except IOError:
		np.savetxt( postdata_path + 'elip.txt', np.array( [y_mesh] + list(compute_ELIP()) ).T )
		return get_ELIP() # elip_U: [4,Ny+1];  elip_V: [4,Ny+1]
	

def get_CORTS_ELIP(layer_id):
	corTSs = []
	d_x, d_t = list( np.meshgrid(delta_x, delta_t) )
	cor1Ds = phys( get_ES2D(layer_id) ) [:,0,np.argsort(k_x)]
	elip_U, elip_V = list( np.array(get_ELIP())[:,:,layer_id] )

	for elpU, elpV, cor1D in zip( elip_U, elip_V, cor1Ds ):
		delta_xc = ( (d_x - elpU*d_t)**2 + (elpV*d_t)**2 )**0.5
		corTSs.append( np.interp(delta_xc, delta_x, cor1D, left=0.0, right=0.0) )

	return np.array(corTSs) # [4, Nt, Nx]

def get_ESTS_ELIP(layer_id):
	kaiser_beta = 10
	window = np.kaiser(Nx,kaiser_beta) * np.kaiser(Nt,kaiser_beta).reshape([Nt,1])
	window /= ( 1.0/Nx/Nt * np.sum(window**2) )**0.5
	corTSs = get_CORTS_ELIP(layer_id) # [4, Nt, Nx]
	esTSs = abs( spec(corTSs * window) )
	esTSs = esTSs * ( np.max(corTSs, axis=(-1,-2)) / (np.sum(esTSs, axis=(-1,-2)) + 1e-20) ).reshape([4,1,1])
	return esTSs # [4, Nt, Nx]



## reconstruction: LAMW

def compute_LAMW():
	lamw_omega_c = np.zeros([4, Ny+1, Nx])
	lamw_B = np.zeros([4, Ny+1, Nx])

	for j in range(1, Ny/2+1):
		print 'Computing LAMW parameters of slice %i ...'%j
		qq = np.sum( get_ES2D(j), axis=-2 ) # [4, Nx]
		qdqdt, dqdtdqdt = [ np.sum(temp, axis=-2) for temp in get_ES2D_DT(j) ] # [4, Nx], [4, Nx]

		lamw_omega_c[:,j] = lamw_omega_c[:,Ny-j] = - qdqdt.imag / ( qq + 1e-20 )
		lamw_B[:,j] = lamw_B[:,Ny-j] = dqdtdqdt / ( qq + 1e-20 ) - lamw_omega_c[:,j]**2

	return lamw_omega_c, lamw_B # [4, Ny+1,Nx], [4, Ny+1,Nx]

def get_LAMW():
	try:
		file_types = ('UU','VV','WW','PP')
		lamw_omega_c = np.array([ np.loadtxt(open( postdata_path + 'lamw_omega_c_%s.txt'%ft )) for ft in file_types ])
		lamw_B = np.array([ np.loadtxt(open( postdata_path + 'lamw_B_%s.txt'%ft )) for ft in file_types ])
		return lamw_omega_c, lamw_B # [4, Ny+1,Nx], [4, Ny+1,Nx]
	except IOError:
		lamw_omega_c, lamw_B = compute_LAMW()
		np.savetxt(postdata_path + 'lamw_omega_c_UU.txt', lamw_omega_c[0])
		np.savetxt(postdata_path + 'lamw_omega_c_VV.txt', lamw_omega_c[1])
		np.savetxt(postdata_path + 'lamw_omega_c_WW.txt', lamw_omega_c[2])
		np.savetxt(postdata_path + 'lamw_omega_c_PP.txt', lamw_omega_c[3])
		np.savetxt(postdata_path + 'lamw_B_UU.txt', lamw_B[0])
		np.savetxt(postdata_path + 'lamw_B_VV.txt', lamw_B[1])
		np.savetxt(postdata_path + 'lamw_B_WW.txt', lamw_B[2])
		np.savetxt(postdata_path + 'lamw_B_PP.txt', lamw_B[3])
		return get_LAMW() # [4, Ny+1,Nx], [4, Ny+1,Nx]
	

def get_ESTS_LAMW(layer_id):
	esTSs = np.zeros([4, Nt, Nx])
	delta_omega = 2*np.pi/Lt

	for tstep in tsteps:
		print 'Computing time-space energy spectra of slice %i using LAMW at time step %i ...' %(layer_id, tstep)
		for j in (layer_id, Ny-layer_id):
			u, v, w, p = [ ifft( read_channel_layer(get_path_name(tstep,ft), j), axis=-2 ) for ft in ('U','V','W','P') ]
			dudt, dvdt, dwdt, dpdt = [ ifft( read_channel_layer(get_path_name(tstep,ft), j), axis=-2 ) for ft in ('UT','VT','WT','PT') ]

			u[:,0] -= U_mean[j]
			p[:,0] -= P_mean[j]

			for q, dqdt, esTS in zip( (u,v,w,p), (dudt,dvdt,dwdt,dpdt), esTSs ):
				dqdt_q = dqdt / ( q + 1e-20 )

				local_omega = - dqdt_q.imag + dqdt_q.real
				for t in range(Nt):
					delta = 1.0/delta_omega * ( abs(local_omega-k_t[t]) <= delta_omega/2 )
					esTS[t] += np.mean(delta * abs(q)**2, axis=-2)

				local_omega = - dqdt_q.imag - dqdt_q.real
				for t in range(Nt):
					delta = 1.0/delta_omega * ( abs(local_omega-k_t[t]) <= delta_omega/2 )
					esTS[t] += np.mean(delta * abs(q)**2, axis=-2)

	esTSs /= 2 * len(tsteps) * 2 # average between: time steps, upper/lower half of channel, two local frequencies
	esTSs *= 2*np.pi/Lt # continuous space to discrete space

	return esTSs # [4, Nt, Nx]


## 3D energy spectra

def write_ESTS_3D(figure_path):
	for j in range(len(j_probes)):
		fpns = [ postdata_path+'ESTS/ESTS_'+ft+str(j_probes[j]).zfill(8)+'.BIN' for ft in ('UU','VV','WW','PP') ]
		esTSs = np.zeros([4,Nt,Nz,Nx])
		for t in range(Nt):
			esTS_xz = np.array([ abs( read_channel_layer(fpn, t) ) for fpn in fpns ])
			t1, t2 = t, (Nt-t)%Nt # the complemented half of Nx should be positioned at the symetric t, except for t==0
			esTSs[:, t1, :, :Nx/2] = esTS_xz[:,:,:Nx/2]
			esTSs[:, t2, :, Nx/2:] = esTS_xz[:,:,Nx/2:]

		esTSs[:3] /= u_tau**2 * (2*np.pi/(Lx/delta_nu))* (2*np.pi/(Lz/delta_nu)) * (2*np.pi/(Lt/t_nu))
		esTSs[-1] /= tau_w**2 * (2*np.pi/(Lx/delta_nu))* (2*np.pi/(Lz/delta_nu)) * (2*np.pi/(Lt/t_nu))
		esTSs = np.log10( esTSs + 1e-20 )

		file_name = 'ESTS_xzt_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "3D time-space energy spectra at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$k_x^+$", r"$k_z^+$", r"$\omega^+$", r"$E_{uu}^+$", r"$E_{vv}^+$", r"$E_{ww}^+$", r"$E_{pp}^+$") )
			fp.write( 'zone i = %i, j = %i, k = %i\n' %(Nt, Nz, Nx) )
			for i in np.argsort(k_x_plus):
				print '\tESTS 3D writing progress %.2f%% ...'%(100.0*(k_x_plus[i]-min(k_x_plus))/(max(k_x_plus)-min(k_x_plus)))
				for k in np.argsort(k_z_plus):
					for t in np.argsort(k_t_plus):
						for data in [ k_x_plus[i], k_z_plus[k], k_t_plus[t] ] + list( esTSs[:,t,k,i] ):
							fp.write('%.18e\t'%data)
						fp.write('\n')

## for debug

def demo_LAMW():
	import matplotlib.pyplot as plt

	tstep = tsteps[0]
	j = 21
	i = 0
	n = 0

	delta_omega = delta_omega = 2*np.pi/Lt * 2
	u, v, w, p = [ ifft( read_channel_layer(get_path_name(tstep,ft), j), axis=-2 ) for ft in ('U','V','W','P') ]
	dudt, dvdt, dwdt, dpdt = [ ifft( read_channel_layer(get_path_name(tstep,ft), j), axis=-2 ) for ft in ('UT','VT','WT','PT') ]
	
	u[:,0] -= U_mean[j]
	p[:,0] -= P_mean[j]

	q = (u,v,w,p)[n]
	dqdt = (dudt,dvdt,dwdt,dpdt)[n]
	dqdt_q = dqdt / (q + 1e-20)

	local_omega = - dqdt_q.imag + dqdt_q.real
	plt.scatter(range(Nz), local_omega[:,i], s=5)
	local_omega = - dqdt_q.imag - dqdt_q.real
	plt.scatter(range(Nz), local_omega[:,i], s=5)

	plt.plot([[kt+delta_omega/2 for kt in k_t] for cnt in range(Nz)])
	plt.plot([[kt-delta_omega/2 for kt in k_t] for cnt in range(Nz)])
	plt.ylim(min(local_omega[:,i]), max(local_omega[:,i]))
	plt.ylabel(r"$\omega'$")
	plt.title(r"$k_x$ = %.2f"%k_x[i])
	plt.show()


	esTS1 = get_ESTS(j)[0][n]
	esTS2 = get_ESTS_LAMW(j)[n]

	plt.contourf(np.sort(k_x), np.sort(k_t), np.log10(esTS1[np.argsort(k_t)][:,np.argsort(k_x)]), levels=(-8,-7,-6,-5,-4))
	plt.contour (np.sort(k_x), np.sort(k_t), np.log10(esTS2[np.argsort(k_t)][:,np.argsort(k_x)]), levels=(-8,-7,-6,-5,-4))
	plt.xlabel(r"$k_x$")
	plt.ylabel(r"$\omega$")
	plt.show()

	plt.semilogy(np.sort(k_t), esTS1[np.argsort(k_t),i])
	plt.semilogy(np.sort(k_t), esTS2[np.argsort(k_t),i])
	plt.xlabel(r"$\omega$")
	plt.ylabel(r"$E_{qq}$")
	plt.show()

	plt.semilogy(np.sort(k_t), np.sum(esTS1, axis=-1)[np.argsort(k_t)])
	plt.semilogy(np.sort(k_t), np.sum(esTS2, axis=-1)[np.argsort(k_t)])
	plt.xlabel(r"$\omega$")
	plt.ylabel(r"$E_{qq}$")
	plt.show()




# demo_LAMW()




