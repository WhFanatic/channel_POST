#!/work1/cuigx2_work/whn/anaconda_install/anaconda2/bin/python
from post_channel import *
from post_statis import *
from post_spacetime import *


figure_path = postdata_path + 'figures/'


# figure parameters
plot_x = range(1, Nx/2)
plot_z = range(1, Nz/2)
plot_y = range(1, Ny/2+1)


# coordinates
np.savetxt(figure_path+'y_plus.dat', y_plus)
np.savetxt(figure_path+'j_probes.dat', j_probes)
np.savetxt(figure_path+'Re_tau.dat', [Re_tau])
np.savetxt(figure_path+'u_tau.dat', [u_tau])


### basic format ###
# set file_name
# get data to be ploted
# normalize data
# open file
# write file head including title, variables and zone
# write data of every line
####################



file_name = 'means_plot.dat'
with open(figure_path+file_name, 'w') as fp:
	fp.write( 'Title = "mean values of basic variables"\n' )
	fp.write( 'variables = "%s", "%s", "%s", "%s", "%s"\n' %(r"$y^+$", r"$U^+$", r"$V^+$", r"$W^+$", r"$P^+$") )
	fp.write( 'zone i = %i\n' %len(plot_y) )
	for j in plot_y:
		for data in [ y_plus[j], U_mean[j]/u_tau, V_mean[j]/u_tau, W_mean[j]/u_tau, P_mean[j]/tau_w ]:
			fp.write('%.18e\t'%data)
		fp.write('\n')



file_name = 'statis_plot.dat'
RS_uu, RS_vv, RS_ww, RS_uv, \
CP_pp, CP_up, CP_vp, CP_wp, \
VX_xx, VX_yy, VX_zz, VX_xy = list( get_statis() ) # 12 * [Ny+1]
RS_uu /= u_tau**2
RS_vv /= u_tau**2
RS_ww /= u_tau**2
RS_uv /= u_tau**2
CP_pp /= tau_w**2
CP_up /= (u_tau*tau_w)
CP_vp /= (u_tau*tau_w)
CP_wp /= (u_tau*tau_w)
VX_xx /= (u_tau/delta_nu)**2
VX_yy /= (u_tau/delta_nu)**2
VX_zz /= (u_tau/delta_nu)**2
VX_xy /= (u_tau/delta_nu)**2
with open(figure_path+file_name, 'w') as fp:
	fp.write( 'Title = "2nd order statistics of basic variables"\n' )
	fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %(
		r"$y^+$",
		r"$<u'u'>^+$", r"$<v'v'>^+$", r"$<w'w'>^+$", r"$<u'v'>^+$",
		r"$<p'p'>^+$", r"$<u'p'>^+$", r"$<v'p'>^+$", r"$<w'p'>^+$",
		r"$<\omega_x'\omega_x'>^+$", r"$<\omega_y'\omega_y'>^+$", r"$<\omega_z'\omega_z'>^+$", r"$<\omega_x'\omega_y'>^+$" ) )
	fp.write( 'zone i = %i\n' %len(plot_y) )
	for j in plot_y:
		for data in [ y_plus[j],
		RS_uu[j], RS_vv[j], RS_ww[j], RS_uv[j],
		CP_pp[j], CP_up[j], CP_vp[j], CP_wp[j],
		VX_xx[j], VX_yy[j], VX_zz[j], VX_xy[j] ]:
			fp.write('%.18e\t'%data)
		fp.write('\n')



file_name = 'budgets_plot.dat'
bgts = get_BGT()[0] # bgts:[3, 7, Ny+1]
bgts = np.array([ bgts[0], bgts[1], bgts[2], 0.5*np.sum(bgts, axis=0) ]) / (u_tau**3/delta_nu)
with open(figure_path+file_name, 'w') as fp:
	fp.write( 'Title = "budget terms of 3 Reynolds stress & turbulent energy transportation equations"\n' )
	# fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %(
	fp.write( ('variables = "%s"' + 28*', "%s"' + '\n') %(
		r"$y^+$",
		r"$P_{uu}^+$", r"$T_{uu}^{d+}$", r"$T_{uu}^{s+}$", r"$\Pi_{uu}^{d+}$", r"$\Pi_{uu}^{s+}$", r"$\nu_{uu}^{d+}$", r"$\epsilon_{uu}^+$",
		r"$P_{vv}^+$", r"$T_{vv}^{d+}$", r"$T_{vv}^{s+}$", r"$\Pi_{vv}^{d+}$", r"$\Pi_{vv}^{s+}$", r"$\nu_{vv}^{d+}$", r"$\epsilon_{vv}^+$",
		r"$P_{ww}^+$", r"$T_{ww}^{d+}$", r"$T_{ww}^{s+}$", r"$\Pi_{ww}^{d+}$", r"$\Pi_{ww}^{s+}$", r"$\nu_{ww}^{d+}$", r"$\epsilon_{ww}^+$",
		r"$P_{k}^+$" , r"$T_{k}^{d+}$" , r"$T_{k}^{s+}$" , r"$\Pi_{k}^{d+}$" , r"$\Pi_{k}^{s+}$" , r"$\nu_{k}^{d+}$" , r"$\epsilon_{k}^+$" )	)
	fp.write( 'zone i = %i\n' %len(plot_y) )
	for j in plot_y:
		for data in [ y_plus[j] ] + list(np.ravel( bgts[:,:,j] )):
			fp.write('%.18e\t'%data)
		fp.write('\n')



es2Ds = np.array([ get_ES2D(j) for j in j_probes ]) # [len(j_probes, 4, Nz, Nx)]
es2Ds[:,:3] *= k_x/(2*np.pi/Lx) * k_z.reshape([Nz,1])/(2*np.pi/Lz) / u_tau**2
es2Ds[:,-1] *= k_x/(2*np.pi/Lx) * k_z.reshape([Nz,1])/(2*np.pi/Lz) / tau_w**2
for j in range(len(j_probes)):
	file_name = 'ES2D_xz_jprb%i.dat' %j
	with open(figure_path+file_name, 'w') as fp:
		fp.write( 'Title = "pre-multiplied 2D energy spectra at y_plus %.2f"\n' %y_plus[j_probes[j]] )
		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$\lambda_x^+$", r"$\lambda_z^+$", r"$k_xk_zE_{uu}^+$", r"$k_xk_zE_{vv}^+$", r"$k_xk_zE_{ww}^+$", r"$k_xk_zE_{pp}^+$") )
		fp.write( 'zone i = %i, j = %i\n' %(len(plot_z), len(plot_x)) )
		for i in plot_x:
			for k in plot_z:
				for data in [ lambda_x_plus[i], lambda_z_plus[k] ] + list( es2Ds[j,:,k,i] ):
					fp.write('%.18e\t'%data)
				fp.write('\n')



es1Ds_xy, es1Ds_zy = get_ES1D() # [4, Ny+1, Nx], [4, Ny+1, Nz]
es1Ds_xy[:3] *= k_x/(2*np.pi/Lx) / u_tau**2
es1Ds_xy[-1] *= k_x/(2*np.pi/Lx) / tau_w**2
es1Ds_zy[:3] *= k_z/(2*np.pi/Lz) / u_tau**2
es1Ds_zy[-1] *= k_z/(2*np.pi/Lz) / tau_w**2

file_name = 'ES1D_xy.dat'
with open(figure_path+file_name, 'w') as fp:
	fp.write( 'Title = "pre-multiplied 1D energy spectra of X & Y"\n' )
	fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$\lambda_x^+$", r"$y^+$", r"$k_xE_{uu}^+$", r"$k_xE_{vv}^+$", r"$k_xE_{ww}^+$", r"$k_xE_{pp}^+$") )
	fp.write( 'zone i = %i, j = %i\n' %(len(plot_y), len(plot_x)) )
	for i in plot_x:
		for j in plot_y:
			for data in [ lambda_x_plus[i], y_plus[j] ] + list( es1Ds_xy[:,j,i] ):
				fp.write('%.18e\t'%data)
			fp.write('\n')

file_name = 'ES1D_zy.dat'
with open(figure_path+file_name, 'w') as fp:
	fp.write( 'Title = "pre-multiplied 1D energy spectra of Z & Y"\n' )
	fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$\lambda_z^+$", r"$y^+$", r"$k_zE_{uu}^+$", r"$k_zE_{vv}^+$", r"$k_zE_{ww}^+$", r"$k_zE_{pp}^+$") )
	fp.write( 'zone i = %i, j = %i\n' %(len(plot_y), len(plot_z)) )
	for k in plot_z:
		for j in plot_y:
			for data in [ lambda_z_plus[k], y_plus[j] ] + list( es1Ds_zy[:,j,k] ):
				fp.write('%.18e\t'%data)
			fp.write('\n')




## time-space DNS data required

try:

	# space-time spectra & correlation from time-resolved DNS data
	for j in range(len(j_probes)):
		esTSs_xt_zt = get_ESTS(j_probes[j]) # [4, Nt, Nx], [4, Nt, Nz]
		corTSs_xt_zt = get_CORTS(j_probes[j]) # [4, Nt, Nx], [4, Nt, Nz]

		# XT spectra & correlation from time-resolved DNS data
		esTSs = esTSs_xt_zt[0] # [4, Nt, Nx]
		# corTSs = phys(esTSs) # [4, Nt, Nx]
		esTSs[:3] /= u_tau**2 * (2*np.pi/(Lx/delta_nu)) * (2*np.pi/(Lt/t_nu))
		esTSs[-1] /= tau_w**2 * (2*np.pi/(Lx/delta_nu)) * (2*np.pi/(Lt/t_nu))
		esTSs = np.log10( esTSs + 1e-20 )
		# corTSs = np.array([ normalize(corTS)[np.argsort(k_t)][:,np.argsort(k_x)] for corTS in corTSs ])
		corTSs = corTSs_xt_zt[0]


		file_name = 'ESTS_xt_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "time-space energy spectra at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$k_x^+$", r"$\omega^+$", r"$E_{uu}^+$", r"$E_{vv}^+$", r"$E_{ww}^+$", r"$E_{pp}^+$") )
			fp.write( 'zone i = %i, j = %i\n' %(Nt, Nx) )
			for i in np.argsort(k_x_plus):
				for t in np.argsort(k_t_plus):
					for data in [ k_x_plus[i], k_t_plus[t] ] + list( esTSs[:,t,i] ):
						fp.write('%.18e\t'%data)
					fp.write('\n')

		file_name = 'CORTS_xt_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "time-space correlation at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$\Delta x^+$", r"$\Delta t^+$", r"$R_{uu}$", r"$R_{vv}$", r"$R_{ww}$", r"$R_{pp}$") )
			fp.write( 'zone i = %i, j = %i\n' %(Nt, Nx) )
			for i in range(Nx):
				for t in range(Nt):
					for data in [ delta_x_plus[i], delta_t_plus[t] ] + list( corTSs[:,t,i] ):
						fp.write('%.18e\t'%data)
					fp.write('\n')


		# ZT spectra & correlation from time-resolved DNS data
		esTSs = esTSs_xt_zt[1] # [4, Nt, Nz]
		# corTSs = phys(esTSs) # [4, Nt, Nz]
		esTSs[:3] /= u_tau**2 * (2*np.pi/(Lz/delta_nu)) * (2*np.pi/(Lt/t_nu))
		esTSs[-1] /= tau_w**2 * (2*np.pi/(Lz/delta_nu)) * (2*np.pi/(Lt/t_nu))
		esTSs = np.log10( esTSs + 1e-20 )
		# corTSs = np.array([ normalize(corTS)[np.argsort(k_t)][:,np.argsort(k_z)] for corTS in corTSs ])
		corTSs = corTSs_xt_zt[1]

		file_name = 'ESTS_zt_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "time-space energy spectra at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$k_z^+$", r"$\omega^+$", r"$E_{uu}^+$", r"$E_{vv}^+$", r"$E_{ww}^+$", r"$E_{pp}^+$") )
			fp.write( 'zone i = %i, j = %i\n' %(Nt, Nz) )
			for k in np.argsort(k_z_plus):
				for t in np.argsort(k_t_plus):
					for data in [ k_z_plus[k], k_t_plus[t] ] + list( esTSs[:,t,k] ):
						fp.write('%.18e\t'%data)
					fp.write('\n')

		file_name = 'CORTS_zt_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "time-space correlation at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$\Delta z^+$", r"$\Delta t^+$", r"$R_{uu}$", r"$R_{vv}$", r"$R_{ww}$", r"$R_{pp}$") )
			fp.write( 'zone i = %i, j = %i\n' %(Nt, Nz) )
			for k in range(Nz):
				for t in range(Nt):
					for data in [ delta_z_plus[k], delta_t_plus[t] ] + list( corTSs[:,t,k] ):
						fp.write('%.18e\t'%data)
					fp.write('\n')


	# LAMW parameters extracted from time-resolved DNS data
	elip_U, elip_V = [], []
	for j in range(len(j_probes)):
		temp1, temp2, lamw_omega_c, lamw_B = get_ESTS_PARA(j_probes[j]) # [4], [4], [4,Nx], [4,Nx]
		elip_U.append(temp1)
		elip_V.append(temp2)
		lamw_omega_c *= -t_nu / (k_x_plus + 1e-20)
		lamw_B_sqrt = lamw_B**0.5 * t_nu / (k_x_plus + 1e-20)

		file_name = 'ESTS_moments_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "first and second order moments of time-space energy spectra at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %(
				r"$k_x^+$",
				r"$\frac{-\omega_{c,uu}^+}{k_x^+}$", r"$\frac{-\omega_{c,vv}^+}{k_x^+}$", r"$\frac{-\omega_{c,ww}^+}{k_x^+}$", r"$\frac{-\omega_{c,pp}^+}{k_x^+}$",
				r"$\frac{\sqrt{B_{uu}^+}}{k_x^+}$", r"$\frac{\sqrt{B_{vv}^+}}{k_x^+}$", r"$\frac{\sqrt{B_{ww}^+}}{k_x^+}$", r"$\frac{\sqrt{B_{pp}^+}}{k_x^+}$"	) )
			fp.write( 'zone i = %i\n' %len(plot_x) )
			for i in plot_x:
				for data in [ k_x_plus[i] ] + list( lamw_omega_c[:,i] ) + list( lamw_B_sqrt[:,i] ):
					fp.write('%.18e\t'%data)
				fp.write('\n')

	# elliptic parameters extracted from time-resolved DNS data
	file_name = 'ESTS_convswep.dat'
	elip_U = np.array(elip_U).T / u_tau
	elip_V = np.array(elip_V).T / u_tau
	with open(figure_path+file_name, 'w') as fp:
		fp.write( 'Title = "elliptic model parameters from time-space energy spectra"\n' )
		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %(
			r"$y^+$",
			r"$U_{uu}^+$", r"$U_{vv}^+$", r"$U_{ww}^+$", r"$U_{pp}^+$",
			r"$V_{uu}^+$", r"$V_{vv}^+$", r"$V_{ww}^+$", r"$V_{pp}^+$"	)	)
		fp.write( 'zone i = %i\n' %len(j_probes) )
		for j in range(len(j_probes)):
			for data in [ y_plus[j_probes[j]] ] + list(elip_U[:,j]) + list(elip_V[:,j]):
				fp.write('%.18e\t'%data)
			fp.write('\n')

except IOError:
	print 'ESTS data not avaliable!\n'





## time derivative data required

try:

	# convection velocities reconstructed from space data as in Jimenez 2009
	for j in range(len(j_probes)):
		conv2Ds = get_CONV2D(j_probes[j]) # [4, Nz, Nx]
		conv1Ds_x, conv1Ds_z = get_CONV1D(j_probes[j]) # [4, Nx], [4, Nz]
		conv2Ds /= u_tau
		conv1Ds_x /= u_tau
		conv1Ds_z /= u_tau

		file_name = 'CONV2D_xz_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "2D scale dependent convection velocity at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$\lambda_x^+$", r"$\lambda_z^+$", r"$c_{u}^+$", r"$c_{v}^+$", r"$c_{w}^+$", r"$c_{p}^+$") )
			fp.write( 'zone i = %i, j = %i\n' %(len(plot_z), len(plot_x)) )
			for i in plot_x:
				for k in plot_z:
					for data in [ lambda_x_plus[i], lambda_z_plus[k] ] + list( conv2Ds[:,k,i] ):
						fp.write('%.18e\t'%data)
					fp.write('\n')

		file_name = 'CONV1D_x_jprb%i.dat' % j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "x scale dependent convection velocity at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s"\n' %(r"$\lambda_x^+$", r"$c_{u}^+$", r"$c_{v}^+$", r"$c_{w}^+$", r"$c_{p}^+$") )
			fp.write( 'zone i = %i\n' %len(plot_x) )
			for i in plot_x:
				for data in [ lambda_x_plus[i] ] + list( conv1Ds_x[:,i] ):
					fp.write('%.18e\t'%data)
				fp.write('\n')

		file_name = 'CONV1D_z_jprb%i.dat' % j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "z scale dependent convection velocity at y_plus %.2f"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s"\n' %(r"$\lambda_z^+$", r"$c_{u}^+$", r"$c_{v}^+$", r"$c_{w}^+$", r"$c_{p}^+$") )
			fp.write( 'zone i = %i\n' %len(plot_z) )
			for k in plot_z:
				for data in [ lambda_z_plus[k] ] + list( conv1Ds_z[:,k] ):
					fp.write('%.18e\t'%data)
				fp.write('\n')



	# elliptic parameters reconstructed from space data
	file_name = 'elip_plot.dat'
	elip_U, elip_V = get_ELIP() # [4, Ny+1], [4, Ny+1]
	elip_U /= u_tau
	elip_V /= u_tau
	with open(figure_path+file_name, 'w') as fp:
		fp.write( 'Title = "elliptic model parameters for space-time correlation reconstructed from space data"\n' )
		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %(
			r"$y^+$",
			r"$U_{uu}^+$", r"$U_{vv}^+$", r"$U_{ww}^+$", r"$U_{pp}^+$",
			r"$V_{uu}^+$", r"$V_{vv}^+$", r"$V_{ww}^+$", r"$V_{pp}^+$"	)	)
		fp.write( 'zone i = %i\n' %len(plot_y) )
		for j in plot_y:
			for data in [ y_plus[j] ] + list(elip_U[:,j]) + list(elip_V[:,j]):
				fp.write('%.18e\t'%data)
			fp.write('\n')

	# space-time correlation reconstructed from space data using elliptic model
	corTSs = np.array([ get_CORTS_ELIP(j) for j in j_probes ]) # [len(j_probes), 4, Nt, Nx]
	corTSs = np.array([ [ normalize(cor) for cor in corTS ] for corTS in corTSs ])
	for j in range(len(j_probes)):
		file_name = 'CORTS_ELIP_xt_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "time-space correlation at y_plus %.2f reconstructed from space data"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$\Delta x^+$", r"$\Delta t^+$", r"$R_{uu}$", r"$R_{vv}$", r"$R_{ww}$", r"$R_{pp}$") )
			fp.write( 'zone i = %i, j = %i\n' %(Nt, Nx) )
			for i in range(Nx):
				for t in range(Nt):
					for data in [ delta_x_plus[i], delta_t_plus[t] ] + list( corTSs[j,:,t,i] ):
						fp.write('%.18e\t'%data)
					fp.write('\n')



	# LAMW parameters reconstructed from space data
	lamw_omega_c, lamw_B = get_LAMW() # [4, Ny+1, Nx], [4, Ny+1, Nx]
	lamw_omega_c *= -t_nu / (k_x_plus + 1e-20)
	lamw_B_sqrt = lamw_B**0.5 * t_nu / (k_x_plus + 1e-20)
	for j in range(len(j_probes)):
		file_name = 'lamw_plot_jprb%i.dat' %j
		with open(figure_path+file_name, 'w') as fp:
			fp.write( 'Title = "LAMW parameters of time-space energy spectra at y_plus %.2f reconstructed from space data"\n' %y_plus[j_probes[j]] )
			fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %(
				r"$k_x^+$",
				r"$\frac{-\omega_{c,uu}^+}{k_x^+}$", r"$\frac{-\omega_{c,vv}^+}{k_x^+}$", r"$\frac{-\omega_{c,ww}^+}{k_x^+}$", r"$\frac{-\omega_{c,pp}^+}{k_x^+}$",
				r"$\frac{\sqrt{B_{uu}^+}}{k_x^+}$", r"$\frac{\sqrt{B_{vv}^+}}{k_x^+}$", r"$\frac{\sqrt{B_{ww}^+}}{k_x^+}$", r"$\frac{\sqrt{B_{pp}^+}}{k_x^+}$"	) )
			fp.write( 'zone i = %i\n' %len(plot_x) )
			for i in plot_x:
				for data in [ k_x_plus[i] ] + list( lamw_omega_c[:,j_probes[j],i] ) + list( lamw_B_sqrt[:,j_probes[j],i] ):
					fp.write('%.18e\t'%data)
				fp.write('\n')

	# # space-time spectra reconstructed from space data using LAMW
	# esTSs = np.array([ get_ESTS_LAMW(j) for j in j_probes ]) # [len(j_probes), 4, Nt, Nx]
	# esTSs[:,:3] /= u_tau**2 * (2*np.pi/(Lx/delta_nu)) * (2*np.pi/(Lt/t_nu))
	# esTSs[:,-1] /= tau_w**2 * (2*np.pi/(Lx/delta_nu)) * (2*np.pi/(Lt/t_nu))
	# esTSs = np.log10( esTSs + 1e-20 )
	# for j in range(len(j_probes)):
	# 	file_name = 'ESTS_LAMW_xt_jprb%i.dat' %j
	# 	with open(figure_path+file_name, 'w') as fp:
	# 		fp.write( 'Title = "time-space energy spectra at y_plus %.2f"\n' %y_plus[j_probes[j]] )
	# 		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s"\n' %(r"$k_x^+$", r"$\omega^+$", r"$E_{uu}^+$", r"$E_{vv}^+$", r"$E_{ww}^+$", r"$E_{pp}^+$") )
	# 		fp.write( 'zone i = %i, j = %i\n' %(Nt, Nx) )
	# 		for i in np.argsort(k_x_plus):
	# 			for t in np.argsort(k_t_plus):
	# 				for data in [ k_x_plus[i], k_t_plus[t] ] + list( esTSs[j,:,t,i] ):
	# 					fp.write('%.18e\t'%data)
	# 				fp.write('\n')

except IOError:
	print 'Time derivative data not avaliable!'











