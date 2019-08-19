from post_pressure import *
import matplotlib.pyplot as plt


def get_rhs(tstep):
	fpn_rhs = pressure_path + 'RHS' + str(tstep).zfill(8) + '.BIN'
	try:
		rhs = np.array([ read_channel_layer(fpn_rhs, j)[:,:Nx/2] for j in range(Ny+1) ])
		print 'Attention: read RHS from former computation !'
	except IOError:
		rhs = rhs_generate(tstep)
		print 'Writing new RHS ...'
		temp = np.zeros([Nz,Nx], dtype=complex)
		write_channel_infosec(fpn_rhs, [0])
		for j in range(Ny+1):
			temp[:,:Nx/2] = rhs[j]
			write_channel_layer(fpn_rhs, j, temp)
	return rhs

def get_bp(tstep):
	fpn_bp = pressure_path + 'BP' + str(tstep).zfill(8) + '.BIN'
	try:
		bp = np.array([ read_channel_layer(fpn_bp, j)[:,:Nx/2] for j in range(2) ])
		print 'Attention: read BP from former computation !'
	except IOError:
		bp = bp_generate(tstep, True)
		print 'Writing new BP ...'
		temp = np.zeros([Nz,Nx], dtype=complex)
		write_channel_infosec(fpn_bp, [0])
		for j in range(2):
			temp[:,:Nx/2] = bp[j]
			write_channel_layer(fpn_bp, j, temp)
	return bp

def get_poisson(tstep, bp, rhs, file_type):
	file_name = file_type + str(tstep).zfill(8) + '.BIN'
	fpn = pressure_path + file_name
	try:
		read_channel_infosec(fpn)
		print 'Attention: use %s from former computation !'%file_name
	except IOError:
		fpn = poisson(tstep=tstep, bp=bp, rhs=rhs, file_type=file_type)
	return fpn

##############################
##### solve functions, apply to both full & off-wall channel
##############################

##### solve pressure based on velocity at all time steps #####
def solve_all():
	fpns = []
	for tstep in tsteps:
		fpns.append( poisson( tstep, bp_generate(tstep), rhs_generate(tstep) ) )
	return fpns

##### solve several modes respectively #####
def solve_modes():
	tstep = tsteps[0]
	kxids = [0,0,0,0,0,0,0,0,0]
	kzids = [0,0,0,0,0,0,0,0,0]

	rhs = rhs_generate(tstep)
	bp = bp_generate(tstep)
	fpn = poisson(tstep, bp=bp, rhs=rhs)

	p = np.array([ read_channel_layer(fpn, j) for j in range(Ny+1) ])
	dpdy = diff1[:,0].reshape([Ny+1,1,1]) * p[[0]+range(Ny-1)+[Ny-2]] + diff1[:,1].reshape([Ny+1,1,1]) * p[[1]+range(1,Ny)+[Ny-1]] + diff1[:,2].reshape([Ny+1,1,1]) * p[[2]+range(2,Ny+1)+[Ny]]

	file_name = 'rhs_%s.dat' %data_path.strip('/').split('/')[-1][5:]
	with open('data/'+file_name, 'w') as fp:
		fp.write( 'Title = "abs( rhs )^+"\n' )
		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %tuple( [ r"$y^+$" ] + [ r"$\lambda_x^+=%i, \lambda_z^+=%i$"%(lambda_x_plus[i],lambda_z_plus[k]) for i,k in zip(kxids,kzids) ] ) )
		fp.write( 'zone i = %i\n' %(Ny+1) )
		for j in range(Ny+1):
			for data in [ y_plus[j] ] + [ abs(rhs[j,k,i])/(tau_w/delta_nu**2) for i,k in zip(kxids,kzids) ]:
				fp.write('%.18e\t'%data)
			fp.write('\n')

	file_name = 'bp_%s.dat' %data_path.strip('/').split('/')[-1][5:]
	with open('data/'+file_name, 'w') as fp:
		fp.write( 'Title = "abs( bp )^+"\n' )
		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %tuple( [ r"$y^+$" ] + [ r"$\lambda_x^+=%i, \lambda_z^+=%i$"%(lambda_x_plus[i],lambda_z_plus[k]) for i,k in zip(kxids,kzids) ] ) )
		fp.write( 'zone i = 2\n' )
		for j in (0,-1):
			for data in [ y_plus[j] ] + [ abs(bp[j,k,i])/(tau_w/delta_nu) for i,k in zip(kxids,kzids) ]:
				fp.write('%.18e\t'%data)
			fp.write('\n')

	file_name = 'p_%s.dat' %data_path.strip('/').split('/')[-1][5:]
	with open('data/'+file_name, 'w') as fp:
		fp.write( 'Title = "abs( p )^+"\n' )
		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %tuple( [ r"$y^+$" ] + [ r"$\lambda_x^+=%i, \lambda_z^+=%i$"%(lambda_x_plus[i],lambda_z_plus[k]) for i,k in zip(kxids,kzids) ] ) )
		fp.write( 'zone i = %i\n' %(Ny+1) )
		for j in range(Ny+1):
			for data in [ y_plus[j] ] + [ abs(p[j,k,i])/(tau_w) for i,k in zip(kxids,kzids) ]:
				fp.write('%.18e\t'%data)
			fp.write('\n')

	file_name = 'dpdy_%s.dat' %data_path.strip('/').split('/')[-1][5:]
	with open('data/'+file_name, 'w') as fp:
		fp.write( 'Title = "abs( dpdy )^+"\n' )
		fp.write( 'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' %tuple( [ r"$y^+$" ] + [ r"$\lambda_x^+=%i, \lambda_z^+=%i$"%(lambda_x_plus[i],lambda_z_plus[k]) for i,k in zip(kxids,kzids) ] ) )
		fp.write( 'zone i = %i\n' %(Ny+1) )
		for j in range(Ny+1):
			for data in [ y_plus[j] ] + [ abs(dpdy[j,k,i])/(tau_w/delta_nu) for i,k in zip(kxids,kzids) ]:
				fp.write('%.18e\t'%data)
			fp.write('\n')


##############################
##### decomp functions, apply to full channel only
##############################

##### scale analysis #####
def decomp_scales():
	mfu_scale = np.ones([Nz,Nx], dtype=bool)
	mfu_scale[abs(lambda_z_plus) > 100*np.pi] = False
	mfu_scale[:,abs(lambda_x_plus) > 1000*np.pi] = False # should streamwise average be retained?
	big_scale = ( mfu_scale != True )

	pp1 = [ np.zeros(Ny+1) for n in range(3) ]
	pp2= [ np.zeros(Ny+1) for n in range(3) ]

	# for tstep in tsteps:
	# 	bp = np.zeros([2,Nz,Nx/2], dtype=complex)
	# 	rhs1 = np.array([ -np.sum(rhs_generate_decomp(tstep, j, mfu_scale),axis=0) for j in range(Ny+1) ])
	# 	rhs2 = np.array([ -np.sum(rhs_generate_decomp(tstep, j, big_scale),axis=0) for j in range(Ny+1) ])

	# 	fpn = get_path_name(tstep, 'P')
	# 	fpn1 = poisson( tstep, bp, rhs1 )
	# 	fpn2 = poisson( tstep, bp, rhs2 )

	# 	for j in range(Ny+1):
	# 		p = read_channel_layer_fluc( , j, Qm='default0' )
	# 		pp1[0] += np.sum( abs(p)**2 )
	# 		pp1[1] += np.sum( abs(p[mfu_scale])**2 )
	# 		pp1[2] += np.sum( abs(p[big_scale])**2 )

	# pp = ( pp + pp[::-1] ) / 2 / len(tsteps)
	# pp1 = ( pp1 + pp1[::-1] ) / 2 / len(tsteps)
	# pp2 = ( pp2 + pp2[::-1] ) / 2 / len(tsteps)

	# plt.plot(y_plus[:Ny/2+1], pp[:Ny/2+1])
	# plt.plot(y_plus[:Ny/2+1], pp1[:Ny/2+1])
	# plt.plot(y_plus[:Ny/2+1], pp2[:Ny/2+1])
	# plt.show()


##### rapid & slow term decomposition, for statistics #####
def decomp_rs_stat():
	# pp1,pp2,pp3,pp4,pp5,pp6 = [np.zeros(Ny+1) for n in range(6)]

	# for tstep in tsteps:
	# 	fpn1, fpn2, fpn4, fpn3 = pressure_decompose(tstep)
	# 	fpn4 = get_path_name(tstep, 'P')

	# 	for j in range(Ny+1):
	# 		p1,p2,p3,p4 = [read_channel_layer(fpn, j) for fpn in (fpn1,fpn2,fpn3,fpn4)]
	# 		p1[0,0] = p2[0,0] = p3[0,0] = p4[0,0] = 0

	# 		pp1[j] += np.sum(abs(p1)**2)
	# 		pp2[j] += np.sum(abs(p2)**2)
	# 		pp3[j] += np.sum(abs(p3)**2)
	# 		pp4[j] += np.sum(abs(p4)**2)
	# 		pp5[j] += np.sum(abs(p4-p3)**2)
	# 		pp6[j] += 2 * np.sum(np.conj(p4)*p3).real

	# for pp in [pp1,pp2,pp3,pp4,pp5,pp6]:
	# 	pp[:] = (pp + pp[::-1]) / 2 / len(tsteps)

	# np.savetxt('data/decomp_rs_stat.txt', np.array([y_plus,pp1,pp2,pp3,pp4,pp5,pp6]).T)
	pp1,pp2,pp3,pp4,pp5,pp6 = np.loadtxt('data/decomp_rs_stat.txt').T[1:]

	fig_name = 'P_rs_decomposition.png'
	plt.figure(num=fig_name, figsize=(4,4))
	fig, ax = plt.gcf(), plt.gca()

	ax.plot(y_plus, pp1/tau_w**2, label=r"$<p'_r^2>^+$")
	ax.plot(y_plus, pp2/tau_w**2, label=r"$<p'_s^2>^+$")
	ax.plot(y_plus, pp3/tau_w**2, label=r"$<p'_{st}^2>^+$")
	ax.plot(y_plus, pp4/tau_w**2, label=r"$<p'^2>^+$")
	ax.plot(y_plus, pp5/tau_w**2, label=r"$<(p'-p'_{st})^2>^+$")
	ax.plot(y_plus, pp6/tau_w**2, label=r"$2<p'p'_{st}>^+$")

	ax.legend(loc='best')
	ax.set_xlim([0,Re_tau])
	ax.set_ylim([0,9])
	ax.set_xlabel(r"$y^+$")

	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()

##### rapid & slow term decomposition, for instantaneous field #####
def decomp_rs_inst(tstep = tsteps[0]):
	fpn1, fpn2, fpn4, fpn3 = pressure_decompose(tstep)
	fpn4 = get_path_name(tstep, 'P')
	p1, p2, p3, p4 = [ read_channel_fluc_easy(fpn,z=[0]) for fpn in (fpn1, fpn2, fpn3, fpn4) ]

	fig_name = 'P_rs_decomp_field.png'
	fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, squeeze=False, num=fig_name, figsize=(6,4))

	ax = axs[0,0]
	ax.contour(x_mesh, y_mesh, p1/tau_w, levels=range(-4,5))
	ax.set_title(r"$p'_r$")

	ax = axs[0,1]
	ax.contour(x_mesh, y_mesh, p2/tau_w, levels=range(-4,5))
	ax.set_title(r"$p'_s$")

	ax = axs[1,0]
	ax.contour(x_mesh, y_mesh, p4/tau_w, levels=range(-4,5), colors='black', linestyles='--')
	ax.contour(x_mesh, y_mesh, (p1+p2)/tau_w, levels=range(-4,5))
	ax.set_title(r"($p'_r+p'_s$) versus $p'$")

	ax = axs[1,1]
	ax.contour(x_mesh, y_mesh, (p4-p3)/tau_w, levels=range(-4,5), colors='black', linestyles='--')
	ax.contour(x_mesh, y_mesh, (p1+p2)/tau_w, levels=range(-4,5))
	ax.set_title(r"($p'_r+p'_s$) versus ($p'-p'_{st}$)")

	for ax in axs[:,0]:
		ax.set_ylabel('y')
	for ax in axs[-1,:]:
		ax.set_xlabel('x')

	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()


##### near-wall & far-from-wall decomposition, for statistics #####
def decomp_nfw_stat(y_indexes = [10,17,24,31]):
	intersects_filename = 'intersects'+data_path.strip('/').split('/')[-1][4:]+'.txt'
	try:
		intersects_record, yids, intersects = False, y_indexes, np.loadtxt('data/'+intersects_filename)
	except IOError:
		intersects_record, yids, intersects = True, range(5,Ny/2,10), []
		for j in y_indexes:
			if not (j in yids):
				yids.append(j)
		yids.sort()

	fig_name = 'P_nfw_decomposition.png'
	fig, axs = plt.subplots(2, len(y_indexes), sharex='col', sharey=False, squeeze=False, num=fig_name, figsize=(3*len(y_indexes),6))

	for jprb0 in yids:
		print 'jprb0 = ',jprb0

		# jprb1 = Ny - jprb0
		# H = np.zeros([Ny+1,1,1])
		# H[jprb0:jprb1+1] = 1

		# pp1,pp2,pp3,pp4,pp5,pp14,pp24,pp12 = [np.zeros(Ny+1) for n in range(8)]

		# for tstep in tsteps:
		# 	rhs, bp = get_rhs(tstep), get_bp(tstep) # Attention: from former computation !!!
		# 	fpn1 = poisson(tstep, bp=bp*0, rhs=rhs*H, file_type='PFW')
		# 	fpn2 = poisson(tstep, bp=bp*0, rhs=rhs*(1-H), file_type='PNW')
		# 	fpn3 = get_poisson(tstep, bp=bp, rhs=rhs*0, file_type='PST') # Attention: from former computation !!!
		# 	fpn4 = get_path_name(tstep, 'P')

		# 	for j in range(Ny+1):
		# 		p1,p2,p3,p4 = [read_channel_layer(fpn, j) for fpn in (fpn1,fpn2,fpn3,fpn4)]

		# 		p1[0,0] = p2[0,0] = p3[0,0] = p4[0,0] = 0

		# 		pp1[j] += np.sum(abs(p1)**2)
		# 		pp2[j] += np.sum(abs(p2)**2)
		# 		pp3[j] += np.sum(abs(p1+p2)**2)
		# 		pp4[j] += np.sum(abs(p4-p3)**2)
		# 		pp5[j] += np.sum(abs(p4)**2)

		# 		pp14[j] += np.sum(np.conj(p1)*p4).real
		# 		pp24[j] += np.sum(np.conj(p2)*p4).real
		# 		pp12[j] += np.sum(np.conj(p1)*p2).real

		# for pp in [pp1,pp2,pp3,pp4,pp5]:
		# 	pp[:] = (pp + pp[::-1]) / 2 / len(tsteps)
		# pp14 = (pp14 + pp14[::-1]) / 2 / len(tsteps) / (pp1 * pp5)**0.5
		# pp24 = (pp24 + pp24[::-1]) / 2 / len(tsteps) / (pp2 * pp5)**0.5
		# pp12 = (pp12 + pp12[::-1]) / 2 / len(tsteps) / (pp1 * pp2)**0.5

		# if intersects_record:
		# 	intersects.append(
		# 		[ y_plus[jprb0], pp1[jprb0], pp2[jprb0], pp3[jprb0], pp14[jprb0], pp24[jprb0], pp12[jprb0], pp1[0], pp2[0], pp3[0], pp14[0], pp24[0], pp12[0] ]
		# 	)
		
		if jprb0 in y_indexes:
			# np.savetxt('data/decomp_nfw_stat_jprb%i.txt'%jprb0, np.array([y_plus,pp1,pp2,pp3,pp4,pp5,pp14,pp24,pp12]).T)
			pp1,pp2,pp3,pp4,pp5,pp14,pp24,pp12 = np.loadtxt('data/decomp_nfw_stat_jprb%i.txt'%jprb0).T[1:]

			ax = axs[0,y_indexes.index(jprb0)]
			ax.semilogx(y_plus, pp1/tau_w**2, label=r"$<p'_{fw}^2>^+$")
			ax.semilogx(y_plus, pp2/tau_w**2, label=r"$<p'_{nw}^2>^+$")
			ax.semilogx(y_plus, pp3/tau_w**2, label=r"$<(p'_{fw}+p'_{nw})^2>^+$")
			ax.semilogx(y_plus, pp4/tau_w**2, label=r"$<(p'-p'_{st})^2>^+$", linestyle='--')
			if ax == axs[0,0]:
				ax.legend(loc='lower right', fontsize=6.5)
			ax.set_xlim([1,Re_tau])
			# ax.set_ylim([0,20])
			ax.set_title(r"$y_b^+$ = %i" %int(np.round(y_plus[jprb0])))
			ax.semilogx([y_plus[jprb0],y_plus[jprb0]], ax.get_ylim(), color='black', linestyle='-.')

			ax = axs[1,y_indexes.index(jprb0)]
			ax.semilogx(y_plus, pp14, label=r"$Cor(p'_{fw},p')$")
			ax.semilogx(y_plus, pp24, label=r"$Cor(p'_{nw},p')$")
			ax.semilogx(y_plus, pp12, label=r"$Cor(p'_{fw},p'_{nw})$")
			if ax == axs[1,0]:
				ax.legend(loc='lower right', fontsize=6.5)
			ax.set_xlim([1,Re_tau])
			ax.set_ylim([-1,1])
			ax.set_xlabel(r"$y^+$")
			ax.semilogx([y_plus[jprb0],y_plus[jprb0]], ax.get_ylim(), ax.get_xlim(), [0,0], color='black', linestyle='-.')

	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()

	intersects = np.array(intersects)
	if intersects_record:
		np.savetxt('data/'+intersects_filename, intersects)


	fig_name = 'P_nfw_decomp_varying.png'
	fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', squeeze=False, num=fig_name, figsize=(6,6))

	for n in range(2):
		ypls = intersects.T[0]
		pp1,pp2,pp3,pp14,pp24,pp12 = intersects.T[1+6*n:7+6*n]
		
		ax = axs[0,n]
		ax.plot(ypls, pp1/tau_w**2, label=r"$<p'_{fw}^2>^+$")
		ax.plot(ypls, pp2/tau_w**2, label=r"$<p'_{nw}^2>^+$")
		ax.plot(ypls, pp3/tau_w**2, label=r"$<(p'_{fw}+p'_{nw})^2>^+$")
		ax.set_xlim([-1,Re_tau])
		if n == 0:
			ax.legend(loc='best')
			ax.set_title(r"intercept at $y = y_b$")
		else:
			ax.set_title(r"intercept at $y = 0$")

		ax = axs[1,n]
		ax.plot(ypls, pp14, label=r"$Cor(p'_{fw},p')$")
		ax.plot(ypls, pp24, label=r"$Cor(p'_{nw},p')$")
		ax.plot(ypls, pp12, label=r"$Cor(p'_{fw},p'_{nw})$")
		if n == 0:
			ax.legend(loc='best')
		ax.set_xlim([-10,Re_tau])
		ax.set_ylim([-1,1])
		ax.set_xlabel(r"$y_b^+$")
		ax.plot(ax.get_xlim(), [0,0], color='black', linestyle='-.')

	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()


##### near-wall & far-from-wall decomposition, for instantaneous fields #####
def decomp_nfw_inst(tstep = tsteps[0], y_indexes = [10,17,24,31]):
	fig_name = 'P_fnw_decomp_field.png'
	fig, axs = plt.subplots(3, len(y_indexes), sharex=True, sharey=True, squeeze=False, num=fig_name, figsize=(3*len(y_indexes),6))

	rhs = rhs_generate(tstep)
	bp = bp_generate(tstep, True)

	for jprb0 in y_indexes:

		jprb1 = Ny - jprb0
		H = np.zeros([Ny+1,1,1])
		H[jprb0:jprb1+1] = 1
			
		fpn1 = poisson(tstep, bp=bp*0, rhs=rhs*H, file_type='PFW')
		fpn2 = poisson(tstep, bp=bp*0, rhs=rhs*(1-H), file_type='PNW')
		# fpn3 = poisson(tstep, bp=bp, rhs=rhs*0, file_type='PST')
		fpn4 = get_path_name(tstep, 'P')

		p1 = read_channel_fluc_easy(fpn1, z=[0])
		p2 = read_channel_fluc_easy(fpn2, z=[0])
		# p3 = read_channel_fluc_easy(fpn3, z=[0])
		p4 = read_channel_fluc_easy(fpn4, z=[0])

		ax = axs[0,y_indexes.index(jprb0)]
		ax.contour(x_mesh, y_mesh, p1/tau_w, levels=np.arange(-4,5))
		ax.set_title(r"$y_b^+$ = %i"%int(np.round(y_plus[jprb0])) + "\n" + r"$p'_{fw}$")

		ax = axs[1,y_indexes.index(jprb0)]
		ax.contour(x_mesh, y_mesh, p2/tau_w, levels=np.arange(-4,5))
		ax.set_title(r"$p'_{nw}$")

		ax = axs[2,y_indexes.index(jprb0)]
		ax.contour(x_mesh, y_mesh, p4/tau_w, levels=np.arange(-4,5), colors='black', linestyles='--')
		ax.contour(x_mesh, y_mesh, (p1+p2)/tau_w, levels=np.arange(-4,5))
		ax.set_title(r"($p'_{fw}+p'_{nw}$) versus $p'$")

	for ax in axs[:,0]:
		ax.set_ylabel('y')
	for ax in axs[-1,:]:
		ax.set_xlabel('x')

	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()


##############################
##### demo functions, apply to full channel only
##############################

##### plot green function #####
def demo_green_function():
	kxids = [1,2,10,63]
	kzids = [0,2,10,63]
	yids = [0,10,17,24,31,64]

	fig_name = 'Green_function.png'
	fig, axs = plt.subplots(2, len(kxids), sharex='col', sharey=False, squeeze=False, num=fig_name, figsize=(3*len(kxids),6))

	for n in range(len(kxids)):
		kx, kz = k_x[kxids[n]], k_z[kzids[n]]
		G, dGdy = green_function((kx**2+kz**2)**0.5)

		ax = axs[0,n]
		ax.contourf(y_mesh, y_mesh, -G, cmap=plt.cm.rainbow, extend='both')
		ax.plot(y_mesh, np.ones([Ny+1,len(yids)])*y_mesh[yids], lw=1, ls='--')
		ax.set_title(r"k = %i, $\lambda^+$ = %i" %((kx**2+kz**2)**0.5, 2*np.pi/(kx**2+kz**2)**0.5/delta_nu))

		ax = axs[1,n]
		ax.plot(y_mesh, -G[yids,:].T)
		ax.set_xlabel(r"$\eta$")

	axs[0,0].set_ylabel('y')
	axs[1,0].set_ylabel('-G')
	axs[1,-1].legend([ r"$y^+$ = %i"%y_plus[j] for j in yids ])

	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()


	fig_name = 'Green_function_dy.png'
	fig, axs = plt.subplots(2, len(kxids), sharex='col', sharey=False, squeeze=False, num=fig_name, figsize=(3*len(kxids),6))

	for n in range(len(kxids)):
		kx, kz = k_x[kxids[n]], k_z[kzids[n]]
		G, dGdy = green_function((kx**2+kz**2)**0.5)

		ax = axs[0,n]
		ax.contourf(y_mesh, y_mesh, dGdy, cmap=plt.cm.rainbow, extend='both')
		ax.plot(y_mesh, np.ones([Ny+1,len(yids)])*y_mesh[yids], lw=1, ls='--')
		ax.set_title(r"k = %i, $\lambda^+$ = %i" %((kx**2+kz**2)**0.5, 2*np.pi/(kx**2+kz**2)**0.5/delta_nu))

		ax = axs[1,n]
		ax.plot(y_mesh, dGdy[yids,:].T)
		ax.set_xlabel(r"$\eta$")

	axs[0,0].set_ylabel('y')
	axs[1,0].set_ylabel(r"$\partial G / \partial y$")
	axs[1,-1].legend([ r"$y^+$ = %i"%y_plus[j] for j in yids ])

	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()

##### RHS test #####
def demo_rhs():
	rhs = np.zeros([Ny+1,3])
	for tstep in tsteps:
		print tstep
		temp1 = -rhs_generate(tstep)[:,0,0].real
		temp2 = rhs2_generate(tstep)[:,0,0].real
		temp3 = temp1 - temp2
		rhs += np.array([temp1,temp2,temp3]).T

	rhs = (rhs + rhs[::-1,:]) / 2 / len(tsteps)

	np.savetxt('rhs.txt', np.array([y_mesh]+list(rhs.T)).T)

	plt.plot(y_mesh, rhs)
	plt.legend(['total', '3 terms', '6 terms'])
	plt.ylim([-0.08,0.08])
	plt.show()


##############################
##### test functions, apply to full channel only
##############################

##### error analysis of post solved pressure #####
def test_postcodes():
	try:
		p_statis = np.loadtxt( 'data/p_statis.txt' ).T[1:]
		p1_statis = np.loadtxt( 'data/p1_statis.txt' ).T[1:]
		dp_statis = np.loadtxt( 'data/dp_statis.txt' ).T[1:]
		uvw_rms = np.loadtxt( 'data/uvw_rms.txt' ).T[1:]
	except IOError:
		p_statis = np.zeros([6,Ny+1])
		p1_statis = np.zeros([6,Ny+1])
		dp_statis = np.zeros([6,Ny+1])
		uvw_rms = np.zeros([3,Ny+1])

		for tstep in tsteps:
			print tstep
			fpn1 = poisson( tstep, bp_generate(tstep, True), rhs_generate(tstep) )
			fpn2,fpn3,fpn4,fpn5 = [ get_path_name(tstep, ft) for ft in ('P','U','V','W') ]

			for j in range(Ny+1):
				p1,p,u,v,w = [ read_channel_layer(fpn, j) for fpn in (fpn1,fpn2,fpn3,fpn4,fpn5) ]
				p1[0,0] = p[0,0] = u[0,0] = v[0,0] = w[0,0] = 0
				dp = p - p1

				p_statis[0,j] += np.sum( abs(p)**2 )
				p_statis[1,j] += np.sum( np.conj(u) * p ).real
				p_statis[2,j] += np.sum( np.conj(v) * p ).real
				p_statis[3,j] += np.sum( np.conj(w) * p ).real
				p_statis[4,j] += convol(convol(p),p)[0,0].real
				p_statis[5,j] += convol(convol(p))[0,0].real

				p1_statis[0,j] += np.sum( abs(p1)**2 )
				p1_statis[1,j] += np.sum( np.conj(u) * p1 ).real
				p1_statis[2,j] += np.sum( np.conj(v) * p1 ).real
				p1_statis[3,j] += np.sum( np.conj(w) * p1 ).real
				p1_statis[4,j] += convol(convol(p1),p1)[0,0].real
				p1_statis[5,j] += convol(convol(p1))[0,0].real

				dp_statis[0,j] += np.sum( abs(dp)**2 )
				dp_statis[1,j] += np.sum( np.conj(u) * dp ).real
				dp_statis[2,j] += np.sum( np.conj(v) * dp ).real
				dp_statis[3,j] += np.sum( np.conj(w) * dp ).real
				dp_statis[4,j] += np.sum( np.conj(p) * dp ).real
				dp_statis[5,j] += np.sum( np.conj(p) * p1 ).real

				uvw_rms[0,j] += np.sum( abs(u)**2 )
				uvw_rms[1,j] += np.sum( abs(v)**2 )
				uvw_rms[2,j] += np.sum( abs(w)**2 )

		p_statis /= len(tsteps)
		p1_statis /= len(tsteps)
		dp_statis /= len(tsteps)
		uvw_rms /= len(tsteps)

		np.savetxt( 'data/p_statis.txt', np.array([y_mesh]+list(p_statis)).T )
		np.savetxt( 'data/p1_statis.txt', np.array([y_mesh]+list(p1_statis)).T )
		np.savetxt( 'data/dp_statis.txt', np.array([y_mesh]+list(dp_statis)).T )
		np.savetxt( 'data/uvw_rms.txt', np.array([y_mesh]+list(uvw_rms)).T )
	# got p_statis, p1_statis, dp_statis, uvw_rms

	ypls = y_plus[:Ny/2+1]
	p_statis = ( p_statis + (p_statis.T[::-1] * [1,1,-1,1,1,1]).T )[:,:Ny/2+1] / 2
	p1_statis = ( p1_statis + (p1_statis.T[::-1] * [1,1,-1,1,1,1]).T )[:,:Ny/2+1] / 2
	dp_statis = ( dp_statis + (dp_statis.T[::-1] * [1,1,-1,1,1,1]).T )[:,:Ny/2+1] / 2
	uvw_rms = ( uvw_rms + uvw_rms[:,::-1] )[:,:Ny/2+1] / 2

	p_statis_plus = ( p_statis.T /[tau_w**2, u_tau*tau_w, u_tau*tau_w, u_tau*tau_w, tau_w**3, tau_w**4] ).T
	p1_statis_plus= ( p1_statis.T/[tau_w**2, u_tau*tau_w, u_tau*tau_w, u_tau*tau_w, tau_w**3, tau_w**4] ).T
	dp_statis_plus= ( dp_statis.T/[tau_w**2, u_tau*tau_w, u_tau*tau_w, u_tau*tau_w, tau_w**2, tau_w**2] ).T
	uvw_rms_plus = uvw_rms / u_tau**2

	fig_name = 'p_statis.png'
	fig, axs = plt.subplots(2, 3, sharex=True, num=fig_name, figsize=(10,6))
	for row in range(2):
		for col in range(3):
			ax = axs[row,col]
			ax.plot( ypls, p_statis_plus[(0,4,5,1,2,3)[3*row+col]] )
			ax.plot( ypls, p1_statis_plus[(0,4,5,1,2,3)[3*row+col]] )
			ax.legend([
				[r"$<p'^2>^+$", r"$<p'_{post}^2>^+$"],
				[r"$<p'^3>^+$", r"$<p'_{post}^3>^+$"],
				[r"$<p'^4>^+$", r"$<p'_{post}^4>^+$"],
				[r"$<u'p'>^+$", r"$<u'p'_{post}>^+$"],
				[r"$<v'p'>^+$", r"$<v'p'_{post}>^+$"],
				[r"$<w'p'>^+$", r"$<w'p'_{post}>^+$"]	][3*row+col], fontsize=8, loc='best')
	for ax in axs[1,:]:
		ax.set_xlabel(r"$y^+$")
	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()

	fig_name = 'p_statis_error.png'
	fig, ax = plt.figure(num=fig_name, figsize=(5,5)), plt.gca()
	ax.plot(ypls, dp_statis[5]/(p_statis[0]*p1_statis[0])**0.5)
	ax.plot(ypls, (p_statis/p1_statis)[[1,2,3,0,4,5]].T)
	ax.legend([r"$Cor(p',p'_{post})$", r"$<u'p'_{post}>/<u'p'>$", r"$<v'p'_{post}>/<v'p'>$", r"$<w'p'_{post}>/<w'p'>$", r"$<p'_{post}^2>/<p'^2>$", r"$<p'_{post}^3>/<p'^3>$", r"$<p'_{post}^4>/<p'^4>$"])
	ax.set_ylim([0.85,1.25])
	ax.set_xlabel(r"$y^+$")
	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()

	fig_name = 'p_cor.png'
	fig, axs = plt.subplots(1, 2, sharex=True, num=fig_name, figsize=(8,4))
	ax = axs[0]
	ax.plot(ypls, (dp_statis[0]/p_statis[0])**0.5)
	ax.plot(ypls, p_statis[1]/(p_statis[0]*uvw_rms[0])**0.5)
	ax.plot(ypls, p_statis[2]/(p_statis[0]*uvw_rms[1])**0.5)
	ax.plot(ypls, p_statis[3]/(p_statis[0]*uvw_rms[2])**0.5)
	ax.legend([r"$dp_{rms}/p_{rms}$", r"$Cor(u',p')$", r"$Cor(v',p')$", r"$Cor(w',p')$"])
	ax = axs[1]
	ax.plot(ypls, dp_statis[4]/(p_statis[0]*dp_statis[0])**0.5)
	ax.plot(ypls, dp_statis[1]/(uvw_rms[0]*dp_statis[0])**0.5)
	ax.plot(ypls, dp_statis[2]/(uvw_rms[1]*dp_statis[0])**0.5)
	ax.plot(ypls, dp_statis[3]/(uvw_rms[2]*dp_statis[0])**0.5)
	ax.legend([r"$Cor(p',dp')$", r"$Cor(u',dp')$", r"$Cor(v',dp')$", r"$Cor(w',dp')$"])
	for ax in axs:
		ax.plot(ypls, np.zeros(len(ypls)), color='black', ls='--')
		ax.set_xlabel(r"$y^+$")
	fig.tight_layout()
	fig.savefig(fig_name, dpi=250)
	plt.close()





decomp_nfw_stat(y_indexes=[10,18,25,47])
# decomp_rs_stat()