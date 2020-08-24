#!/work1/cuigx2_work/whn/anaconda_install/anaconda3/bin/python
from basic import *
from statis import Statis
from spacetime import SpaceTime


workpath = 'postdata_M4000/'

para = DataSetInfo("/back1/cuigx2_back1/whn/data/DNS4000M/")
para.inner_scale()

stas = Statis(para)
sptm = SpaceTime(para)


################################################
stas.calc_statis(tsteps=para.tsteps[:])
stas.write_es2d(workpath + 'ES2D/')

sptm.calc_es2d_dt(tsteps=para.tsteps[:])
sptm.write_ES2D_DT(workpath + 'ESDT/')

for j in para.j_probes[::-1]:
	print('for y+ = %.2f:'%(para.ys[j]/para.lc))

	sptm.calc_corts_elip(j)
	sptm.write_CORTS(workpath + 'RECONS/', j)

	sptm.calc_ests_lamw(j, tsteps=para.tsteps[:])
	sptm.write_ESTS(workpath + 'RECONS/', j)

	sptm.calc_ests(j)
	sptm.write_ESTS(workpath + 'ESTS/', j)

	sptm.calc_corts(j)
	sptm.write_CORTS(workpath + 'ESTS/', j)
################################################
sptm.read_ES2D_DT(workpath + 'ESDT/')
################################################




casename = para.datapath.split('/')[-2]
jrange = range(1, para.Ny//2+1) #range(1, para.Ny) #
krange = range(1, para.Nzc)
irange = range(1, para.Nxc)

para.inner_scale()







header = \
	'Title = "elliptic model parameters from time-space energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	%(	r"$y^+$",
		r"$U_{uu}^+$", r"$U_{vv}^+$", r"$U_{ww}^+$", r"$U_{pp}^+$",
		r"$V_{uu}^+$", r"$V_{vv}^+$", r"$V_{ww}^+$", r"$V_{pp}^+$"	) + \
	'zone t = "%s", i = %i\n' %(casename, len(jrange))

data = []

for j in para.j_probes:
	sptm.read_ESTS(workpath + 'ESTS/', j)
	sptm.calc_para_direct()

	temp = np.vstack([
		para.ys[j]/para.lc,
		sptm.elip[0][0]/para.uc,
		sptm.elip[1][0]/para.uc,
		sptm.elip[2][0]/para.uc,
		sptm.elip[3][0]/para.uc,
		sptm.elip[0][1]/para.uc,
		sptm.elip[1][1]/para.uc,
		sptm.elip[2][1]/para.uc,
		sptm.elip[3][1]/para.uc ])

	data.append(temp.T)

np.savetxt(workpath + 'Para_elip_direct.dat', data, header=header, comments='')






header = \
	'Title = "first and second order moments of streamwise space-time spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	%(	r"$k_x^+$",
		r"$\frac{-\omega_{c,uu}^+}{k_x^+}$",
		r"$\frac{-\omega_{c,vv}^+}{k_x^+}$",
		r"$\frac{-\omega_{c,ww}^+}{k_x^+}$",
		r"$\frac{-\omega_{c,pp}^+}{k_x^+}$",
		r"$\frac{\sqrt{B_{uu}^+}}{k_x^+}$",
		r"$\frac{\sqrt{B_{vv}^+}}{k_x^+}$",
		r"$\frac{\sqrt{B_{ww}^+}}{k_x^+}$",
		r"$\frac{\sqrt{B_{pp}^+}}{k_x^+}$"	)

with open(workpath + 'Para_lamw_direct.dat', 'w') as fp:
	fp.write(header)
	for j in para.j_probes:
		header = 'zone t = "%s j = %i", i = %i\n' %(casename, j, len(irange))

		sptm.read_ESTS(workpath + 'ESTS/', j)
		sptm.calc_para_direct()

		data = np.vstack([
			para.kx[:para.Nxc]*para.lc,
			sptm.lamw[0][0][j]/para.uc,
			sptm.lamw[1][0][j]/para.uc,
			sptm.lamw[2][0][j]/para.uc,
			sptm.lamw[3][0][j]/para.uc,
			sptm.lamw[0][1][j]/para.uc,
			sptm.lamw[1][1][j]/para.uc,
			sptm.lamw[2][1][j]/para.uc,
			sptm.lamw[3][1][j]/para.uc ])
		data = data.T[irange]
		np.savetxt('tempfile', data, header=header, comments='')

		with open('tempfile') as temp:
			for line in temp: fp.write(line)










sptm.calc_para_model()


header = \
	'Title = "elliptic model parameters for space-time correlation reconstructed from space data"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	%(	r"$y^+$",
		r"$U_{uu}^+$", r"$U_{vv}^+$", r"$U_{ww}^+$", r"$U_{pp}^+$",
		r"$V_{uu}^+$", r"$V_{vv}^+$", r"$V_{ww}^+$", r"$V_{pp}^+$"	) + \
	'zone t = "%s", i = %i\n' %(casename, len(jrange))
data = np.vstack([
	para.ys/para.lc,
	sptm.elip[0][0]/para.uc,
	sptm.elip[1][0]/para.uc,
	sptm.elip[2][0]/para.uc,
	sptm.elip[3][0]/para.uc,
	sptm.elip[0][1]/para.uc,
	sptm.elip[1][1]/para.uc,
	sptm.elip[2][1]/para.uc,
	sptm.elip[3][1]/para.uc ])
data = data.T[jrange]
np.savetxt(workpath + 'Para_elip_recons.dat', data, header=header, comments='')




header = \
	'Title = "LAMW parameters reconstructed from space data"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	%(	r"$k_x^+$",
		r"$\frac{-\omega_{c,uu}^+}{k_x^+}$",
		r"$\frac{-\omega_{c,vv}^+}{k_x^+}$",
		r"$\frac{-\omega_{c,ww}^+}{k_x^+}$",
		r"$\frac{-\omega_{c,pp}^+}{k_x^+}$",
		r"$\frac{\sqrt{B_{uu}^+}}{k_x^+}$",
		r"$\frac{\sqrt{B_{vv}^+}}{k_x^+}$",
		r"$\frac{\sqrt{B_{ww}^+}}{k_x^+}$",
		r"$\frac{\sqrt{B_{pp}^+}}{k_x^+}$"	)

with open(workpath + 'Para_lamw_recons.dat', 'w') as fp:
	fp.write(header)
	for j in para.j_probes:
		header = 'zone t = "%s j = %i", i = %i\n' %(casename, j, len(irange))
		data = np.vstack([
			para.kx[:para.Nxc]*para.lc,
			sptm.lamw[0][0][j]/para.uc,
			sptm.lamw[1][0][j]/para.uc,
			sptm.lamw[2][0][j]/para.uc,
			sptm.lamw[3][0][j]/para.uc,
			sptm.lamw[0][1][j]/para.uc,
			sptm.lamw[1][1][j]/para.uc,
			sptm.lamw[2][1][j]/para.uc,
			sptm.lamw[3][1][j]/para.uc ])
		data = data.T[irange]
		np.savetxt('tempfile', data, header=header, comments='')

		with open('tempfile') as temp:
			for line in temp: fp.write(line)



exit()



with open(workpath + 'wallscale.txt', 'w') as fp:
	para.inner_scale()
	fp.write("Re_tau = %.18e\n"%para.Ret)
	fp.write("u_tau = %.18e\n"%para.utau)
	fp.write("tau_w = %.18e\n"%para.tauw)
	fp.write("delta_nu = %.18e\n"%para.dnu)
	fp.write("t_nu = %.18e\n"%para.tnu)
	fp.write("dy_min_plus = %.18e\n"%((para.ys[1]-para.ys[0])/para.dnu))
	fp.write("dy_max_plus = %.18e\n"%((para.ys[para.Ny//2+1]-para.ys[para.Ny//2])/para.dnu))





stas.flipy()


header = \
	'Title = "mean values of basic variables"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s"\n' %(r"$y^+$", r"$U^+$", r"$V^+$", r"$W^+$", r"$P^+$") + \
	'zone t = "%s", i = %i\n' %(casename, len(jrange))
data = np.vstack([
	para.ys/para.lc,
	stas.Um/para.uc,
	stas.Vm/para.uc,
	stas.Wm/para.uc,
	stas.Pm/para.pc ])
data = data.T[jrange]
np.savetxt(workpath + 'means_plot.dat', data, header=header, comments='')



header = \
	'Title = "2nd order statistics of basic variables"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	%(	r"$y^+$",
		r"$<u'u'>^+$", r"$<v'v'>^+$", r"$<w'w'>^+$", r"$<u'v'>^+$",
		r"$<p'p'>^+$", r"$<u'p'>^+$", r"$<v'p'>^+$", r"$<w'p'>^+$" ) + \
	'zone t = "%s", i = %i\n' %(casename, len(jrange))
data = np.vstack([
	para.ys/para.lc,
	stas.R11/para.uc**2,
	stas.R22/para.uc**2,
	stas.R33/para.uc**2,
	stas.R12/para.uc**2,
	stas.Rpp/para.pc**2,
	stas.Rpu/para.uc/para.pc,
	stas.Rpv/para.uc/para.pc,
	stas.Rpw/para.uc/para.pc ])
data = data.T[jrange]
np.savetxt(workpath + 'statis_plot.dat', data, header=header, comments='')



exit()




# header = \
# 	'Title = "profiles of budgets"\n' + \
# 	'variables = "%s", "%s"\n' % ( "y<sup>+</sup>", "<greek>e</greek><sup>+</sup>" ) + \
# 	'zone t = "%s", i = %i' %( casename, len(jrange) )
# data = np.vstack([ para.yc/stas.lc, bgts.epsl/(stas.uc**3/stas.lc) ])
# data = data.T[jrange]
# np.savetxt(para.postpath+"budgets.dat", data, header=header, comments='')





header = \
	'Title = "2D energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
		"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
		"y<sup>+</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
	'zone t = "%s", i = %i, j = %i, k = %i' %( casename, len(jrange), len(krange), len(irange) )

data = np.empty([10, len(irange), len(krange), len(jrange)])
for i in irange:
	for k in krange:
		for j in jrange:
			data[:, irange.index(i), krange.index(k), jrange.index(j)] = [
				para.kx[i], para.kz[k], para.yc[j],
				stas.Euu[j,k,i],
				stas.Evv[j,k,i],
				stas.Eww[j,k,i],
				stas.Epp[j,k,i],
				stas.Euv[j,k,i],
				stas.Evw[j,k,i],
				stas.Euw[j,k,i]	]

data[3:] *= data[0] * data[1] / (4*np.pi**2 / para.Lx / para.Lz) / stas.uc**2
data[6] *= stas.uc**2 / stas.pc**2
data[:2] = np.log10(2*np.pi / data[:2] / stas.lc)
data[2] /= stas.lc
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES2D.dat"
np.savetxt(pame, data, header=header, comments='')
if not system("preplot " + pame):
	system("rm -f " + pame)





header = \
	'Title = "1D streamwise energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
		"log<sub>10</sub>(y<sup>+</sup>)",
		"k<sub>x</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
	'zone t = "%s", i = %i, j = %i' %( casename, len(jrange), len(irange) )

data = np.empty([9, len(irange), len(jrange)])
for i in irange:
	for j in jrange:
		data[:, irange.index(i), jrange.index(j)] = [
			para.kx[i], para.yc[j],
			np.sum(stas.Euu[j,:,i]),
			np.sum(stas.Evv[j,:,i]),
			np.sum(stas.Eww[j,:,i]),
			np.sum(stas.Epp[j,:,i]),
			np.sum(stas.Euv[j,:,i]),
			np.sum(stas.Evw[j,:,i]),
			np.sum(stas.Euw[j,:,i])	]

data[2:] *= data[0] / (2*np.pi / para.Lx) / stas.uc**2
data[5] *= stas.uc**2 / stas.pc**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / stas.lc)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_xy.dat"
np.savetxt(pame, data, header=header, comments='')
if not system("preplot " + pame):
	system("rm -f " + pame)





header = \
	'Title = "1D spanwise energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
		"log<sub>10</sub>(y<sup>+</sup>)",
		"k<sub>z</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
	'zone t = "%s", i = %i, j = %i' %( casename, len(jrange), len(krange) )

data = np.empty([9, len(krange), len(jrange)])
for k in krange:
	for j in jrange:
		data[:, krange.index(k), jrange.index(j)] = [
			para.kz[k], para.yc[j],
			np.sum(stas.Euu[j,k]),
			np.sum(stas.Evv[j,k]),
			np.sum(stas.Eww[j,k]),
			np.sum(stas.Epp[j,k]),
			np.sum(stas.Euv[j,k]),
			np.sum(stas.Evw[j,k]),
			np.sum(stas.Euw[j,k])	]

data[2:] *= data[0] / (2*np.pi / para.Lz) / stas.uc**2
data[5] *= stas.uc**2 / stas.pc**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / stas.lc)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_zy.dat"
np.savetxt(pame, data, header=header, comments='')
if not system("preplot " + pame):
	system("rm -f " + pame)





