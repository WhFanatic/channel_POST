from basic import *
from operators import Operators


class SpaceTime:
	def __init__(self, para):
		self.para = para
		self.feld = Field(para)

		self.calc_Nt(12)

		self.dltx = np.arange(-para.Nx//2, para.Nx//2) * (para.Lx/para.Nx)
		self.dltz = np.arange(-para.Nz//2, para.Nz//2) * (para.Lz/para.Nz)
		self.dltt = np.arange(-self.Nt//2, self.Nt//2) * (self.Lt/self.Nt)

	def calc_Nt(self, Cpls):
		self.para.inner_scale()
		dt = self.para.nsprb * self.para.dt
		nt = (self.para.Lx / (Cpls * self.para.uc)) / dt
		self.Nt = int(2**np.ceil(np.log2(nt))) # nextpow2
		self.Lt = self.Nt * dt
		self.ts = np.arange(self.Nt) * dt
		self.kt = fftfreq(self.Nt, d=dt) * 2*np.pi



	def calc_es2d_dt(self, tsteps=None):
		Ny = self.para.Ny
		Nz = self.para.Nz
		Nxc= self.para.Nxc
		if tsteps is None: tsteps = self.para.tsteps

		self.uu,  self.vv,  self.ww,  self.pp  = (np.zeros([Ny+1,Nz,Nxc]) for _ in range(4))
		self.uut, self.vvt, self.wwt, self.ppt = (np.zeros([Ny+1,Nz,Nxc], dtype=complex) for _ in range(4))
		self.utut,self.vtvt,self.wtwt,self.ptpt= (np.zeros([Ny+1,Nz,Nxc]) for _ in range(4))

		for tstep in tsteps:
			print("Calculating spacetime statis: tstep", tstep)

			u, um = self.feld.read_fluc_mean("U%08i.BIN"%tstep)
			v, vm = self.feld.read_fluc_mean("V%08i.BIN"%tstep)
			w, wm = self.feld.read_fluc_mean("W%08i.BIN"%tstep)
			p, pm = self.feld.read_fluc_mean("P%08i.BIN"%tstep)

			ut, utm = self.feld.read_fluc_mean("UT%08i.BIN"%tstep)
			vt, vtm = self.feld.read_fluc_mean("VT%08i.BIN"%tstep)
			wt, wtm = self.feld.read_fluc_mean("WT%08i.BIN"%tstep)
			pt, ptm = self.feld.read_fluc_mean("PT%08i.BIN"%tstep)

			self.uu += np.abs(u)**2 / len(tsteps) # energy spectra here is equal to that in Statis, but without flipk or flipy
			self.vv += np.abs(v)**2 / len(tsteps)
			self.ww += np.abs(w)**2 / len(tsteps)
			self.pp += np.abs(p)**2 / len(tsteps)

			self.uut += u.conj() * ut / len(tsteps)
			self.vvt += v.conj() * vt / len(tsteps)
			self.wwt += w.conj() * wt / len(tsteps)
			self.ppt += p.conj() * pt / len(tsteps)

			self.utut += np.abs(ut)**2 / len(tsteps)
			self.vtvt += np.abs(vt)**2 / len(tsteps)
			self.wtwt += np.abs(wt)**2 / len(tsteps)
			self.ptpt += np.abs(pt)**2 / len(tsteps)

	def calc_para_model(self):
		kx = self.para.kx[:self.para.Nxc]

		U1, V1 = self.__calc_elip_para(self.uu, self.uut, self.utut, kx)
		U2, V2 = self.__calc_elip_para(self.vv, self.vvt, self.vtvt, kx)
		U3, V3 = self.__calc_elip_para(self.ww, self.wwt, self.wtwt, kx)
		U4, V4 = self.__calc_elip_para(self.pp, self.ppt, self.ptpt, kx)

		omgc1, B1 = self.__calc_lamw_para(self.uu, self.uut, self.utut)
		omgc2, B2 = self.__calc_lamw_para(self.vv, self.vvt, self.vtvt)
		omgc3, B3 = self.__calc_lamw_para(self.ww, self.wwt, self.wtwt)
		omgc4, B4 = self.__calc_lamw_para(self.pp, self.ppt, self.ptpt)

		u1, v1 = - omgc1 / (kx + 1e-20), B1**.5 / (kx + 1e-20)
		u2, v2 = - omgc2 / (kx + 1e-20), B2**.5 / (kx + 1e-20)
		u3, v3 = - omgc3 / (kx + 1e-20), B3**.5 / (kx + 1e-20)
		u4, v4 = - omgc4 / (kx + 1e-20), B4**.5 / (kx + 1e-20)

		self.elip = [(U1,V1), (U2,V2), (U3,V3), (U4,V4)]
		self.lamw = [(u1,v1), (u2,v2), (u3,v3), (u4,v4)]

	def calc_corts_elip(self, j):
		Ny = self.para.Ny
		kx = self.para.kx[:self.para.Nxc]

		uu = .5 * (self.uu[j] + self.uu[Ny-j])
		vv = .5 * (self.vv[j] + self.vv[Ny-j])
		ww = .5 * (self.ww[j] + self.ww[Ny-j])
		pp = .5 * (self.pp[j] + self.pp[Ny-j])
		uut = .5 * (self.uut[j] + self.uut[Ny-j])
		vvt = .5 * (self.vvt[j] + self.vvt[Ny-j])
		wwt = .5 * (self.wwt[j] + self.wwt[Ny-j])
		ppt = .5 * (self.ppt[j] + self.ppt[Ny-j])
		utut = .5 * (self.utut[j] + self.utut[Ny-j])
		vtvt = .5 * (self.vtvt[j] + self.vtvt[Ny-j])
		wtwt = .5 * (self.wtwt[j] + self.wtwt[Ny-j])
		ptpt = .5 * (self.ptpt[j] + self.ptpt[Ny-j])

		U1, V1 = self.__calc_elip_para(uu, uut, utut, kx)
		U2, V2 = self.__calc_elip_para(vv, vvt, vtvt, kx)
		U3, V3 = self.__calc_elip_para(ww, wwt, wtwt, kx)
		U4, V4 = self.__calc_elip_para(pp, ppt, ptpt, kx)

		self.Ruu = self.__calc_corts_elip(uu, U1, V1, self.dltx, self.dltt)
		self.Rvv = self.__calc_corts_elip(vv, U2, V2, self.dltx, self.dltt)
		self.Rww = self.__calc_corts_elip(ww, U3, V3, self.dltx, self.dltt)
		self.Rpp = self.__calc_corts_elip(pp, U4, V4, self.dltx, self.dltt)

	def calc_ests_lamw(self, j, tsteps=None):
		Ny = self.para.Ny

		if tsteps is None: tsteps = self.para.tsteps

		u = np.vstack([[self.feld.read_fluc_layer('U%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])
		v = np.vstack([[self.feld.read_fluc_layer('V%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])
		w = np.vstack([[self.feld.read_fluc_layer('W%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])
		p = np.vstack([[self.feld.read_fluc_layer('P%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])

		ut = np.vstack([[self.feld.read_fluc_layer('UT%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])
		vt = np.vstack([[self.feld.read_fluc_layer('VT%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])
		wt = np.vstack([[self.feld.read_fluc_layer('WT%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])
		pt = np.vstack([[self.feld.read_fluc_layer('PT%08i.BIN'%t, _) for _ in (j,Ny-j)] for t in tsteps])

		# u = np.vstack([self.feld.read_fluc('U%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])
		# v = np.vstack([self.feld.read_fluc('V%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])
		# w = np.vstack([self.feld.read_fluc('W%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])
		# p = np.vstack([self.feld.read_fluc('P%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])

		# ut = np.vstack([self.feld.read_fluc('UT%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])
		# vt = np.vstack([self.feld.read_fluc('VT%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])
		# wt = np.vstack([self.feld.read_fluc('WT%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])
		# pt = np.vstack([self.feld.read_fluc('PT%08i.BIN'%tstep)[[j,Ny-j]] for tstep in tsteps])

		self.Euu = self.__calc_ests_lamw(u, ut, self.Nt, self.Lt)
		self.Evv = self.__calc_ests_lamw(v, vt, self.Nt, self.Lt)
		self.Eww = self.__calc_ests_lamw(w, wt, self.Nt, self.Lt)
		self.Epp = self.__calc_ests_lamw(p, pt, self.Nt, self.Lt)



	def calc_ests(self, j):
		self.Euu = self.__calc_ests(self.feld.read_probe('U%08i.BIN'%j), self.Nt)
		self.Evv = self.__calc_ests(self.feld.read_probe('V%08i.BIN'%j), self.Nt)
		self.Eww = self.__calc_ests(self.feld.read_probe('W%08i.BIN'%j), self.Nt)
		self.Epp = self.__calc_ests(self.feld.read_probe('P%08i.BIN'%j), self.Nt)

	def calc_corts(self, j):
		self.Ruu = self.__calc_corts(self.feld.read_probe('U%08i.BIN'%j), self.Nt)
		self.Rvv = self.__calc_corts(self.feld.read_probe('V%08i.BIN'%j), self.Nt)
		self.Rww = self.__calc_corts(self.feld.read_probe('W%08i.BIN'%j), self.Nt)
		self.Rpp = self.__calc_corts(self.feld.read_probe('P%08i.BIN'%j), self.Nt)

	def calc_para_direct(self):
		kx = self.para.kx[:self.para.Nxc]

		U1, V1, omgc1, B1 = self.__calc_ests_para(self.Euu, kx, self.kt)
		U2, V2, omgc2, B2 = self.__calc_ests_para(self.Evv, kx, self.kt)
		U3, V3, omgc3, B3 = self.__calc_ests_para(self.Eww, kx, self.kt)
		U4, V4, omgc4, B4 = self.__calc_ests_para(self.Epp, kx, self.kt)

		u1, v1 = - omgc1 / (kx + 1e-20), B1**.5 / (kx + 1e-20)
		u2, v2 = - omgc2 / (kx + 1e-20), B2**.5 / (kx + 1e-20)
		u3, v3 = - omgc3 / (kx + 1e-20), B3**.5 / (kx + 1e-20)
		u4, v4 = - omgc4 / (kx + 1e-20), B4**.5 / (kx + 1e-20)

		self.elip = [(U1,V1), (U2,V2), (U3,V3), (U4,V4)]
		self.lamw = [(u1,v1), (u2,v2), (u3,v3), (u4,v4)]



	def write_ES2D_DT(self, path):
		write_bin(path + 'ES2D_UU.BIN', self.uu)
		write_bin(path + 'ES2D_VV.BIN', self.vv)
		write_bin(path + 'ES2D_WW.BIN', self.ww)
		write_bin(path + 'ES2D_PP.BIN', self.pp)

		write_bin(path + 'ES2D_DUDTDUDT.BIN', self.utut)
		write_bin(path + 'ES2D_DVDTDVDT.BIN', self.vtvt)
		write_bin(path + 'ES2D_DWDTDWDT.BIN', self.wtwt)
		write_bin(path + 'ES2D_DPDTDPDT.BIN', self.ptpt)

		write_channel(path + 'ES2D_UDUDT.BIN', self.uut)
		write_channel(path + 'ES2D_VDVDT.BIN', self.vvt)
		write_channel(path + 'ES2D_WDWDT.BIN', self.wwt)
		write_channel(path + 'ES2D_PDPDT.BIN', self.ppt)

	def read_ES2D_DT(self, path):
		self.uu = read_bin(path + 'ES2D_UU.BIN')
		self.vv = read_bin(path + 'ES2D_VV.BIN')
		self.ww = read_bin(path + 'ES2D_WW.BIN')
		self.pp = read_bin(path + 'ES2D_PP.BIN')

		self.utut = read_bin(path + 'ES2D_DUDTDUDT.BIN')
		self.vtvt = read_bin(path + 'ES2D_DVDTDVDT.BIN')
		self.wtwt = read_bin(path + 'ES2D_DWDTDWDT.BIN')
		self.ptpt = read_bin(path + 'ES2D_DPDTDPDT.BIN')

		self.uut = read_channel(path + 'ES2D_UDUDT.BIN', self.para.Ny+1, self.para.Nz, self.para.Nxc)
		self.vvt = read_channel(path + 'ES2D_VDVDT.BIN', self.para.Ny+1, self.para.Nz, self.para.Nxc)
		self.wwt = read_channel(path + 'ES2D_WDWDT.BIN', self.para.Ny+1, self.para.Nz, self.para.Nxc)
		self.ppt = read_channel(path + 'ES2D_PDPDT.BIN', self.para.Ny+1, self.para.Nz, self.para.Nxc)

	def write_ESTS(self, path, j):
		write_bin(path + 'ESTS_UU%08i.BIN'%j, self.Euu)
		write_bin(path + 'ESTS_VV%08i.BIN'%j, self.Evv)
		write_bin(path + 'ESTS_WW%08i.BIN'%j, self.Eww)
		write_bin(path + 'ESTS_PP%08i.BIN'%j, self.Epp)

	def read_ESTS(self, path, j):
		self.Euu = read_bin(path + 'ESTS_UU%08i.BIN'%j)
		self.Evv = read_bin(path + 'ESTS_VV%08i.BIN'%j)
		self.Eww = read_bin(path + 'ESTS_WW%08i.BIN'%j)
		self.Epp = read_bin(path + 'ESTS_PP%08i.BIN'%j)

	def write_CORTS(self, path, j):
		write_bin(path + 'CORTS_UU%08i.BIN'%j, self.Ruu)
		write_bin(path + 'CORTS_VV%08i.BIN'%j, self.Rvv)
		write_bin(path + 'CORTS_WW%08i.BIN'%j, self.Rww)
		write_bin(path + 'CORTS_PP%08i.BIN'%j, self.Rpp)

	def read_CORTS(self, path, j):
		self.Ruu = write_bin(path + 'CORTS_UU%08i.BIN'%j)
		self.Rvv = write_bin(path + 'CORTS_VV%08i.BIN'%j)
		self.Rww = write_bin(path + 'CORTS_WW%08i.BIN'%j)
		self.Rpp = write_bin(path + 'CORTS_PP%08i.BIN'%j)

	def flipy(self):
		self.uu += self.uu[::-1]
		self.vv += self.vv[::-1]
		self.ww += self.ww[::-1]
		self.pp += self.pp[::-1]
		self.uu /= 2.
		self.vv /= 2.
		self.ww /= 2.
		self.pp /= 2.

		self.uut += self.uut[::-1]
		self.vvt += self.vvt[::-1]
		self.wwt += self.wwt[::-1]
		self.ppt += self.ppt[::-1]
		self.uut /= 2.
		self.vvt /= 2.
		self.wwt /= 2.
		self.ppt /= 2.

		self.utut += self.utut[::-1]
		self.vtvt += self.vtvt[::-1]
		self.wtwt += self.wtwt[::-1]
		self.ptpt += self.ptpt[::-1]
		self.utut /= 2.
		self.vtvt /= 2.
		self.wtwt /= 2.
		self.ptpt /= 2.

	@staticmethod
	def __calc_ests(q, Nt):
		''' direct method to compute space-time spectrum over (omg,kz,kx) '''
		print('computing space-time spectrum using direct method ...')

		Ns = len(q) // (Nt//2) - 1

		window = np.hanning(Nt) / (np.sum(np.hanning(Nt)**2)/Nt)**.5 # the normalized window defined in de Kat 2015

		qq = np.zeros([Nt] + list(q.shape[1:]))

		for s in range(Ns):
			qq += np.abs(ifft(q[Nt//2*s:Nt//2*s+Nt].T * window).T)**2 / Ns

		qq *= np.mean(np.abs(q[:Nt//2*(Ns+1)])**2, axis=0) / (np.sum(qq, axis=0) + 1e-20) # rescale to recover the total energy altered by window function

		return qq

	@staticmethod
	def __calc_corts(q, Nt):
		''' direct method to compute space-time spectrum over (dt,dz,dx), time complexity O( Nx*Nz*Nt * (Nt*Ns + log(Nx*Nz) ) '''
		print('computing space-time correlation using direct method ...')

		Ns = len(q) // (Nt//2) - 1

		qq = np.zeros([Nt//2] + list(q.shape[1:]), dtype=np.complex64)

		for s in range(Ns):
			q2 = q[Nt//2*s:Nt//2*(s+1)]
			q1c = q2.conj()
			for t in range(Nt//2):
				qq[t] += np.mean(q1c*q2, axis=0)
				q2 = q[Nt//2*s+t:Nt//2*(s+1)+t]

		Nxc, Nzc = q.shape[-1], q.shape[-2]//2

		return np.roll(np.roll(phys(qq) / Ns, Nxc, axis=-1), Nzc, axis=-2)

	@staticmethod
	def __calc_ests_para(Eqq, kx, kt):
		''' extract EA & LAMW parameters from directly computed space-time spectrum '''
		Ext = np.sum(Eqq, axis=-2)
		Ex  = np.sum(Ext, axis=-2) + 1e-20

		dRdxx = np.sum(SpaceTime.__flipk(-kx**2 * Ext)) + 1e-20
		dRdtt = np.sum(SpaceTime.__flipk(-kt**2 * Ext.T))
		dRdxt = np.sum(SpaceTime.__flipk(-kt*(kx*Ext).T))

		U = - dRdxt / dRdxx
		V = ( dRdtt / dRdxx - U**2 )**.5

		# calculate 1st and 2nd order moments around the ridge of the energy spectrum
		# to avoid the effect of sidelobe due to limited time resolution
		ridge = - U * kx

		kts = np.array([omg-ridge for omg in kt]) # [Nt, Nxc]
		kts[kts<np.min(kt)] += (kt[1] - kt[0]) * len(kt)
		kts[kts>np.max(kt)] -= (kt[1] - kt[0]) * len(kt)

		omgc = np.sum(kts * Ext, axis=-2) / Ex
		B = np.sum((kts-omgc)**2 * Ext, axis=-2) / Ex
		omgc += ridge

		return U, V, omgc, B

	@staticmethod
	def __calc_elip_para(qq, qqt, qtqt, kx):
		''' calculate EA parameters from space data over y '''
		dRdxx = np.sum(SpaceTime.__flipk(-kx**2 * qq), axis=(-1,-2))
		dRdxt = np.sum(SpaceTime.__flipk(-kx*1j * qqt), axis=(-1,-2)).real
		dRdtt = np.sum(SpaceTime.__flipk(-        qtqt), axis=(-1,-2))

		Us = - dRdxt / (dRdxx + 1e-20)
		Vs = ( dRdtt / (dRdxx + 1e-20) - Us**2 )**.5

		return Us, Vs

	@staticmethod
	def __calc_lamw_para(qq, qqt, qtqt):
		''' calculate LAMW parameters from space data over (kx,y) '''
		Eqq = np.sum(qq, axis=-2) + 1e-20
		Eqqt = np.sum(qqt, axis=-2)
		Eqtqt = np.sum(qtqt, axis=-2)

		omgcs = - Eqqt.imag / Eqq
		Bs = Eqtqt / Eqq - omgcs**2

		return omgcs, Bs

	@staticmethod
	def __calc_corts_elip(qq, U, V, dltx, dltt):
		''' reconstruct space-time correlation of a single layer over (dx,dt) using EA '''
		Rxz = phys(qq) # qq must not be flipk
		Rx = np.roll(Rxz.T[:,0].T, Rxz.shape[-1]//2, axis=-1)
		dltxc = [((dltx - U * t)**2 + (V * t)**2)**.5 for t in dltt]
		return np.expand_dims(np.interp(dltxc, dltx, Rx), axis=-2)

	@staticmethod
	def __calc_ests_lamw(q, qt, Nt, Lt):
		''' reconstruct space-time energy spectrum of a single layer over (kx,omg) using LAMW '''
		print('reconstructing space-time spectrum using LAMW ...')

		Nxc = q.shape[-1]
		domg = 2*np.pi / Lt

		qkx = fft(q, axis=-2) # from 2D spec to 1D spec
		qtq = fft(qt, axis=-2) / (qkx + 1e-20) # dqdt / q
		qq = np.abs(qkx)**2 # 1D energy spec

		Ext = np.zeros([Nt, Nxc])
		
		for omg_prime in (- qtq.imag + qtq.real, - qtq.imag - qtq.real):
			# relate local frequency to discrete frequency grids, negative ones happen to land on the trailing half
			omgid = np.floor(omg_prime/domg + .5).astype(int)

			for ek,ts,es in zip(Ext.T, omgid.T, qq.T): # for every kx
				# eliminate out-of-range local-frequencies
				es[ts<-Nt//2] = ts[ts<-Nt//2] = 0 # the former is evaluated first, then the later
				es[ts>=Nt//2] = ts[ts>=Nt//2] = 0
				# distribute energy to discrete frequencies
				for t,e in zip(ts.ravel(), es.ravel()):
					ek[t] += .5 * e # average between 2 local frequencies

		# note: not divided by spatial energy spectrum because we aim to solve non-normalized space-time spectrum
		# note: not divided by delta omega because we aim to compute discrete spectrum, not the continuous one
		return np.expand_dims(Ext, axis=-2)

	@staticmethod
	def __calc_CONV2D(qqt, qq, kx):
		return qqt.imag / (kx * qq + 1e-20)

	@staticmethod
	def __calc_CONV1D(qqt, qq, kx):
		# kx = self.para.kx[:self.para.Nxc]
		cq = SpaceTime.calc_CONV2D(qqt, qq, kx)
		cqx = np.sum(cq * qq, axis=-2) / (np.sum(qq, axis=-2) + 1e-20)
		cqz = np.sum(cq * qq * kx**2, axis=-1) / (np.sum(qq * kx**2, axis=-1) + 1e-20)
		return cqx, cqz

	@staticmethod
	def __flipk(q):
		p = q.copy()
		p.T[1:] *= 2
		return p









