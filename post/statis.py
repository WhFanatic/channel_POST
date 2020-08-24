from basic import *
from operators import Operators


class Statis:
	def __init__(self, para):
		self.para = para
		self.feld = Field(para)

	def calc_statis(self, tsteps=None):
		Ny = self.para.Ny
		Nz = self.para.Nz
		Nxc= self.para.Nxc
		if tsteps is None: tsteps = self.para.tsteps

		self.R11, self.R22, self.R33          = (np.zeros(Ny+1) for _ in range(3))
		self.R12, self.R23, self.R13          = (np.zeros(Ny+1) for _ in range(3))
		self.Rpu, self.Rpv, self.Rpw, self.Rpp= (np.zeros(Ny+1) for _ in range(4))
		self.Um , self.Vm , self.Wm , self.Pm = (np.zeros(Ny+1) for _ in range(4))
		self.Euu, self.Evv, self.Eww, self.Epp= (np.zeros([Ny+1, Nz, Nxc]) for _ in range(4))
		self.Euv, self.Evw, self.Euw          = (np.zeros([Ny+1, Nz, Nxc]) for _ in range(3))

		for tstep in tsteps:
			print("Reading statis: tstep", tstep)
			u, um = self.feld.read_fluc_mean("U%08i.BIN"%tstep)
			v, vm = self.feld.read_fluc_mean("V%08i.BIN"%tstep)
			w, wm = self.feld.read_fluc_mean("W%08i.BIN"%tstep)
			p, pm = self.feld.read_fluc_mean("P%08i.BIN"%tstep)

			self.Um += um / len(tsteps)
			self.Vm += vm / len(tsteps)
			self.Wm += wm / len(tsteps)
			self.Pm += pm / len(tsteps)

			self.Euu += np.abs(u)**2 / len(tsteps)
			self.Evv += np.abs(v)**2 / len(tsteps)
			self.Eww += np.abs(w)**2 / len(tsteps)
			self.Epp += np.abs(p)**2 / len(tsteps)
			self.Euv +=(np.conj(u)*v).real / len(tsteps)
			self.Evw +=(np.conj(v)*w).real / len(tsteps)
			self.Euw +=(np.conj(u)*w).real / len(tsteps)

			u = phys(u)
			v = phys(v)
			w = phys(w)
			p = phys(p)

			self.R11 += np.mean(u**2,axis=(-1,-2)) / len(tsteps)
			self.R22 += np.mean(v**2,axis=(-1,-2)) / len(tsteps)
			self.R33 += np.mean(w**2,axis=(-1,-2)) / len(tsteps)
			self.R12 += np.mean(u*v, axis=(-1,-2)) / len(tsteps)
			self.R23 += np.mean(v*w, axis=(-1,-2)) / len(tsteps)
			self.R13 += np.mean(u*w, axis=(-1,-2)) / len(tsteps)

			self.Rpu += np.mean(p*u, axis=(-1,-2)) / len(tsteps)
			self.Rpv += np.mean(p*v, axis=(-1,-2)) / len(tsteps)
			self.Rpw += np.mean(p*w, axis=(-1,-2)) / len(tsteps)
			self.Rpp += np.mean(p**2,axis=(-1,-2)) / len(tsteps)

	def write_es2d(self, path):
		write_bin(path + 'ES2D_UU.BIN', self.Euu)
		write_bin(path + 'ES2D_VV.BIN', self.Evv)
		write_bin(path + 'ES2D_WW.BIN', self.Eww)
		write_bin(path + 'ES2D_PP.BIN', self.Epp)
		write_bin(path + 'ES2D_UV.BIN', self.Euv)

	def read_es2d(self, path):
		self.Euu = read_bin(path + 'ES2D_UU.BIN')
		self.Evv = read_bin(path + 'ES2D_VV.BIN')
		self.Eww = read_bin(path + 'ES2D_WW.BIN')
		self.Epp = read_bin(path + 'ES2D_PP.BIN')
		self.Euv = read_bin(path + 'ES2D_UV.BIN')

	def flipy(self):
		self.Um += self.Um[::-1]
		self.Vm -= self.Vm[::-1]
		self.Wm += self.Wm[::-1]
		self.Pm += self.Pm[::-1]
		self.Um /= 2.
		self.Vm /= 2.
		self.Wm /= 2.
		self.Pm /= 2.

		self.R11 += self.R11[::-1]
		self.R22 += self.R22[::-1]
		self.R33 += self.R33[::-1]
		self.R12 -= self.R12[::-1]
		self.R23 -= self.R23[::-1]
		self.R13 += self.R13[::-1]
		self.R11 /= 2.
		self.R22 /= 2.
		self.R33 /= 2.
		self.R12 /= 2.
		self.R23 /= 2.
		self.R13 /= 2.

		self.Rpu += self.Rpu[::-1]
		self.Rpv -= self.Rpv[::-1]
		self.Rpw += self.Rpw[::-1]
		self.Rpp += self.Rpp[::-1]
		self.Rpu /= 2.
		self.Rpv /= 2.
		self.Rpw /= 2.
		self.Rpp /= 2.

		self.Euu += self.Euu[::-1]
		self.Evv += self.Evv[::-1]
		self.Eww += self.Eww[::-1]
		self.Epp += self.Epp[::-1]
		self.Euv -= self.Euv[::-1]
		self.Evw -= self.Evw[::-1]
		self.Euw += self.Euw[::-1]
		self.Euu /= 2.
		self.Evv /= 2.
		self.Eww /= 2.
		self.Epp /= 2.
		self.Euv /= 2.
		self.Evw /= 2.
		self.Euw /= 2.

	@staticmethod
	def __flipk(q):
		q.T[1:] *= 2
		return q

	@staticmethod
	def __flipk2(q):
		''' fold all energy to the [:nzc,:nxc] range  '''
		nzcu = q.shape[-2]//2 + q.shape[-2]%2
		nzcd = q.shape[-2]//2
		q.T[1:] *= 2 # .T usually returns a view which can act as left value
		q.T[:,1:nzcu] += q.T[:,:nzcd:-1]
		return (q.T[:,:nzcd+1]).T





