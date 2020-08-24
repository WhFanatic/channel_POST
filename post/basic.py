#!/work1/cuigx2_work/whn/anaconda_install/anaconda3/bin/python
import numpy as np
from numpy.fft import fft,ifft, fft2,ifft2, hfft, ihfft, fftfreq
from scipy.integrate import trapz
from os import listdir, system


##### common functions #####

def spec(q): return ifft(ihfft(q).T[:-1].T, axis=-2)
def phys(q): return hfft(fft(q, axis=-2), n=q.shape[-1]*2)

def write_bin(pame, q):
	ny, nz, nx = q.shape
	info = np.array([ny, nz, nx]+[0]*(nx*nz-3))
	np.hstack([
		np.frombuffer(info.astype(np.int32).tobytes(), np.float32),
		np.ravel(q).astype(np.float32) ]).tofile(pame)

def read_bin(pame):
	with open(pame, 'rb') as fp:
		ny, nz, nx = np.frombuffer(fp.read(12), np.int32)
	q = np.fromfile(pame, np.float32).reshape([ny+1, nz, nx])
	return q[1:]

def write_channel(pame, q):
	ny, nz, nx = q.shape
	info = np.array([ny, nz, nx]+[0]*(nx*nz*2-3))
	data = np.transpose([q.real.T, q.imag.T])
	np.hstack([
		np.frombuffer(info.astype(np.int32).tobytes(), np.float32),
		np.ravel(data).astype(np.float32) ]).tofile(pame)

def read_channel(pame, ny, nz, nx):
	q = np.fromfile(pame, np.float32).reshape([ny+1, nz, nx, 2])
	return q[1:,:,:,0] + 1j * q[1:,:,:,1]

def read_channel_layer(pame, nz, nx, j):
	with open(pame, 'rb') as fp:
		fp.seek((j+1)*nz*nx*8)
		q = np.frombuffer(fp.read(nx*nz*8), np.float32).reshape([nz, nx, 2])
	return q[:,:,0] + 1j * q[:,:,1]

###########################


class DataSetInfo:
	def __init__(self, path):
		self.datapath = path

		self.read_XIN(self.datapath)
		self.read_GRD(self.runtimepath)

		self.Nxc, self.kx = self.__get_wavenumber(self.Nx, self.Lx)
		self.Nzc, self.kz = self.__get_wavenumber(self.Nz, self.Lz)

		self.tsteps = self.__get_tsteps(self.fieldpath)

		self.calc_wallscale()

	def read_XIN(self, path):

		self.fieldpath = path + 'fielddata/'
		self.runtimepath = path + 'runtimedata/'
		self.probepath = path + 'probedata/'

		with open(path + 'XINDAT') as fp:
			lines = fp.readlines()

			info = lines[1].split()
			self.nu = 1. / float(info[0])
			self.dt = float(info[3])
			self.Nx = int(info[4])
			self.Ny = int(info[5])
			self.Nz = int(info[6])
			self.Lx = 2.*np.pi / float(info[1])
			self.Ly = 2.
			self.Lz = 2.*np.pi / float(info[2])

			info = lines[3].split()
			self.nsp = int(info[1])
			self.nsprb = int(info[4])

			info = lines[17].split()
			self.j_probes = [int(n) for n in info]

	def read_GRD(self, path):
		
		try:
			self.ys = 1. + np.loadtxt(path + 'INITFIELD.DAT').T[0]
		except IOError:
			self.ys = 1. - np.cos( np.arange(self.Ny+1) * np.pi / self.Ny )
		
		xs = self.Lx * np.arange(self.Nx) / self.Nx
		zs = self.Lz * np.arange(self.Nz) / self.Nz

	def __get_wavenumber(self, Nx, Lx):
		kx = fftfreq(self.Nx, d=self.Lx/self.Nx) * 2*np.pi
		Nxc = Nx // 2
		return Nxc, kx

	def __get_tsteps(self, path):
		s, e = np.loadtxt(path + 'datainfo').astype(int)
		return range(s, e+1, self.nsp)


	def calc_wallscale(self, tsteps=None):
		if tsteps is None: tsteps = self.tsteps

		Re = 1. / self.nu
		logs = np.loadtxt(self.runtimepath+'XUTP.DAT')
		
		self.tauw = np.mean([log[1] for log in logs if np.min(np.abs(log[0]/self.dt-np.array(tsteps))) < .5])
		self.utau = self.tauw**.5
		self.dnu = 1./Re / self.utau
		self.tnu = self.dnu / self.utau
		self.Ret = 1./self.dnu # channel height is taken for 2.0

	def inner_scale(self):
		self.lc = self.dnu
		self.tc = self.tnu
		self.uc = self.utau
		self.pc = self.tauw

	def outer_scale(self):
		self.lc = 1.
		self.tc = 1.
		self.uc = 1.
		self.pc = 1.



class Field:
	def __init__(self, para):
		self.para = para

	def read(self, name):
		return read_channel(self.para.fieldpath + self.__infer_path(name), self.para.Ny+1, self.para.Nz, self.para.Nxc)

	def read_layer(self, name, j):
		return read_channel_layer(self.para.fieldpath + self.__infer_path(name), self.para.Nz, self.para.Nxc, j)

	def read_probe(self, name):
		q = np.reshape(np.fromfile(
			self.para.probepath + self.__infer_path(name), np.float32),
			[-1, self.para.Nz, self.para.Nxc, 2])
		q.T[0,0] = 0
		return q[1:,:,:,0] + 1j * q[1:,:,:,1]

	def read_fluc_mean(self, name):
		q = self.read(name)
		qm = q.T[0,0].real.copy()
		q.T[0,0] = 0
		return q, qm

	def read_mean(self, name):
		return self.read_fluc_mean(name)[1]

	def read_fluc(self, name):
		return self.read_fluc_mean(name)[0]

	def read_fluc_layer(self, name, j):
		q = self.read_layer(name, j)
		q[0,0] = 0
		return q

	@staticmethod
	def __infer_path(name):
		prefix = name[:name.find(''.join(filter(str.isdigit, name)))]
		return prefix + '/' + name




