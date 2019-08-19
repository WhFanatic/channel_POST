#!/work1/cuigx2_work/whn/anaconda_install/anaconda3/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Tecplot:

	def __parse_head(self, fpn):
		with open(fpn) as fp:
			line = fp.readline()

			line = fp.readline().strip()
			records = []
			record = False
			for c in line:
				if c == '"':
					if record == False:
						record = ''
					else:
						records.append(record)
						record = False
				elif record != False:
					record += c
			labels = records

			line = fp.readline().strip()
			records = []
			record = False
			for c in line:
				if c == '=':
					record = ''
				elif c == ',':
					records.append(int(record))
					record = False
				elif record != False:
					record += c
			coords = records + [int(record)]

		return labels, coords



	def label(self, fpn, n):
		labels, coords = self.__parse_head(fpn)
		return labels[n-1+len(coords)]

	def xlabel(self, fpn):
		return self.__parse_head(fpn)[0][0]

	def ylabel(self, fpn):
		return self.__parse_head(fpn)[0][1]

	def coord(self, fpn, n):
		return self.__parse_head(fpn)[1][n]

	def getData(fpn, n, rows=[]):
		labels, coords = self.__parse_head(fpn)
		ncol = len(labels)
		nrow = np.prod(coords)
		rowrange = range(nrow) if rows == [] else rows

		with open(fpn) as fp:
			lines = fp.readlines()
		data = np.array([ [ float(a) for a in lines[row+3].strip().split() ] for row in rowrange ])
		
		return data[:, list(range(len(coords)))+[n-1+len(coords)]]



	def curve(self, ax, fpn, n, **kwarg):
		ylen = self.coord(fpn, 0)
		x, y = [ np.zeros(ylen) for i in range(2) ]
		with open(fpn) as fp:
			for i in range(3): fp.readline()
			for j in range(ylen):
				line = fp.readline().strip().split()
				x[j] = float( line[0] )
				y[j] = float( line[n] )
		c = ax.semilogx(x, y, label=self.label(fpn, n), lw=2, **kwarg)
		ax.set_xlabel(self.xlabel(fpn))
		ax.set_ylabel(self.label(fpn, n))
		return c

	def contour(self, ax, fpn, n, filled=0, **kwarg):
		ylen, xlen = self.coord(fpn, 0), self.coord(fpn, 1)
		x, y, z = [ np.zeros([ylen, xlen]) for i in range(3) ]
		with open(fpn) as fp:
			for i in range(3): fp.readline()
			for i in range(xlen):
				for j in range(ylen):
					line = fp.readline().strip().split()
					x[j,i] = float( line[0] )
					y[j,i] = float( line[1] )
					z[j,i] = float( line[n+1] )

		if filled == 0: c = ax.contour(x, y, z, linewidths=2, **kwarg)
		elif filled == 1: c = ax.contourf(x, y, z, cmap=plt.cm.rainbow, extend='both', **kwarg)
		elif filled == 2:
			ax.contourf(x, y, z, cmap=plt.cm.rainbow, extend='both', **kwarg)
			c = ax.contour(x, y, z, linewidths=1, **kwarg)

		ax.set_xlabel(self.xlabel(fpn))
		ax.set_ylabel(self.ylabel(fpn))
		ax.set_title (self.label(fpn, n))
		return c



	def mark_out(self, ax, xs=[], ys=[], **kwarg):
		xmin, xmax, ymin, ymax = ax.axis()
		for x in xs: ax.plot([x,x], [ymin, ymax], '--k', **kwarg)
		for y in ys: ax.plot([xmin, xmax], [y,y], '--k', **kwarg)

	def mark_out_domain(self, ax, xmin, xmax, ymin, ymax):
		ax.plot( # mark out the computational domain size
			[xmin, xmax], [ymin, ymin], '--k',
			[xmin, xmax], [ymax, ymax], '--k',
			[xmin, xmin], [ymin, ymax], '--k',
			[xmax, xmax], [ymin, ymax], '--k', lw=0.5	)




class Case:
	def __init__(self, path, name="case0", color=None, style=None, marker=None, markevery=None):

		self.path = path
		self.name = name
		self.color = color
		self.style = style
		self.marker = marker
		self.markevery = markevery

		self.plot = Tecplot()

		self.y_plus = np.loadtxt(open(self.path+'y_plus.dat'))
		self.j_probes = np.array([ int(n) for n in np.ravel(np.loadtxt(open(self.path+'j_probes.dat'))) ])
		self.Re_tau = np.loadtxt(open(self.path+'Re_tau.dat'))
		self.u_tau = np.loadtxt(open(self.path+'u_tau.dat'))
		self.delta_nu = 1.0 / self.Re_tau
		self.t_nu = self.delta_nu / self.u_tau
		self.tau_w = self.u_tau**2

	def copy(self, path=None, name=None, color=None, style=None, marker=None, markevery=None):
		return Case(
			path = path if path else self.path,
			name = name if name else self.name,
			color = color if color else self.color,
			style = style if style else self.style,
			marker = marker if marker else self.marker,
			markevery = markevery if markevery else self.markevery	)

	def curve(self, ax, filename, n, **kwarg):
		self.plot.curve(ax, self.path+filename, n, color=self.color, ls=self.style, marker=self.marker, markevery=self.markevery, **kwarg)

	def contour(self, ax, filename, n, filled=0, **kwarg):
		self.plot.contour(ax, self.path+filename, n, filled, colors=self.color, linestyles=self.style, **kwarg)


