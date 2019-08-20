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


	def coord(self, fpn, n):
		return self.__parse_head(fpn)[1][n]

	def dim(self, fpn):
		return len(self.__parse_head(fpn)[1])

	def label(self, fpn, n):
		labels, coords = self.__parse_head(fpn)
		return labels[n-1+self.dim(fpn)]

	def xlabel(self, fpn):
		return self.__parse_head(fpn)[0][0]

	def ylabel(self, fpn):
		return self.__parse_head(fpn)[0][1]


	def getData(self, fpn, n, rows=[]):
		labels, coords = self.__parse_head(fpn)
		ncol = len(labels)
		nrow = np.prod(coords)
		rowrange = range(nrow) if rows==[] else rows

		with open(fpn) as fp:	lines = fp.readlines()
		data = np.array([ [float(a) for a in lines[row+3].strip().split()] for row in rowrange ])
		
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
		c = ax.semilogx(x, y, label=self.label(fpn, n), **kwarg)
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

		if filled == 0: c = ax.contour(x, y, z, **kwarg)
		elif filled == 1: c = ax.contourf(x, y, z, cmap=plt.cm.rainbow, extend='both', **kwarg)
		elif filled == 2:
			ax.contourf(x, y, z, cmap=plt.cm.rainbow, extend='both', **kwarg)
			c = ax.contour(x, y, z, **kwarg)

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
	def __init__(self, path, name="case0", color=None, style=None):

		self.path = path
		self.name = name

		self.color = color
		self.style = style
		self.width = 2

		self.plot = Tecplot()

		self.y_plus = np.loadtxt(open(self.path+'y_plus.dat'))
		self.j_probes = np.array([ int(n) for n in np.ravel(np.loadtxt(open(self.path+'j_probes.dat'))) ])
		self.Re_tau = np.loadtxt(open(self.path+'Re_tau.dat'))
		self.u_tau = np.loadtxt(open(self.path+'u_tau.dat'))
		self.delta_nu = 1.0 / self.Re_tau
		self.t_nu = self.delta_nu / self.u_tau
		self.tau_w = self.u_tau**2

	def copy(self, path=-1, name=-1, color=-1, style=-1):
		return Case(
			path = self.path if path==-1 else path,
			name = self.name if name==-1 else name,
			color = self.color if color==-1 else color,
			style = self.style if style==-1 else style	)

	def curve(self, ax, filename, n, **kwarg):
		ka = {key:kwarg[key] for key in kwarg.keys()}
		if "color" not in ka.keys(): ka["color"] = self.color
		if "linestyle" not in ka.keys() and "ls" not in ka.keys(): ka["linestyle"] = self.style
		if "linewidth" not in ka.keys() and "lw" not in ka.keys(): ka["linewidth"] = self.width

		c = self.plot.curve(ax, self.path+filename, n, **ka)

	def contour(self, ax, filename, n, filled=0, **kwarg):
		ka = {key:kwarg[key] for key in kwarg.keys()}
		if "colors" not in ka.keys(): ka["colors"] = self.color
		if "linestyles" not in ka.keys(): ka["linestyles"] = self.style
		if "linewidths" not in ka.keys(): ka["linewidths"] = self.width

		c = self.plot.contour(ax, self.path+filename, n, filled, **ka)



# class Cases:
# 	def __init__(pathfile = "workpath.txt"):

# 		with open(pathfile) as fp:
# 			for line in fp.readlines():
# 				line = line.strip().split().strip()
# 				if line and line[0] != '#':


# test
if __name__ == "__main__":
	case = Case('')
	case2 = case.copy()
	print(case2)





