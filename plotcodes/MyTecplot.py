#!/work1/cuigx2_work/whn/anaconda_install/anaconda3/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Tecplot:

	def parse_head(self, fpn):
		with open(fpn) as fp:
			line = fp.readline()


			line = fp.readline().strip()
			records = line.split('=')[-1].strip().split(',')
			labels = [ record.strip().strip('"') for record in records ]

			line = fp.readline().strip()
			line = line[line.index("zone")+4:].strip()
			records = line.split(',')
			dic = { record.split('=')[0].strip(): record.split('=')[1].strip() for record in records }
			coords = []
			if 'i' in dic.keys(): coords.append( int(dic['i']) )
			if 'j' in dic.keys(): coords.append( int(dic['j']) )
			if 'k' in dic.keys(): coords.append( int(dic['k']) )

			# line = fp.readline().strip()
			# records = []
			# record = False
			# for c in line:
			# 	if c == '"':
			# 		if record == False:
			# 			record = ''
			# 		else:
			# 			records.append(record)
			# 			record = False
			# 	elif record != False:
			# 		record += c
			# labels = records

			# line = fp.readline().strip()
			# records = []
			# record = False
			# for c in line:
			# 	if c == '=':
			# 		record = ''
			# 	elif c == ',':
			# 		records.append(int(record))
			# 		record = False
			# 	elif record != False:
			# 		record += c
			# coords = records + [int(record)]

		return labels, coords

	def write_head(self, fpn, title, labels, coords):
		with open(fpn, 'w') as fp:
			fp.write( 'Title = "%s"\n' %title )
			fp.write( 'variables = "%s"\n' %( '", "'.join(labels) ) )
			if len(coords) == 1: fp.write( 'zone i = %i\n' %tuple(coords) )
			if len(coords) == 2: fp.write( 'zone i = %i, j = %i\n' %tuple(coords) )
			if len(coords) == 3: fp.write( 'zone i = %i, j = %i, k = %i\n' %tuple(coords) )



	def coord(self, fpn, n):
		return self.parse_head(fpn)[1][n]

	def dim(self, fpn):
		return len(self.parse_head(fpn)[1])

	def label(self, fpn, n):
		labels, coords = self.parse_head(fpn)
		return labels[n-1+self.dim(fpn)]

	def xlabel(self, fpn):
		return self.parse_head(fpn)[0][0]

	def ylabel(self, fpn):
		return self.parse_head(fpn)[0][1]


	def getData(self, fpn, n, rows=[]):
		labels, coords = self.parse_head(fpn)
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






# class Cases:
# 	def __init__(pathfile = "workpath.txt"):

# 		with open(pathfile) as fp:
# 			for line in fp.readlines():
# 				line = line.strip().split().strip()
# 				if line and line[0] != '#':


# test
if __name__ == "__main__":
	tcplt = Tecplot()
	labels, coords = tcplt.parse_head("E:\\POST\\data\\post_M1000\\CORTS_xt_jprb2.dat")
	print(labels)
	print(coords)




