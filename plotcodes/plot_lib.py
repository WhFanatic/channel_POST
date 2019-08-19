#!/work1/cuigx2_work/whn/anaconda_install/anaconda2/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def parse_head(file_path_name):
	with open(file_path_name) as fp:
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

def plot_line(ax, file_path_name, n, **kwarg):
	labels, coords = parse_head(file_path_name)
	ylen, = coords
	x = np.zeros(ylen)
	y = np.zeros(ylen)
	with open(file_path_name) as fp:
		line = fp.readline()
		line = fp.readline()
		line = fp.readline()
		for j in range(ylen):
			line = fp.readline().strip().split()
			x[j] = float( line[0] )
			y[j] = float( line[n] )

	ax.semilogx(x, y, label=labels[n], lw=2, **kwarg)

	return [ labels[0], labels[n] ]

def plot_contour(ax, file_path_name, n, filled=0, **kwarg):
	labels, coords = parse_head(file_path_name)
	ylen, xlen = coords
	x = np.zeros([ylen, xlen])
	y = np.zeros([ylen, xlen])
	z = np.zeros([ylen, xlen])
	with open(file_path_name) as fp:
		line = fp.readline()
		line = fp.readline()
		line = fp.readline()
		for i in range(xlen):
			for j in range(ylen):
				line = fp.readline().strip().split()
				x[j,i] = float( line[0] )
				y[j,i] = float( line[1] )
				z[j,i] = float( line[n+1] )

	if filled == 0:
		cs = ax.contour(x, y, z, linewidths=2, **kwarg)
	elif filled == 1:
		cs = ax.contourf(x, y, z, cmap=plt.cm.rainbow, extend='both', **kwarg)
	elif filled == 2:
		ax.contourf(x, y, z, cmap=plt.cm.rainbow, extend='both', **kwarg)
		cs = ax.contour(x, y, z, linewidths=1, **kwarg)

	return [ labels[0], labels[1], labels[n+1] ], cs

def mark_out(ax, xs=[], ys=[], **kwarg):
	xmin, xmax, ymin, ymax = ax.axis()
	for x in xs:
		ax.plot([x,x], [ymin, ymax], '--k', **kwarg)
	for y in ys:
		ax.plot([xmin, xmax], [y,y], '--k', **kwarg)


def mark_out_domain(ax, xmin, xmax, ymin, ymax):
	ax.plot( # mark out the computational domain size
		[xmin, xmax], [ymin, ymin], '--k',
		[xmin, xmax], [ymax, ymax], '--k',
		[xmin, xmin], [ymin, ymax], '--k',
		[xmax, xmax], [ymin, ymax], '--k', lw=0.5	)



