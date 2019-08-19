#!/work1/cuigx2_work/whn/anaconda_install/anaconda2/bin/python
from plot_lib import *

case_names = ['b0', 'b5', 'b15', 'b30', 'b50']
plot_colors = {'b0':'black', 'b5':'blue', 'b15':'red', 'b30':'green', 'b50':'purple'}
plot_styles = {'b0':'-', 'b5':'-', 'b15':'-', 'b30':'-', 'b50':'-'}
plot_markers = {'b0':'o', 'b5':'s', 'b15':'^', 'b30':'v', 'b50':'d'}
file_paths = {}
file_paths[case_names[0]] = r'E:\POST\data\data_F180\postdata\figures\\'
file_paths[case_names[1]] = r'E:\POST\data\data_F180_ofw_b5_arti\postdata\figures\\'
file_paths[case_names[2]] = r'E:\POST\data\data_F180_ofw_b15_arti\postdata\figures\\'
file_paths[case_names[3]] = r'E:\POST\data\data_F180_ofw_b30_arti\postdata\figures\\'
file_paths[case_names[4]] = r'E:\POST\data\data_F180_ofw_b50_arti\postdata\figures\\'
figure_path = 'figures/'
js = [0,1,2,3]
ns = [1,2,3,4]

# figure parameters
lambda_x_plus_lim = [20, 6e3]
lambda_z_plus_lim = [10, 3e3]
k_x_plus_lim = [-0.35,0.35]
k_t_plus_lim = [-5.1,5.1]
y_plus_lim = [1, 500]

# data parameters
y_plus = np.zeros(0)
j_probes = np.zeros(0, dtype=int)
Re_tau = u_tau = delta_nu = t_nu = tau_w = 0

# basic functions
def set_path(case_name):
	global y_plus, j_probes, Re_tau, u_tau, delta_nu, t_nu, tau_w
	path = file_paths[case_name]
	y_plus = np.loadtxt(open(path+'y_plus.dat'))
	j_probes = np.array([ int(n) for n in np.ravel(np.loadtxt(open(path+'j_probes.dat'))) ])
	Re_tau = np.loadtxt(open(path+'Re_tau.dat'))
	u_tau = np.loadtxt(open(path+'u_tau.dat'))
	delta_nu = 1.0 / Re_tau
	t_nu = delta_nu / u_tau
	tau_w = u_tau**2
	return path





### basic statistics


fig_name = 'MeanProfile'
fig, axs = plt.subplots(2, 2, sharex=True, num=fig_name, figsize=(6, 4))
ylims = ([0,25], [-1,1], [-1,1], [-1.5,0])
for n in range(4):
	ax = axs.ravel()[n]
	for case_name in case_names:
		path = set_path(case_name)
		labels = plot_line(ax, path+'means_plot.dat', n+1, color=plot_colors[case_name], linestyle=plot_styles[case_name])
	if n == 0:
		ax.legend(case_names, fontsize=7.5, loc='best', handlelength=5, frameon=False)
	ax.set_xlim(y_plus_lim)
	ax.set_ylim(ylims[n])
	if n in [2, 3]:
		ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
# fig.align_labels()
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()

fig_name = 'ReynoldsStress'
fig, axs = plt.subplots(2, 2, sharex=True, num=fig_name, figsize=(6, 4))
ylims = ([0,10], [0,1.5], [0,2.4], [-1,0])
for n in range(4):
	ax = axs.ravel()[n]
	for case_name in case_names:
		path = set_path(case_name)
		labels = plot_line(ax, path+'statis_plot.dat', n+1, color=plot_colors[case_name], linestyle=plot_styles[case_name])
	if n == 0:
		ax.legend(case_names, fontsize=7.5, loc='best', handlelength=5, frameon=False)
	ax.set_xlim(y_plus_lim)
	ax.set_ylim(ylims[n])
	if n in [2, 3]:
		ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()

fig_name = 'PressureCovariance'
fig, axs = plt.subplots(2, 2, sharex=True, num=fig_name, figsize=(6, 4))
ylims = ([0,10], [-0.2,1], [-0.2,0.05], [-0.5,0.5])
for n in range(4):
	ax = axs.ravel()[n]
	for case_name in case_names:
		path = set_path(case_name)
		labels = plot_line(ax, path+'statis_plot.dat', n+5, color=plot_colors[case_name], linestyle=plot_styles[case_name])
	if n == 0:
		ax.legend(case_names, fontsize=7.5, loc='best', handlelength=5, frameon=False)
	ax.set_xlim(y_plus_lim)
	ax.set_ylim(ylims[n])
	if n in [2, 3]:
		ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()

fig_name = 'VorticityCovariance'
fig, axs = plt.subplots(2, 2, sharex=True, num=fig_name, figsize=(6, 4))
ylims = ([0,0.06], [0,0.05], [0,0.2], [-0.01,0.02])
for n in range(4):
	ax = axs.ravel()[n]
	for case_name in case_names:
		path = set_path(case_name)
		labels = plot_line(ax, path+'statis_plot.dat', n+9, color=plot_colors[case_name], linestyle=plot_styles[case_name])
	if n == 0:
		ax.legend(case_names, fontsize=7.5, loc='best', handlelength=5, frameon=False)
	ax.set_xlim(y_plus_lim)
	ax.set_ylim(ylims[n])
	if n in [2, 3]:
		ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()





fig_name = 'EnergySpectra1D'
fig, axs = plt.subplots(len(ns), 2, sharex='col', sharey='row', squeeze=False, num=fig_name, figsize=(6, 3*len(ns)))
contour_levels = [
	[(0.125, 0.25, 0.5, 0.9), (0.2, 0.4, 0.8, 1.6)],
	[(0.025, 0.05, 0.1 ,0.15), (0.03, 0.06, 0.12, 0.2)],
	[(0.05, 0.1, 0.2, 0.28), (0.045, 0.09, 0.18, 0.26)],
	[(0.2, 0.5, 0.7, 0.9), (0.1,0.2,0.4,0.7)]	]
for col in range(2):
	for row in range(len(ns)):
		ax = axs[row, col]
		for case_name in case_names:
			path = set_path(case_name)
			file_name = ('ES1D_xy.dat', 'ES1D_zy.dat') [col]
			labels,cs = plot_contour(ax, path+file_name, ns[row], levels=contour_levels[row][col], colors=plot_colors[case_name], linestyles=plot_styles[case_name])
		ax.set_xlim( (lambda_x_plus_lim, lambda_z_plus_lim) [col] )
		ax.set_ylim([1,1000])
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_title(labels[2])
		if row == len(ns)-1:
			ax.set_xlabel(labels[0])
		if col == 0:
			ax.set_ylabel(labels[1])
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()


fig_name = 'EnergySpectra2D'
fig, axs = plt.subplots(len(ns), len(js), sharex=True, sharey=True, squeeze=False, num=fig_name, figsize=(3*len(js), 3*len(ns)))
for col in range(len(js)):
	for row in range(len(ns)):
		ax = axs[row, col]

		case_name = case_names[0]
		path = set_path(case_name)
		file_name = 'ES2D_xz_jprb%i.dat' %js[col]
		labels,cs1 = plot_contour(ax, path+file_name, ns[row], filled=True)

		case_name = case_names[col+1]
		path = set_path(case_name)
		file_name = 'ES2D_xz_jprb0.dat'
		labels,cs2 = plot_contour(ax, path+file_name, ns[row], colors='black', levels=cs1.levels)

		case_name = case_names[0]
		path = set_path(case_name)

		ax.set_xlim(lambda_x_plus_lim)
		ax.set_ylim(lambda_z_plus_lim)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_title(labels[2]+r", $y^+$ = %i"%int(y_plus[j_probes[js[col]]]))
		if row == len(ns)-1:
			ax.set_xlabel(labels[0])
		if col == 0:
			ax.set_ylabel(labels[1])
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()


exit()




for case_name in case_names:
	path = set_path(case_name)

	fig_name = 'EnergySpectra2D_%s'%case_name
	fig, axs = plt.subplots(len(ns), len(js), sharex=True, sharey=True, squeeze=False, num=fig_name, figsize=(3*len(js), 3*len(ns)))
	for col in range(len(js)):
		for row in range(len(ns)):
			ax = axs[row, col]
			file_name = 'ES2D_xz_jprb%i.dat' %js[col]
			labels,cs = plot_contour(ax, path+file_name, ns[row])
			ax.set_xlim(lambda_x_plus_lim)
			ax.set_ylim(lambda_z_plus_lim)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_title(labels[2]+r", $y^+$ = %.0f" %y_plus[j_probes[js[col]]])
			if row == len(ns)-1:
				ax.set_xlabel(labels[0])
			if col == 0:
				ax.set_ylabel(labels[1])
	fig.tight_layout()
	fig.savefig(figure_path+fig_name+'.png', dpi=200)
	plt.close()


fig_name = 'Mean_U'
fig, ax = plt.figure(num=fig_name, figsize=(3, 3)), plt.gca()
for case_name in case_names:
	path = set_path(case_name)
	labels = plot_line(ax, path+'means_plot.dat', 1, color=plot_colors[case_name], linestyle=plot_styles[case_name])
ax.legend(case_names, fontsize=7.5, loc='upper left', handlelength=5, frameon=False)
ax.set_xlim(y_plus_lim)
ax.set_ylim([0,25])
ax.set_xlabel(labels[0])
ax.set_ylabel(labels[1])
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()

fig_name = 'FlucIntens'
fig, axs = plt.subplots(2, 2, sharex=True, num=fig_name, figsize=(6, 4))
ylims = ([0,10], [0,1.5], [0,2.4], [0,10])
for n in range(4):
	ax = axs.ravel()[n]
	for case_name in case_names:
		path = set_path(case_name)
		labels = plot_line(ax, path+'statis_plot.dat', (1,2,3,5)[n], color=plot_colors[case_name], linestyle=plot_styles[case_name])
	if n == 1:
		ax.legend(case_names, fontsize=7.5, loc='upper left', handlelength=5, frameon=False)
	ax.set_xlim(y_plus_lim)
	ax.set_ylim(ylims[n])
	if n in [2, 3]:
		ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
fig.tight_layout()
fig.savefig(figure_path+fig_name+'.png', dpi=200)
plt.close()





