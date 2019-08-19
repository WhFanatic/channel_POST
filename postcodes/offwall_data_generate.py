## Artificially generate offwall data from full channel data
from post_channel import *

data_path2 = '../data/data_F180_ofw_b100_arti/'

jprb0 = 45
jprb1 = Ny - jprb0
jmids = np.arange(jprb0, jprb1+1)


file_name = 'XINDAT'
with open(data_path + file_name) as fp:
	lines = fp.readlines()
lines[1] = ' '.join(lines[1].split()[:-2] + [str(len(jmids)-1)] + [lines[1].split()[-1]]) + '\n'
lines[17] = str(jprb0) + ' ' + str(jprb1) + '\n'
with open(data_path2 + file_name, 'w') as fp:
	for line in lines:
		fp.write(line)

file_name = 'datainfo'
with open(data_path + 'fielddata/' + file_name) as fp:
	lines = fp.readlines()
with open(data_path2 + 'fielddata/' + file_name, 'w') as fp:
	for line in lines:
		fp.write(line)

file_name = 'INITFIELD.DAT'
with open(data_path + 'runtimedata/' + file_name) as fp:
	lines = fp.readlines()
with open(data_path2 + 'runtimedata/' + file_name, 'w') as fp:
	for line in lines[jprb0:jprb1+1]:
		fp.write(line)


file_names = ('XFR.DAT', 'XUTP.DAT')
for file_name in file_names:
	with open(data_path + 'runtimedata/' + file_name) as fp:
		lines = fp.readlines()
	with open(data_path2 + 'runtimedata/' + file_name, 'w') as fp:
		for line in lines:
			fp.write(line)


for ft in ('U','V','W','P','UU','VV','WW','UV','VW','WU','UT','VT','WT','PT'):
	for tstep in tsteps:
		fpn1 = get_path_name(tstep, ft)
		fpn2 = data_path2 + 'fielddata/' + ft + '/' + ft + str(tstep).zfill(8) + '.BIN'

		write_channel_infosec(fpn2, read_channel_infosec(fpn1))

		for j in range(len(jmids)):
			write_channel_layer(fpn2, j, read_channel_layer(fpn1, jmids[j]))




