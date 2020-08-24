from basic import *


class Operators:
	def __init__(self, para):
		self.para = para

	@staticmethod
	def convol(u, v=None):
		return spec(phys(u) * phys(v)) if v is not None else spec(phys(u)**2)

	@staticmethod
	def normalize(q):
		return q / np.max(np.abs(q))

	def diff(self):
		## centeral difference coefficients for 1st & 2nd order derivative in y direction (upwind difference at wall)
		pass
		# diff1 = np.zeros([Ny+1, 3])
		# diff2 = np.zeros([Ny+1, 3])
		# for j in range(Ny+1):
		# 	if 0 < j < Ny:
		# 		dy1 = y_mesh[j-1] - y_mesh[j]
		# 		dy2 = y_mesh[j+1] - y_mesh[j]
		# 		diff1[j] = np.array([dy2/dy1/(dy2-dy1), -(1/dy1+1/dy2), -dy1/dy2/(dy2-dy1)])
		# 		diff2[j] = 2/(dy1*dy2) * np.array([-dy2/(dy2-dy1), 1, dy1/(dy2-dy1)])
		# 	elif j == 0:
		# 		dy1 = y_mesh[1] - y_mesh[0]
		# 		dy2 = y_mesh[2] - y_mesh[0]
		# 		diff1[j] = np.array([-(1/dy1+1/dy2), dy2/dy1/(dy2-dy1), -dy1/dy2/(dy2-dy1)])
		# 		diff2[j] = 2/(dy1*dy2) * np.array([1, -dy2/(dy2-dy1), dy1/(dy2-dy1)])
		# 	elif j == Ny:
		# 		dy1 = y_mesh[-2] - y_mesh[-1]
		# 		dy2 = y_mesh[-3] - y_mesh[-1]
		# 		diff1[j] = np.array([-dy1/dy2/(dy2-dy1), dy2/dy1/(dy2-dy1), -(1/dy1+1/dy2)])
		# 		diff2[j] = 2/(dy1*dy2) * np.array([dy1/(dy2-dy1), -dy2/(dy2-dy1), 1])



