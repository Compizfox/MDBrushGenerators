"""
Exports the PolydisperseKGBrushGenerator and PoissonKGBrushGenerator classes
"""
from typing import Optional

import numpy as np
from scipy.stats._distn_infrastructure import rv_frozen
from scipy.stats import poisson

from KremerGrestBrushGenerator import KremerGrestBrushGenerator


class PolydisperseKGBrushGenerator(KremerGrestBrushGenerator):
	"""
	Generate a LAMMPS data file containing a polydisperse Kremer-Grest polymer brush grafted to a planar wall in a
	rectangular box.
	"""
	def __init__(self, box_size: tuple[float, float, Optional[float]], rng_seed: Optional[int], cld: rv_frozen,
	             cg_factor: float = 1, graft: bool = True):
		"""
		:param box_size:  3-tuple of floats describing the dimensions of the rectangular box. If the third
		                  (z) value is None, it will be automatically sized to contain the longest chain.
		:param rng_seed:  Seed used to initialize the PRNG. May be None, in which case a random seed will be used.
		:param cld:       Chain length distribution.
		:param cg_factor: Coarse-graining factor (number of real monomers per coarse-grained bead) that the chain
		                  lengths drawn from the CLD will be divided by.
		:param graft:     Generates grafted brushes when True, and non-grafted films when False.
		"""
		self.cld = cld
		self.cld_factor = cg_factor
		self.rng = np.random.RandomState(seed=rng_seed)

		super().__init__(box_size, rng_seed, 0, graft)

	def _build_chain(self) -> float:
		z_max = 0
		# Loop over chains
		for mol_id, i in enumerate(self.coordinates):
			# Randomly draw a number from the chain length distribution
			n = round(self.cld.rvs(random_state=self.rng) / self.cld_factor)
			# Loop over successive beads in chain
			for j in range(0, n + 1):
				z = self._build_bead(mol_id, i, j)
				if z > z_max:
					z_max = z

		return z_max


class PoissonKGBrushGenerator(PolydisperseKGBrushGenerator):
	"""
	Generate a LAMMPS data file containing a Kremer-Grest polymer brush, with chain lengths obeying a Poisson
	distribution, grafted to a planar wall in a	rectangular box.
	"""
	def __init__(self, box_size: tuple[float, float, Optional[float]], rng_seed: Optional[int], n_mean: int,
	             cg_factor: float = 1, graft: bool = True):
		"""
		:param box_size:  3-tuple of floats describing the dimensions of the rectangular box.
		:param rng_seed:  Seed used to initialize the PRNG. May be None, in which case a random seed will be used.
		:param n_mean:    Number-average chain length.
		:param cg_factor: Coarse-graining factor (number of real monomers per coarse-grained bead) that will be taken
		                  into account for the Poisson distribution.
		:param graft:     Generates grafted brushes when True, and non-grafted films when False.
		"""
		cld = poisson(n_mean * cg_factor)
		super().__init__(box_size, rng_seed, cld, cg_factor, graft)
