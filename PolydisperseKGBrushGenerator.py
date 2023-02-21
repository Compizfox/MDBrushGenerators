"""
Exports the PolydisperseKGBrushGenerator and PoissonKGBrushGenerator classes
"""
from typing import Optional

import numpy as np
from scipy.stats._distn_infrastructure import rv_discrete, rv_frozen
from scipy.stats import poisson
from scipy.special import gamma

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
				z = self._build_bead(mol_id, i, j, n)
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


class SchulzZimmKGBrushGenerator(PolydisperseKGBrushGenerator):
	"""
	Generate a LAMMPS data file containing a Kremer-Grest polymer brush, with chain lengths obeying a Schulz-Zimm
	distribution, grafted to a planar wall in a	rectangular box.
	"""
	class SchulzZimm_gen(rv_discrete):
		"""
		Schulz-Zimm distribution

		(1) de Vos, W. M.; Leermakers, F. A. M. Modeling the Structure of a Polydisperse Polymer Brush.
		Polymer 2009, 50 (1), 305–316. https://doi.org/10.1016/j.polymer.2008.10.025.
		"""

		def _pmf(self, N_i, N_n, x):
			return x ** (x + 1) / gamma(x + 1) * N_i ** (x - 1) / N_n ** x * np.exp(-x * N_i / N_n)

	SchulzZimm = SchulzZimm_gen(name="Schulz-Zimm")

	def __init__(self, box_size: tuple[float, float, Optional[float]], rng_seed: Optional[int], n_mean: int, x: float,
	             cg_factor: float = 1, graft: bool = True):
		"""
		:param box_size:  3-tuple of floats describing the dimensions of the rectangular box.
		:param rng_seed:  Seed used to initialize the PRNG. May be None, in which case a random seed will be used.
		:param n_mean:    Number-average chain length.
		:param x:         Broadness of the Schulz-Zimm distribution
		:param cg_factor: Coarse-graining factor (number of real monomers per coarse-grained bead) that will be taken
		                  into account for the distribution.
		:param graft:     Generates grafted brushes when True, and non-grafted films when False.
		"""
		cld = self.SchulzZimm(n_mean * cg_factor, x)
		super().__init__(box_size, rng_seed, cld, cg_factor, graft)
