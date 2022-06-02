"""
Exports the PAAmMartiniBG class
"""

from enum import Enum
from typing import Optional

import numpy as np

from BrushGenerator import BrushGenerator


class PAAmMartiniBG(BrushGenerator):
	"""
	Generate a LAMMPS data file containing a brush of Martini CG poly(acrylamide) grafted to a planar wall in a
	rectangular box.

	Perrin, E.; Schoen, M.; Coudert, F.-X.; Boutin, A.
	Structure and Dynamics of Solvated Polymers near a Silica Surface: On the Different Roles Played by Solvent.
	J. Phys. Chem. B 2018, 122 (16), 4573–4582. https://doi.org/10.1021/acs.jpcb.7b11753.
	"""

	AtomTypes = Enum('AtomTypes', ['graft', 'C', 'A', 'water'])
	BondTypes = Enum('BondTypes', ['CC', 'CA'])
	AngleTypes = Enum('AngleTypes', ['CCC', 'CCA'])

	masses = {
		AtomTypes.graft: 2*12 + 3*1,
		AtomTypes.C:     2*12 + 3*1,
		AtomTypes.A:     1*12 + 1*16 + 1*14 + 2*1,
		AtomTypes.water: 4*(1*16+2*1),
	}

	styles = {
		'pair': 'lj/gromacs',
		'bond': 'harmonic',
		'angle': 'cosine/squared',
	}

	pair_ij_coeffs = {
		#                                    ε (kcal/mol)   # σ (Å)
		(AtomTypes.graft, AtomTypes.graft): [0,             0],
		(AtomTypes.graft, AtomTypes.C):     [0.239 * 2.625, 4.3],
		(AtomTypes.graft, AtomTypes.A):     [0.239 * 2.042, 4.3],
		(AtomTypes.graft, AtomTypes.water): [0.239 * 2.300, 4.5],
		(AtomTypes.C, AtomTypes.C):         [0.239 * 2.625, 4.3],
		(AtomTypes.C, AtomTypes.A):         [0.239 * 2.042, 4.3],
		(AtomTypes.C, AtomTypes.water):     [0.239 * 2.300, 4.5],
		(AtomTypes.A, AtomTypes.A):         [0.239 * 3.375, 4.3],
		(AtomTypes.A, AtomTypes.water):     [0.239 * 4.500, 4.5],
		(AtomTypes.water, AtomTypes.water): [0.239 * 4.500, 4.7],
		# conversion from kJ/mol
	}

	bond_coeffs = {
		#              K (kcal/mol/A^2)           R0 (Å)
		BondTypes.CC: [0.239 * 566.7 / 10**2 / 2, 0.249 * 10],
		BondTypes.CA: [0.239 * 666.7 / 10**2 / 2, 0.237 * 10],
		# conversion from kJ/mol/nm^2 and including 1/2 factor
	}

	angle_coeffs = {
		#                K (kcal/mol)       θ_0 (°)
		AngleTypes.CCC: [0.239 * 116.7 / 2, 127.5],
		AngleTypes.CCA: [0.239 * 233.3 / 2, 85.5],
		# conversion from kJ/mol and including 1/2 factor
	}

	def __init__(self, box_size: tuple[float, float, Optional[float]], rng_seed: Optional[int], n_beads: int,
	             graft: bool = True):
		"""
		:param box_size: 3-tuple of floats describing the dimensions of the rectangular box. If the third (z) value
		                 is None, it will be automatically sized to contain the longest chain.
		:param rng_seed: Seed used to initialize the PRNG. May be None, in which case a random seed will be used.
		:param n_beads:  Chain length.
		:param graft:    Generates grafted brushes when True, and non-grafted films when False
		"""
		bead_size = 5  # (sigma)
		bottom_padding = 1  # (sigma)
		self.graft = graft
		self.z_spacing = self.bead_size / 2
		super().__init__(box_size, rng_seed, bead_size, n_beads, bottom_padding)

	def _build_bead(self, mol_id: int, graft_coord: np.ndarray, bead_id: int) -> float:
		if bead_id == 0:
			# Omit grafting bead if graft=False
			if self.graft:
				self._atoms_list.append({'mol_id'   : mol_id + 1,
				                         'atom_type': self.AtomTypes.graft.value,
				                         'q'        : 0,
				                         'x'        : graft_coord[0],
				                         'y'        : graft_coord[1],
				                         'z'        : 0
				                         })
		else:
			# Atoms
			# Backbone bead
			self._atoms_list.append({'mol_id'   : mol_id + 1,
			                         'atom_type': self.AtomTypes.C.value,
			                         'q'        : 0,
			                         'x'        : graft_coord[0],
			                         'y'        : graft_coord[1],
			                         'z'        : float(bead_id * self.z_spacing)
			                         })
			# Side group bead
			self._atoms_list.append({'mol_id'   : mol_id + 1,
			                         'atom_type': self.AtomTypes.A.value,
			                         'q'        : 0,
			                         'x'        : graft_coord[0] + self.bead_size,
			                         'y'        : graft_coord[1],
			                         'z'        : float(bead_id * self.z_spacing)
			                         })

			# Bonds
			# atom_id is set to the id of the last atom (the side group)
			atom_id = len(self._atoms_list)
			# Backbone bond
			if bead_id == 1 and self.graft:
				# Grafting bond
				self._bonds_list.append({'bond_type': self.BondTypes.CC.value,
				                         'atom1'    : atom_id - 2,  # grafting atom
				                         'atom2'    : atom_id - 1   # current backbone atom
				                         })
			elif bead_id > 1:
				# Normal backbone bond
				self._bonds_list.append({'bond_type': self.BondTypes.CC.value,
				                         'atom1'    : atom_id - 3,  # previous backbone atom
				                         'atom2'    : atom_id - 1   # current backbone atom
				                         })

			# Side group bond
			self._bonds_list.append({'bond_type': self.BondTypes.CA.value,
			                         'atom1'    : atom_id,     # current side group atom
			                         'atom2'    : atom_id - 1  # current backbone atom
			                         })

			# Angles
			# Backbone (not including grafting atom)
			if bead_id > 2:
				self._angles_list.append({'angle_type': self.AngleTypes.CCC.value,
				                          'atom1'     : atom_id - 5,
				                          'atom2'     : atom_id - 3,
				                          'atom3'     : atom_id - 1
				                          })

			if bead_id > 1:
				# Side group left-down (not including grafting atom)
				self._angles_list.append({'angle_type': self.AngleTypes.CCA.value,
				                          'atom1'     : atom_id - 3,  # previous backbone atom
				                          'atom2'     : atom_id - 1,  # current backbone atom
				                          'atom3'     : atom_id       # current side group atom
				                          })
				# Side group down-right
				self._angles_list.append({'angle_type': self.AngleTypes.CCA.value,
				                          'atom1'     : atom_id - 2,  # previous side group atom
				                          'atom2'     : atom_id - 3,  # previous backbone atom
				                          'atom3'     : atom_id - 1   # current backbone atom
				                          })

		return float(bead_id*self.z_spacing)
