"""
Exports the BrushGenerator class
"""

import gzip
from abc import ABC, abstractmethod
from enum import Enum
from typing import Optional

import numpy as np
import pandas as pd

from PoissonDiskGenerator import PoissonDiskGenerator


class BrushGenerator(ABC):
	"""
	Generate a LAMMPS data file containing a coarse-grained polymer brush grafted to a planar wall in a rectangular box.
	See https://lammps.sandia.gov/doc/read_data.html
	"""

	AtomTypes: Enum = Enum('AtomTypes', [])
	BondTypes: Enum = Enum('BondTypes', [])
	AngleTypes: Enum = Enum('AngleTypes', [])
	DihedralTypes: Enum = Enum('DihedralTypes', [])

	masses: dict = {}
	pair_ij_coeffs: dict = {}
	bond_coeffs: dict = {}
	angle_coeffs: dict = {}
	dihedral_coeffs: dict = {}

	styles: dict[str] = {
		'pair': '',
		'bond': '',
		'angle': '',
		'dihedral': '',
	}

	def __init__(self, box_size: tuple[float, float, Optional[float]], rng_seed: Optional[int], bead_size: float,
	             n_beads: int, bottom_padding: float = 0):
		"""
		:param box_size:       3-tuple of floats describing the dimensions of the rectangular box. If the third (z)
		                       value is None, it will be automatically sized to contain the longest chain.
		:param rng_seed:       Seed used to initialize the PRNG. May be None, in which case a random seed will be used.
		:param bead_size:      Size of the 'grafting beads': used as minimum distance for the Poisson-disk point set
		                       generator.
		:param n_beads:        Chain length.
		:param bottom_padding: Distance between the bottom edge of the box and the grafting layer. Must be positive.
		"""
		self.box_size = list(box_size)
		self.rng_seed = rng_seed
		self.bead_size = bead_size
		self.n_beads = n_beads
		self.bottom_padding = bottom_padding

		self.coordinates: np.ndarray = np.array([])

		# Non-final lists
		self._atoms_list = []
		self._bonds_list = []
		self._angles_list = []
		self._dihedrals_list = []

		# Final dataframes
		self.atoms: pd.DataFrame = pd.DataFrame()
		self.bonds: pd.DataFrame = pd.DataFrame()
		self.angles: pd.DataFrame = pd.DataFrame()
		self.dihedrals: pd.DataFrame = pd.DataFrame()

	@abstractmethod
	def _build_bead(self, mol_id: int, graft_coord: np.ndarray, bead_id: int) -> float:
		"""
		Adds a bead to the instance's atom/bond/angle/dihedral lists.
		Override this and implement according to the polymer model used. Note that LAMMPS expects ids to be 1-indexed.
		:param mol_id:      0-indexed molecule (chain) id
		:param graft_coord: 2-element ndarray containing xy coordinate of the grafting point
		:param bead_id:     0-indexed bead id
		:return Maximum z height
		"""
		pass

	def generate_grafting_layer(self, n_chains: int, max_overlap_iter: int = 10**3) -> int:
		"""
		Generate coordinates of the grafting layer using a Poisson-disk point set generator.
		:param n_chains:         Number of grafting points (chains).
		:param max_overlap_iter: Iteration limit for the Poisson-disk point set generator.
		:return Number of grafting points generated.
		"""
		# Generate grafting point coordinates
		pdg = PoissonDiskGenerator(self.rng_seed)
		self.coordinates = pdg.generate(n_chains, self.bead_size, self.box_size[0:2], max_overlap_iter)

		return len(self.coordinates)

	def _build_chain(self) -> float:
		"""
		Create chains by looping over chains and beads in each chain, calling _build_bead() for each bead.
		:return Maximum z height
		"""
		z_max = 0
		# Loop over chains
		for mol_id, i in enumerate(self.coordinates):
			# Loop over successive beads in chain
			for j in range(0, self.n_beads + 1):
				z = self._build_bead(mol_id, i, j)
				if z > z_max:
					z_max = z

		return z_max

	def build(self) -> None:
		"""
		Create atom positions and molecular topology for a randomly-grafted monodisperse AdG-brush.
		"""
		# Set z size to z_max from _build_chain() if not set
		z_max = self._build_chain()
		if not self.box_size[2]:
			self.box_size[2] = z_max + 1

		# Make dataframes from the non-final lists created by _build_chain()
		self.atoms = pd.DataFrame(self._atoms_list, columns=['mol_id', 'atom_type', 'q', 'x', 'y', 'z'])
		self.bonds = pd.DataFrame(self._bonds_list, columns=['bond_type', 'atom1', 'atom2'])
		self.angles = pd.DataFrame(self._angles_list, columns=['angle_type', 'atom1', 'atom2', 'atom3'])
		self.dihedrals = pd.DataFrame(self._dihedrals_list, columns=['dihedral_type', 'atom1', 'atom2', 'atom3',
		                                                             'atom4'])

		# LAMMPS ids start at 1
		self.atoms.index += 1
		self.bonds.index += 1
		self.angles.index += 1
		self.dihedrals.index += 1

	def write(self, filename: str) -> None:
		"""
		Write the LAMMPS data file.
		:param filename: Filename for the output file.
		"""
		num_atoms = len(self.atoms)
		num_bonds = len(self.bonds)
		num_angles = len(self.angles)
		num_dihedrals = len(self.dihedrals)

		num_atom_types = len(self.AtomTypes)
		num_bond_types = len(self.BondTypes) if not self.bonds.empty else 0
		num_angle_types = len(self.AngleTypes) if not self.angles.empty else 0
		num_dihedral_types = len(self.DihedralTypes) if not self.dihedrals.empty else 0

		# Auto-detect gzipped files
		o = gzip.open if filename.endswith('.gz') else open

		with o(filename, 'xt', newline='\n') as f:
			# Header
			f.write("#Header\n")
			f.write(f"{num_atoms} atoms\n")
			if num_bonds > 0:     f.write(f"{num_bonds} bonds\n")
			if num_angles > 0:    f.write(f"{num_angles} angles\n")
			if num_dihedrals > 0: f.write(f"{num_dihedrals} dihedrals\n\n")

			f.write(f"{num_atom_types} atom types\n")
			if num_bond_types > 0:     f.write(f"{num_bond_types} bond types\n")
			if num_angle_types > 0:    f.write(f"{num_angle_types} angle types\n")
			if num_dihedral_types > 0: f.write(f"{num_dihedral_types} dihedral types\n\n")

			# Box geometry
			f.write(f"0 {self.box_size[0]} xlo xhi\n")
			f.write(f"0 {self.box_size[1]} ylo yhi\n")
			f.write(f"-{self.bottom_padding} {self.box_size[2]} zlo zhi\n\n")

			# Force field coeffs
			f.write("Masses\n\n")
			for k, v in self.masses.items():
				f.write(f"{k.value} {v}\n")
			f.write("\n")

			if len(self.pair_ij_coeffs) > 0:
				f.write("PairIJ Coeffs" + (f" # {self.styles['pair']}" if self.styles['pair'] else '') + "\n\n")
				for k, vs in self.pair_ij_coeffs.items():
					f.write(f"{k[0].value} {k[1].value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			if len(self.bond_coeffs) > 0:
				f.write("Bond Coeffs" + (f" # {self.styles['bond']}" if self.styles['bond'] else '') + "\n\n")
				for k, vs in self.bond_coeffs.items():
					f.write(f"{k.value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			if len(self.angle_coeffs) > 0:
				f.write("Angle Coeffs" + (f" # {self.styles['angle']}" if self.styles['angle'] else '') + "\n\n")
				for k, vs in self.angle_coeffs.items():
					f.write(f"{k.value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			if len(self.dihedral_coeffs) > 0:
				f.write("Dihedral Coeffs" + (f" # {self.styles['dihedral']}" if self.styles['dihedral'] else '') + "\n\n")
				for k, vs in self.dihedral_coeffs.items():
					f.write(f"{k.value} " + ' '.join([f'{v:.3g}' for v in vs]) + "\n")
				f.write("\n")

			# Atom properties
			f.write("Atoms # full\n\n")
			self.atoms.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.3g')
			f.write("\n")

			# Molecular topology
			if len(self.bonds) > 0:
				f.write(f"Bonds\n\n")
				self.bonds.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.3g')
				f.write("\n")
			if len(self.angles) > 0:
				f.write("Angles\n\n")
				self.angles.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.3g')
				f.write("\n")
			if len(self.dihedrals) > 0:
				f.write("Dihedrals\n\n")
				self.dihedrals.to_csv(f, sep=' ', header=False, index=True, line_terminator='\n', float_format='%.3g')
				f.write("\n")
