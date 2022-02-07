"""
Batch generates initial data comprising a Kremer-Grest brush/film, by specifying the chain length, box size,
and grafting density from a list.
Checks for over-density.
"""

from KremerGrestBrushGenerator import KremerGrestBrushGenerator

# List of tuples of  (n, gd)
entries = [
	(50, 0.2),
]


def generate(size, n, graft, chain_length, filename):
	kgbg = KremerGrestBrushGenerator(size, 34882, chain_length, graft=graft)
	n_actual = kgbg.generate_grafting_layer(n, 10**8)

	if n_actual != n:
		print(f'Warning: Grafting layer too dense. {n_actual} grafting points instead of {n}.')
	else:
		print(f'{n_actual} points')

	kgbg.build()
	kgbg.write(filename)


for chain_length, density in entries:
	x = 100
	y = 100
	z = 60
	size = (x, y, z)  # sigma
	n = int(density*size[0]*size[1])
	
	generate(size, n, True, chain_length, f'output/brush_gd{density}_N{chain_length}_{x}x{y}x{z}.pos.gz')
	generate(size, n, False, chain_length, f'output/film_gd{density}_N{chain_length}_{x}x{y}x{z}.pos.gz')
