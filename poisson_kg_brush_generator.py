import matplotlib.pyplot as plt

from PolydisperseKGBrushGenerator import PoissonKGBrushGenerator

chain_length = 100
density = 0.25  # sigma^-2
x = chain_length
y = chain_length
z = chain_length + 1

size = (x, y, z)  # sigma
n = int(density * size[0] * size[1])

bg = PoissonKGBrushGenerator(size, 34882, chain_length, cg_factor=3, graft=True)

n_actual = bg.generate_grafting_layer(n, 10**6)

if n_actual != n:
	print(f'Warning: Grafting layer too dense. {n_actual} grafting points instead of {n}.')
else:
	print(f'{n_actual} points')

bg.build()
bg.write(f'brush_pdPoisson_gd{density}_N{chain_length}_{x}x{y}x{z}.pos.gz')
