"""
Generates initial data comprising a Kremer-Grest brush, by specifying the chain length, box size, and grafting density.
Plots the grafting surface using Matplotlib, and checks for over-density.
"""
import matplotlib.pyplot as plt

from KremerGrestBrushGenerator import KremerGrestBrushGenerator

chain_length = 100
density = 0.25  # sigma^-2
x = chain_length
y = chain_length
z = chain_length + 1

size = (x, y, z)  # sigma
n = int(density * size[0] * size[1])

kgbg = KremerGrestBrushGenerator(size, 34882, chain_length, graft=True)

n_actual = kgbg.generate_grafting_layer(n, 10**6)

if n_actual != n:
	print(f'Warning: Grafting layer too dense. {n_actual} grafting points instead of {n}.')
else:
	print(f'{n_actual} points')

fig, ax = plt.subplots()
ax.scatter(kgbg.coordinates[:, 0], kgbg.coordinates[:, 1])
plt.show()

kgbg.build()
kgbg.write(f'brush_gd{density}_N{chain_length}_{x}x{y}x{z}.pos.gz')
