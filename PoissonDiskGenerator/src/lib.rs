use std::f64::consts::SQRT_2;

use itertools::Itertools;
use pyo3::prelude::*;
use rand::prelude::*;


/// Generate a 2D Poisson-disk point set using cell list accelerated dart throwing.
/// * `n` - Number of points
/// * `bead_size` - Minimum distance between points
/// * `size` - 2-tuple corresponding to the domain size
/// * `seed` - Optional seed number for PRNG
/// * `max_iter` - Optional iteration limit
/// Returns a list of tuples corresponding to the point coordinates.
#[pyfunction]
fn gen(n: usize, bead_size: f64, size: (f64, f64), seed: Option<u64>, max_iter: Option<i32>) -> PyResult<Vec<(f64, f64)>> {
    // Create a grid of square background cells. The cells should be as large as possible but while still being fully
    // covered by the size of a point within it.
    let cell_size = bead_size / SQRT_2;

    // Round desired domain size up to nearest multiple of cell size
    let max = (
        (size.0 / cell_size).ceil() as usize,
        (size.1 / cell_size).ceil() as usize
    );

    // Flat list of indices of cells. Coordinate of bottom left corner = (i*cell_size, j*cell_size)
    let mut active_cells: Vec<(usize, usize)> = (0..max.0).cartesian_product(0..max.1).collect();
    // 2D cell list
    let mut cell_list: Vec<Vec<Option<(f64, f64)>>> = vec![vec![None; max.1]; max.0];
    // List of point coordinates
    let mut coordinates: Vec<(f64, f64)> = Vec::new();

    // Initialize PRNG
    let mut rng = match seed {
        Some(x) => SmallRng::seed_from_u64(x),
        None => SmallRng::from_entropy(),
    };

    /// Check for overlap in the cells neighboring the current cell.
    /// * `cell_list` - The cell list
    /// * `central_cell` - Indices of the central cell
    /// * `a` - Coordinates of the point
    /// Returns true if point overlaps, false if not.
    let check_overlap = |cell_list: &Vec<Vec<Option<(f64, f64)>>>, central_cell: (usize, usize), a: (f64, f64)| -> bool {
        // Check 20 neighbors (5x5 excluding corners and center)
        for (i, j) in [(-2, -1), (-2, 0), (-2, 1), (-1, -2), (-1, -1), (-1, 0), (-1, 1), (-1, 2), (0, -2), (0, -1), (0, 1), (0, 2), (1, -2), (1, -1), (1, 0), (1, 1), (1, 2), (2, -1), (2, 0), (2, 1)] {
            let neigh_cell = (
                central_cell.0 as i32 + i,
                central_cell.1 as i32 + j,
            );
            if (neigh_cell.0 > 0) & ((neigh_cell.0 as usize) < max.0 - 1) &&
               (neigh_cell.1 > 0) & ((neigh_cell.1 as usize) < max.1 - 1) {
                if cell_list[neigh_cell.0 as usize][neigh_cell.1 as usize]
                    .is_some_and(|b| (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) < bead_size.powi(2)) {
                    return true
                }
            }
        }
        return false
    };

    /// Check if the point is out-of-bounds, which can happen because the domain the algorithm runs on is slightly
    /// larger than the requested domain because of rasterising errors.
    /// * `point` - Coordinates of the point
    /// Returs true if point is out-of-bounds, false if not.
    let check_oob = |point: (f64, f64)| -> bool {
        point.0 > size.0 || point.1 > size.1
    };

    let max_iter = max_iter.unwrap_or(1000);
    for _ in 0..max_iter {
        if active_cells.len() == 0 {
            println!("Maximal point set reached.");
            break
        }

        // Choose random active cell
        let cell_id = rng.gen_range(0..active_cells.len());
        let cell = active_cells[cell_id];

        // Throw a dart
        let point = (
            (cell.0 as f64 + rng.gen::<f64>()) * cell_size,
            (cell.1 as f64 + rng.gen::<f64>()) * cell_size,
        );

        // Check whether point overlaps in neighbouring cells or is OOB
        if !check_overlap(&cell_list, cell, point) && !check_oob(point) {
            coordinates.push(point);
            cell_list[cell.0][cell.1] = Some(point);
            active_cells.swap_remove(cell_id);
        }

        if coordinates.len() >= n {
            break;
        }
    }

    println!("Generated {} points", coordinates.len());
    Ok(coordinates)
}

#[pymodule]
fn PoissonDiskGenerator(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gen, m)?)?;
    Ok(())
}
