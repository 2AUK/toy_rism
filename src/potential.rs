use ndarray::Array1;
use crate::grid::Grid;
use ndarray_stats::QuantileExt;
use approx::assert_relative_eq;

// Only implementing LJ potential at the moment
pub struct LJPotential<'a> {
    temp: f64,
    epsilon: f64,
    sigma: f64,
    pub grid: &'a Grid,
    pub func: Array1<f64>,
}

impl<'a> LJPotential<'a>{
    pub fn new(temp: f64, epsilon: f64, sigma: f64, grid: &'a Grid) -> LJPotential {

        let result = (1.0 / temp) * 4.0 * epsilon *
            ( (sigma / &grid.ri).mapv(|a| a.powf(12.0)) -
               ( (sigma / &grid.ri).mapv(|a| a.powf(6.0))  ));

        LJPotential {
            temp: temp,
            epsilon: epsilon,
            sigma: sigma,
            grid: grid,
            func: result,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_argon() {
        // Initialise potential struct
        // and calculate LJ over grid
        let grid = Grid::new(163840, 0.0001);
        // parameters for liquid argon at 85 K
        let pot = LJPotential::new(85.0, 120.0, 3.4, &grid);

        // Compute analytical r_min of LJ
        let rmin = 3.4 * f64::powf(2.0, 1.0 / 6.0);
        let min_pot = (1.0 / 85.0) * 4.0 * 120.0 *
            ( (3.4 / rmin).powf(12.0) -
               ( (3.4 / rmin).powf(6.0)) );

        // find minimum of tabulated LJ
        let min = &pot.func.min().unwrap();

        // Check analytical potential minimum == grid minimum
        // for a relatively fine grid (0.00025 dr might be too
        // small a spacing for an actual calculation)
        assert_relative_eq!(min_pot, min, epsilon=1E-8);
    }
}
