use ndarray::{Array1, Array, arr1};
use std::f64::consts::PI;
use approx::assert_relative_eq;

pub struct Grid {
    pub npts: i32,
    pub radius: f64,
    dr: f64,
    dk: f64,
    pub ri: Array1<f64>,
    pub ki: Array1<f64>,
}

impl Grid {
    pub fn new(npts: i32, dr: f64) -> Grid {

        let radius: f64 = npts as f64 * dr;
        let dk: f64 = 2.0 * PI / (2.0 * npts as f64 * dr);
        let ri: Array1<f64> = Array::range(0.5, npts as f64, 1.0) * dr;
        let ki: Array1<f64> = Array::range(0.5, npts as f64, 1.0) * dk;

        Grid {
            npts: npts,
            radius: radius,
            dr: dr,
            dk: dk,
            ri: ri,
            ki: ki,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_gen() {
        let grid = Grid::new(10, 1.0);
        //check length of grid is as expected
        assert_eq!(grid.ri.len(), 10);
        //radius calculation check
        assert_relative_eq!(grid.radius, 10.0, epsilon=1E-8);
        // dk calculation check
        assert_relative_eq!(grid.dk, 0.3141592653589793, epsilon=1E-8);
        // r-grid check - values from pyRISM
        assert_relative_eq!(grid.ri, arr1(&[0.5, 1.5, 2.5,
                                    3.5, 4.5, 5.5,
                                    6.5, 7.5, 8.5,
                                    9.5]), epsilon=1E-8);
        // k-grid check - values from pyRISM
        assert_relative_eq!(grid.ki, arr1(&[0.15707963, 0.4712389, 0.78539816,
                                            1.09955743, 1.41371669, 1.72787596,
                                            2.04203522, 2.35619449, 2.67035376,
                                            2.98451302]), epsilon=1E-8);
    }
}
