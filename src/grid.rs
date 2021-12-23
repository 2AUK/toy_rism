use ndarray::{Array1, Array, arr1};
use std::f64::consts::PI;
use std::sync::Arc;
use approx::assert_relative_eq;
use rustdct::{DctPlanner, TransformType4};


pub struct Grid {
    pub npts: u32,
    pub radius: f64,
    pub dr: f64,
    pub dk: f64,
    pub ri: Array1<f64>,
    pub ki: Array1<f64>,
}

impl Grid {
    pub fn new(npts: u32, dr: f64) -> Grid {

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

    pub fn dfbt(grid1: &Array1<f64>,
                grid2: &Array1<f64>,
                func: &Array1<f64>,
                fac: f64,
                dst4: Arc<dyn TransformType4<f64>>) -> Array1<f64> {

        let mut buffer: Vec<f64> = (grid1 * func).to_vec();

        dst4.process_dst4(&mut buffer);
        fac * Array1::from_vec(buffer) / grid2
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

    #[test]
    fn transforms() {
        let grid = Grid::new(16384, 0.0001);
        let func: Array1<f64> =  (1.0 / 85.0) * 4.0 * 120.0 *
            ( (3.4 / &grid.ri).mapv(|a| a.powf(12.0)) -
               ( (3.4 / &grid.ri).mapv(|a| a.powf(6.0))  ));
        let mayer_f = func.mapv(|a| f64::exp(-a)) - 1.0;
        let mut planner = DctPlanner::new();
        let dst4 = planner.plan_dst4(grid.npts as usize);
        let r2k_fac = 2.0 * PI * grid.dr;
        let k2r_fac = grid.dk / (PI * PI);
        let kspace = Grid::dfbt(&grid.ri, &grid.ki, &mayer_f, r2k_fac, Arc::clone(&dst4));
        let rspace = Grid::dfbt(&grid.ki, &grid.ri,  &kspace, k2r_fac, Arc::clone(&dst4));
       

        assert_relative_eq!(mayer_f, rspace, epsilon=1E-5);
    }
}
