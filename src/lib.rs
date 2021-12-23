use ndarray::{Array3, Array2, Array1, Array};
use toml;

pub mod grid;
pub mod potential;
pub mod util;

// A fun way to implement:
// RISM_planner - similar to FFTW planner??
// Pass enums as flags to determine options
// Options derived from toml
// hmmmm
//
// lets not fret about that now and just hardcode Ar and N2

pub struct RISM{
    pub u: Array3<f64>,
    pub ω: Array3<f64>,
    pub h: Array3<f64>,
    pub c: Array3<f64>,
    pub t: Array3<f64>
}

impl RISM{
    pub fn new(npts: usize, ns1: usize, ns2: usize) -> RISM {
        RISM {
            u: Array::zeros((npts, ns1, ns2)),
            ω: Array::zeros((npts, ns1, ns2)),
            h: Array::zeros((npts, ns1, ns2)),
            c: Array::zeros((npts, ns1, ns2)),
            t: Array::zeros((npts, ns1, ns2))
        }
    }
}

pub fn initialise(input_toml: toml::Value) {
    let params = input_toml["params"].as_table().unwrap();
    let system = input_toml["system"].as_table().unwrap();
    let solvent = input_toml["solvent"].as_table().unwrap();
    //let solute = input_toml["solute"].as_table().expect("Could not find solute system");
    println!("{:?}", solvent);
    //let nsv = input_toml["solvent"].as_table().unwrap().remove("nsv").unwrap().as_integer().unwrap();
    //let nspv = input_toml["solvent"]["nspv"].as_integer().unwrap();
    //for entry in input_toml["solvent"].as_table().unwrap().iter() {
    //    println!("{:?}", entry);
    //}
}

pub fn solve(input_sys: RISM) -> RISM {
    input_sys
}

fn intramolecular_corr(grid: grid::Grid) { //-> Array3<f64> {

}

fn distance_matrix(coords: Array2<f64>) { //-> Array2<f64> {

}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
