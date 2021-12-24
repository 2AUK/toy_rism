use ndarray::{Array3, Array2, Array1, Array};
use serde::Deserialize;
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

#[derive(Deserialize, Debug)]
struct problem {
    params: Parameters,
    system: System,
    solvent: Solvent,
    solute: Option<Solute>,
}

#[derive(Deserialize, Debug)]
struct Parameters {
    potential: String,
    closure: String,
    IE: String,
    solver: String,
    picard_damping: f64,
    itermax: i32,
    tol: f64,
}

#[derive(Deserialize, Debug)]
struct System {
    temp: f64,
    kT: f64,
    kU: f64,
    charge_coeff: f64,
    npts: i32,
    dr: f64,
    lam: i32,
}

#[derive(Deserialize, Debug)]
struct Solvent {
    nsv: i32,
    nspv: i32,
    species: toml::value::Table,
}

#[derive(Deserialize, Debug)]
struct Solute {
    nsv: i32,
    nspv: i32,
    species: toml::value::Table,
}

pub fn initialise(input_file: &str) {
    let inp: problem = toml::from_str(input_file).unwrap();
    //let solute = input_toml["solute"].as_table().expect("Could not find solute system");
    println!("{:?}\n", inp);
    let site_info = inp.solvent.species.iter().take(inp.solvent.nspv as usize);
    for i in site_info.into_iter() {
        let species_name = &i.0;
        let mut params = i.1.as_table().unwrap().clone();
        let ns = params.remove("ns").unwrap().as_integer().unwrap();
        let density = params.remove("dens").unwrap().as_float().unwrap();
        //let density = i.1["dens"].as_float().unwrap();
        //let ns = i.1["ns"].as_integer().unwrap();
        println!("{}\n", species_name);
        println!("{:?}\n", density);
        println!("{:?}\n", ns);
        for j in params.iter() {
            let site_name = j.0;
            let site_info: Vec<Vec<f64>> = j.1
                .as_array()
                .unwrap()
                .iter()
                .map(|a| a.as_array()
                     .unwrap()
                     .iter()
                     .map(|a| a.as_float()
                          .unwrap())
                     .collect())
                .collect();
            let param = &site_info[0];
            let coord = &site_info[1];
            println!("{:?}", site_name);
            println!("{:?}", param);
            println!("{:?}\n", coord);
        }
    }
}
    //let nsv = input_toml["solvent"].as_table().unwrap().remove("nsv").unwrap().as_integer().unwrap();
    //let nspv = input_toml["solvent"]["nspv"].as_integer().unwrap();
    //for entry in input_toml["solvent"].as_table().unwrap().iter() {
    //    println!("{:?}", entry);
    //}

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
