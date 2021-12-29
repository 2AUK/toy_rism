use ndarray::{Array3, Array2, Array1, Array, s};
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
//
//
//
// The data structures do not make any sense at the moment but its a very rough implementation
// so it doesn't really matter at the moment

#[derive(Debug)]
pub struct RISM{
    pub u: Array3<f64>,
    pub ω: Array3<f64>,
    pub h: Array3<f64>,
    pub c: Array3<f64>,
    pub t: Array3<f64>,
    pub params: Parameters,
    pub sys: System,
    pub vv: Vec<SolventInfo>,
    pub uv: Option<Vec<SoluteInfo>>,
}

impl RISM{
    pub fn new(npts: usize,
               ns1: usize,
               ns2: usize,
               params: Parameters,
               sys: System,
               vv: Vec<SolventInfo>,
               uv: Option<Vec<SoluteInfo>>
    ) -> RISM {
        RISM {
            u: Array::zeros((npts, ns1, ns2)),
            ω: Array::zeros((npts, ns1, ns2)),
            h: Array::zeros((npts, ns1, ns2)),
            c: Array::zeros((npts, ns1, ns2)),
            t: Array::zeros((npts, ns1, ns2)),
            params: params,
            sys: sys,
            vv: vv,
            uv: uv,
        }
    }
}

#[derive(Deserialize, Debug)]
pub struct problem {
    params: Parameters,
    system: System,
    solvent: Solvent,
    solute: Option<Solute>,
}

#[derive(Deserialize, Debug)]
pub struct Parameters {
    potential: String,
    closure: String,
    IE: String,
    solver: String,
    picard_damping: f64,
    itermax: i32,
    tol: f64,
}

#[derive(Deserialize, Debug)]
pub struct System {
    temp: f64,
    kT: f64,
    kU: f64,
    charge_coeff: f64,
    npts: i32,
    dr: f64,
    lam: i32,
}

#[derive(Deserialize, Debug)]
pub struct Solvent {
    nsv: i32,
    nspv: i32,
    species: toml::value::Table,
}

#[derive(Deserialize, Debug)]
pub struct Solute {
    nsv: i32,
    nspv: i32,
    species: toml::value::Table,
}

#[derive(Debug)]
pub struct SolventInfo {
    pub ns: usize,
    pub dens: f64,
    pub sites: Vec<Site>
}

#[derive(Debug)]
pub struct SoluteInfo {
    pub ns: usize,
    pub dens: f64,
    pub sites: Vec<Site>
}

#[derive(Debug)]
pub struct Site {
    pub epsilon: f64,
    pub sigma: f64,
    pub charge: f64,
    pub coord: Array1<f64>,
}

impl Site {
    pub fn new(eps: f64, sig: f64, charge: f64, coord: Array1<f64>) -> Site {
        Site {
            epsilon: eps,
            sigma: sig,
            charge: charge,
            coord: coord,
        }
    }
}

impl SolventInfo {
    pub fn new(ns: usize, dens: f64, sites: Vec<Site>) -> SolventInfo {
        SolventInfo {
            ns: ns,
            dens: dens,
            sites: sites,
        }
    }
}

pub fn initialise(input_file: &str) -> RISM {
    let inp: problem = toml::from_str(input_file).unwrap();
    let site_info = inp.solvent.species.iter().take(inp.solvent.nspv as usize);
    let mut species = Vec::new();
    for i in site_info.into_iter() {
        let species_name = &i.0;
        let mut params = i.1.as_table().unwrap().clone();
        let ns = params.remove("ns").unwrap().as_integer().unwrap() as usize;
        let density = params.remove("dens").unwrap().as_float().unwrap();
        let mut sites = Vec::new();
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
            let curr_site = Site::new(param[0], param[1], param[2], Array1::from_vec(coord.clone()));
            sites.push(curr_site);
        }
        let species_info = SolventInfo::new(ns, density, sites);
        species.push(species_info);
    }
    let rism = RISM::new(inp.system.npts as usize,
                         inp.solvent.nspv as usize,
                         inp.solvent.nspv as usize, inp.params, inp.system, species, None);

    rism
}

pub fn solve(input_sys: RISM) -> RISM {
    let mut coords = Vec::new();
    for species in input_sys.vv.iter(){
        println!("{:?}", species);
        for site_info in species.sites.iter(){
            println!("{:?}", site_info.coord.to_vec());
            coords.push(site_info.coord.to_vec());
        }
    }
    let w = intramolecular_corr(input_sys.grid, input_sys.vv[0].ns, coords);
    input_sys
}


fn intramolecular_corr(grid: grid::Grid, ns: usize, coords: Vec<Vec<f64>>) -> Array3<f64> {
    let mut I: Array1<f64> = Array1::ones(grid.npts as usize);
    let dists: Array2<f64> = distance_matrix(ns, &coords);
    let mut w: Array3<f64> = Array::zeros((grid.npts as usize, ns, ns));
    for (i, is) in coords.iter().enumerate() {
        for (j, js) in coords.iter().enumerate() {
            let mut w_slice = w.slice_mut(s![.., i, j]);
            if dists[[i, j]] == 0.0 {
                w_slice = I.view_mut();
            } else {
                w_slice = ((&grid.ki * dists[[i, j]]).mapv(|a| a.sin()) / (&grid.ki * dists[[i, j]])).view_mut();
            }
        }
    }
    w
}

fn distance_matrix(ns: usize, coords: &Vec<Vec<f64>>) -> Array2<f64> {
    let mut dist_mat = Array2::zeros((ns, ns));
    for (i, is) in coords.iter().enumerate() {
        for (j, js) in coords.iter().enumerate() {
            dist_mat[[i, j]] = dist(Array::from_vec(is.clone()), Array::from_vec(js.clone()));
        }
    }
    dist_mat
}

fn dist(a: Array1<f64>, b: Array1<f64>) -> f64 {
    let diff = a - b;
    let result = (&diff * &diff).sum();
    result.sqrt()
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
