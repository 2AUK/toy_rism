use ndarray::{Array3, Array2, Array1, Array, s, Axis, ArrayView};
use ndarray_stats::QuantileExt;
use ndarray_linalg::Inverse;
use serde::Deserialize;
use toml;
use rustdct::{DctPlanner};
use std::f64::consts::PI;
use std::sync::Arc;

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
    pub p: Array2<f64>,
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
            p: Array::zeros((ns1, ns2)),
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
    npts: u32,
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
    pub params: Array1<f64>,
    pub coord: Array1<f64>,
}

impl Site {
    pub fn new(params: Array1<f64>, coord: Array1<f64>) -> Site {
        Site {
            params: params,
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
            let curr_site = Site::new(Array1::from_vec(param.clone()), Array1::from_vec(coord.clone()));
            sites.push(curr_site);
        }
        let species_info = SolventInfo::new(ns, density, sites);
        species.push(species_info);
    }
    let rism = RISM::new(inp.system.npts as usize,
                         inp.solvent.nsv as usize,
                         inp.solvent.nsv as usize, inp.params, inp.system, species, None);

    rism
}

pub fn solve(input_sys: &mut RISM) {
    let grid = grid::Grid::new(input_sys.sys.npts, input_sys.sys.dr);
    let mut coords = Vec::new();
    let mut params = Vec::new();
    let mut dens = Vec::new();
    for species in input_sys.vv.iter(){
        for site_info in species.sites.iter(){
            println!("{:?}", site_info.coord.to_vec());
            println!("{:?}", species.dens);
            dens.push(species.dens);
            coords.push(site_info.coord.to_vec());
            params.push(site_info.params.to_vec());
        }
    }
    input_sys.ω = intramolecular_corr(&grid, input_sys.vv[0].ns, coords);
    input_sys.p = Array2::from_diag(&Array::from_vec(dens));
    println!("{}", input_sys.ω);
    println!("{}", input_sys.p);
    for j in 1..=(input_sys.sys.lam) as usize {
        let lam = 1.0 * j as f64 / input_sys.sys.lam as f64;
        input_sys.u = potential(&grid, input_sys.vv[0].ns, &params, &input_sys.sys, lam);
        let mut i = 0;

        while i < input_sys.params.itermax {
            let prev_c = input_sys.c.clone();
            IE(input_sys, &grid, input_sys.vv[0].ns);
            let c_A = HNC(&input_sys.u, &input_sys.t);
            input_sys.c = &prev_c + (input_sys.params.picard_damping * (c_A - &prev_c));
            let diff = &input_sys.c - &prev_c;
            let rms: f64 = (grid.dr *
                ((diff).mapv(|a| a.powf(2.0)).sum()) /
                            (ArrayView::from(input_sys.c.shape()).product() as f64)).sqrt();
            let diff_max = diff.max().unwrap();
            if rms < input_sys.params.tol {
                println!("lambda: {}\niteration: {}\nRMS: {:10.3e}\ndiff: {:10.3e}", lam, i, rms, diff_max);
                break;
            }
            if i % 100 == 0 {
                println!("lambda: {}\niteration: {}\nRMS: {:10.3e}\ndiff: {:10.3e}", lam, i, rms, diff_max);
            }
            i += 1;
        }
    }
    let gr = 1.0 + &input_sys.c + &input_sys.t;
    println!("{}", gr);
    util::write_csv(&grid.ri, &gr.slice(s![.., 0, 0]).to_owned(), "rgr.csv").unwrap();

}

fn HNC(pot: &Array3<f64>, t: &Array3<f64>) -> Array3<f64> {
    (-pot + t).mapv(f64::exp) - 1.0 - t
}

fn IE(input_sys: &mut RISM, grid: &grid::Grid, ns: usize) {
    let mut planner = DctPlanner::new();
    let dst4 = planner.plan_dst4(grid.npts as usize);
    let r2k_fac = 2.0 * PI * grid.dr;
    let k2r_fac = grid.dk / (4.0 * PI * PI);
    let I: Array2<f64> = Array::eye(ns);
    let mut ck: Array3<f64> = Array::zeros((grid.npts as usize, ns, ns));
    let mut hk: Array3<f64> = Array::zeros((grid.npts as usize, ns, ns));
    for i in 0..ns {
        for j in 0..ns {
            let c_slice = input_sys.c.slice(s![.., i, j]);
            let mut ck_slice = ck.slice_mut(s![.., i, j]);
            ck_slice.assign(&grid::Grid::dfbt(&grid.ri, &grid.ki, &c_slice.to_owned(), r2k_fac, Arc::clone(&dst4)));
        }
    }
    for ri in 0..(grid.npts as usize) {
        let w_slice = input_sys.ω.slice(s![ri, .., ..]);
        let ck_slice = ck.slice(s![ri, .., ..]);
        let mut h_slice = hk.slice_mut(s![ri, .., ..]);
        let iwcp: Array2<f64> = (&I - w_slice.dot(&ck_slice.dot(&input_sys.p))).inv().unwrap();
        let wcw: Array2<f64> = w_slice.dot(&ck_slice.dot(&w_slice));
        h_slice.assign(&iwcp.dot(&wcw));
    }
    for i in 0..ns {
        for j in 0..ns {
            let mut c_slice = input_sys.c.slice_mut(s![.., i, j]);
            let mut h_slice = input_sys.h.slice_mut(s![.., i, j]);
            let ck_slice = ck.slice(s![.., i, j]);
            let hk_slice = hk.slice(s![.., i, j]);
            c_slice.assign(&grid::Grid::dfbt(&grid.ki, &grid.ri, &ck_slice.to_owned(), k2r_fac, Arc::clone(&dst4)));
            h_slice.assign(&grid::Grid::dfbt(&grid.ki, &grid.ri, &hk_slice.to_owned(), k2r_fac, Arc::clone(&dst4)));
        }
    }
    input_sys.t = &input_sys.h - &input_sys.c;
}

fn intramolecular_corr(grid: &grid::Grid, ns: usize, coords: Vec<Vec<f64>>) -> Array3<f64> {
    let I: Array1<f64> = Array1::ones(grid.npts as usize);
    let dists: Array2<f64> = distance_matrix(ns, &coords);
    let mut w: Array3<f64> = Array::zeros((grid.npts as usize, ns, ns));
    for (i, is) in coords.iter().enumerate() {
        for (j, js) in coords.iter().enumerate() {
            let mut w_slice = w.slice_mut(s![.., i, j]);
            if dists[[i, j]] == 0.0 {
                w_slice.assign(&I);
            } else {
                w_slice.assign(&((&grid.ki * dists[[i, j]]).mapv(|a| a.sin()) / (&grid.ki * dists[[i, j]])));
            }
        }
    }
    w
}

fn potential(grid: &grid::Grid, ns: usize, LJ_param: &Vec<Vec<f64>>, sys: &System, lam: f64) -> Array3<f64> {
    let mut pot: Array3<f64> = Array::zeros((grid.npts as usize, ns, ns));
    for (i, is) in LJ_param.iter().enumerate() {
        for (j, js) in LJ_param.iter().enumerate() {
            let mut pot_slice = pot.slice_mut(s![.., i, j]);
            let tabulate_pot = potential::LJPotential::new(sys.temp, sys.kT, js[0], js[1], &grid);
            pot_slice.assign(&(tabulate_pot.func * lam));
        }
    }
    pot
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
