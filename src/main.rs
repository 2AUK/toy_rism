use toml;
use rustdct::DctPlanner;
use ndarray::{Array1};
use std::io;
use std::sync::Arc;
use std::f64::consts::PI;
use csv;
use toy_rism::grid::Grid;
use toy_rism::potential::LJPotential;
use toy_rism::util::write_csv;
use toy_rism::{initialise, solve};

static ARGON: &str = r#"[system]
temp = 85
kT = 1.0
kU = 0.00198720414667
charge_coeff = 167101.0
npts = 512
dr = 0.01
lam = 10

[params]
potential = "LJ"
closure = "HNC"
IE = "XRISM"
solver = "Ng"
picard_damping = 0.001
itermax = 10000
tol = 1E-12
diel = 1.5053
adbcor = 0.5

[solvent]
nsv = 1
nspv = 1

[solvent.species.argon]
dens = 0.021017479720736955
ns = 1
"Ar" = [
    [120.0, 3.4, 0.0],
    [0.00000000e+00, 0.00000000e+00, 0.00000000e+00]
]
"#;

static N2: &str = r#"[system]
temp = 72
kT = 1.0
kU = 0.00198720414667
charge_coeff = 167101.0
npts = 2048
dr = 20.48
lam = 10

[params]
potential = "LJ"
closure = "PY"
IE = "XRISM"
solver = "Ng"
picard_damping = 0.001
itermax = 10000
tol = 1E-7

[solvent]
nsv = 2
nspv = 1

[solvent.species.N2]
dens = 0.01867
ns = 2
"N1" = [
    [44.0, 3.341, 0.0],
    [0.00000000e+00, 0.00000000e+00, 0.00000000e+00]
]

"N2" = [
    [44.0, 3.341, 0.0],
    [1.10000000e+00, 0.00000000e+00, 0.00000000e+00]
]"#;

fn main() {
    //how do use toml lole
    //println!("{}", argon_toml);
    //let param_array: Array1<f64> = argon_toml["solvent"]["argon"]["Ar"][0].as_array()
    //    .unwrap()
    //    .iter()
    //    .map(|x| x.as_float().unwrap())
    //    .collect();
    //let sum: f64 = param_array.iter().sum();
    //println!("{:?}", param_array);
    //println!("{} + {} + {} = {}", param_array[0], param_array[1], param_array[2], sum);

    // let grid = Grid::new(16384, 0.001);
    // let pot = LJPotential::new(85.0, 120.0, 3.4, &grid);
    // let r2k_fac = 2.0 * PI * grid.dr;
    // let k2r_fac = grid.dk / (PI * PI);
    // let mut planner = DctPlanner::new();
    // let dst4 = planner.plan_dst4(grid.npts as usize);
    // let mayer_f = &pot.func.mapv(|a| f64::exp(-a)) - 1.0;
    // let kspace = Grid::dfbt(&grid.ri, &grid.ki, &mayer_f, r2k_fac, Arc::clone(&dst4));
    // let rspace = Grid::dfbt(&grid.ki, &grid.ri,  &kspace, k2r_fac, Arc::clone(&dst4));


    //println!("{:10.3e}\n\n{:10.3e}\n\n{:10.3e}", &mayer_f, kspace, &rspace);

    //write_csv(&grid.ri, &rspace, "arLJ.csv").expect("Could not create file!");
    let argon = initialise(&ARGON);
    println!("{:?}", argon);
    let n2 = initialise(&N2);
    println!("{:?}", n2);
    solve(n2);
}
