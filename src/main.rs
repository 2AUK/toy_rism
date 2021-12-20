use toml;
use ndarray::{Array1};
use std::io;
use csv;
use toy_rism::grid::Grid;
use toy_rism::potential::LJPotential;
use toy_rism::util::write_csv;

static ARGON: &str = r#"[system]
temp = 85
kT = 1.0
kU = 0.00198720414667
charge_coeff = 167101.0
npts = 512
radius = 20.48
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

[solvent.argon]
dens = 0.021017479720736955
ns = 1
"Ar" = [
    [120.0, 3.4, 0.0],
    [0.00000000e+00, 0.00000000e+00, 0.00000000e+00]
]"#;

static N2: &str = r#"[system]
temp = 72
kT = 1.0
charge_coeff = 167101.0
npts = 2048
radius = 20.48
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

[solvent.N2]
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
    let argon_toml: toml::Value = toml::from_str(&ARGON).unwrap();
    let n2_toml: toml::Value = toml::from_str(&N2).unwrap();

    //how do use toml lole
    //println!("{}", argon_toml);
    let param_array: Array1<f64> = argon_toml["solvent"]["argon"]["Ar"][0].as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_float().unwrap())
        .collect();
    let sum: f64 = param_array.iter().sum();
    //println!("{:?}", param_array);
    //println!("{} + {} + {} = {}", param_array[0], param_array[1], param_array[2], sum);

    let grid = Grid::new(16384*10, 0.0001);
    let pot = LJPotential::new(85.0, 120.0, 3.4, &grid);
    let u_LJ = &pot.func;

    write_csv(&grid.ri, &u_LJ, "arLJ.csv").expect("Could not create file!");
}
