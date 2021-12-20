use std::fs::File;
use std::io::prelude::*;
use csv;

pub fn write_csv<T>(a: &T, b: &T, fname: &str) -> Result<(), std::io::Error>
where
    for<'a> &'a T: IntoIterator<Item = &'a f64>
{
    let ofile = File::create(fname);

    let mut ofile = match ofile {
        Ok(file) => file,
        Err(error) => return Err(error),
    };
    ofile.write(b"#Tabulated potential\n").unwrap();
    let mut wtr = csv::Writer::from_writer(ofile);
    wtr.write_record(&["r", "f(r)"]).unwrap();
    for (&r, &fr) in a.into_iter().zip(b.into_iter()) {
        wtr.serialize((&r, &fr)).unwrap();
    }

    Ok(())
}
