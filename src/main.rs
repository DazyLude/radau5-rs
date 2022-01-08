#![allow(dead_code)]
#![allow(unused)]

pub mod radau;
// y' = x + y, solution is y = exp(x) - x - 1 with y(0) = 0
fn test_1d_function(x: f64, y: &Vec<f64>) -> Vec<f64> {
    let mut vec = Vec::new();
    vec.push(x + y[0]);
    vec
}

fn main() {
    // radau::radau5::integrate_radau5_1d(test_1d_function, &0., &0., &10.).unwrap();
    let radlad: radau::radau5::Radau5Field = radau::radau5::Radau5Field::new();

    radlad.init(
        1,
        test_1d_function,
        0.,
        &vec![0.; 1],
        10.,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    );
}
