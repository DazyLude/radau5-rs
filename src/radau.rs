// This code is a copy(cat) of Fortran code, located at "http://www.unige.ch/~hairer/prog/stiff/radau5.f"
// Since it is a redistribution of source code of sorts, I'm obliged to retain this copyright notice: http://www.unige.ch/~hairer/prog/licence.txt
// It was written for educational purposes
// take care, as there is no warranty that it works correctly :P
use std::f64;

pub mod radau5 {

    pub struct Radau5Field {
        n: usize,
        fcn: fn(f64, &Vec<f64>) -> Vec<f64>,
        x: f64,
        y: Vec<f64>,
        x_end: f64,
        h: f64,
        r_tol: Vec<f64>,
        a_tol: Vec<f64>,
        i_tol: bool,
        // jac, i_jac, ml_jac, mu_jac,
        // mas, i_mas, ml_mas, mu_mas,
        // sol_out, i_out,
        work: Vec<f64>,
        l_work: usize,
        i_work: Vec<i64>,
        l_i_work: usize,
        u_round: f64,
        n_fcn: i64,
        n_jac: i64,
        n_step: i64,
        n_accpt: i64,
        n_rejct: i64,
        n_dec: i64,
        n_sol: i64,
        arret: bool,
        exp_m: f64,
        quot: f64,
        n_max: i64,
        n_it: i64,
        start_n: bool,
        nind1: i64,
        nind2: i64,
        nind3: i64,
    }

    impl Radau5Field {
        pub fn new() -> Radau5Field {
            Radau5Field {
                n: 1,
                fcn: |x, y| vec![x + y[0]],
                x: 0.,
                y: vec![0.],
                x_end: 1.,
                h: 0.,
                r_tol: vec![0.],
                a_tol: vec![0.],
                i_tol: false,
                work: vec![0.; 20],
                l_work: 1,
                i_work: vec![0; 20],
                l_i_work: 1,
                u_round: 0.,
                n_fcn: 0,
                n_jac: 0,
                n_step: 0,
                n_accpt: 0,
                n_rejct: 0,
                n_dec: 0,
                n_sol: 0,
                arret: false,
                exp_m: 0.,
                quot: 0.,
                n_max: 0,
                n_it: 0,
                start_n: false,
                nind1: 0,
                nind2: 0,
                nind3: 0,
            }
        }

        pub fn init(
            mut self,
            n: usize,
            fcn: fn(f64, &Vec<f64>) -> Vec<f64>,
            x: f64,
            y: &Vec<f64>,
            x_end: f64,
            h: Option<f64>,
            r_tol: Option<Vec<f64>>,
            a_tol: Option<Vec<f64>>,
            i_tol: Option<bool>,
            // jac, i_jac, ml_jac, mu_jac,
            // mas, i_mas, ml_mas, mu_mas,
            // sol_out, i_out,
            work: Option<Vec<f64>>,
            l_work: Option<usize>,
            i_work: Option<Vec<i64>>,
            l_i_work: Option<usize>,
            // r_par, i_par, idid
        ) {
            // Declarations and setting of the params
            self.n_fcn = 0;
            self.n_jac = 0;
            self.n_step = 0;
            self.n_accpt = 0;
            self.n_rejct = 0;
            self.n_dec = 0;
            self.n_sol = 0;
            self.arret = false;

            match n == y.len() {
                true => self.n = n,
                _ => println!(
                    "dimensions of initial y vector are not equal to explicitly provided dimensions of system"
                ),
            }
            self.work = match work {
                None => vec![0.; 20],
                Some(x) => x.clone(),
            };
            self.i_work = match i_work {
                None => vec![0; 20],
                Some(x) => x.clone(),
            };
            self.i_tol = match i_tol {
                None => false,
                Some(x) => x,
            };
            self.a_tol = match a_tol {
                None => vec![0.; n],
                Some(x) => x.clone(),
            };
            self.r_tol = match r_tol {
                None => vec![0.; n],
                Some(x) => x.clone(),
            };

            // setting u_round
            if self.work[0] == 0. {
                self.u_round = 1e-16_f64;
            } else {
                self.u_round = self.work[0];
                if (self.u_round <= 1.0e-19) || (self.u_round >= 1.) {
                    println!("Coefficients have 20 digits, u_round = {}", self.u_round);
                    self.arret = true;
                }
            }
            // setting tolerances
            self.exp_m = 2.0_f64 / 3.0_f64;
            // i_tol defines if tolerances are the same for all variables or individual
            if !self.i_tol {
                if (self.a_tol[0] <= 0.) || (self.r_tol[0] <= 10. * self.u_round) {
                    println!("Tolerances are too small");
                    self.arret = true;
                } else {
                    self.quot = self.a_tol[0] / self.r_tol[0];
                    self.r_tol[0] = (0.1 * self.r_tol[0]).powf(self.exp_m);
                    self.a_tol[0] = self.r_tol[0] * self.quot;
                }
            } else {
                for i in 0..self.n {
                    if (self.a_tol[i] <= 0.) || (self.r_tol[i] <= 10. * self.u_round) {
                        println!("Tolerances for {}-th element are too small", i);
                        self.arret = true;
                    } else {
                        self.quot = self.a_tol[i] / self.r_tol[i];
                        self.r_tol[i] = (0.1 * self.r_tol[0]).powf(self.exp_m);
                        self.a_tol[i] = self.r_tol[i] * self.quot;
                    }
                }
            }
            // setting maximal number of steps
            if self.i_work[1] == 0 {
                self.n_max = 100000
            } else if self.i_work[1] < 0 {
                println!("wrong input of i_work[1]: {}", self.i_work[1]);
                self.arret = true;
            } else {
                self.n_max = self.i_work[1];
            }
            // setting maximal number of newthon iterations
            match self.i_work[2].cmp(&0) {
                core::cmp::Ordering::Equal => self.n_it = 7,
                core::cmp::Ordering::Greater => self.n_it = self.i_work[2],
                core::cmp::Ordering::Less => {
                    println!("curious input of i_work[2]: {}", self.i_work[2]);
                    self.arret = true;
                }
            }
            // some weird switch for newthon idk he ded
            self.start_n = self.i_work[3] != 0;
            // parameters for differential-algebraic components
            if self.i_work[4] == 0 {
                self.nind1 = self.n as i64;
            } else {
                self.nind1 = self.i_work[4];
            }
            if (self.nind1 + self.i_work[5] + self.i_work[6]) != self.n as i64 {
                println!(
                    "curious input of i_work[4..6]: {}, {}, {}",
                    self.i_work[4], self.i_work[5], self.i_work[6]
                );
                self.arret = true;
            } else {
                self.nind2 = self.i_work[5];
                self.nind3 = self.i_work[6];
            }
            //
        }

        fn radau5_core() {}

        //sol_out dummy
        fn newthon() {}
    }
}
