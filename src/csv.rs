extern crate nalgebra;
use nalgebra::*;
use koe::*;
use cb::*;
use std::rc::*;

/* Cartesian State Vectors */
#[derive(Clone)]
pub struct CSV {
    pub r: Vec3<f64>,
    pub v: Vec3<f64>,
    pub cb: Rc<CB>,
}

impl CSV {
    pub fn new(r: Vec3<f64>, v: Vec3<f64>, cb: Rc<CB>) -> CSV {
        CSV { r: r, v: v, cb: cb }
    }

    pub fn from_koe(koe: KOE) -> CSV {
        let m0 = koe.m0;
        let iterations = 10;
        let ea = CSV::newton_raphson(&m0, &koe.e, &iterations);
        let ta = 2.0*((1.0+koe.e).sqrt()*(ea/2.0).sin())
            .atan2((1.0-koe.e).sqrt()*(ea/2.0).cos());
        let dist = koe.a*(1.0-koe.e*ea.cos());
        let mut r = (koe.cb.i*ta.cos() + koe.cb.j*ta.sin()) * dist;
        let mut v = (koe.cb.i*(-ea.sin()) +
                    koe.cb.j*((1.0-koe.e.powf(2.0)).sqrt()*ea.cos())) * ((koe.cb.mu*koe.a).sqrt()/dist);
        r = koe.rot.transform(&r);
        v = koe.rot.transform(&v);
        CSV::new(r, v, koe.cb.clone())
    }

    fn newton_raphson(m0: &f64, e: &f64, iterations: &i32) -> f64 {
        let mut ea = m0.clone();
        for _ in 0..*iterations {
            ea  = ea - (ea - e*ea.sin() - m0)/(1.0 - e*ea.cos());
        }
        ea
    }
}
