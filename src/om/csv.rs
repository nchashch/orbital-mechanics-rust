extern crate nalgebra;
use nalgebra::*;
use om::koe::*;
use om::cb::*;
use tick::*;
use std::rc::*;

/// Cartesian State Vectors.
#[derive(Clone)]
pub struct CSV {
    /// Position.
    pub r: Vec3<f64>,
    /// Velocity.
    pub v: Vec3<f64>,
    /// Reference to the central body that this object orbits.
    pub cb: Rc<CB>,
}

impl Tick for CSV {
    fn tick(&self, dt: f64) -> Self {
        CSV {
            r: self.r + self.v * dt,
            cb: self.cb.clone(),
            ..*self
        }
    }
}

impl CSV {
    /// Construct CSV from position and velocity.
    pub fn new(r: Vec3<f64>, v: Vec3<f64>, cb: Rc<CB>) -> CSV {
        CSV {
            r: r,
            v: v,
            cb: cb,
        }
    }

    /// Construct CSV from KOE.
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
