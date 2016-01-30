extern crate nalgebra;
use nalgebra::*;
use om::koe::*;
use om::cb::*;
use tick::*;
use push::*;
use std::rc::*;

/// #Cartesian State Vectors
/// This structure represents an orbit using a
/// radius vector and a velocity vector.
/// It holds a reference to the central body.
#[derive(Clone)]
pub struct CSV {
    /// Radius vector.
    pub r: Vec3<f64>,
    /// Velocity.
    pub v: Vec3<f64>,
    /// Reference to the central body.
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

impl Push for CSV {
    fn push(&self, dv: Vec3<f64>) -> Self {
        CSV {
            v: self.v + dv,
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
        // Mean anomaly
        let m0 = koe.m0;
        // Number of iterations for newton_raphson
        let iterations = 10;
        // Eccentric anomaly
        let ea = CSV::newton_raphson(&m0, &koe.e, &iterations);
        // True anomaly
        let ta = 2.0*((1.0+koe.e).sqrt()*(ea/2.0).sin())
            .atan2((1.0-koe.e).sqrt()*(ea/2.0).cos());
        // Distance to the center of the central body
        let dist = koe.a*(1.0-koe.e*ea.cos());
        // Radius vector in i, j plane
        let mut r = (koe.cb.i*ta.cos() + koe.cb.j*ta.sin()) * dist;
        // Velocity in i, j plane
        let mut v = (koe.cb.i*(-ea.sin()) +
                    koe.cb.j*((1.0-koe.e.powf(2.0)).sqrt()*ea.cos())) * ((koe.cb.mu*koe.a).sqrt()/dist);
        // Radius vector in orbital plane
        r = koe.rot.transform(&r);
        // Velocity in orbital plane
        v = koe.rot.transform(&v);
        CSV::new(r, v, koe.cb.clone())
    }

    // Function that numerically solves Kepler's equation
    fn newton_raphson(m0: &f64, e: &f64, iterations: &i32) -> f64 {
        let mut ea = m0.clone();
        for _ in 0..*iterations {
            ea  = ea - (ea - e*ea.sin() - m0)/(1.0 - e*ea.cos());
        }
        ea
    }
}
