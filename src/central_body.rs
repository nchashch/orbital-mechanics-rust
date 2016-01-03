extern crate nalgebra as na;

use self::na::*;

pub struct CentralBody {
    pub mu: f64,
    pub up: Vec3<f64>,
    pub reference: Vec3<f64>,
    pub right: Vec3<f64>,
}

impl CentralBody {
    pub fn new(mu: f64, up: Vec3<f64>, reference: Vec3<f64>) -> CentralBody {
        assert!(na::approx_eq(&na::dot(&up, &reference), &0.0));
        assert!(na::approx_eq(&na::norm(&up), &1.0));
        assert!(na::approx_eq(&na::norm(&reference), &1.0));
        let right = cross(&up, &reference);
        CentralBody {
            mu: mu,
            up: up,
            right: right,
            reference: reference,
        }
    }
    
    pub fn approx_eq(&self, other: &CentralBody) -> bool {
        na::approx_eq(&self.mu, &other.mu) &&
        self.up.approx_eq(&other.up) &&
        self.reference.approx_eq(&other.reference) &&
        self.right.approx_eq(&other.right)
    }
}
