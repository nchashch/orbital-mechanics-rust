extern crate nalgebra;
use nalgebra::*;
use cb::*;
use csv::*;
use std::rc::*;
use tick::*;
use std::f64::consts::*;

/// Keplerian Orbital Elements.
#[derive(Clone)]
pub struct KOE {
    /// Semi Major Axis.
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination.
    pub inc: f64,
    /// Longitude of Ascending Node.
    pub lan: f64,
    /// Argument of Periapsis.
    pub ap: f64,
    /// Mean anomaly.
    pub m0: f64,
    /// A matrix that transforms a vector lying in i, j plane into a corresponding vector lying in the orbital plane
    /// defined by inc and lan.
    pub rot: Rot3<f64>,
    /// Reference to the central body that this object orbits.
    pub cb: Rc<CB>,
}

impl Tick for KOE {
    fn tick(&self, dt: f64) -> Self {
        let n = (self.cb.mu/self.a.powf(3.0)).sqrt();
        KOE {
            a: self.a,
            e: self.e,
            inc: self.inc,
            lan: self.lan,
            ap: self.ap,
            m0: self.m0 + n * dt,
            rot: self.rot,
            cb: self.cb.clone(),
        }
    }
}

impl KOE {
    /// Construct KOE from orbital elements.
    pub fn new(a: f64, e: f64, inc: f64, lan: f64, ap: f64, m0: f64, cb: Rc<CB>) -> KOE {
        let mut rot = Rot3::new_identity(3);
        if approx_eq(&inc, &0.0) {
            if !approx_eq(&e, &0.0) {
                rot = Rot3::new(cb.k * ap);
            }
        } else {
            let lan_axisangle = cb.k * lan;
            let inc_axisangle = cb.i * inc;
            let ap_axisangle = cb.k * ap;
            rot = Rot3::new(lan_axisangle);
            rot = rot.prepend_rotation(&inc_axisangle);
            if !approx_eq(&e, &0.0) {
                rot = rot.prepend_rotation(&ap_axisangle);
            }
        }
        KOE {
            a: a,
            e: e,
            inc: inc,
            lan: lan,
            ap: ap,
            m0: m0,
            rot: rot,
            cb: cb,
        }
    }

    /// Construct KOE from CSV.
    pub fn from_csv(csv: CSV) -> KOE {
        let r = csv.r;
        let v = csv.v;

        let mu = &csv.cb.mu;

        let h = cross(&r, &v);
        let e = cross(&v, &h)/ *mu - normalize(&r);
        let n = cross(&csv.cb.k, &h);

        let cos_inc = dot(&h, &csv.cb.k)/(norm(&h));
        let inc = if approx_eq(&cos_inc, &1.0) {    
            0.0
        } else {
            cos_inc.acos()
        };

        let es = norm(&e);
        
        let mut lan = if dot(&n, &csv.cb.j) >= 0.0 {
            (dot(&n, &csv.cb.i)/norm(&n)).acos()
        } else {
            2.0*PI - (dot(&n, &csv.cb.i)/norm(&n)).acos()
        };
        
        let right = cross(&h, &n);
        let mut ap = if dot(&e, &right) >= 0.0 {
            (dot(&n, &e)/(norm(&n)*norm(&e))).acos()
        } else {
            2.0*PI - (dot(&n, &e)/(norm(&n)*norm(&e))).acos()
        };
        
        if approx_eq(&es, &0.0) {
            ap = 0.0;
        }
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            if !approx_eq(&es, &0.0) {
                ap = if dot(&e, &csv.cb.j) >= 0.0 {
                    (dot(&csv.cb.i, &e)/norm(&e)).acos()
                } else {
                    2.0*PI - (dot(&csv.cb.i, &e)/norm(&e)).acos()
                };
            } else {
                ap = 0.0;
            }
        }

        let mut ta = if dot(&r, &v) >= 0.0 {
            (dot(&e, &r)/(norm(&e)*norm(&r))).acos()
        } else {
            2.0*PI - (dot(&e, &r)/(norm(&e)*norm(&r))).acos()
        };
        
        if approx_eq(&es, &0.0) {
            // Compute argument of latitude
            if approx_eq(&inc, &0.0) {
                ta = if dot(&csv.cb.i, &v) <= 0.0 {
                    (dot(&csv.cb.i, &r)/(norm(&csv.cb.i)*norm(&r))).acos()
                } else {
                    2.0*PI - (dot(&csv.cb.i, &r)/(norm(&csv.cb.i)*norm(&r))).acos()
                }
            } else {
                ta = if dot(&n, &v) <= 0.0 {
                    (dot(&n, &r)/(norm(&n)*norm(&r))).acos()
                } else {
                    2.0*PI - (dot(&n, &r)/(norm(&n)*norm(&r))).acos()
                }
            }
        }

        let ea = 2.0*((ta/2.0).tan()/((1.0+es)/(1.0-es)).sqrt()).atan();
        let m0 = ea - es * ea.sin();
        let a = 1.0/(2.0/norm(&r) - sqnorm(&v)/mu);
        KOE::new(a, es, inc, lan, ap, m0, csv.cb.clone())
    }
}
