extern crate nalgebra;

use self::nalgebra::*;
use central_body::*;
use koe::*;
use std::f64::consts::*;

// Cartesian State Vectors
#[derive(Clone, Debug)]
pub struct CSV<'a> {
    pub r: Vec3<f64>,
    pub v: Vec3<f64>,
    pub cb: &'a CentralBody,
}

impl<'a> CSV<'a> {
    pub fn new(r: Vec3<f64>, v: Vec3<f64>, cb: &'a CentralBody ) -> CSV {
        CSV { r: r, v: v, cb: cb }
    }

    pub fn approx_eq(&self, other: &CSV) -> bool {
        let eps = 1.0e-4;
        self.r.approx_eq_eps(&other.r, &eps) &&
        self.v.approx_eq_eps(&other.v, &eps) &&
        self.cb.approx_eq(&other.cb)
    }

    pub fn approx_eq_eps(&self, other: &CSV, eps: &f64) -> bool {
        self.r.approx_eq_eps(&other.r, eps) &&
        self.v.approx_eq_eps(&other.v, eps) &&
        self.cb.approx_eq(&other.cb)
    }
    
    pub fn to_koe(&self) -> KOE {
        let r = self.r.clone();
        let v = self.v.clone();

        let h = cross(&r, &v);
        let e = cross(&v, &h)/self.cb.mu - normalize(&r);
        let n = cross(&self.cb.up, &h);

        let cos_inc = dot(&h, &self.cb.up)/(norm(&h));
        let inc = if approx_eq(&cos_inc, &1.0) {    
            0.0
        } else {
            cos_inc.acos()
        };

        let es = norm(&e);
        
        let mut lan = if dot(&n, &self.cb.right) >= 0.0 {
            (dot(&n, &self.cb.reference)/norm(&n)).acos()
        } else {
            2.0*PI - (dot(&n, &self.cb.reference)/norm(&n)).acos()
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
        let mut lp = 0.0;
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            ap = 0.0;
            if !approx_eq(&es, &0.0) {
                lp = if dot(&e, &self.cb.right) >= 0.0 {
                    (dot(&self.cb.reference, &e)/norm(&e)).acos()
                } else {
                    2.0*PI - (dot(&self.cb.reference, &e)/norm(&e)).acos()
                };
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
                ta = if dot(&self.cb.reference, &v) <= 0.0 {
                    (dot(&self.cb.reference, &r)/(norm(&self.cb.reference)*norm(&r))).acos()
                } else {
                    2.0*PI - (dot(&self.cb.reference, &r)/(norm(&self.cb.reference)*norm(&r))).acos()
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
        let a = 1.0/(2.0/norm(&r) - sqnorm(&v)/&self.cb.mu);
        KOE::new(a, es, inc, lan, ap, lp, m0, self.cb)
    }
}