extern crate nalgebra;
use self::nalgebra::*;
use central_body::*;
use csv::*;
use std::rc::*;
use tick::*;

// Keplerian Orbital Elements
#[derive(Clone, Debug)]
pub struct KOE {
    pub a: f64,
    pub e: f64,
    pub inc: f64,
    pub lan: f64,
    pub ap: f64,
    pub m0: f64, // Mean anomaly
    pub cb: Rc<CentralBody>,
    pub to_csv_rot: Rot3<f64>,
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
            cb: self.cb.clone(),
            to_csv_rot: self.to_csv_rot,
        }
    }
}

impl KOE {
    pub fn new(a: f64, e: f64, inc: f64, lan: f64, ap: f64, m0: f64, cb: Rc<CentralBody>) -> KOE {
        let mut rot = Rot3::new_identity(3);
        if approx_eq(&inc, &0.0) {
            if !approx_eq(&e, &0.0) {
                rot = Rot3::new(cb.up * ap);
            }
        } else {
            let lan_axisangle = cb.up * lan;
            let inc_axisangle = cb.reference * inc;
            let ap_axisangle = cb.up * ap;
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
            cb: cb,
            to_csv_rot: rot,
        }
    }
    
    pub fn approx_eq(&self, other: &KOE) -> bool {
        let eps = 1.0e-5;
        approx_eq_eps(&self.a, &other.a, &eps) &&
        approx_eq_eps(&self.e, &other.e, &eps) &&
        approx_eq_eps(&self.inc, &other.inc, &eps) &&
        approx_eq_eps(&self.lan, &other.lan, &eps) &&
        approx_eq_eps(&self.ap, &other.ap, &eps) &&
        approx_eq_eps(&self.m0, &other.m0, &eps) &&
        self.cb.approx_eq(&other.cb)
    }

    pub fn approx_eq_eps(&self, other: &KOE, eps: &f64) -> bool {
        approx_eq_eps(&self.a, &other.a, &eps) &&
        approx_eq_eps(&self.e, &other.e, &eps) &&
        approx_eq_eps(&self.inc, &other.inc, &eps) &&
        approx_eq_eps(&self.lan, &other.lan, &eps) &&
        approx_eq_eps(&self.ap, &other.ap, &eps) &&
        approx_eq_eps(&self.m0, &other.m0, &eps) &&
        self.cb.approx_eq(&other.cb)
    }

    pub fn to_csv(&self) -> CSV {
        let m0 = self.m0;
        let iterations = 10;
        let ea = KOE::newton_raphson(&m0, &self.e, &iterations);
        let ta = 2.0*((1.0+self.e).sqrt()*(ea/2.0).sin())
            .atan2((1.0-self.e).sqrt()*(ea/2.0).cos());
        let dist = self.a*(1.0-self.e*ea.cos());
        let mut r = (self.cb.reference*ta.cos() + self.cb.right*ta.sin()) * dist;
        let mut v = (self.cb.reference*(-ea.sin()) +
                    self.cb.right*((1.0-self.e.powf(2.0)).sqrt()*ea.cos())) * ((self.cb.mu*self.a).sqrt()/dist);
        r = self.to_csv_rot.transform(&r);
        v = self.to_csv_rot.transform(&v);
        CSV::new(r, v, self.cb.clone())
    }

    fn newton_raphson(m0: &f64, e: &f64, iterations: &i32) -> f64 {
        let mut ea = m0.clone();
        for _ in 0..*iterations {
            ea  = ea - (ea - e*ea.sin() - m0)/(1.0 - e*ea.cos());
        }
        ea
    }
}
