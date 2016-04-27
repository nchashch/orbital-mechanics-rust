extern crate nalgebra;
use nalgebra::*;
use om::cb::*;
use om::csv::*;
use std::rc::*;
use tick::*;
use std::f64::consts::*;

/// #Keplerian Orbital Elements
/// This structure represents an orbit using
/// six keplerian element. It also holds mean motion
/// to avoid recomputing it for every tick() call, and
/// rot matrix to avoid recomputing it for every CSV::from_koe() call.
/// Like CSV it holds a reference to the central body.
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
    /// Mean motion.
    pub n: f64,
    /// A matrix that transforms a vector lying in i, j plane into a corresponding vector lying in the orbital plane
    /// defined by inc and lan, it is stored in KOE to avoid recomputing it for every CSV::from_koe() call.
    pub rot: Rot3<f64>,
    /// Reference to the central body.
    pub cb: Rc<CB>,
}

impl Tick for KOE {
    fn tick(&self, dt: f64) -> Self {
        KOE {
            m0: self.m0 + self.n * dt,
            cb: self.cb.clone(),
            ..*self
        }
    }
}

impl KOE {
    /// Construct KOE from orbital elements and CB reference.
    pub fn new(a: f64, e: f64, inc: f64, lan: f64, ap: f64, m0: f64, cb: Rc<CB>) -> KOE {
        // Vector cb.i points towards intersection of 0th meridian and equator
        // Vector cb.k points towards north pole
        // cb.j == cross(&cb.k, &cb.i)

        // rot transformation matrix is not expected to change very often
        // and it is stored in KOE to avoid recomputing it
        // for every CSV::from_koe() call
        let mut rot = Rot3::new_identity(3);
        // Do lan and inc rotations only if
        // the orbit is not equatorial
        if !approx_eq(&inc, &0.0) {
            let lan_axisangle = cb.k * lan;
            let inc_axisangle = cb.i * inc;
            rot = rot.prepend_rotation(&lan_axisangle);
            rot = rot.prepend_rotation(&inc_axisangle);
        }
        // Do ap rotation only if
        // the orbit is not circular
        if !approx_eq(&e, &0.0) {
            let ap_axisangle = cb.k * ap;
            rot = rot.prepend_rotation(&ap_axisangle);
        }
        // Mean motion
        let n = (cb.mu/a.powf(3.0)).sqrt();
        KOE {
            a: a,
            e: e,
            inc: inc,
            lan: lan,
            ap: ap,
            m0: m0,
            n: n,
            rot: rot,
            cb: cb,
        }
    }

    /// Construct KOE from CSV.
    pub fn from_csv(csv: CSV) -> KOE {
        // Vector csv.cb.i points towards intersection of 0th meridian and equator
        // Vector csv.cb.k points towards north pole
        // csv.cb.j == cross(&csv.cb.k, &csv.cb.i)

        // Radius vector
        let r = csv.r;
        // Velocity
        let v = csv.v;

        // Standard gravitational parameter
        let mu = &csv.cb.mu;

        // Specific angular momentum
        let h = cross(&r, &v);

        // Eccentricity vector
        // (It is pointed towards periapsis and it's length
        // is equal to eccentricity of the orbit)
        let e = cross(&v, &h)/ *mu - normalize(&r);

        // Node vector
        // (vector pointing towards the ascending node)
        // Ascending node is a point where satelite
        // is above equator and goes north
        let n = cross(&csv.cb.k, &h);

        let cos_inc = dot(&h, &csv.cb.k)/(norm(&h));
        // cos_inc is sometimes greater than 1.0
        // and without this fix cos_inc.acos() is NaN
        // for cos_inc > 1.0 cases
        let inc = if cos_inc > 1.0 {
            0.0
        } else {
            cos_inc.acos()
        };

        // Eccentricity
        let es = norm(&e);

        // Longitude of Ascending Node
        // (angle between vector csv.cb.i and ascending node)
        let mut lan = if dot(&n, &csv.cb.j) >= 0.0 {
            (dot(&n, &csv.cb.i)/norm(&n)).acos()
        } else {
            2.0*PI - (dot(&n, &csv.cb.i)/norm(&n)).acos()
        };

        let right = cross(&h, &n);
        // Argument of periapsis
        // (angle between ascending node and periapsis)
        let mut ap = if dot(&e, &right) >= 0.0 {
            (dot(&n, &e)/(norm(&n)*norm(&e))).acos()
        } else {
            2.0*PI - (dot(&n, &e)/(norm(&n)*norm(&e))).acos()
        };

        // If the orbit is circular ap is 0.0
        // (ap doesn't make sense for circular orbits)
        if approx_eq(&es, &0.0) {
            ap = 0.0;
        }
        // If the orbit is equatorial lan is 0.0
        // (lan doesn't make sense for equatorial orbits)
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            // If it is equatorial, non circular orbit ap is Longitude of Periapsis
            // (angle between vector csv.cb.i and periapsis)
            if !approx_eq(&es, &0.0) {
                ap = if dot(&e, &csv.cb.j) >= 0.0 {
                    (dot(&csv.cb.i, &e)/norm(&e)).acos()
                } else {
                    2.0*PI - (dot(&csv.cb.i, &e)/norm(&e)).acos()
                };
            }
        }

        // True anomaly
        // (angle between periapsis and radius vector)
        let mut ta = if dot(&r, &v) >= 0.0 {
            (dot(&e, &r)/(norm(&e)*norm(&r))).acos()
        } else {
            2.0*PI - (dot(&e, &r)/(norm(&e)*norm(&r))).acos()
        };

        if approx_eq(&es, &0.0) {
            // For circular equatorial orbit use longitude
            // (angle between vector csv.cb.i and radius vector)
            if approx_eq(&inc, &0.0) {
                ta = if dot(&csv.cb.i, &v) <= 0.0 {
                    (dot(&csv.cb.i, &r)/(norm(&csv.cb.i)*norm(&r))).acos()
                } else {
                    2.0*PI - (dot(&csv.cb.i, &r)/(norm(&csv.cb.i)*norm(&r))).acos()
                }
            // For circular non equatorial orbit use argument of latitude
            // (angle between ascending node and radius vector)
            } else {
                ta = if dot(&n, &v) <= 0.0 {
                    (dot(&n, &r)/(norm(&n)*norm(&r))).acos()
                } else {
                    2.0*PI - (dot(&n, &r)/(norm(&n)*norm(&r))).acos()
                }
            }
        }

        // Eccentric anomaly (intermidiate step to compute mean anomaly)
        let ea = 2.0*((ta/2.0).tan()/((1.0+es)/(1.0-es)).sqrt()).atan();
        // Mean anomaly (it is used because it changes linearly with time,
        // and for that reason is cheap to update)
        let m0 = ea - es * ea.sin();
        // Semi Major Axis
        let a = 1.0/(2.0/norm(&r) - sqnorm(&v)/mu);
        KOE::new(a, es, inc, lan, ap, m0, csv.cb.clone())
    }
}
