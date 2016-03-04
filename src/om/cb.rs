use nalgebra::*;

/// #Central body
/// It holds standard gravitational parameter (product of body's mass and the gravitational constant)
/// and a local coordinate system associated with this body.
pub struct CB {
    /// Standard gravitational parameter.
    pub mu: f64,
    /// A unit vector pointing towards
    /// the intersection of 0th meridian and equator.
    pub i: Vec3<f64>,
    /// j = cross(&k, &i)
    pub j: Vec3<f64>,
    /// A unit vector pointing towards the north pole.
    pub k: Vec3<f64>,
}

impl CB {
    /// Create a new central body asserting that i, j, k compose an orthonormal right-handed vector basis.
    pub fn new(mu: f64, i: Vec3<f64>, j: Vec3<f64>, k: Vec3<f64>) -> CB {
        assert!(approx_eq(&dot(&i, &j), &0.0));
        assert!(approx_eq(&dot(&i, &k), &0.0));
        assert!(approx_eq(&dot(&j, &k), &0.0));
        assert!(approx_eq(&norm(&i), &1.0));
        assert!(approx_eq(&norm(&j), &1.0));
        assert!(approx_eq(&norm(&k), &1.0));
        assert!(approx_eq(&cross(&i, &j), &k));
        CB {
            mu: mu,
            i: i,
            j: j,
            k: k,
        }
    }
}
