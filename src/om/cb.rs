use nalgebra::*;

/// Central body.
pub struct CB {
    /// Standard gravitational parameter.
    pub mu: f64,
    /// Basis vector.
    pub i: Vec3<f64>,
    /// Basis vector.
    pub j: Vec3<f64>,
    /// Basis vector.
    pub k: Vec3<f64>,
}

impl CB {
    /// Create a new central body asserting that i, j, k compose an orthonormal right-handed coordinate system basis.
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

    /// Find what value does a given vector has in the coordinate system associated with a different central body.
    pub fn change_cb(&self, vec: Vec3<f64>, other: CB) -> Vec3<f64> {
        other.i * dot(&vec, &self.i) +
        other.j * dot(&vec, &self.j) +
        other.k * dot(&vec, &self.k)
    }
}
