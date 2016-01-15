use nalgebra::*;

/* Central body */
pub struct CB {
    pub mu: f64,
    pub i: Vec3<f64>,
    pub j: Vec3<f64>,
    pub k: Vec3<f64>,
}

impl CB {
    pub fn new(mu: f64, i: Vec3<f64>, j: Vec3<f64>, k: Vec3<f64>) -> CB {
        assert!(approx_eq(&dot(&i, &j), &0.0));
        assert!(approx_eq(&dot(&i, &k), &0.0));
        assert!(approx_eq(&dot(&j, &k), &0.0));
        assert!(approx_eq(&norm(&i), &1.0));
        assert!(approx_eq(&norm(&j), &1.0));
        assert!(approx_eq(&norm(&k), &1.0));
        assert!(approx_eq(&cross(&i, &j), &k));
        CB {mu: mu, i: i, j: j, k: k }
    }
    pub fn change_cb(&self, vec: Vec3<f64>, other: CB) -> Vec3<f64> {
        other.i * dot(&vec, &self.i) +
        other.j * dot(&vec, &self.j) +
        other.k * dot(&vec, &self.k)
    }
}
