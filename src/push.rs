extern crate nalgebra;
use nalgebra::*;

pub trait Push {
    fn push(&self, dv: Vec3<f64>) -> Self;
}
