pub trait Tickable {
    fn tick(&self, dt: &f64) -> Self;
}
