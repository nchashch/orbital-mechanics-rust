pub trait Tick {
    fn tick(&self, dt: f64) -> Self;
}
