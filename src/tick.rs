/// Objects that implement this trait can tell what their
/// future state will be.
pub trait Tick {
    /// Return state of self in dt seconds.
    fn tick(&self, dt: f64) -> Self;
}
