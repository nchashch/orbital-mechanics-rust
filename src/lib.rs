extern crate nalgebra;
extern crate rand;

pub mod om;
mod tick;
mod push;
pub use tick::*;
pub use push::*;

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;
    use std::rc::*;
    use om::*;
    use rand::*;
    use nalgebra::*;
    
    const ITERATIONS: u32 = 10000;
    const EPS: f64 = 1.0e-3;

    #[test]
    fn koe_invariance() {
        let mut a_ok = 0;
        let mut e_ok = 0;
        let mut inc_ok = 0;
        let mut lan_ok = 0;
        let mut ap_ok = 0;
        let mut m0_ok = 0;
        for _ in 0..ITERATIONS {
            let cb = get_random_cb();
            let koe0 = get_random_koe(Rc::<CB>::new(cb));
            let csv0 = CSV::from_koe(koe0.clone());
            let koe1 = KOE::from_csv(csv0);
            if approx_eq_eps(&koe0.a, &koe1.a, &EPS) {
                a_ok += 1;
            }
            if approx_eq_eps(&koe0.e, &koe1.e, &EPS) {
                e_ok += 1;
            }
            if approx_eq_eps(&koe0.inc, &koe1.inc, &EPS) {
                inc_ok += 1;
            }
            if approx_eq_eps(&koe0.lan, &koe1.lan, &EPS) {
                lan_ok += 1;
            }
            if approx_eq_eps(&koe0.ap, &koe1.ap, &EPS) {
                ap_ok += 1;
            }
            if approx_eq_eps(&koe0.m0, &koe1.m0, &EPS) {
                m0_ok += 1;
            }
        }
        let a_ok_f = a_ok as f32 / ITERATIONS as f32;
        let e_ok_f = e_ok as f32 / ITERATIONS as f32;
        let inc_ok_f = inc_ok as f32 / ITERATIONS as f32;
        let lan_ok_f = lan_ok as f32 / ITERATIONS as f32;
        let ap_ok_f = ap_ok as f32 / ITERATIONS as f32;
        let m0_ok_f = m0_ok as f32 / ITERATIONS as f32;
        let threshold = 0.99;
        let mut failed = false;

        println!("");
        println!("EPS = {}", EPS);
        println!("threshold = {}", threshold);
        if a_ok_f < threshold {
            println!("a_ok_f = {}", a_ok_f);
            if !failed {
                failed = true;
            }
        }
        if e_ok_f < threshold {
            println!("e_ok_f = {}", e_ok_f);
            if !failed {
                failed = true;
            }
        }
        if inc_ok_f < threshold {
            println!("inc_ok_f = {}", inc_ok_f);
            if !failed {
                failed = true;
            }
        }
        if lan_ok_f < threshold {
            println!("lan_ok_f = {}", lan_ok_f);
            if !failed {
                failed = true;
            }
        }
        if ap_ok_f < threshold {
            println!("ap_ok_f = {}", ap_ok_f);
            if !failed {
                failed = true;
            }
        }
        if m0_ok_f < threshold {
            println!("m0_ok_f = {}", m0_ok_f);
            if !failed {
                failed = true;
            }
        }
        if failed {
            panic!("success rate for one of the parameters is below threshold");
        }
    }

    #[test]
    fn csv_invariance() {
        let mut r_ok = 0;
        let mut v_ok = 0;
        for _ in 0..ITERATIONS {
            let cb = get_random_cb();
            let koe0 = get_random_koe(Rc::<CB>::new(cb));
            let csv0 = CSV::from_koe(koe0);
            let koe1 = KOE::from_csv(csv0.clone());
            let csv1 = CSV::from_koe(koe1);
            if csv0.r.approx_eq_eps(&csv1.r, &EPS) {
                r_ok += 1;
            }
            if csv0.v.approx_eq_eps(&csv1.v, &EPS) {
                v_ok += 1;
            }
        }
        let r_ok_f = r_ok as f32 / ITERATIONS as f32;
        let v_ok_f = v_ok as f32 / ITERATIONS as f32;
        let threshold = 0.99;
        let mut failed = false;

        println!("");
        println!("EPS = {}", EPS);
        println!("threshold = {}", threshold);
        if r_ok_f < threshold {
            println!("r_ok_f = {}", r_ok_f);
            if !failed {
                failed = true;
            }
        }
        if v_ok_f < threshold {
            println!("v_ok_f = {}", v_ok_f);
            if !failed {
                failed = true;
            }
        }
        if failed {
            panic!("success rate for some of the parameters is below the threshold");
        }
    }

    fn get_random_cb() -> CB {
        let mut rng = thread_rng();
        let rand_rot_vec =
            normalize(
                &Vec3::<f64>::new( // Up vector
                    rng.gen_range(-1.0, 1.0),
                    rng.gen_range(-1.0, 1.0),
                    rng.gen_range(-1.0, 1.0)
                    )
            )*rng.gen_range(0.0, 2.0*PI);

        let rand_rot = Rot3::<f64>::new(rand_rot_vec);

        let i = Vec3::<f64>::new(1.0, 0.0, 0.0);
        let j = Vec3::<f64>::new(0.0, 1.0, 0.0);
        let k = Vec3::<f64>::new(0.0, 0.0, 1.0);

        let i = rand_rot.transform(&i);
        let j = rand_rot.transform(&j);
        let k = rand_rot.transform(&k);
        
        CB::new(
            rng.gen_range(1.0e9, 1.0e25), // Standard Gravitational Parameter
            i,
            j,
            k
        )
    }

    fn get_random_koe(cb: Rc<CB>) -> KOE {
        let mut rng = thread_rng();
        let ecc = rng.gen_range(0.0, 1.0);
        let inc = rng.gen_range(0.0, PI);
        let mut lan = rng.gen_range(0.0, 2.0*PI);
        let mut ap = rng.gen_range(0.0, 2.0*PI);
        if approx_eq(&ecc, &0.0) {
            ap = 0.0;
        }
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            ap = 0.0;
            if !approx_eq(&ecc, &0.0) {
                ap = rng.gen_range(0.0, 2.0*PI);
            }
        }
        KOE::new(
            rng.gen_range(1.0, 1.0e10), // Semi major axis
            ecc, // Eccentricity
            inc, // Inclination
            lan, // Longitude of Ascending Node
            ap,  // Argument of Periapsis
            rng.gen_range(-PI, PI), // Mean anomaly
            cb
        )
    }
}
