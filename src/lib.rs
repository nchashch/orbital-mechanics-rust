extern crate nalgebra;
extern crate rand;
use rand::*;
use nalgebra::*;

mod csv;
mod koe;
mod central_body;
use std::f64::consts::PI;
use central_body::*;
use csv::*;
use koe::*;

// TODO: This file is too long, find a better way of testing!
#[test]
fn koe_invariance() {
    let iterations = 10000;
    let eps = 1.0e-3;
    let mut rng = rand::thread_rng();
    let mut a_ok = 0;
    let mut e_ok = 0;
    let mut inc_ok = 0;
    let mut lan_ok = 0;
    let mut ap_ok = 0;
    let mut lp_ok = 0;
    let mut m0_ok = 0;
    for _ in 0..iterations {
        let up = normalize(
                &Vec3::new( // Up vector
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0)
                    )
            );
        let x_axis = Vec3::<f64>::new(1.0, 0.0, 0.0);
        let reference = normalize(&cross(&x_axis, &up));
        let cb = CentralBody::new(
            rng.gen_range(1.0e9, 1.0e25), // Standard Gravitational Parameter
            up,
            reference
        );
        let ecc = rng.gen_range(0.0, 1.0);
        let inc = rng.gen_range(0.0, PI);
        let mut lan = rng.gen_range(0.0, 2.0*PI);
        let mut ap = rng.gen_range(0.0, 2.0*PI);
        if approx_eq(&ecc, &0.0) {
            ap = 0.0;
        }
        let mut lp = 0.0;
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            ap = 0.0;
            if !approx_eq(&ecc, &0.0) {
                lp = rng.gen_range(0.0, 2.0*PI);
            }
        }
        let koe0 = KOE::new(
            rng.gen_range(1.0, 1.0e10), // Semi major axis
            ecc, // Eccentricity
            inc, // Inclination
            lan, // Longitude of Ascending Node
            ap,  // Argument of Periapsis
            lp,  // Longitude of Periapsis
            rng.gen_range(-PI, PI), // Mean anomaly
            &cb
            );
        let csv0 = koe0.to_csv();
        let koe1 = csv0.to_koe();
        if approx_eq_eps(&koe0.a, &koe1.a, &eps) {
            a_ok += 1;
        }
        if approx_eq_eps(&koe0.e, &koe1.e, &eps) {
            e_ok += 1;
        }
        if approx_eq_eps(&koe0.inc, &koe1.inc, &eps) {
            inc_ok += 1;
        }
        if approx_eq_eps(&koe0.lan, &koe1.lan, &eps) {
            lan_ok += 1;
        }
        if approx_eq_eps(&koe0.ap, &koe1.ap, &eps) {
            ap_ok += 1;
        }
        if approx_eq_eps(&koe0.lp, &koe1.lp, &eps) {
            lp_ok += 1;
        }
        if approx_eq_eps(&koe0.m0, &koe1.m0, &eps) {
            m0_ok += 1;
        }
    }
    let a_ok_f = a_ok as f32 / iterations as f32;
    let e_ok_f = e_ok as f32 / iterations as f32;
    let inc_ok_f = inc_ok as f32 / iterations as f32;
    let lan_ok_f = lan_ok as f32 / iterations as f32;
    let ap_ok_f = ap_ok as f32 / iterations as f32;
    let lp_ok_f = lp_ok as f32 / iterations as f32;
    let m0_ok_f = m0_ok as f32 / iterations as f32;
    let threshold = 0.99;
    let mut failed = false;

    println!("");
    println!("eps = {}", eps);
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
    if ap_ok_f < threshold {
        println!("ap_ok_f = {}", ap_ok_f);
        if !failed {
            failed = true;
        }
    }
    if lp_ok_f < threshold {
        println!("lp_ok_f = {}", lp_ok_f);
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
    let iterations = 1000;
    let eps = 1.0e-3;
    let mut rng = rand::thread_rng();
    let mut r_ok = 0;
    let mut v_ok = 0;
    for _ in 0..iterations {
        let up = normalize(
                &Vec3::new( // Up vector
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0)
                    )
            );
        let x_axis = Vec3::<f64>::new(1.0, 0.0, 0.0);
        let reference = normalize(&cross(&x_axis, &up));
        let cb = CentralBody::new(
            rng.gen_range(1.0e9, 1.0e25), // Standard Gravitational Parameter
            up,
            reference
        );
        let ecc = rng.gen_range(0.0, 1.0);
        let inc = rng.gen_range(0.0, PI);
        let mut lan = rng.gen_range(0.0, 2.0*PI);
        let mut ap = rng.gen_range(0.0, 2.0*PI);
        let mut lp = lan + ap;
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            ap = 0.0;
        };
        let koe0 = KOE::new(
            rng.gen_range(1.0, 1.0e10), // Semi major axis
            ecc, //rng.gen_range(0.0, 1.0),    // Eccentricity
            inc, // Inclination
            lan, // Longitude of Ascending Node
            ap, // Argument of Periapsis
            lp, // Longitude of Periapsis
            rng.gen_range(-PI, PI), // Mean anomaly
            &cb
            );
        let csv0 = koe0.to_csv();
        let koe1 = csv0.to_koe();
        let csv1 = koe1.to_csv();
        if csv0.r.approx_eq_eps(&csv1.r, &eps) {
            r_ok += 1;
        }
        if csv0.v.approx_eq_eps(&csv1.v, &eps) {
            v_ok += 1;
        }
    }
    let r_ok_f = r_ok as f32 / iterations as f32;
    let v_ok_f = v_ok as f32 / iterations as f32;
    let threshold = 0.99;
    let mut failed = false;

    println!("");
    println!("eps = {}", eps);
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

#[test]
fn koe_invariance_inc_eq_0() {
    let iterations = 10000;
    let eps = 1.0e-3;
    let mut rng = rand::thread_rng();
    let mut a_ok = 0;
    let mut e_ok = 0;
    let mut inc_ok = 0;
    let mut lan_ok = 0;
    let mut ap_ok = 0;
    let mut lp_ok = 0;
    let mut m0_ok = 0;
    for _ in 0..iterations {
        let up = normalize(
                &Vec3::new( // Up vector
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0)
                    )
            );
        let x_axis = Vec3::<f64>::new(1.0, 0.0, 0.0);
        let reference = normalize(&cross(&x_axis, &up));
        let cb = CentralBody::new(
            rng.gen_range(1.0e9, 1.0e25), // Standard Gravitational Parameter
            up,
            reference
        );
        let ecc = rng.gen_range(0.0, 1.0);
        let inc = 0.0;//rng.gen_range(0.0, PI);
        let mut lan = rng.gen_range(0.0, 2.0*PI);
        let mut ap = rng.gen_range(0.0, 2.0*PI);
        if approx_eq(&ecc, &0.0) {
            ap = 0.0;
        }
        let mut lp = 0.0;
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            ap = 0.0;
            if !approx_eq(&ecc, &0.0) {
                lp = rng.gen_range(0.0, 2.0*PI);
            }
        }
        let koe0 = KOE::new(
            rng.gen_range(1.0, 1.0e10), // Semi major axis
            ecc, // Eccentricity
            inc, // Inclination
            lan, // Longitude of Ascending Node
            ap,  // Argument of Periapsis
            lp,  // Longitude of Periapsis
            rng.gen_range(-PI, PI), // Mean anomaly
            &cb
            );
        let csv0 = koe0.to_csv();
        let koe1 = csv0.to_koe();
        if approx_eq_eps(&koe0.a, &koe1.a, &eps) {
            a_ok += 1;
        }
        if approx_eq_eps(&koe0.e, &koe1.e, &eps) {
            e_ok += 1;
        }
        if approx_eq_eps(&koe0.inc, &koe1.inc, &eps) {
            inc_ok += 1;
        } else {
            println!("");
            println!("koe0.inc = {}", koe0.inc);
            println!("koe1.inc = {}", koe1.inc);
        }
        if approx_eq_eps(&koe0.lan, &koe1.lan, &eps) {
            lan_ok += 1;
        }
        if approx_eq_eps(&koe0.ap, &koe1.ap, &eps) {
            ap_ok += 1;
        } else {
            println!("");
            println!("koe0.ap = {}", koe0.ap);
            println!("koe1.ap = {}", koe1.ap);
        }
        if approx_eq_eps(&koe0.lp, &koe1.lp, &eps) {
            lp_ok += 1;
        } else {
            println!("");
            println!("koe0.lp = {}", koe0.lp);
            println!("koe1.lp = {}", koe1.lp);
        }
        if approx_eq_eps(&koe0.m0, &koe1.m0, &eps) {
            m0_ok += 1;
        }
    }
    let a_ok_f = a_ok as f32 / iterations as f32;
    let e_ok_f = e_ok as f32 / iterations as f32;
    let inc_ok_f = inc_ok as f32 / iterations as f32;
    let lan_ok_f = lan_ok as f32 / iterations as f32;
    let ap_ok_f = ap_ok as f32 / iterations as f32;
    let lp_ok_f = lp_ok as f32 / iterations as f32;
    let m0_ok_f = m0_ok as f32 / iterations as f32;
    let threshold = 0.99;
    let mut failed = false;

    println!("");
    println!("eps = {}", eps);
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
    if ap_ok_f < threshold {
        println!("ap_ok_f = {}", ap_ok_f);
        if !failed {
            failed = true;
        }
    }
    if lp_ok_f < threshold {
        println!("lp_ok_f = {}", lp_ok_f);
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
fn csv_invariance_inc_eq_0() {
    let iterations = 1000;
    let eps = 1.0e-3;
    let mut rng = rand::thread_rng();
    let mut r_ok = 0;
    let mut v_ok = 0;
    for _ in 0..iterations {
        let up = normalize(
                &Vec3::new( // Up vector
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0),
                    rng.gen_range(0.1, 1.0)
                    )
            );
        let x_axis = Vec3::<f64>::new(1.0, 0.0, 0.0);
        let reference = normalize(&cross(&x_axis, &up));
        let cb = CentralBody::new(
            rng.gen_range(1.0e9, 1.0e25), // Standard Gravitational Parameter
            up,
            reference
        );
        let ecc = rng.gen_range(0.0, 1.0);
        let inc = 0.0;//rng.gen_range(0.0, PI);
        let mut lan = rng.gen_range(0.0, 2.0*PI);
        let mut ap = rng.gen_range(0.0, 2.0*PI);
        let mut lp = lan + ap;
        if approx_eq(&inc, &0.0) {
            lan = 0.0;
            ap = 0.0;
        };
        let koe0 = KOE::new(
            rng.gen_range(1.0, 1.0e10), // Semi major axis
            ecc, //rng.gen_range(0.0, 1.0),    // Eccentricity
            inc, // Inclination
            lan, // Longitude of Ascending Node
            ap, // Argument of Periapsis
            lp, // Longitude of Periapsis
            rng.gen_range(-PI, PI), // Mean anomaly
            &cb
            );
        let csv0 = koe0.to_csv();
        let koe1 = csv0.to_koe();
        let csv1 = koe1.to_csv();
        if csv0.r.approx_eq_eps(&csv1.r, &eps) {
            r_ok += 1;
        }
        if csv0.v.approx_eq_eps(&csv1.v, &eps) {
            v_ok += 1;
        }
    }
    let r_ok_f = r_ok as f32 / iterations as f32;
    let v_ok_f = v_ok as f32 / iterations as f32;
    let threshold = 0.99;
    let mut failed = false;

    println!("");
    println!("eps = {}", eps);
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
