use elliptic_curve_bsd::EllipticCurveQ;

fn main() {
    println!("=== Audit F_p: Ab Initio Hecke Fourier Coefficients a_p via Point Counting ===");
    // Canonical curve 37a1: in short Weierstrass y^2 = x^3 - x + 1 (model over F_p)
    let curve = EllipticCurveQ::new(-1, 1);
    let primes = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];
    
    println!("{:<6} {:<10} {:<10} {:<15} {:<10}", "Prime p", "#E(F_p)", "a_p", "Hasse Bound 2√p", "Status");
    println!("---------------------------------------------------------------");
    for &p in &primes {
        let a_p = curve.compute_a_p(p);
        let n_points = p + 1 - a_p;
        let bound = 2.0 * (p as f64).sqrt();
        let pass = (a_p.abs() as f64) <= bound;
        println!("{:<6} {:<10} {:<10} {:<15.2} {:<10}", p, n_points, a_p, bound, if pass { "PASS" } else { "FAIL" });
        assert!(pass, "Hasse bound violated!");
    }
    println!("---------------------------------------------------------------");
    println!("[SUCCESS] All Hecke eigenvalues a_p generated ab initio satisfy Hasse bound.");
}
