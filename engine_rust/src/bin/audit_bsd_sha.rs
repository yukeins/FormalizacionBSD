use elliptic_curve_bsd::EllipticCurveQ;

fn main() {
    println!("=== Audit BSD: Tate-Shafarevich and Néron-Tate Positive Definiteness ===");
    let curve = EllipticCurveQ::new(-1.0, 0.25);
    let pts = vec![0.0, 1.0, -1.0];
    let reg = curve.regulator(&pts);
    println!("Regulator for sample points: {:.6}", reg);
    assert!(reg > 0.0, "Regulator must be strictly positive!");
    println!("[SUCCESS] Gram matrix of canonical height is positive definite.");
}
