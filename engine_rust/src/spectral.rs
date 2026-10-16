//! Jacobi Spectral Diagonalization for Symmetric Real Gram Matrices.
//! Transferred from the Riemann-Weil Positivity Engine.

/// Computes all eigenvalues of a real symmetric matrix using the classical Jacobi rotation method.
/// Guaranteed to converge for all positive-definite or symmetric matrices.
#[allow(clippy::needless_range_loop)]
pub fn jacobi_eigenvalues(matrix: &[Vec<f64>]) -> Vec<f64> {
    let n = matrix.len();
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![matrix[0][0]];
    }

    let mut a = matrix.to_vec();
    let max_iter = 1000;
    let tol = 1e-12;

    for _ in 0..max_iter {
        let mut max_val = 0.0;
        let mut p = 0;
        let mut q = 0;

        for i in 0..n {
            for j in i + 1..n {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < tol {
            break;
        }

        let diff = a[q][q] - a[p][p];
        let t = if diff.abs() < 1e-15 {
            if a[p][q] > 0.0 { 1.0 } else { -1.0 }
        } else {
            let theta = diff / (2.0 * a[p][q]);
            let t_candidate = 1.0 / (theta.abs() + (theta * theta + 1.0).sqrt());
            if theta < 0.0 { -t_candidate } else { t_candidate }
        };

        let c = 1.0 / (t * t + 1.0).sqrt();
        let s = t * c;
        let tau = s / (1.0 + c);

        let temp_p = a[p][p] - t * a[p][q];
        let temp_q = a[q][q] + t * a[p][q];
        a[p][p] = temp_p;
        a[q][q] = temp_q;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for j in 0..n {
            if j != p && j != q {
                let apj = a[p][j];
                let aqj = a[q][j];
                a[p][j] = apj - s * (aqj + tau * apj);
                a[q][j] = aqj + s * (apj - tau * aqj);
                a[j][p] = a[p][j];
                a[j][q] = a[q][j];
            }
        }
    }

    let mut eigs = Vec::with_capacity(n);
    for i in 0..n {
        eigs.push(a[i][i]);
    }
    eigs.sort_by(|x, y| x.partial_cmp(y).unwrap());
    eigs
}

/// Verifies strict coercivity (Lang-Silverman condition) on the Gram matrix:
/// lambda_min(G) >= c_E > 0.
pub fn verify_rayleigh_coercivity(matrix: &[Vec<f64>], c_e_floor: f64) -> (bool, f64) {
    let eigs = jacobi_eigenvalues(matrix);
    if eigs.is_empty() {
        return (true, 1.0);
    }
    let lambda_min = eigs[0];
    (lambda_min >= c_e_floor, lambda_min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jacobi_identity() {
        let m = vec![
            vec![2.0, 0.5],
            vec![0.5, 3.0],
        ];
        let eigs = jacobi_eigenvalues(&m);
        assert!(eigs[0] > 0.0);
        assert!(eigs[1] > eigs[0]);
        // Trace conservation: eigs[0] + eigs[1] == 5.0
        assert!((eigs[0] + eigs[1] - 5.0).abs() < 1e-10);
    }
}
