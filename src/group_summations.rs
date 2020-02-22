use std::iter::Sum;

use crate::{Bls12Base, G1AffinePoint, G1ProjectivePoint};

impl Sum<G1AffinePoint> for G1ProjectivePoint {
    fn sum<I: Iterator<Item=G1AffinePoint>>(mut iter: I) -> Self {
        let points: Vec<G1AffinePoint> = iter.collect();
        best_affine_summation(points)
    }
}

impl Sum for G1ProjectivePoint {
    fn sum<I: Iterator<Item=G1ProjectivePoint>>(mut iter: I) -> Self {
        iter.fold(G1ProjectivePoint::ZERO, |acc, x| acc + x)
    }
}

fn best_affine_summation(points: Vec<G1AffinePoint>) -> G1ProjectivePoint {
    // This threshold is chosen based on data from the summation benchmarks.
    if points.len() < 23 {
        pairwise_affine_summation(points)
    } else {
        pairwise_affine_summation_batch_inversion(points)
    }
}

/// Adds each pair of points using an affine + affine = projective formula, then adds up the
/// intermediate sums using a projective formula.
pub fn pairwise_affine_summation(points: Vec<G1AffinePoint>) -> G1ProjectivePoint {
    let mut reduced_points: Vec<G1ProjectivePoint> = Vec::new();
    for chunk in points.chunks(2) {
        match chunk.len() {
            1 => reduced_points.push(chunk[0].to_projective()),
            2 => reduced_points.push(chunk[0] + chunk[1]),
            _ => panic!(),
        }
    }
    // TODO: Avoid copying (deref)
    reduced_points.iter().fold(G1ProjectivePoint::ZERO, |sum, x| sum + *x)
}

pub fn pairwise_affine_summation_batch_inversion(points: Vec<G1AffinePoint>) -> G1ProjectivePoint {
    // Short Weierstrass group laws are interesting in that they require an inversion but only three
    // multiplications. Normally the inversion would be prohibitively expensive, but in the context
    // of batch addition, we can apply Montgomery's trick to amortize the cost. Montgomery's trick
    // adds three multiplications per inversion, bringing the cost to six multiplications per group
    // operation, which is still quite efficient compared to a mixed addition (projective + affine)
    // sum. Note that this approach only allows us to add pairs of points simultaneously, so we
    // still end up with batch_size / 2 points which need to be summed up separately.

    let n = points.len();

    let mut elements_to_invert = Vec::new();
    for i in (0..n - 1).step_by(2) {
        let x1 = points[i].x;
        let x2 = points[i + 1].x;
        let y1 = points[i].y;
        elements_to_invert.push(if x1 == x2 { y1.double() } else { x1 - x2 })
    }

    let inverses: Vec<Option<Bls12Base>> = Bls12Base::batch_multiplicative_inverse_opt(&elements_to_invert);

    let mut reduced_points = Vec::new();
    for i in (0..n - 1).step_by(2) {
        let p1 = points[i];
        let p2 = points[i + 1];
        let G1AffinePoint { x: x1, y: y1 } = p1;
        let G1AffinePoint { x: x2, y: y2 } = p2;

        let inverse = inverses[i / 2];
        let sum = if p1.is_zero() {
            p2
        } else if p2.is_zero() {
            p1
        } else if p1 == -p2 {
            G1AffinePoint::ZERO
        } else if p1 == p2 {
            // This is the doubling case.
            let quotient = x1.square().triple() * inverse.unwrap();
            let x3 = quotient.square() - x1.double();
            let y3 = quotient * (x1 - x3) - y1;
            G1AffinePoint { x: x3, y: y3 }
        } else {
            // This is the general case. We use the incomplete addition formulas 4.3 and 4.4.
            let quotient = (y1 - y2) * inverse.unwrap();
            let x3 = quotient.square() - x1 - x2;
            let y3 = quotient * (x1 - x3) - y1;
            G1AffinePoint { x: x3, y: y3 }
        };
        reduced_points.push(sum);
    }

    // If n is odd, the last point was not part of a pair.
    if n % 2 == 1 {
        reduced_points.push(points[n - 1]);
    }

    // Recurse with our smaller set of points.
    best_affine_summation(reduced_points)
}

#[cfg(test)]
mod tests {
    use crate::group_summations::{pairwise_affine_summation, pairwise_affine_summation_batch_inversion};
    use crate::{G1AffinePoint, G1_GENERATOR_AFFINE};

    #[test]
    fn test_pairwise_affine_summation() {
        let g_affine = G1_GENERATOR_AFFINE;
        let g2_affine = (g_affine + g_affine).to_affine();
        let g3_affine = (g_affine + g_affine + g_affine).to_affine();
        let g_proj = g_affine.to_projective();
        let g2_proj = g2_affine.to_projective();
        let g3_proj = g3_affine.to_projective();
        assert_eq!(pairwise_affine_summation(vec![g_affine, g_affine]), g2_proj);
        assert_eq!(pairwise_affine_summation(vec![g_affine, g2_affine]), g3_proj);
        assert_eq!(pairwise_affine_summation(vec![g_affine, g_affine, g_affine]), g3_proj);
    }

    #[test]
    fn test_pairwise_affine_summation_batch_inversion() {
        let g = G1_GENERATOR_AFFINE;
        let g_proj = g.to_projective();
        assert_eq!(pairwise_affine_summation_batch_inversion(vec![g, g]), g_proj + g_proj);
        assert_eq!(pairwise_affine_summation_batch_inversion(vec![g, g, g]), g_proj + g_proj + g_proj);
    }
}
