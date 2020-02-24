use chashmap::CHashMap;

use lazy_static::lazy_static;

use crate::Bls12Scalar;

lazy_static! {
    static ref SUBGROUPS_BY_ORDER_POWER: CHashMap<usize, Vec<Bls12Scalar>> = CHashMap::new();
}

/// Permutes `arr` such that each index is mapped to its reverse in binary.
fn reverse_index_bits<T: Copy>(arr: Vec<T>) -> Vec<T> {
    let n = arr.len();
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        // TODO: We can reverse either index; is one faster?
        result[i.reverse_bits()] = arr[i];
    }
    result
}

fn get_subgroup(order_power: usize) -> Vec<Bls12Scalar> {
    match SUBGROUPS_BY_ORDER_POWER.get(&order_power) {
        Some(subgroup) => subgroup.clone(),
        None => {
            let subgroup = Bls12Scalar::cyclic_subgroup_unknown_order(
                Bls12Scalar::primitive_root_of_unity(order_power));
            SUBGROUPS_BY_ORDER_POWER.insert(order_power, subgroup.clone());
            subgroup
        }
    }
}

// TODO: Make coefficients mutable and do it in-place?
fn fft(subgroup: Vec<Bls12Scalar>, coefficients: Vec<Bls12Scalar>) -> Vec<Bls12Scalar> {
    let mut degree_pow = 0;
    while 1 << degree_pow < coefficients.len() {
        degree_pow += 1;
    }
    let degree = 1 << degree_pow;
    assert_eq!(coefficients.len(), degree,
               "Coefficients must be padded to a power of 2 size");
    let half_degree = degree >> 1;

//    let coefficients = reverse_index_bits(coefficients);
    let mut buffer = Vec::with_capacity(degree);

    // In the first phase, we have a size of 2, so we're evaluating degree 1 polynomials.
    let base_x_0 = subgroup[0];
    let base_x_1 = subgroup[degree >> 1];
    for i_0 in 0..half_degree {
        let i_1 = i_0 + half_degree;
        buffer[i_0] = coefficients[i_0] * base_x_0 + coefficients[i_1];
        buffer[i_0 + 1] = coefficients[i_0] * base_x_1 + coefficients[i_1];
    }

    let size = 2;
    loop {}

    buffer
}
