# Plonky

Plonky is a prototype implementation of recursive arguments. It is loosely based on [PLONK](https://eprint.iacr.org/2019/953), with a few customizations:

* While PLONK uses [KZG](https://www.iacr.org/cryptodb/data/paper.php?pubkey=23846)'s pairing-based polynomial commitment scheme, we use a batched variant of the [Halo](https://eprint.iacr.org/2019/1021) technique to recursively verify discrete log based polynomial commitments.
* The standard PLONK model was designed for arithmetic circuits; it uses a single constraint to verify additive and multiplicative relationships. We use a variety of custom gates, such as a gate which performs a step of a [Rescue](https://eprint.iacr.org/2019/426) permutation. The maximum degree of our constraints is 8, compared to 3 in standard PLONK.
* In standard PLONK, each gate interacts with three wires, which are typically thought of as two input wires and one output wire. We use a much higher arity -- 9 wires per gate -- although only 6 of them are involved in the permutation argument. The other 3 can be thought of as "advice" wires.
* The zero-knowledge technique in the Plonk paper would cause each witness polynomial's degree to slightly exceed a power of two, which doesn't work well with Halo. We use a [different blinding method](https://mirprotocol.org/blog/Adding-zero-knowledge-to-Plonk-Halo).

For more details, see [Fast recursive arguments based on Plonk and Halo](https://mirprotocol.org/blog/Fast-recursive-arguments-based-on-Plonk-and-Halo).


## Disclaimer

This code has not been thoroughly reviewed or tested, and should not be used in any production systems.
