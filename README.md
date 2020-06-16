# Plonky

Plonky is a prototype implementation of recursive arguments. It is loosely based on [PLONK](https://eprint.iacr.org/2019/953), with various modifications:

* While PLONK uses [KZG](https://www.iacr.org/cryptodb/data/paper.php?pubkey=23846)'s pairing-based polynomial commitment scheme, we use the [Halo](https://eprint.iacr.org/2019/1021) technique to recursively verify discrete log based polynomial commitments.
* The standard PLONK model was designed for arithmetic circuits; it uses a single constraint to verify additive and multiplicative relationships. We use a variety of custom gates, such as a gate which performs a full round of a [Rescue](https://eprint.iacr.org/2019/426) permutation. The maximum degree of our constraints is 7.
* In the standard version of PLONK, each gate interacts with three wires, which are typically thought of as two input wires and one output wire. We use a much higher arity -- 11 wires per gate -- although only 6 of them are involved in the permutation argument. The other 5 can be thought of as "advice" wires.
* In PLONK, the verifier generates a challenge point `x` and polynomials are opened at both `x` and `g x` (where `g` is some multiplicative group generator). This makes it possible for a gate to access the wire values of the "following" gate at no additional cost, a la TurboPLONK, which can be helpful in eliminating redundant state. We add a third opening at `g^65 x` which, if we imagine gates arranged on a grid with a width of 65, allows each gate to access neighboring gates along both dimensions.

## Disclaimer

This code has not been thoroughly reviewed or tested, and should not be used in any production systems.
