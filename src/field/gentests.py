#
# Generate input/output pairs for testing modular arithmetic
#
# Based on: https://github.com/unzvfu/cuda-fixnum/blob/master/tests/gentests.py
#

from itertools import chain, product
from collections import deque
from timeit import default_timer as timer
from gmpy2 import gcdext

def write_int(dest, sz, n):
    dest.write(n.to_bytes(sz, byteorder = 'little'))

def write_vector(dest, elt_sz, v):
    for n in v:
        write_int(dest, elt_sz, n)

def mktests(op, xs, nargs, modulus):
    # TODO: Refactor this.
    if nargs == 1:
        yield zip(*[op(x, modulus) for x in xs])
    elif nargs == 2:
        ys = deque(xs)
        for i in range(len(xs)):
            yield zip(*[op(x, y, modulus) for x, y in zip(xs, ys)])
            ys.rotate(1)
    else:
        raise NotImplementedError()

def write_tests(fname, arg):
    op, xs, nargs, nres, modulus = arg
    vec_len = len(xs)
    ntests = vec_len**nargs
    t = timer()
    print('Writing {} tests into "{}"... '.format(ntests, fname), end='', flush=True)
    with open(fname, 'wb') as f:
        fixnum_bytes = (modulus.bit_length() + 7) // 8
        write_int(f, 4, fixnum_bytes)
        write_int(f, 4, vec_len)
        write_int(f, 4, nres)
        write_vector(f, fixnum_bytes, xs)
        for v in mktests(op, xs, nargs, modulus):
            v = list(v)
            assert len(v) == nres, 'bad result length; expected {}, got {}'.format(nres, len(v))
            for res in v:
                write_vector(f, fixnum_bytes, res)
    t = timer() - t
    print('done ({:.2f}s).'.format(t))
    return fname

# TODO: Not used yet
def reduce_mod(x, modulus):
    return [x % modulus]

def add_mod(x, y, modulus):
    return [(x + y) % modulus]

def sub_mod(x, y, modulus):
    return [(x - y) % modulus]

def neg_mod(x, modulus):
    return [modulus - x]

def mul_mod(x, y, modulus):
    return [(x * y) % modulus]

def sqr_mod(x, modulus):
    return [(x * x) % modulus]

def div_mod(x, y, modulus):
    if y == 0:
        # Minor hack: Indicate that we were given a zero divisor by
        # returning modulus. (Don't return 0, as that's a valid answer
        # for 0/y.)
        return [modulus]
    g, yinv, _ = gcdext(y, modulus)
    assert g == 1, 'bad modulus'
    yinv = int(yinv) # convert mpz -> int
    return [(x * yinv) % modulus]

def test_inputs(modulus):
    wordbits = 64
    modbits = modulus.bit_length()
    modwords = (modbits + wordbits - 1) // wordbits
    totalbits = modwords * wordbits

    print(f'   modulus length: {modbits} <= {wordbits}*{modwords} = {totalbits}')

    res = []

    # Numbers close to word boundaries
    nums = [1, 2, 3, 4, 5];
    nums.extend([2**wordbits - n for n in nums])
    for i in range(modwords):
        res.extend(n << wordbits*i for n in nums)
    res += [modulus - r for r in res] + [0]

    ## There are too many inputs here; find a way to reduce. Maybe
    ## just use mpz_rrandomb() to generate a fixed number of inputs
    # for i in range(2, totalbits - 1):
    #     # b = 0xF, 0xFF, 0xFFFF, 0xFFFFFFFF, ...
    #     e = 1 << i
    #     b = (1 << e) - 1
    #     c = sum(b << 2*e*j for j in range(totalbits // (2*e)))
    #     res.extend([c, (1 << totalbits) - c - 1])

    # Return all elements in the standard range.
    return [r % modulus for r in res]

def generate_tests(modulusname, modulus, tests):
    print(f'Generating input arguments for {modulusname}...')

    t = timer()
    xs = test_inputs(modulus)
    t = timer() - t
    print('done ({:.2f}s). Created {} arguments.'.format(t, len(xs)))

    # ops maps a test name to a tuple of (function name, generated
    # inputs, number of inputs, number of outputs, modulus)
    ops = {
        'add': (add_mod, xs, 2, 1, modulus),
        'sub': (sub_mod, xs, 2, 1, modulus),
        'neg': (neg_mod, xs, 1, 1, modulus),
        'mul': (mul_mod, xs, 2, 1, modulus),
        'sqr': (sqr_mod, xs, 1, 1, modulus),
        'div': (div_mod, xs, 2, 1, modulus)
    }
    test_names = ops.keys() & tests if len(tests) > 0 else ops.keys()
    fnames = map(lambda fn: modulusname + '_' + fn, test_names)
    return list(map(write_tests, fnames, [ops[fn] for fn in test_names]))


def print_usage(progname):
    print(f"""Please specify the functions for which you want to generate test cases:

    $ python3 {progname} <func_1> ... <func_n>

where each <fn_i> is one of 'add', 'sub', 'neg', 'mul', 'sqr', 'div'.
Specifying no functions will generate all of them (you will need to do this
at least once).""")


if __name__ == '__main__':
    import sys
    if len(sys.argv[1:]) > 0 and sys.argv[1] == '-h':
        print_usage(sys.argv[0])
    else:
        moduli = { # '25519': 2**255 - 19,
                  'tweedledum': 2**254 + 4707489545178046908921067385359695873,
                  'tweedledee': 2**254 + 4707489544292117082687961190295928833}
        for name, modulus in moduli.items():
            generate_tests(name, modulus, sys.argv[1:])
