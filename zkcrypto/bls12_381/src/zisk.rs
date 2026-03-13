//! Ziskos precompile syscall bindings for BLS12-381 field and curve operations.
//!
//! This module provides low-level access to ziskos precompile syscalls for accelerating
//! BLS12-381 cryptographic operations when running inside the zisk zkVM.

#![allow(unsafe_code)]

// ============================================================================
// Constants
// ============================================================================

/// BLS12-381 base field modulus p
pub const MODULUS: [u64; 6] = [
    0xb9fe_ffff_ffff_aaab,
    0x1eab_fffe_b153_ffff,
    0x6730_d2a0_f6b0_f624,
    0x6477_4b84_f385_12bf,
    0x4b1b_a7b6_434b_acd7,
    0x1a01_11ea_397f_e69a,
];

/// p - 1 (used for negation: -x = x * (p-1) mod p)
pub const P_MINUS_ONE: [u64; 6] = [
    0xb9fe_ffff_ffff_aaaa,
    0x1eab_fffe_b153_ffff,
    0x6730_d2a0_f6b0_f624,
    0x6477_4b84_f385_12bf,
    0x4b1b_a7b6_434b_acd7,
    0x1a01_11ea_397f_e69a,
];

/// R = 2^384 mod p (Montgomery parameter, internal representation of 1)
pub const R_RAW: [u64; 6] = [
    0x7609_0000_0002_fffd,
    0xebf4_000b_c40c_0002,
    0x5f48_9857_53c7_58ba,
    0x77ce_5853_7052_5745,
    0x5c07_1a97_a256_ec6d,
    0x15f6_5ec3_fa80_e493,
];

/// R^-1 mod p (modular inverse of Montgomery parameter)
pub const R_INV_RAW: [u64; 6] = [
    0xf4d3_8259_380b_4820,
    0x7fe1_1274_d898_fafb,
    0x343e_a979_1495_6dc8,
    0x1797_ab14_58a8_8de9,
    0xed5e_6427_3c4f_538b,
    0x14fe_c701_e8fb_0ce9,
];

pub const ZERO_6: [u64; 6] = [0u64; 6];
pub const ONE_6: [u64; 6] = [1, 0, 0, 0, 0, 0];

// ============================================================================
// Syscall parameter structures (must match ziskos ABI exactly)
// ============================================================================

/// Parameters for `syscall_arith384_mod`: computes d = (a * b + c) mod module
#[repr(C)]
pub struct Arith384ModParams<'a> {
    pub a: &'a [u64; 6],
    pub b: &'a [u64; 6],
    pub c: &'a [u64; 6],
    pub module: &'a [u64; 6],
    pub d: &'a mut [u64; 6],
}

/// Complex field element (Fp2) for BLS12-381 complex syscalls
#[repr(C)]
pub struct Complex384 {
    pub x: [u64; 6], // c0 (real part)
    pub y: [u64; 6], // c1 (imaginary part)
}

/// Parameters for `syscall_bls12_381_complex_add`
#[repr(C)]
pub struct ComplexAddParams<'a> {
    pub f1: &'a mut Complex384,
    pub f2: &'a Complex384,
}

/// Parameters for `syscall_bls12_381_complex_sub`
#[repr(C)]
pub struct ComplexSubParams<'a> {
    pub f1: &'a mut Complex384,
    pub f2: &'a Complex384,
}

/// Parameters for `syscall_bls12_381_complex_mul`
#[repr(C)]
pub struct ComplexMulParams<'a> {
    pub f1: &'a mut Complex384,
    pub f2: &'a Complex384,
}

/// Curve point (G1 affine) for BLS12-381 curve syscalls
#[repr(C)]
pub struct Point384 {
    pub x: [u64; 6],
    pub y: [u64; 6],
}

/// Parameters for `syscall_bls12_381_curve_add`
#[repr(C)]
pub struct CurveAddParams<'a> {
    pub p1: &'a mut Point384,
    pub p2: &'a Point384,
}

// ============================================================================
// Extern syscall declarations
// ============================================================================

extern "C" {
    /// Generic 384-bit modular arithmetic: d = (a * b + c) mod module
    pub fn syscall_arith384_mod(params: &mut Arith384ModParams);

    /// Fp2 complex addition: f1 = f1 + f2
    pub fn syscall_bls12_381_complex_add(params: &mut ComplexAddParams);

    /// Fp2 complex subtraction: f1 = f1 - f2
    pub fn syscall_bls12_381_complex_sub(params: &mut ComplexSubParams);

    /// Fp2 complex multiplication: f1 = f1 * f2
    pub fn syscall_bls12_381_complex_mul(params: &mut ComplexMulParams);

    /// G1 curve point addition: p1 = p1 + p2 (affine coordinates, non-Montgomery)
    pub fn syscall_bls12_381_curve_add(params: &mut CurveAddParams);

    /// G1 curve point doubling: p1 = 2*p1 (affine coordinates, non-Montgomery)
    pub fn syscall_bls12_381_curve_dbl(p1: &mut Point384);
}

// ============================================================================
// Fp-level helpers (Montgomery form arithmetic via syscall_arith384_mod)
// ============================================================================

/// Raw modular multiply: out = a * b mod p
/// No Montgomery correction — operates on raw u64 limbs.
#[inline]
pub fn fp_mul_raw(a: &[u64; 6], b: &[u64; 6]) -> [u64; 6] {
    let mut out = [0u64; 6];
    let mut params = Arith384ModParams {
        a,
        b,
        c: &ZERO_6,
        module: &MODULUS,
        d: &mut out,
    };
    unsafe { syscall_arith384_mod(&mut params) };
    out
}

/// Montgomery multiplication: computes abR mod p from inputs aR, bR
/// (aR * bR) mod p = abR², then multiply by R^-1 to get abR
#[inline]
pub fn fp_mul_mont(a: &[u64; 6], b: &[u64; 6]) -> [u64; 6] {
    let tmp = fp_mul_raw(a, b); // abR²
    fp_mul_raw(&tmp, &R_INV_RAW) // abR² * R^-1 = abR
}

/// Montgomery squaring: computes a²R mod p from input aR
#[inline]
pub fn fp_square_mont(a: &[u64; 6]) -> [u64; 6] {
    fp_mul_mont(a, a)
}

/// Montgomery addition: computes (a+b)R mod p from inputs aR, bR
/// Uses arith384: d = aR * 1 + bR mod p = (a+b)R
#[inline]
pub fn fp_add_mont(a: &[u64; 6], b: &[u64; 6]) -> [u64; 6] {
    let mut out = [0u64; 6];
    let mut params = Arith384ModParams {
        a,
        b: &ONE_6,
        c: b,
        module: &MODULUS,
        d: &mut out,
    };
    unsafe { syscall_arith384_mod(&mut params) };
    out
}

/// Montgomery subtraction: computes (a-b)R mod p from inputs aR, bR
/// Uses arith384: d = bR * (p-1) + aR mod p = -bR + aR = (a-b)R
#[inline]
pub fn fp_sub_mont(a: &[u64; 6], b: &[u64; 6]) -> [u64; 6] {
    let mut out = [0u64; 6];
    let mut params = Arith384ModParams {
        a: b,
        b: &P_MINUS_ONE,
        c: a,
        module: &MODULUS,
        d: &mut out,
    };
    unsafe { syscall_arith384_mod(&mut params) };
    out
}

/// Montgomery negation: computes (-a)R mod p from input aR
/// Uses arith384: d = aR * (p-1) + 0 mod p = -aR
#[inline]
pub fn fp_neg_mont(a: &[u64; 6]) -> [u64; 6] {
    let mut out = [0u64; 6];
    let mut params = Arith384ModParams {
        a,
        b: &P_MINUS_ONE,
        c: &ZERO_6,
        module: &MODULUS,
        d: &mut out,
    };
    unsafe { syscall_arith384_mod(&mut params) };
    out
}

/// Convert from Montgomery form (aR) to raw form (a)
/// Computes: a = aR * R^-1 mod p
#[inline]
pub fn to_raw(x: &[u64; 6]) -> [u64; 6] {
    fp_mul_raw(x, &R_INV_RAW)
}

/// Convert from raw form (a) to Montgomery form (aR)
/// Computes: aR = a * R mod p
#[inline]
pub fn to_mont(x: &[u64; 6]) -> [u64; 6] {
    fp_mul_raw(x, &R_RAW)
}

/// Sum of products: computes (a[0]*b[0] + a[1]*b[1] + ... + a[T-1]*b[T-1]) in Montgomery form.
///
/// Inputs are in Montgomery form (a_i*R, b_i*R).
/// Uses chained arith384_mod calls:
///   t = a[0]R * b[0]R + 0       = a[0]*b[0]*R²
///   t = a[1]R * b[1]R + t       = (a[0]*b[0] + a[1]*b[1])*R²
///   ...
///   t = a[T-1]R * b[T-1]R + t   = (sum of all products)*R²
///   result = t * R_INV + 0       = (sum)*R  (Montgomery form)
///
/// Total: T+1 syscalls for T products.
#[inline]
pub fn fp_sum_of_products<const T: usize>(a: &[[u64; 6]; T], b: &[[u64; 6]; T]) -> [u64; 6] {
    let mut acc = ZERO_6;
    for i in 0..T {
        // acc = a[i]*R * b[i]*R + acc  (mod p)
        let mut out = [0u64; 6];
        let mut params = Arith384ModParams {
            a: &a[i],
            b: &b[i],
            c: &acc,
            module: &MODULUS,
            d: &mut out,
        };
        unsafe { syscall_arith384_mod(&mut params) };
        acc = out;
    }
    // acc = (sum of products) * R².  Multiply by R_INV to get (sum)*R.
    fp_mul_raw(&acc, &R_INV_RAW)
}

// ============================================================================
// Scalar field (256-bit) operations via arith256_mod syscall
// ============================================================================

/// Parameters for `syscall_arith256_mod`: computes d = (a * b + c) mod module
#[repr(C)]
pub struct Arith256ModParams<'a> {
    pub a: &'a [u64; 4],
    pub b: &'a [u64; 4],
    pub c: &'a [u64; 4],
    pub module: &'a [u64; 4],
    pub d: &'a mut [u64; 4],
}

extern "C" {
    pub fn syscall_arith256_mod(params: &mut Arith256ModParams);
}

/// BLS12-381 scalar field modulus q
pub const SCALAR_MODULUS: [u64; 4] = [
    0xffff_ffff_0000_0001,
    0x53bd_a402_fffe_5bfe,
    0x3339_d808_09a1_d805,
    0x73ed_a753_299d_7d48,
];

/// R^-1 mod q where R = 2^256 (for Montgomery correction)
pub const SCALAR_R_INV_RAW: [u64; 4] = [
    0x13f7_5b69_fe75_c040,
    0xab6f_ca8f_09dc_705f,
    0x7204_078a_4f77_266a,
    0x1bbe_8693_3000_9d57,
];

const SCALAR_ZERO_4: [u64; 4] = [0, 0, 0, 0];

/// Raw modular multiply: d = (a * b) mod q
#[inline]
pub fn scalar_mul_raw(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    let mut out = [0u64; 4];
    let mut params = Arith256ModParams {
        a,
        b,
        c: &SCALAR_ZERO_4,
        module: &SCALAR_MODULUS,
        d: &mut out,
    };
    unsafe { syscall_arith256_mod(&mut params) };
    out
}

/// Montgomery multiply: given aR, bR, returns abR mod q (2 syscalls)
#[inline]
pub fn scalar_mul_mont(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    // Step 1: (aR * bR) mod q = abR² mod q
    let ab_r2 = scalar_mul_raw(a, b);
    // Step 2: abR² * R⁻¹ mod q = abR mod q
    scalar_mul_raw(&ab_r2, &SCALAR_R_INV_RAW)
}

/// Montgomery square: given aR, returns a²R mod q (2 syscalls)
#[inline]
pub fn scalar_square_mont(a: &[u64; 4]) -> [u64; 4] {
    scalar_mul_mont(a, a)
}
