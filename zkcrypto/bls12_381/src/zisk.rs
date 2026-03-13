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
    0x7908_5da1_6a58_f0e1,
    0xb3a6_3d2d_1da8_a4e6,
    0x52ae_4a89_9dae_4dcf,
    0x61a0_e07c_b5e3_f16a,
    0x5ed2_96f0_7cb1_80ac,
    0x0d97_3e7a_96c2_ac64,
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
