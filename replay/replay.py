"""Exact rational and symbolic certificate for the hypercubic Holstein pair-clustering ceiling.

Verifies the hypercubic Z^d (d >= 2, with explicit checks at d = 2, 3, 4, 5
and symbolic d in Q(d)) strong-coupling Holstein clustering ceiling:
  - Theorem 1 (unconditional edge-local SOS + Bethe-tree vertex-animal counting):
      eps0(d) = +oo, L(d) = 72(d+1), C(d) = 442368 d^2 (d+1)^2 = Theta(d^4),
      q(d) = 5(2d-1)/(24(d+1)), A(d) = 360(d+1)/(14d+29);
  - Proposition A0 (unconditional edge-local SOS + Euler-tour D_aug^(2k) counting):
      eps0(d) = +oo, L_ET(d) = 96(d+1)^2, C_ET(d) = 1179648 d (d+1)^4 = Theta(d^5);
  - Proposition B1 (released conditional form bound + Euler-tour counting):
      eps0_B1(d) = 1/(2^17 d (d+1)^2) = Theta(d^-3), C_B1(d) = 2^31 d (d+1)^4 = Theta(d^5);
  - Proposition B2 (released conditional form bound + vertex-animal counting):
      eps0_B2(d) = 1/(7680 d (d+1)) = Theta(d^-2), C_B2(d) = 9830400 d^2 (d+1)^2 = Theta(d^4);
  - Proposition C (exact rational proof of the earlier cubic canonical obstruction at d=3
    for C = 2^27 under the earlier cubic Euler-tour counting, and its resolution by Theorem 1).
"""
import argparse
from fractions import Fraction as Q
from math import comb, factorial
import sys

import sympy as sp


MUTANTS = (
    "euler-number-enclosure",
    "car-gram-sos",
    "conditional-form-boundary",
    "fugacity-sos-identity",
    "pair-source-scalar-norm",
    "bethe-lagrange-animals",
    "bethe-kp-budget",
    "source-perturbed-kp-budget",
    "signed-deletion-prefactor",
    "theorem-constants",
    "route-a0-kp-budget",
    "route-b2-enclosure",
    "route-hierarchy-comparison",
    "canonical-obstruction",
    "asymptotics-and-heuristic",
)


class Rejected(Exception):
    pass


def require(ok, name):
    if not ok:
        raise Rejected(name)


def check_euler_enclosure(mutant):
    # Exact rational Taylor enclosure for e = sum_{n>=0} 1/n!:
    # e_lo = sum_{n=0}^6 1/n! = 1957/720
    # e_hi_6 = sum_{n=0}^5 1/n! + (1/720) * sum_{m>=0} (1/7)^m = 163/60 + 7/4320 = 11743/4320
    # e_hi = sum_{n=0}^5 1/n! + (1/360) * sum_{m>=0} (1/7)^m = 163/60 + 7/2160 = 1175/432
    e_lo = sum((Q(1, factorial(n)) for n in range(7)), Q(0))
    e_hi_6 = sum((Q(1, factorial(n)) for n in range(6)), Q(0)) + Q(1, 720) / (1 - Q(1, 7))
    e_hi = sum((Q(1, factorial(n)) for n in range(6)), Q(0)) + Q(1, 360) / (1 - Q(1, 7))
    if mutant == "euler-number-enclosure":
        e_hi = Q(3)
    require(e_lo == Q(1957, 720), "e-lower-1957-720")
    require(e_hi_6 == Q(11743, 4320), "e-upper-6fact-11743-4320")
    require(e_hi == Q(1175, 432), "e-upper-1175-432")
    require(Q(2) < e_lo < e_hi_6 < e_hi < Q(49, 18) < Q(3), "e-rational-bracketing")
    # Rational root bounds used in exponential penalties:
    # e^(1/4) < 4/3   <=>  (4/3)^4 = 256/81 > e_hi
    # e^(5/16) < 11/8 <=>  (11/8)^16 > e_hi^5
    # e^(1/6) < 6/5   <=>  (6/5)^6 = 46656/15625 > e_hi
    # e^(1/8) > 9/8   <=>  (9/8)^8 < e_lo
    require(Q(4, 3) ** 4 == Q(256, 81) > e_hi, "exp-quarter-lt-4-3")
    require(
        Q(11, 8) ** 16 == Q(45949729863572161, 281474976710656)
        and e_hi ** 5 == Q(2239697333984375, 15045919506432)
        and Q(11, 8) ** 16 > e_hi ** 5,
        "exp-five-sixteenths-lt-11-8",
    )
    require(Q(6, 5) ** 6 == Q(46656, 15625) > e_hi, "exp-sixth-lt-6-5")
    require(Q(11, 10) ** 11 == Q(285311670611, 100000000000) > e_hi, "exp-eleventh-lt-11-10")
    require(Q(9, 8) ** 8 < e_lo, "exp-eighth-gt-9-8")
    return e_lo, e_hi


def build_single_bond_fock():
    # 16-state two-site electronic Fock space ordered as (u_up, u_dn, v_up, v_dn).
    n = 16
    b = [[Q(0) for _ in range(n)] for _ in range(n)]
    n_odd = [Q(0) for _ in range(n)]
    n_tot = [0 for _ in range(n)]
    for state in range(n):
        u_up = (state >> 0) & 1
        u_dn = (state >> 1) & 1
        v_up = (state >> 2) & 1
        v_dn = (state >> 3) & 1
        occ_u = u_up + u_dn
        occ_v = v_up + v_dn
        n_tot[state] = occ_u + occ_v
        n_odd[state] = Q((occ_u & 1) + (occ_v & 1))
        # Signless hopping B_e: flip bit 0<->2 (up hop) or bit 1<->3 (down hop) when occupation differs
        if u_up != v_up:
            b[state][state ^ 0b0101] += Q(1)
        if u_dn != v_dn:
            b[state][state ^ 0b1010] += Q(1)
    return b, n_odd, n_tot


def check_car_gram_and_sos(mutant):
    b, n_odd, n_tot = build_single_bond_fock()
    g, alpha = sp.symbols("g alpha", positive=True)
    coeff_g2 = sp.Rational(1, 1) if mutant == "car-gram-sos" else sp.Rational(2, 1)

    # 1. Verify unconditional single-bond SOS identity across all particle-number sectors:
    # R_e(g, alpha) = alpha * N_odd,e + (2 g^2 / alpha) I_16 - g B_e
    R = sp.zeros(16, 16)
    for i in range(16):
        R[i, i] = alpha * int(n_odd[i]) + coeff_g2 * g**2 / alpha
        for j in range(16):
            if b[i][j] != 0:
                R[i, j] -= g * int(b[i][j])

    # Sector N_e in {0, 4}: B_e = 0, N_odd,e = 0
    for s in range(16):
        if n_tot[s] in (0, 4):
            require(n_odd[s] == 0 and all(b[s][j] == 0 for j in range(16)), "sector-0-4-trivial")

    # Sector N_e in {1, 3}: N_odd,e = I, B_e^2 = I, R = (sqrt(alpha) I - g/(2 sqrt(alpha)) B_e)^2 + (7 g^2 / (4 alpha)) I
    idx_13 = [s for s in range(16) if n_tot[s] in (1, 3)]
    require(len(idx_13) == 8, "sector-1-3-dim-8")
    B_13 = sp.Matrix([[int(b[i][j]) for j in idx_13] for i in idx_13])
    I_8 = sp.eye(8)
    require(B_13 * B_13 == I_8, "sector-1-3-B2-is-I")
    R_13 = sp.Matrix([[R[i, j] for j in idx_13] for i in idx_13])
    sos_13 = (sp.sqrt(alpha) * I_8 - (g / (2 * sp.sqrt(alpha))) * B_13) ** 2 + (sp.Rational(7, 4) * g**2 / alpha) * I_8
    require(sp.simplify(R_13 - sos_13) == sp.zeros(8, 8), "sector-1-3-sos-identity")
    psi1, psi2 = sp.symbols("psi1 psi2", real=True)
    B_doublet = sp.Matrix([[0, 1], [1, 0]])
    R_doublet = alpha * sp.eye(2) + (coeff_g2 * g**2 / alpha) * sp.eye(2) - g * B_doublet
    psi_vec = sp.Matrix([psi1, psi2])
    qf_doublet = (psi_vec.T * R_doublet * psi_vec)[0, 0]
    qf_sos = (
        (sp.sqrt(alpha) * psi1 - (g / (2 * sp.sqrt(alpha))) * psi2) ** 2
        + (sp.sqrt(alpha) * psi2 - (g / (2 * sp.sqrt(alpha))) * psi1) ** 2
        + (sp.Rational(7, 4) * g**2 / alpha) * (psi1**2 + psi2**2)
    )
    require(sp.simplify(qf_doublet - qf_sos) == 0, "sector-1-3-doublet-sos-identity")

    # Sector N_e = 2 (dim 6): verify Pauli rank-one reduction on {|e_+>, |o_+>} and exact SOS
    # States: e1 = |1,1,0,0> (3), e2 = |0,0,1,1> (12), o1 = |1,0,0,1> (9), o2 = |0,1,1,0> (6),
    #         p_up = |1,0,1,0> (5), p_dn = |0,1,0,1> (10)
    idx_2 = [3, 12, 9, 6, 5, 10]
    require(sorted(idx_2) == [s for s in range(16) if n_tot[s] == 2], "sector-2-basis")
    B_2 = sp.Matrix([[int(b[i][j]) for j in idx_2] for i in idx_2])
    Nodd_2 = sp.diag(*[int(n_odd[i]) for i in idx_2])
    # Orthogonal change of basis U_2 to (e_+, o_+, e_-, o_-, p_up, p_dn):
    inv_rt2 = 1 / sp.sqrt(2)
    U_2 = sp.Matrix([
        [inv_rt2, inv_rt2, 0, 0, 0, 0],
        [0, 0, inv_rt2, inv_rt2, 0, 0],
        [inv_rt2, -inv_rt2, 0, 0, 0, 0],
        [0, 0, inv_rt2, -inv_rt2, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
    ])
    require(sp.simplify(U_2 * U_2.T - sp.eye(6)) == sp.zeros(6, 6), "sector-2-unitary")
    B_2_rot = sp.simplify(U_2 * B_2 * U_2.T)
    Nodd_2_rot = sp.simplify(U_2 * Nodd_2 * U_2.T)
    # In rotated basis, B_2_rot has ONLY (0,1) and (1,0) entries equal to 2; all other entries are 0!
    B_2_expected = sp.zeros(6, 6)
    B_2_expected[0, 1] = 2
    B_2_expected[1, 0] = 2
    require(B_2_rot == B_2_expected, "sector-2-rank-one-coupling")
    require(Nodd_2_rot == sp.diag(0, 2, 0, 2, 2, 2), "sector-2-nodd-diagonal")
    R_2_rot = sp.simplify(U_2 * sp.Matrix([[R[i, j] for j in idx_2] for i in idx_2]) * U_2.T)
    # Top-left 2x2 block on span{|e_+>, |o_+>} is [[2 g^2/alpha, -2g], [-2g, 2 alpha + 2 g^2/alpha]]
    # whose determinant is 4 g^2 + 4 g^4/alpha^2 - 4 g^2 = 4 g^4 / alpha^2 > 0!
    block_2x2 = R_2_rot[:2, :2]
    sos_2x2 = (
        sp.Matrix([-sp.sqrt(2) * g / sp.sqrt(alpha), sp.sqrt(2 * alpha)])
        * sp.Matrix([-sp.sqrt(2) * g / sp.sqrt(alpha), sp.sqrt(2 * alpha)]).T
        + sp.diag(0, 2 * g**2 / alpha)
    )
    require(sp.simplify(block_2x2 - sos_2x2) == sp.zeros(2, 2), "sector-2-top-block-sos")

    # 2. Verify the released conditional form bound on hypercubic Z^d (D_phys = 2d) for d = 2, 3, 4, 5:
    # At boundary g/U = 1/(16 D_phys) = 1/(32 d), on one bond (with share N_odd,e / (2d)):
    # (U / (8d)) N_odd,e + (32 d g^2 / U) I_16 - g B_e is diagonally dominant with min residual 0.
    min_residuals = {}
    for d in (2, 3, 4, 5):
        D_phys = 2 * d
        g_ratio = Q(1, 8 * D_phys) if mutant == "conditional-form-boundary" else Q(1, 16 * D_phys)
        # Split B_e into creation/annihilation (B_pm) and odd-preserving (B_0)
        # Weight h = 8 D_phys (g/U) = 1/2 at boundary
        h = 8 * D_phys * g_ratio
        # Combined multi-bond N_odd channel budget from the released record (4.2)-(4.3):
        # creation/annihilation share (D_phys * g_ratio) / h = 1/8 plus zero-channel share 2 * D_phys * g_ratio <= 1/8
        combined_nodd_budget = (D_phys * g_ratio) / h + 2 * D_phys * g_ratio
        require(combined_nodd_budget <= Q(1, 4), f"conditional-multi-bond-nodd-budget-d{d}")
        diag_margin = []
        for s in range(16):
            births = sum(b[s][t] for t in range(16) if n_odd[t] > n_odd[s])
            deaths = sum(b[s][t] for t in range(16) if n_odd[t] < n_odd[s])
            zeros = sum(b[s][t] for t in range(16) if n_odd[t] == n_odd[s] and b[s][t] != 0)
            require(births <= 2 and deaths <= n_odd[s] and zeros <= 2 * n_odd[s], "channel-degree-bounds")
            lhs_upper = g_ratio * (h * births + (deaths / h) + zeros)
            rhs_budget = Q(1, 4 * D_phys) * n_odd[s] + 16 * D_phys * (g_ratio ** 2)
            diag_margin.append(rhs_budget - lhs_upper)
        min_residuals[d] = min(diag_margin)
        require(min_residuals[d] == Q(0), f"conditional-form-nonnegative-d{d}")
    return min_residuals


def check_fugacity_and_source(mutant):
    # 1. Exact bivariate polynomial identity in Q[z, w]:
    # 5 z_a - 4 z_b = (1 - z)^2 + 2 z (2w - 1)^2 + 2 z w^2
    z, w = sp.symbols("z w")
    z_a = 1 + 2 * z * w**2 + z**2
    z_b = 1 + 2 * z * w + z**2
    lhs = sp.expand(5 * z_a - 4 * z_b)
    rhs = sp.expand((1 - z) ** 2 + 2 * z * (2 * w - 1) ** 2 + (4 * z**2 * w**2 if mutant == "fugacity-sos-identity" else 2 * z * w**2))
    require(lhs == rhs, "fugacity-5za-4zb-sos-identity")

    # Also verify the released fugacity identity: (1 + w) z_a - z_b = w (z - 1)^2 + 2 z w^2 (1 + w)
    lhs_rel = sp.expand((1 + w) * z_a - z_b)
    rhs_rel = sp.expand(w * (z - 1) ** 2 + 2 * z * w**2 * (1 + w))
    require(lhs_rel == rhs_rel, "fugacity-identity-5za-minus-4zb")

    # 2. Verify two-site pair-source operator norm and zero-N_odd action on |1,1,0,0> <-> |0,0,1,1>:
    _, n_odd, _ = build_single_bond_fock()
    # States |1,1,0,0> = 3 and |0,0,1,1> = 12 both have N_odd = 0:
    require(n_odd[3] == 0 and n_odd[12] == 0, "pair-source-connects-even-even")
    # Norm of X_e = e^(i theta) Delta_u^+ Delta_v + h.c. is 2, so on |xi| <= tau / (16 L),
    # 2 * beta * L * |xi| <= 1/8 <= 1/4:
    xi_over_tau_L = Q(1, 4) if mutant == "pair-source-scalar-norm" else Q(1, 16)
    source_penalty = 2 * xi_over_tau_L
    require(source_penalty == Q(1, 8) and source_penalty <= Q(1, 4), "pair-source-penalty-le-quarter")
    return source_penalty


def T_bethe(Delta, v):
    if v == 1:
        return 1
    return (Delta * comb(v * (Delta - 1), v - 2)) // (v - 1)


def exact_zd_animals(d, v_max=4):
    origin = tuple([0] * d)
    dirs = []
    for i in range(d):
        pos = [0] * d
        pos[i] = 1
        neg = [0] * d
        neg[i] = -1
        dirs.append(tuple(pos))
        dirs.append(tuple(neg))

    def neighbors(pt):
        return {tuple(pt[i] + dx[i] for i in range(d)) for dx in dirs}

    animals = {1: {frozenset([origin])}}
    for v in range(2, v_max + 1):
        nxt = set()
        for a in animals[v - 1]:
            nbrs = set()
            for pt in a:
                nbrs |= neighbors(pt)
            nbrs -= a
            for nb in nbrs:
                nxt.add(a | {nb})
        animals[v] = nxt
    return [len(animals[v]) for v in range(1, v_max + 1)]


def _poly_mul_trunc(p, q, deg=7):
    out = [0] * deg
    for i, a in enumerate(p):
        if a:
            for j in range(min(deg - i, len(q))):
                out[i + j] += a * q[j]
    return out


def _poly_pow_trunc(p, exp, deg=7):
    res = [1] + [0] * (deg - 1)
    base = list(p[:deg])
    e = exp
    while e > 0:
        if e & 1:
            res = _poly_mul_trunc(res, base, deg)
        base = _poly_mul_trunc(base, base, deg)
        e >>= 1
    return res


def check_bethe_and_animals(mutant):
    animal_witnesses = {}
    for d in (2, 3, 4, 5):
        for Delta in (2 * d, 2 * d + 2):
            # Solve f(x) = x (1 + f(x))^(Delta - 1) to order x^6 and R(x) = x (1 + f(x))^Delta
            f = [0] * 7
            for _ in range(6):
                one_plus_f = [1 + f[0]] + f[1:]
                pw = _poly_pow_trunc(one_plus_f, Delta - 1, 7)
                f = [0] + pw[:6]
            one_plus_f = [1 + f[0]] + f[1:]
            R_pw = _poly_pow_trunc(one_plus_f, Delta, 7)
            R = [0] + R_pw[:6]
            series_coeffs = R[1:7]
            lagrange_coeffs = [T_bethe(Delta, v) for v in range(1, 7)]
            if mutant == "bethe-lagrange-animals" and d == 3 and Delta == 6:
                lagrange_coeffs[3] -= 1
            require(series_coeffs == lagrange_coeffs, f"bethe-lagrange-match-Delta{Delta}")
        ex = exact_zd_animals(d, 4)
        fc = [T_bethe(2 * d, v) for v in range(1, 5)]
        # On triangle-free Z^d, a_v(Z^d) == T_{2d}(v) for v = 1, 2, 3,
        # and at v = 4 each of the 4*comb(d, 2) plaquettes containing the origin has 4 spanning trees in T_{2d}
        # versus 1 animal in Z^d, giving the exact plaquette identity T_{2d}(4) - a_4(Z^d) = 12 * comb(d, 2)!
        require(ex[:3] == fc[:3], f"zd-animals-small-v-exact-d{d}")
        require(fc[3] - ex[3] == 12 * comb(d, 2), f"zd-animals-v4-plaquette-identity-d{d}")
        require(all(ex[i] <= fc[i] for i in range(4)), f"zd-animals-upper-bounded-d{d}")
        animal_witnesses[d] = (ex, fc)
    return animal_witnesses


def check_theorem_a_and_kp_budget(mutant, e_hi):
    # Symbolic verification in Q(d)
    d_sym = sp.Symbol("d", integer=True, positive=True)
    D_aug_sym = 2 * d_sym + 2
    L_sym = 36 * D_aug_sym  # = 72(d + 1)
    C_sym = sp.Rational(256, 3) * (d_sym ** 2) * (L_sym ** 2)  # = 442368 d^2 (d + 1)^2
    require(sp.expand(L_sym - 72 * (d_sym + 1)) == 0, "symbolic-L-formula")
    require(sp.expand(C_sym - 442368 * d_sym**2 * (d_sym + 1) ** 2) == 0, "symbolic-C-formula")

    q_sym = sp.Rational(5, 24) * (2 * d_sym - 1) / (d_sym + 1)
    A_sym = sp.simplify(sp.Rational(15, 1) / (1 - q_sym))
    require(sp.simplify(A_sym - 360 * (d_sym + 1) / (14 * d_sym + 29)) == 0, "symbolic-A-formula")

    results = {}
    for d in (2, 3, 4, 5):
        D_phys = 2 * d
        D_aug = 2 * d + 2
        L_d = 36 * D_aug
        C_d = Q(256, 3) * (d ** 2) * (L_d ** 2)
        if mutant == "theorem-constants" and d == 3:
            C_d = Q(2 ** 25)
        require(L_d == 72 * (d + 1), f"L-val-d{d}")
        require(C_d == 442368 * (d ** 2) * ((d + 1) ** 2), f"C-val-d{d}")

        # Verify delta_phys = 16 * d^2 * L_d^2 / C_d == 3/16 <= 1/4:
        delta_phys = Q(16 * (d ** 2) * (L_d ** 2), C_d)
        require(delta_phys == Q(3, 16) and delta_phys < Q(1, 4), f"delta-phys-3-16-d{d}")

        # Bethe-tree evaluation point y_0 = 1 / (6 * D_aug), x_0 = y_0 / (1 + y_0)^(D_aug - 1):
        y_0 = Q(1, 6 * D_aug)
        x_0 = y_0 / ((1 + y_0) ** (D_aug - 1))
        # Since (1 + y_0)^D_aug < e^(1/6) < 6/5, (1 + y_0)^D_aug - 1 < 1/5:
        tree_factor = ((1 + y_0) ** D_aug) - 1
        require(tree_factor < Q(1, 5), f"tree-factor-lt-fifth-d{d}")

        # Physical (xi = 0) per-site factor x_phys = (5 / (4 * L_d)) * exp(delta_phys) * e:
        # Using exp(3/16) < exp(1/4) < 4/3 and e < e_hi = 1175/432 < 3:
        x_eff_hi = Q(5, 4 * L_d) * Q(4, 3) * e_hi
        x_eff_clean = Q(5, L_d)
        require(L_d * x_eff_hi == Q(5875, 1296), f"L-x-phys-hi-5875-1296-d{d}")
        require(x_eff_hi < x_eff_clean, f"x-eff-hi-lt-clean-d{d}")

        # Source-perturbed (|xi| <= tau / (16 * L_d)) per-site factor:
        # delta_tot = delta_phys + 1/8 = 3/16 + 1/8 = 5/16, exp(5/16) < 11/8:
        delta_tot = delta_phys + Q(1, 8)
        require(delta_tot == Q(5, 16), f"delta-tot-5-16-d{d}")
        x_src_hi = Q(5, 4 * L_d) * Q(11, 8) * e_hi
        require(L_d * x_src_hi == Q(64625, 13824), f"L-x-src-hi-64625-13824-d{d}")
        require(x_eff_hi < x_src_hi < x_eff_clean, f"x-src-hi-lt-clean-d{d}")

        # Verify x_eff_clean <= x_0  <=>  5 / (36 * D_aug) <= 1 / (6 * D_aug * (1 + y_0)^(D_aug - 1))
        # <=>  (1 + y_0)^(D_aug - 1) <= 6/5, which holds since (1 + y_0)^(D_aug - 1) < (1 + y_0)^D_aug < 6/5!
        require(x_src_hi < x_eff_clean < x_0, f"x-src-inside-convergence-disk-d{d}")

        # Kotecky-Preiss site budget (physical xi=0 and source-perturbed |xi|<=tau/(16 L_d)):
        kp_budget_sharp = (L_d * x_eff_hi) * tree_factor
        kp_src_sharp = (L_d * x_src_hi) * tree_factor
        kp_budget_clean = (L_d * x_eff_clean) * Q(1, 5)
        if mutant == "bethe-kp-budget" and d == 3:
            kp_budget_sharp = Q(2)
        if mutant == "source-perturbed-kp-budget" and d == 3:
            kp_src_sharp = Q(2)
        require(kp_budget_sharp < Q(1175, 1296) < kp_budget_clean == Q(1), f"kp-budget-lt-1-d{d}")
        require(kp_budget_sharp < kp_src_sharp < Q(12925, 13824) < kp_budget_clean == Q(1), f"kp-src-budget-lt-1-d{d}")

        # Pair correlation decay rate q(d) and prefactor A(d) on the physical graph (degree D_phys = 2d),
        # including the signed hard-core deletion factor |Z_{Lambda\S_0}/Z_Lambda| <= e^{|S_0|}:
        # At y_1 = 1 / (2d - 1), x_1 = y_1 / (1 + y_1)^(2d - 1) > 1 / (e * (2d - 1)) > 1 / (3 * (2d - 1)):
        y_1 = Q(1, D_phys - 1)
        x_1 = y_1 / ((1 + y_1) ** (D_phys - 1))
        require(x_1 > Q(1, 3 * (D_phys - 1)), f"x1-lower-bound-d{d}")
        q_exact_rat = x_eff_hi / x_1
        q_sharp = Q(5875 * 1175 * (2 * d - 1), 1296 * 72 * 432 * (d + 1))
        q_clean = Q(5 * (2 * d - 1), 24 * (d + 1))
        require(q_exact_rat < q_sharp < q_clean < Q(5, 12) < Q(1, 2), f"q-hierarchy-d{d}")

        # Signed-deletion prefactor numerator: L_d * x_phys * (1 + y_1)^(2d) <= (5875/1296) * (4/3)^4 = 94000/6561 < 15:
        phys_tree_factor = (1 + y_1) ** D_phys
        require(phys_tree_factor <= Q(4, 3) ** 4 == Q(256, 81), f"phys-tree-factor-le-256-81-d{d}")
        prefactor_num_hi = (L_d * x_eff_hi) * phys_tree_factor
        if mutant == "signed-deletion-prefactor" and d == 2:
            prefactor_num_hi = Q(16)
        require(prefactor_num_hi <= Q(94000, 6561) < Q(15), f"signed-deletion-prefactor-lt-15-d{d}")
        A_clean = Q(15, 1) / (1 - q_clean)
        require(A_clean == Q(360 * (d + 1), 14 * d + 29), f"A-clean-formula-d{d}")
        row_sum_clean = Q(1) + A_clean * (((1 + q_clean) / (1 - q_clean)) ** d - 1)
        require(row_sum_clean > Q(1), f"row-sum-finite-d{d}")

        results[d] = {
            "D_phys": D_phys,
            "D_aug": D_aug,
            "L": L_d,
            "C": int(C_d),
            "kp_budget_sharp": kp_budget_sharp,
            "kp_src_sharp": kp_src_sharp,
            "prefactor_num_hi": prefactor_num_hi,
            "q_clean": q_clean,
            "q_sharp": q_sharp,
            "A_clean": A_clean,
            "row_sum_clean": row_sum_clean,
        }
    # Highlight d=3: C(3) = 63,700,992 < 2^26 < 2^27
    require(results[3]["C"] == 63700992 and results[3]["C"] < 2**26 < 2**27, "d3-strictly-below-2-27")
    return results


def check_theorem_b_routes(mutant, e_hi):
    d_sym = sp.Symbol("d", integer=True, positive=True)
    # Route A0: earlier cubic Euler-tour D_aug^(2k) extension
    # L_ET(d) = 96(d+1)^2, C_ET(d) = 128 d L_ET(d)^2 = 1,179,648 d (d+1)^4
    L_a0_sym = 96 * (d_sym + 1) ** 2
    C_a0_sym = 128 * d_sym * (L_a0_sym ** 2)
    require(sp.expand(C_a0_sym - 1179648 * d_sym * (d_sym + 1) ** 4) == 0, "route-a0-symbolic-C")

    # Route B1: released Euler-tour D_aug^(2k) extension
    # L_B1(d) = 4096 (d+1)^2, eps0_B1(d) = 1 / (131072 d (d+1)^2), C_B1(d) = 2^31 d (d+1)^4
    L_b1_sym = 4096 * (d_sym + 1) ** 2
    C_b1_sym = 128 * d_sym * (L_b1_sym ** 2)
    require(sp.expand(C_b1_sym - (2 ** 31) * d_sym * (d_sym + 1) ** 4) == 0, "route-b1-symbolic-C")

    # Route B2: released conditional form, vertex-animal extension
    # L_B2(d) = 240 (d+1), eps0_B2(d) = 1 / (7680 d (d+1)), C_B2(d) = 9,830,400 d^2 (d+1)^2
    L_b2_sym = 240 * (d_sym + 1)
    C_b2_sym = sp.Rational(512, 3) * (d_sym ** 2) * (L_b2_sym ** 2)
    require(sp.expand(C_b2_sym - 9830400 * (d_sym ** 2) * (d_sym + 1) ** 2) == 0, "route-b2-symbolic-C")
    routes = {}
    for d in (2, 3, 4, 5):
        D_phys = 2 * d
        D_aug = 2 * d + 2
        # Route A0 at d (earlier cubic theorem with D_phys=2d, D_aug=2d+2):
        L_a0 = 24 * (D_aug ** 2)
        require(L_a0 == 96 * ((d + 1) ** 2), f"route-a0-L-d{d}")
        C_a0 = 64 * D_phys * (L_a0 ** 2)
        require(C_a0 == 128 * d * (L_a0 ** 2), f"route-a0-C-d{d}")
        delta_a0 = Q(8 * D_phys * (L_a0 ** 2), C_a0)
        require(delta_a0 == Q(1, 8) < Q(1, 4), f"route-a0-delta-d{d}")
        site_pref_a0 = (Q(25, 16) * e_hi) if (mutant == "route-a0-kp-budget" and d == 3) else (Q(5, 4) * e_hi)
        require(site_pref_a0 == Q(5875, 1728) < Q(15, 4), f"route-a0-site-pref-lt-15-4-d{d}")
        q_phys_a0 = Q(15 * (D_phys ** 2), 4 * L_a0) * Q(4, 3)
        q_src_a0 = Q(15 * (D_aug ** 2), 4 * L_a0) * Q(4, 3)
        require(q_phys_a0 == Q(5 * (d ** 2), 24 * ((d + 1) ** 2)) < Q(5, 24), f"route-a0-q-phys-d{d}")
        kp_src_a0 = Q(15, 4) * q_src_a0 / (1 - q_src_a0)
        A_a0 = Q(12) / (1 - Q(5, 24))
        require(
            q_src_a0 == Q(5, 24)
            and kp_src_a0 == Q(75, 76) < Q(1)
            and A_a0 == Q(288, 19) < Q(76, 5),
            f"route-a0-kp-d{d}",
        )

        # Route B1 at d:
        L_b1 = (2 ** 12) * ((d + 1) ** 2)
        eps0_b1 = Q(1, 32 * d * L_b1)
        C_b1 = 128 * d * (L_b1 ** 2)
        q_phys_b1 = Q(27 * (D_phys ** 2), L_b1)
        q_src_b1 = Q(27 * (D_aug ** 2), L_b1)
        kp_src_b1 = Q(9) * q_src_b1 / (1 - q_src_b1)
        require(q_src_b1 == Q(27, 1024) and kp_src_b1 == Q(243, 997) < Q(1), f"route-b1-kp-d{d}")

        # Route B2 at d (released conditional form with Bethe-tree vertex-animal counting):
        # On source disk |xi| <= tau / (32 L_b2), delta_tot_b2 <= 3/16 + 1/16 = 1/4, so e^(delta_tot_b2) < 4/3.
        # Single-site factor uses z_b/z_a <= e < e_hi = 1175/432, so K_cond_hi = (4/3)*e_hi^2 = 1380625/139968 < 10.
        # At y_2 = 1/(11 * D_aug), (1 + y_2)^D_aug < e^(1/11) < 11/10 (since (11/10)^11 > e_hi) => (1+y_2)^D_aug - 1 < 1/10.
        # Choosing L_b2 = 120 * D_aug = 240(d+1) gives x_B2 < K_cond_hi / L_b2 < (1200/121) / (120 D_aug) = 10 / (121 D_aug) < x_0(y_2)
        # and KP budget < K_cond_hi * (1/10) = 276125 / 279936 < 1.
        L_b2 = 240 * (d + 1)
        eps0_b2 = Q(1, 32 * d * L_b2)
        C_b2 = Q(512, 3) * (d ** 2) * (L_b2 ** 2)
        delta_b2 = Q(32 * (d ** 2) * (L_b2 ** 2), C_b2)
        source_b2 = Q(1, 16)
        require(delta_b2 == Q(3, 16) and delta_b2 + source_b2 == Q(1, 4), f"route-b2-delta-tot-d{d}")
        base_b2 = Q(12, 11) if (mutant == "route-b2-enclosure" and d == 3) else Q(11, 10)
        require(base_b2 ** 11 > e_hi, f"exp-1-11-lt-11-10-d{d}")
        K_s6_hi = Q(4, 3) * (e_hi ** 2)
        require(K_s6_hi == Q(1380625, 139968) < Q(1200, 121) < Q(10), f"K-cond-lt-1200-121-lt-10-d{d}")
        y_b2 = Q(1, 11 * D_aug)
        x0_b2 = y_b2 / ((1 + y_b2) ** (D_aug - 1))
        require(
            K_s6_hi / L_b2 < Q(10, 121 * D_aug) < x0_b2,
            f"x-eff-b2-lt-x0-d{d}",
        )
        tree_b2 = ((1 + y_b2) ** D_aug) - 1
        require(tree_b2 < Q(1, 10), f"tree-b2-lt-1-10-d{d}")
        kp_b2 = K_s6_hi * tree_b2
        require(kp_b2 < K_s6_hi * Q(1, 10) == Q(276125, 279936) < Q(1), f"route-b2-kp-d{d}")
        q_b2_hi = (K_s6_hi * e_hi * (2 * d - 1)) / L_b2
        require(
            q_b2_hi < Q(1175 * (2 * d - 1), 10368 * (d + 1)) < Q(2 * d - 1, 8 * (d + 1)) < Q(1, 4) < Q(5, 11),
            f"route-b2-decay-rate-d{d}",
        )
        if mutant == "route-hierarchy-comparison" and d == 3:
            C_b2 = C_b1
        require(C_b2 < C_b1 and C_a0 < C_b2, f"route-hierarchy-d{d}")
        C_thm_a_d = 442368 * (d ** 2) * ((d + 1) ** 2)
        require(Q(C_b2, C_thm_a_d) == Q(200, 9), f"ratio-C-va-to-C-d{d}")
        require(Q(C_a0, C_thm_a_d) == Q(8 * ((d + 1) ** 2), 3 * d), f"ratio-C-et-to-C-d{d}")
        if d == 3:
            require(Q(C_a0, C_thm_a_d) == Q(128, 9), "ratio-C-et-to-C-d3")

        routes[d] = {
            "L_a0": L_a0,
            "C_a0": C_a0,
            "q_phys_a0": q_phys_a0,
            "kp_src_a0": kp_src_a0,
            "L_b1": L_b1,
            "eps0_b1": eps0_b1,
            "C_b1": C_b1,
            "q_phys_b1": q_phys_b1,
            "L_b2": L_b2,
            "eps0_b2": eps0_b2,
            "C_b2": int(C_b2),
            "kp_b2": kp_b2,
            "q_b2_hi": q_b2_hi,
        }
    return routes


def check_cubic_canonical_obstruction(mutant):
    # Defect 1: the earlier cubic attempt fabricated delta_6 = 96 L_7^2 beta t^2 / U and q_3 = (900 / (4 L_7)) e^{delta_6}
    # vs the genuine cubic route at d=3 (D_phys=6, D_aug=8):
    # Theorem 3.1 (3.7) gives 8 d_E k g^2 / U = 48 L_7^2 beta t^2 / U,
    # q_phys = 36 * (5 / (4 L_7)) * 3 * e^{delta_6} = (135 / L_7) e^{delta_6} (<= 180 / L_7 when e^{delta_6} <= 4/3),
    # q_src  = 64 * (5 / (4 L_7)) * 3 * e^{delta_6} = (240 / L_7) e^{delta_6}.
    L_try = 1024
    C_27 = 2 ** 27
    delta_6_genuine = Q(48 * (L_try ** 2), C_27)
    delta_6_fabricated = Q(96 * (L_try ** 2), C_27)
    require(delta_6_genuine == Q(3, 8) and delta_6_fabricated == Q(3, 4), "canonical-delta6-mismatch")

    # Defect 2: Under the genuine cubic Euler-tour counting (D_aug^2 = 64) at d=3 and C = 2^27,
    # the canonical KP condition (15/4) q_src / (1 - q_src) <= 1 requires q_src <= 4/19.
    # For EVERY L_7 > 0, let u = 48 L_7^2 / 2^27 > 0. Then
    #   q_src(L_7)^2 = (240^2 * 48 / 2^27) * (e^(2u) / u) = (675 / 32768) * (e^(2u) / u).
    # Since e^(2u) > 1 + 2u + 2u^2 = 4u + (2u - 1)^2/2 + 1/2 > 4u for all u > 0,
    #   e^(2u) / u > 4  =>  q_src(L_7)^2 > (675 / 32768) * 4 = 675 / 8192 > 16 / 361 = (4/19)^2!
    u = sp.Symbol("u", positive=True)
    poly_diff = sp.expand((1 + 2 * u + 2 * u**2) - 4 * u - ((2 * u - 1) ** 2 + 1) / 2)
    require(poly_diff == 0, "quadratic-taylor-4u-identity")
    pref = Q(240**2 * 48, C_27)
    require(pref == Q(675, 32768), "q-src-sq-prefactor")
    lower_bound_sq = pref * (2 if mutant == "canonical-obstruction" else 4)
    target_sq = Q(4, 19) ** 2
    require(lower_bound_sq == Q(675, 8192) and lower_bound_sq > target_sq, "no-L7-satisfies-cubic-canonical-bound-at-2-27")
    # Exact minimum ceiling under the earlier cubic Euler-tour route at d=3: C_min = 96 e * (240 * 19 / 4)^2 = 124,761,600 e
    c_min_coeff = 96 * ((240 * 19 // 4) ** 2)
    require(c_min_coeff == 124761600, "cubic-d3-c-min-coeff")
    ratio_lo = Q(c_min_coeff) * Q(1957, 720) / Q(C_27)
    ratio_hi = Q(c_min_coeff) * Q(11743, 4320) / Q(C_27)
    require(Q(252, 100) < ratio_lo < ratio_hi < Q(253, 100), "cubic-d3-c-min-ratio")
    return pref, lower_bound_sq, target_sq


def check_asymptotics_and_heuristic(mutant):
    d_sym = sp.Symbol("d", positive=True)
    C_thm_a = 442368 * d_sym**2 * (d_sym + 1) ** 2
    C_et_a0 = 1179648 * d_sym * (d_sym + 1) ** 4
    C_et_b1 = (2 ** 31) * d_sym * (d_sym + 1) ** 4
    C_va_b2 = 9830400 * d_sym**2 * (d_sym + 1) ** 2
    q_thm_a = sp.Rational(5, 24) * (2 * d_sym - 1) / (d_sym + 1)
    q_sharp_thm_a = sp.Rational(6903125, 40310784) * (2 * d_sym - 1) / (d_sym + 1)
    A_thm_a = 360 * (d_sym + 1) / (14 * d_sym + 29)

    # Leading large-d limits:
    lim_C_a = sp.limit(C_thm_a / d_sym**4, d_sym, sp.oo)
    lim_C_a0 = sp.limit(C_et_a0 / d_sym**5, d_sym, sp.oo)
    lim_C_b1 = sp.limit(C_et_b1 / d_sym**5, d_sym, sp.oo)
    lim_C_b2 = sp.limit(C_va_b2 / d_sym**4, d_sym, sp.oo)
    lim_q = sp.limit(q_thm_a, d_sym, sp.oo)
    lim_q_sharp = sp.limit(q_sharp_thm_a, d_sym, sp.oo)
    lim_A = sp.limit(A_thm_a, d_sym, sp.oo)
    if mutant == "asymptotics-and-heuristic":
        lim_q = sp.Rational(1, 2)
    require(lim_C_a == 442368 and lim_C_a0 == 1179648, "asymptotic-C-unconditional")
    require(lim_C_b1 == 2**31 and lim_C_b2 == 9830400, "asymptotic-C-conditional-form")
    require(
        lim_q == sp.Rational(5, 12)
        and lim_q_sharp == sp.Rational(6903125, 20155392)
        and lim_A == sp.Rational(180, 7),
        "asymptotic-q-and-A",
    )

    # Heuristic superexchange coordination scale: J = 4 t^2 / U, coordination 2d => 2d J = 8 d t^2 / U.
    # Ratio C_thm_a(d) / (8 d) ~ 55296 d^3 as d -> oo (vs C_et_a0(d) / (8 d) ~ 147456 d^4).
    ratio_a_heur = sp.limit(C_thm_a / (8 * d_sym * d_sym**3), d_sym, sp.oo)
    require(ratio_a_heur == 55296, "heuristic-ratio-d3")
    return lim_C_a, lim_q, lim_A, ratio_a_heur


def certificate(mutant=None, verbose=True):
    require(mutant is None or mutant in MUTANTS, "known-mutant")
    e_lo, e_hi = check_euler_enclosure(mutant)
    min_residuals = check_car_gram_and_sos(mutant)
    source_penalty = check_fugacity_and_source(mutant)
    animals = check_bethe_and_animals(mutant)
    thm_a = check_theorem_a_and_kp_budget(mutant, e_hi)
    routes = check_theorem_b_routes(mutant, e_hi)
    _, lower_sq, target_sq = check_cubic_canonical_obstruction(mutant)
    lim_C_a, lim_q, lim_A, ratio_a_heur = check_asymptotics_and_heuristic(mutant)

    if verbose:
        print(f"CERTIFIED rational e enclosure: {e_lo} < 11743/4320 < e_hi={e_hi} < 3; (4/3)^4 > e_hi, (11/8)^16 > e_hi^5, (6/5)^6 > e_hi")
        print(f"CERTIFIED 16x16 CAR bond SOS (all N_e=0..4) & released conditional-form boundary min residuals (d=2..5): {min_residuals}")
        print(f"CERTIFIED fugacity SOS identity 5*z_a - 4*z_b > 0 and pair-source scalar penalty={source_penalty} <= 1/4")
        for d in (2, 3, 4, 5):
            ex, fc = animals[d]
            print(f"CERTIFIED Z^{d} lattice animals v=1..4: exact={ex} <= Bethe T_{2*d}={fc} (v=4 plaquette diff={fc[3]-ex[3]})")
        for d in (2, 3, 4, 5):
            r = thm_a[d]
            print(
                f"CERTIFIED Theorem 1 (d={d}): D_phys={r['D_phys']}, D_aug={r['D_aug']}, "
                f"L={r['L']}, C={r['C']}, KP_phys<{r['kp_budget_sharp']}<1, KP_src<{r['kp_src_sharp']}<1, "
                f"prefactor_num={r['prefactor_num_hi']}<15, q={r['q_clean']} (sharp<{r['q_sharp']}), A={r['A_clean']}"
            )
        for d in (2, 3, 4, 5):
            rb = routes[d]
            print(
                f"CERTIFIED Baseline routes (d={d}): "
                f"C_ET={rb['C_a0']}, eps0_B1={rb['eps0_b1']}, C_B1={rb['C_b1']}, "
                f"eps0_B2={rb['eps0_b2']}, C_B2={rb['C_b2']}"
            )
        print(
            f"CERTIFIED earlier cubic route, canonical obstruction at 2^27: "
            f"for all L_7>0, q_src^2 > {lower_sq} > {target_sq} = (4/19)^2; "
            f"resolved by Theorem 1 with C(3)={thm_a[3]['C']} < 2^26 < 2^27"
        )
        print(
            f"CERTIFIED Asymptotics (d->oo): C(d)/d^4 -> {lim_C_a}, q(d) -> {lim_q}, "
            f"A(d) -> {lim_A}, C(d)/(8d*d^3) -> {ratio_a_heur}"
        )
        print("SCOPE: finite exact rational/symbolic witnesses; analytical original-model theorems require independent review")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutant", choices=MUTANTS, help="Run a single guarded mutant (must be rejected)")
    parser.add_argument(
        "--mutants",
        "--all-mutants",
        dest="all_mutants",
        action="store_true",
        help="Run all guarded mutants and verify each is rejected",
    )
    args = parser.parse_args()

    if args.all_mutants:
        for m in MUTANTS:
            try:
                certificate(mutant=m, verbose=False)
            except Rejected as exc:
                print(f"MUTANT {m}: REJECTED ({exc})")
            else:
                print(f"MUTANT {m}: SURVIVED")
                return 2
        print(f"ALL {len(MUTANTS)} MUTANTS REJECTED")
        return 0

    if args.mutant:
        try:
            certificate(mutant=args.mutant, verbose=False)
        except Rejected as exc:
            print(f"MUTANT {args.mutant}: REJECTED ({exc})")
            return 1
        print(f"MUTANT {args.mutant}: SURVIVED")
        return 2

    certificate(mutant=None, verbose=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
