"""
etab1BE1F_CaseS2 — Case S2 leptogenesis: full phase-space Boltzmann equations
with Delta L=1 scattering (s- and t-channel, Nl -> qt).

This module extends Case 4 (Fermi-Dirac/Bose-Einstein statistics, no kinetic
equilibrium assumption) by adding the Delta L=1 scattering collision integrals
C_{S,s} and C_{S,t} from arXiv:0907.0205 (Hahn-Woernle et al.).

Computational strategy
----------------------
The 9-dimensional scattering integrals reduce analytically to ten sets of
2-dimensional integrals (Appendix A of arXiv:0907.0205).  These are precomputed
over a dense (z, y_N) grid and stored in ulysses/data/AB_grid_CaseS2.npz.
At runtime the solver loads the grid and evaluates A(z,y_N) and B(z,y_N) by
fast bilinear interpolation rather than repeated quadrature.

The RHN distribution function is solved in the rescaled variable
    g(z, y_N) = f_N / f_N^eq,   where f_N^eq = 1/(exp(E_N/T)+1),
which improves numerical stability near equilibrium.  The per-momentum-mode
ODEs are then integrated with scipy.integrate.solve_ivp (Radau for strong
washout, RK45 for weak washout).  The lepton asymmetry ODE is subsequently
solved with the stored f_N grid as input.

Limitations: unflavoured (1F) only; no off-diagonal lepton-flavour
oscillations; gauge-boson processes and Delta L=1 CP-violation not included.
"""

import math
import os
import numpy as np
import ulysses
from scipy.integrate import solve_ivp
from scipy.sparse import diags as sp_diags
from numba import njit
from NumbaQuadpack import quadpack_sig, dqags
import numba as nb

# Top-quark Yukawa coupling squared (enters s- and t-channel prefactors)
_HT = 0.95



def load_AB(infile):
    """Load the precomputed AB grid.  Returns (A_grid, B_grid, z_grid, yn_vals)."""
    data = np.load(infile)
    return (data["A_grid"], data["B_grid"],
            data["z_grid"], data["yn_vals"])




def make_AB_eval(A_grid, B_grid, z_grid, yn_vals):
    """
    Build numba-compatible bilinear interpolators for A(z,yn) and B(z,yn).

    The interpolators are plain @njit functions that do a fast 2D lookup on
    the regular (z, log-yn) grid stored in A_grid / B_grid.

    Returns (eval_A_vec, eval_B_vec) — two @njit callables that accept a yn array.
    """
    # Pre-compute the grid spacing for fast index arithmetic
    dlogz  = np.log(z_grid[1]) - np.log(z_grid[0])  # uniform in log(z)
    logz0  = np.log(z_grid[0])
    n_z    = len(z_grid)

    log_yn  = np.log10(yn_vals)
    dlog_yn = log_yn[1] - log_yn[0]          # uniform in log10(yn)
    log_yn0 = log_yn[0]
    n_yn    = len(yn_vals)

    # Make contiguous C-order copies so numba can stride them
    A = np.ascontiguousarray(A_grid)
    B = np.ascontiguousarray(B_grid)

    @njit
    def _bilinear(grid, z, yn):
        """Bilinear interpolation on the (z, log-yn) grid."""
        # z index
        log_z = math.log(z)
        # Out-of-range: A/B are Boltzmann-suppressed beyond the precomputed grid
        if log_z > logz0 + (n_z - 1) * dlogz:
            return 0.0
        iz = int((log_z - logz0) / dlogz)
        if iz < 0:
            iz = 0
        if iz >= n_z - 1:
            iz = n_z - 2
        tz = (log_z - logz0 - iz * dlogz) / dlogz

        # yn index (log-uniform)
        log_y = math.log10(yn)
        iy = int((log_y - log_yn0) / dlog_yn)
        if iy < 0:
            iy = 0
        if iy >= n_yn - 1:
            iy = n_yn - 2
        ty = (log_y - log_yn0 - iy * dlog_yn) / dlog_yn

        # Bilinear interpolation
        v00 = grid[iz,     iy    ]
        v10 = grid[iz + 1, iy    ]
        v01 = grid[iz,     iy + 1]
        v11 = grid[iz + 1, iy + 1]

        return ((1 - tz) * (1 - ty) * v00 +
                      tz * (1 - ty) * v10 +
                (1 - tz) *       ty * v01 +
                      tz *       ty * v11)

    @njit
    def eval_A_vec(z, yn_array):
        """Vectorised A(z, yn) over an array of yn values at fixed z."""
        out = np.empty(len(yn_array))
        for i in range(len(yn_array)):
            out[i] = _bilinear(A, z, yn_array[i])
        return out

    @njit
    def eval_B_vec(z, yn_array):
        """Vectorised B(z, yn) over an array of yn values at fixed z."""
        out = np.empty(len(yn_array))
        for i in range(len(yn_array)):
            out[i] = _bilinear(B, z, yn_array[i])
        return out

    return eval_A_vec, eval_B_vec


def solve_fN_all_modes(yn_vals, K, zs, eval_A_vec, eval_B_vec, g0=0.0):
    """Solve the RHN phase-space ODE for all momentum modes in one call.

    All modes are independent, so they are stacked into a single ODE system
    of size n_yn and solved with one call to scipy.integrate.solve_ivp.
    The RHS and Jacobian are fully vectorised using numpy operations and
    the precomputed AB grid interpolators.

    Parameters
    ----------
    yn_vals : ndarray, shape (n_yn,)
        Dimensionless momentum grid y_N = |p_N| / T_ref.
    K : float
        Washout parameter k1 = Gamma_1 / H(M_1).
    zs : ndarray, shape (n_z,)
        Evolution variable z = M_1 / T.
    eval_A_vec, eval_B_vec : njit callables
        Vectorised interpolators returning A(z, yn_vals) and B(z, yn_vals)
        as arrays for a given scalar z.
    g0 : float
        Initial value of g = f_N / f_N^eq at zs[0].
        0.0 → vanishing IC;  1.0 → thermal IC.

    Returns
    -------
    fn_store : ndarray, shape (n_yn, n_z)
        f_N(z, y_N) for each momentum mode, ready for Nl_lstore.
    """
    n_yn = len(yn_vals)

    # Pre-mask modes that are fully Boltzmann-suppressed at z_min
    en_min   = np.sqrt(zs[0]**2 + yn_vals**2)
    active   = en_min < 50.0          # bool mask, shape (n_yn,)
    yn_act   = yn_vals[active]        # active modes only

    def rhs(z, g):
        # g has shape (n_yn_active,); each entry is g_i = f_N_i / f_Neq_i
        en      = np.sqrt(z**2 + yn_act**2)
        fNeq    = 1.0 / (1.0 + np.exp(en))
        dlnfNeq = -(1.0 - fNeq) * (z / en)
        A_arr   = eval_A_vec(z, yn_act)
        B_arr   = eval_B_vec(z, yn_act)
        return K * (A_arr * g + B_arr / fNeq) - g * dlnfNeq

    def jac(z, _):
        # Jacobian is diagonal (modes are independent); linear ODE so independent of g
        en      = np.sqrt(z**2 + yn_act**2)
        fNeq    = 1.0 / (1.0 + np.exp(en))
        dlnfNeq = -(1.0 - fNeq) * (z / en)
        A_arr   = eval_A_vec(z, yn_act)
        diag    = K * A_arr - dlnfNeq
        return sp_diags(diag, format='csc')

    y0 = np.full(len(yn_act), g0)

    if K > 3.0:
        sol = solve_ivp(rhs, [zs[0], zs[-1]], y0, method="Radau",
                        t_eval=zs, jac=jac,
                        rtol=1e-7, atol=1e-7)
    else:
        sol = solve_ivp(rhs, [zs[0], zs[-1]], y0, method="RK45",
                        t_eval=zs, rtol=1e-7, atol=1e-7)

    # Convert g → f_N = g * f_Neq for active modes
    fn_store = np.zeros((n_yn, len(zs)), dtype=np.float64)
    if sol.success and sol.y.shape[1] == len(zs):
        en_grid  = np.sqrt(zs[np.newaxis, :]**2 + yn_act[:, np.newaxis]**2)
        fNeq_grid = 1.0 / (1.0 + np.exp(en_grid))          # (n_act, n_z)
        fn_store[active] = sol.y * fNeq_grid
    else:
        # Fallback: fill with equilibrium distribution
        en_grid  = np.sqrt(zs[np.newaxis, :]**2 + yn_act[:, np.newaxis]**2)
        fn_store[active] = 1.0 / (1.0 + np.exp(en_grid))

    return fn_store   # shape (n_yn, n_z)



def Nl_lstore(fn_store,yn_vals,zs):

     # ── interpolator for fN(z, yn) ────────────────────────────────────────
    dlogy  = np.log10(yn_vals[1]) - np.log10(yn_vals[0])
    dlogy0 = np.log10(yn_vals[0])
    nz     = len(zs)
    ny     = len(yn_vals)
    dlogz  = np.log(zs[1]) - np.log(zs[0])
    logz0  = np.log(zs[0])  

    @njit
    def eval(z, yn):
        log_y  = math.log10(yn) if yn > 0 else dlogy0
        yindex = int(math.floor((log_y  - dlogy0) / dlogy))
        zindex = int(math.floor((math.log(z) - logz0) / dlogz))

        # z bounds — ODE should never go outside, but clamp safely
        if zindex < 0:
            zindex = 0
        if zindex >= nz - 1:
            zindex = nz - 2

        # yn upper bound — physically zero beyond stored grid
        if yindex >= ny - 1:
            return 0.0

        # yn lower bound — nearest edge, NOT zero
        if yindex < 0:
            yindex = 0

        # z interpolation — correct for log spacing
        tz = (z - zs[zindex]) / (zs[zindex+1] - zs[zindex])

        fn_lo = (fn_store[yindex    ][zindex  ] * (1 - tz)
            + fn_store[yindex    ][zindex+1] * tz)
        fn_hi = (fn_store[yindex + 1][zindex  ] * (1 - tz)
            + fn_store[yindex + 1][zindex+1] * tz)

        # yn interpolation — linear in yn (grid is log-spaced but fn varies smoothly)
        ty = (yn - yn_vals[yindex]) / (yn_vals[yindex+1] - yn_vals[yindex])
        return fn_lo * (1 - ty) + fn_hi * ty

    # ── safe Bose-Einstein for Higgs ──────────────────────────────────────
    @njit
    def bose_einstein(x):
        # fphi = 1/(exp(x)-1),  stable for x→0 and x→large
        if x > 40.0:
            return math.exp(-x)        # Maxwell-Boltzmann limit
        elif x < 1e-10:
            return 1.0 / 1e-10         # IR cutoff, avoids 1/0
        else:
            return 1.0 / (math.exp(x) - 1.0)




    #Setting up the s channel integral for lepton asymmetry
    #Et integrand by multiplyinng A.39 and A.54
    @nb.cfunc(quadpack_sig)
    def l_Cs1etintegrand(et,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = yl
        yn = math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        #Phase space factor from A.54; rewritten to avoid exp(et) overflow: exp(et)/(1+exp(et)) = 1/(1+exp(-et))
        lmbs = -fl*(math.exp(-en-el)+(1-math.exp(-en-el))*fn)/((1.0+math.exp(-et))*(1+math.exp(et-el-en)))
        alph = 4*yn*(el+en)
        beta = ((en-yn)*(2*el+en-yn))/((en+yn)*(2*el+en+yn))
        if beta <= 0.0:
            return 0.0
        #q integral result from A.39
        Is1  = (alph+z*z*math.log(beta))/(2*(el+en))
        return lmbs*Is1
    l_Cs1etintegrand_ptr = l_Cs1etintegrand.address

    #Solving Et integral from A.56
    @nb.cfunc(quadpack_sig)
    def l_ets1integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = yl
        yn=math.sqrt(en*en-z*z)
        etlowerlim = 1e-10
        ethigherlim = el+en
        integral,abserr,success_et = dqags(l_Cs1etintegrand_ptr,etlowerlim,ethigherlim,data1)
        return integral
    l_ets1integral_ptr = l_ets1integral.address

    #Solving En integral from A.56
    @njit
    def l_ens1integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        el   = yl
        enlowerlim = z
        enhigherlim = math.sqrt(el*el+z*z)
        integral,abserr,success = dqags(l_ets1integral_ptr,enlowerlim,enhigherlim,data)  
        return integral


    #Et integrand by multiplyinng A.41 and A.54
    @nb.cfunc(quadpack_sig)
    def l_Cs2etintegrand(et,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = yl
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        #Phase space factor from A.54; rewritten to avoid exp(et) overflow
        lmbs = -fl*(math.exp(-en-el)+(1-math.exp(-en-el))*fn)/((1.0+math.exp(-et))*(1+math.exp(et-el-en)))
        alph = 4*el*(el+en)
        beta = (en*en-yn*yn)/((2*el+en)*(2*el+en)-yn*yn)
        if beta <= 0.0:
            return 0.0
        #q integral result from A.41
        Is2  = (alph+z*z*math.log(beta))/(2*(el+en))
        return lmbs*Is2
    l_Cs2etintegrand_ptr = l_Cs2etintegrand.address

    #Solving Et integral from A.57
    @nb.cfunc(quadpack_sig)
    def l_ets2integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = yl
        yn=math.sqrt(en*en-z*z)
        etlowerlim = 1e-10 
        ethigherlim = el+en
        integral,abserr,success_et = dqags(l_Cs2etintegrand_ptr,etlowerlim,ethigherlim,data1)  
        return integral
    l_ets2integral_ptr = l_ets2integral.address

    #Solving En integral from A.57
    @njit
    def l_ens2integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        el   = yl
        enlowerlim = math.sqrt(el*el+z*z)
        enhigherlim = 300
        if enlowerlim >= enhigherlim:
            return 0.0
        integral,abserr,success = dqags(l_ets2integral_ptr,enlowerlim,enhigherlim,data)
        return integral


    #Et integrand by multiplyinng A.43 and A.54
    @nb.cfunc(quadpack_sig)
    def l_Cs3etintegrand(et,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = yl
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        #Phase space factor from A.54; rewritten to avoid exp(et) overflow
        lmbs = -fl*(math.exp(-en-el)+(1-math.exp(-en-el))*fn)/((1.0+math.exp(-et))*(1+math.exp(et-el-en)))
        alph = 2*(el+en)*(en-2*et+yn)
        beta = (et*(2*el+en-yn))/((el+en-et)*(en+yn))
        if beta <= 0.0:
            return 0.0
        #q integral result from A.43
        Is3  = -(alph+z*z*math.log(beta))/(2*(el+en))
        return lmbs*Is3
    l_Cs3etintegrand_ptr = l_Cs3etintegrand.address

    #Solving Et integral from A.58
    @nb.cfunc(quadpack_sig)
    def l_ets3integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = yl
        yn=math.sqrt(en*en-z*z)
        etlowerlim = 1e-10 
        ethigherlim = 0.5*(en+yn)
        integral,abserr,success_et = dqags(l_Cs3etintegrand_ptr,etlowerlim,ethigherlim,data1)       
        return integral
    l_ets3integral_ptr = l_ets3integral.address

    #Solving En integral from A.58
    @njit
    def l_ens3integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        el   = yl
        enlowerlim = z
        enhigherlim = math.sqrt(el*el+z*z)
        integral,abserr,success = dqags(l_ets3integral_ptr,enlowerlim,enhigherlim,data)       
        return integral


    #Et integrand by multiplyinng A.45 and A.54
    @nb.cfunc(quadpack_sig)
    def l_Cs4etintegrand(et,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = yl
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        #Phase space factor from A.54; rewritten to avoid exp(et) overflow
        lmbs = -fl*(math.exp(-en-el)+(1-math.exp(-en-el))*fn)/((1.0+math.exp(-et))*(1+math.exp(et-el-en)))
        alph = 2*(el+en)*(2*el+en-2*et-yn)
        beta = ((el+en-et)*(2*el+en-yn))/(et*(en+yn))
        if beta <= 0.0:
            return 0.0
        #q integral result from A.45
        Is4  = (alph-z*z*math.log(beta))/(2*(el+en))
        return lmbs*Is4
    l_Cs4etintegrand_ptr = l_Cs4etintegrand.address

    #Solving Et integral from A.59
    @nb.cfunc(quadpack_sig)
    def l_ets4integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = yl
        yn=math.sqrt(en*en-z*z)
        etlowerlim = 0.5*(2*el+en-yn)
        ethigherlim = el+en
        integral,abserr,success_et = dqags(l_Cs4etintegrand_ptr,etlowerlim,ethigherlim,data1)        
        return integral
    l_ets4integral_ptr = l_ets4integral.address

    #Solving En integral from A.59
    @njit
    def l_ens4integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        el   = yl
        enlowerlim = z
        enhigherlim = math.sqrt(el*el+z*z)
        integral,abserr,success = dqags(l_ets4integral_ptr,enlowerlim,enhigherlim,data) 
        return integral


    #Et integrand by multiplyinng A.47 and A.54
    @nb.cfunc(quadpack_sig)
    def l_Cs5etintegrand(et,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = yl
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        #Phase space factor from A.54; rewritten to avoid exp(et) overflow
        lmbs = -fl*(math.exp(-en-el)+(1-math.exp(-en-el))*fn)/((1.0+math.exp(-et))*(1+math.exp(et-el-en)))
        alph = 2*(el+en)*(2*el+en-2*et-yn)
        beta = ((el+en-et)*(2*el+en-yn))/(et*(en+yn))
        if beta <= 0.0:
            return 0.0
        #q integral result from A.47
        Is5  = -(alph-z*z*math.log(beta))/(2*(el+en))
        return lmbs*Is5
    l_Cs5etintegrand_ptr = l_Cs5etintegrand.address

    #Solving Et integral from A.50
    @nb.cfunc(quadpack_sig)
    def l_ets5integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = yl
        yn=math.sqrt(en*en-z*z)
        etlowerlim = 1e-10
        ethigherlim = 0.5*(2*el+en-yn)
        integral,abserr,success_et = dqags(l_Cs5etintegrand_ptr,etlowerlim,ethigherlim,data1)
        return integral
    l_ets5integral_ptr = l_ets5integral.address

    #Solving En integral from A.60
    @njit
    def l_ens5integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        el   = yl
        enlowerlim = math.sqrt(el*el+z*z)
        enhigherlim = 300
        if enlowerlim >= enhigherlim:
            return 0.0
        integral,abserr,success = dqags(l_ets5integral_ptr,enlowerlim,enhigherlim,data)
        return integral


    #Et integrand by multiplyinng A.49 and A.54
    @nb.cfunc(quadpack_sig)
    def l_Cs6etintegrand(et,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = yl
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        #Phase space factor from A.54; rewritten to avoid exp(et) overflow
        lmbs = -fl*(math.exp(-en-el)+(1-math.exp(-en-el))*fn)/((1.0+math.exp(-et))*(1+math.exp(et-el-en)))
        alph = 2*(el+en)*(en-2*et+yn)
        beta = ((el+en-et)*(en+yn))/(et*(2*el+en-yn))
        if beta <= 0.0:
            return 0.0
        #q integral result from A.49
        Is6  = (alph-z*z*math.log(beta))/(2*(el+en))
        return lmbs*Is6
    l_Cs6etintegrand_ptr = l_Cs6etintegrand.address

    #Solving Et integral from A.61
    @nb.cfunc(quadpack_sig)
    def l_ets6integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = yl
        yn=math.sqrt(en*en-z*z)
        etlowerlim = 0.5*(en+yn)
        ethigherlim = el+en
        integral,abserr,success_et = dqags(l_Cs6etintegrand_ptr,etlowerlim,ethigherlim,data1)
        return integral
    l_ets6integral_ptr = l_ets6integral.address

    #Solving En integral from A.61
    @njit
    def l_ens6integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        el   = yl
        enlowerlim = math.sqrt(el*el+z*z)
        enhigherlim = 300
        if enlowerlim >= enhigherlim:
            return 0.0
        integral,abserr,success = dqags(l_ets6integral_ptr,enlowerlim,enhigherlim,data)
        return integral


    #Setting up the t channel integral for lepton asymmetry
    #Et integrand by multiplyinng It1 (A.82) and A.80
    @nb.cfunc(quadpack_sig)
    def l_Ct1eqintegrand(eq,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = math.sqrt(yl*yl)
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        ah = 1e-5
        #Phase space factor from A.80; guard against exp(en-el) overflow
        _enel = en - el
        if _enel > 709.0:
            return 0.0
        lmbt = fl*(((-1+fn)-math.exp(_enel)*fn))/((1+math.exp(eq))*(1+math.exp(_enel-eq)))
        alph = 2*_enel*(2*(en-eq-el)+ah*z)
        beta = ((eq)*(ah*z))/((eq-en+el)*(2*_enel+ah*z))
        #q integral result from It1 (A.82)
        if en>el:
            if beta <= 0.0:
                return 0.0
            It1  = -(alph+z*z*math.log(beta))/(2*_enel)
        else:
            return 0
        return lmbt*It1
    l_Ct1eqintegrand_ptr = l_Ct1eqintegrand.address


    #Solving Eq integral from A.82
    @nb.cfunc(quadpack_sig)
    def l_eqt1integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = math.sqrt(yl*yl)
        ah = 1e-5
        yn=math.sqrt(en*en-z*z)
        eqlowerlim = (en-el+0.5*ah*z)
        eqhigherlim = 0.5*(en+yn)
        # Inverted limits: kinematic region where beta<0 throughout
        if eqlowerlim >= eqhigherlim:
            return 0.0
        integral,abserr,success_eq = dqags(l_Ct1eqintegrand_ptr,eqlowerlim,eqhigherlim,data1)
        return integral
    l_eqt1integral_ptr = l_eqt1integral.address


    #Solving En integral from A.82
    @njit
    def l_ent1integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        ah = 1e-5
        el   = math.sqrt(yl*yl)
        enlowerlim = (((2*el-ah*z)*(2*el-ah*z)+z*z)/(2*(2*el-ah*z)))
        enhigherlim = 300
        if enlowerlim<enhigherlim:
            integral,abserr,success = dqags(l_eqt1integral_ptr,enlowerlim,enhigherlim,data)
            return integral
        else:
            return 0


    #Et integrand by multiplyinng It2 (A.83) and A.80
    @nb.cfunc(quadpack_sig)
    def l_Ct2eqintegrand(eq,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = math.sqrt(yl*yl)
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        ah = 1e-5
        #Phase space factor from A.80; guard against exp(en-el) overflow
        _enel = en - el
        if _enel > 709.0:
            return 0.0
        lmbt = fl*(((-1+fn)-math.exp(_enel)*fn))/((1+math.exp(eq))*(1+math.exp(_enel-eq)))
        alph = 2*_enel*(2*el-en+yn-ah*z)
        beta = ((en+yn)*(ah*z))/((-en+2*el+yn)*(2*_enel+ah*z))
        #q integral result from It2 (A.83)
        if en>el:
            if beta <= 0.0:
                return 0.0
            It2  = (alph-z*z*math.log(beta))/(2*_enel)
        else:
            return 0
        return lmbt*It2
    l_Ct2eqintegrand_ptr = l_Ct2eqintegrand.address

    #Solving Eq integral from A.83
    @nb.cfunc(quadpack_sig)
    def l_eqt2integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = math.sqrt(yl*yl)
        ah = 1e-5
        yn=math.sqrt(en*en-z*z)
        eqlowerlim = 0.5*(en+yn)
        eqhigherlim = 300
        integral,abserr,success_eq = dqags(l_Ct2eqintegrand_ptr,eqlowerlim,eqhigherlim,data1)
        return integral
    l_eqt2integral_ptr = l_eqt2integral.address

    #Solving En integral from A.83
    @njit
    def l_ent2integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        ah = 1e-5
        el   = math.sqrt(yl*yl)
        enlowerlim = ((2*el-ah*z)*(2*el-ah*z)+z*z)/(2*(2*el-ah*z))
        enhigherlim = 300
        if enlowerlim<enhigherlim:
            integral,abserr,success = dqags(l_eqt2integral_ptr,enlowerlim,enhigherlim,data)
            return integral
        else:
            return 0

    #Et integrand by multiplyinng It3 (A.84) and A.80
    @nb.cfunc(quadpack_sig)
    def l_Ct3eqintegrand(eq,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = math.sqrt(yl*yl)
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        ah = 1e-5
        #Phase space factor from A.80; guard against exp(en-el) overflow
        _enel = en - el
        if _enel > 709.0:
            return 0.0
        lmbt = fl*(((-1+fn)-math.exp(_enel)*fn))/((1+math.exp(eq))*(1+math.exp(_enel-eq)))
        alph = 2*_enel*(2*eq-en+yn)
        beta = ((en-yn)*(en-eq-el))/((en-2*el+yn)*(eq))
        if beta <= 0.0:
            return 0.0
        #q integral result from It3 (A.84)
        It3  = (alph+z*z*math.log(beta))/(2*_enel)
        return lmbt*It3
    l_Ct3eqintegrand_ptr = l_Ct3eqintegrand.address

    #Solving Eq integral from A.84
    @nb.cfunc(quadpack_sig)
    def l_eqt3integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = math.sqrt(yl*yl)
        ah = 1e-5
        yn=math.sqrt(en*en-z*z)
        eqlowerlim = 0.5*(en-yn)
        eqhigherlim = 0.5*(en+yn)
        integral,abserr,success_eq = dqags(l_Ct3eqintegrand_ptr,eqlowerlim,eqhigherlim,data1)
        return integral
    l_eqt3integral_ptr = l_eqt3integral.address

    #Solving En integral from A.84
    @njit
    def l_ent3integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        ah = 1e-5
        el   = math.sqrt(yl*yl)
        enlowerlim = z
        enhigherlim = (4*el*el+z*z)/(4*el)
        integral,abserr,success = dqags(l_eqt3integral_ptr,enlowerlim,enhigherlim,data)
        return integral

    #Et integrand by multiplyinng It4 (A.85) and A.80
    @nb.cfunc(quadpack_sig)
    def l_Ct4eqintegrand(eq,data_):
        data=nb.carray(data_, (4,))
        en,yl,z,Nl = data
        el   = math.sqrt(yl*yl)
        yn=math.sqrt(en*en-z*z)
        fn = eval(z,yn)
        fl= Nl/(1+math.exp(el))*4/3
        ah = 1e-5
        #Phase space factor from A.80; guard against exp(en-el) overflow
        _enel = en - el
        if _enel > 709.0:
            return 0.0
        lmbt = fl*(((-1+fn)-math.exp(_enel)*fn))/((1+math.exp(eq))*(1+math.exp(_enel-eq)))
        alph = 4*_enel*yn
        beta = ((en-yn)*(en-yn-2*el))/((en+yn)*(en+yn-2*el))
        if beta <= 0.0:
            return 0.0
        #q integral result from It4 (A.85)
        It4  = (alph+z*z*math.log(beta))/(2*_enel)
        return lmbt*It4
    l_Ct4eqintegrand_ptr = l_Ct4eqintegrand.address

    #Solving Eq integral from A.85
    @nb.cfunc(quadpack_sig)
    def l_eqt4integral(en,data_):
        data = nb.carray(data_,(3,))
        yl,z,Nl = data
        data1 = np.array([en,yl,z,Nl],np.float64)
        el   = math.sqrt(yl*yl)
        ah = 1e-5
        yn=math.sqrt(en*en-z*z)
        eqlowerlim = 0.5*(en+yn)
        eqhigherlim = 300
        integral,abserr,success_eq = dqags(l_Ct4eqintegrand_ptr,eqlowerlim,eqhigherlim,data1)
        return integral
    l_eqt4integral_ptr = l_eqt4integral.address

    #Solving En integral from A.85
    @njit
    def l_ent4integral(yl,z,Nl):
        data=np.array([yl,z,Nl],np.float64)
        ah = 1e-5
        el   = math.sqrt(yl*yl)
        enlowerlim = z
        enhigherlim = (4*el*el+z*z)/(4*el)
        integral,abserr,success = dqags(l_eqt4integral_ptr,enlowerlim,enhigherlim,data)
        return integral


    #yn integrand for lepton asymmetry evolution as in case D4
    @nb.cfunc(quadpack_sig)
    def D4Lintegrand(yn, data_):
        data = nb.carray(data_,(4,))
        yl, z, Nl, eps = data
        en   = math.sqrt(z * z + yn * yn)
        fn = eval(z,yn)
        fleq = math.exp(-yl)/(1. + math.exp(-yl))
        fphi = bose_einstein(en - yl) 
        p1   =  (yn / en) * ((fphi + fn) * ((4/3) * Nl + 2 * eps) * fleq)
        p2   = (yn/en) * (-2 * eps * fn * (1 + fphi))
        return p1 + p2
    D4Lintegrand_ptr = D4Lintegrand.address

    #Solving the yn integral for lepton aymmetry evolution
    @njit
    def ynintegral(yl, z, Nl, eps,K, highlim = 100, epsrel = 1e-10, epsabs = 1e-10):
        data=np.array([yl,z,Nl,eps],np.float64)
        nlowerlim = np.abs((-z * z + 4. * yl * yl) / (4. * yl))
        integral,abserr,success = dqags(D4Lintegrand_ptr, nlowerlim, highlim, data)
        #4.80822761 comes from 4*ζ(3)
        return -z*z*K/4.80822761*integral

    '''Adding up all the parts of the l_Csintegral and multiplying them with 
    the pre factor as well as a factor of 2 (processes involving antiparticle)'''
    @njit
    def l_Csintegrand(yl,z,Nl,K):
        integsum = l_ens1integral(yl,z,Nl)+l_ens2integral(yl,z,Nl)+l_ens3integral(yl,z,Nl)+\
            l_ens4integral(yl,z,Nl)+l_ens5integral(yl,z,Nl)+l_ens6integral(yl,z,Nl)
        # 0.03160869 = 3/(8*pi^2*ζ(3)); _HT^2 = top Yukawa squared
        return 0.03160869*_HT*_HT*K*integsum

    '''Adding up all the parts of the l_Ctintegral and multiplying them with
    the pre factor as well as a factor of 4 (processes involving antiparticle and u channel)'''
    @njit
    def l_Ctintegrand(yl,z,Nl,K):
        integsum = l_ent1integral(yl,z,Nl)+l_ent2integral(yl,z,Nl)+l_ent3integral(yl,z,Nl)+l_ent4integral(yl,z,Nl)
        # 0.06321738 = 3/(4*pi^2*ζ(3)); _HT^2 = top Yukawa squared
        return 0.06321738*_HT*_HT*K*integsum

    #yl iintegrand for lepton asymmetry evolution
    @nb.cfunc(quadpack_sig)
    def l_ylintegral(yl,data_):
        data = nb.carray(data_,(4,))
        z,Nl,eps,K = data
        integsum =  ynintegral(yl,z,Nl,eps,K)+l_Csintegrand(yl,z,Nl,K)+l_Ctintegrand(yl,z,Nl,K)
        return integsum
    l_ylintegral_ptr = l_ylintegral.address

    #Solving the yl integral for lepton asymmetry evolution
    @njit
    def Nlrhs(z,Nl,eps,K):
        data = np.array([z,Nl[0],eps,K],np.float64)
        llowerlim = 1e-9
        lhigherlim = 100
        integral,abserr,success = dqags(l_ylintegral_ptr,llowerlim,lhigherlim,data)
        return np.array([integral],np.float64)
    
    return Nlrhs

        
class EtaB_1BE1F_CaseS2(ulysses.ULSBase):
    """
    Boltzmann equation (BE) with one decaying sterile (Case S2): drops kinetic equilibrium
    assumption and uses full quantum statistics. Includes Delta L=1 scatterings (s- and
    t-channel) via precomputed AB grids. See arXiv:0907.0205 Eqns. 3.32 and 3.34.
    Flavour off-diagonal oscillations are not included.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def shortname(self): return "1BE1F_CaseS2"

    def flavourindices(self): return [1]

    def flavourlabels(self): return ["$NBL$"]

    @property
    def EtaB(self):
        #Define fixed quantities for BEs
        nevalsy              =  500
        nevalsz              =  200
        yn_vals              =  np.logspace(-2., np.log10(105.), nevalsy).astype(np.float64)
        epstt                =  np.real(self.epsilon1ab(2,2))
        epsmm                =  np.real(self.epsilon1ab(1,1))
        epsee                =  np.real(self.epsilon1ab(0,0))
        self.eps             =  epsee + epsmm + epstt
        self._K              =  np.real(self.k1)
        K                    =  self._K
        eps                  =  self.eps
        self.zs              =  np.logspace(np.log10(self.zmin), np.log10(self.zmax), nevalsz).astype(np.float64)
        zs                   =  self.zs

        data_dir = os.path.dirname(ulysses.__file__)
        PATH_AB_integral = os.path.join(data_dir, './data/AB_grid_CaseS2.npz')
        AB_integral = load_AB(PATH_AB_integral)
        
        eval_A_vec, eval_B_vec = make_AB_eval(AB_integral[0],AB_integral[1],AB_integral[2],AB_integral[3])

        fn_store = solve_fN_all_modes(yn_vals, K=K, zs=zs,
                               eval_A_vec=eval_A_vec, eval_B_vec=eval_B_vec,
                               g0=float(self.initial_abundance))

        print("RHN phase space distribution done")

        Nlrhs = (Nl_lstore(fn_store,yn_vals,zs))

        #Solving differential equation to get N{l-l}(z) for Case S2
        if K > 3.0:

            solLS2 = solve_ivp(
                Nlrhs, [zs[0], zs[-1]], y0=[0.0],
                args=(eps, K),
                method="Radau",
                t_eval=zs,
                rtol=1e-7, atol=1e-7,
            )
        else:
            solLS2 = solve_ivp(
                Nlrhs, [zs[0], zs[-1]], y0=[0.0],
                args=(eps, K),
                method="RK45",
                max_step=1.0/10,
                t_eval=zs,
                rtol=1e-7, atol=1e-7,
            )

        Nl_final = solLS2.y[0]

        ys = np.transpose((Nl_final, Nl_final))
        self.setEvolData(ys)
        print(f"K={K}, EtaB={self.ys[-1][-1]:.4e}")
        return self.ys[-1][-1]