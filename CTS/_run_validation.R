###############################################################################
# _run_validation.R
# Elagolix Calcium Homeostasis & Bone Remodeling QSP Model
# -----------------------------------------------------------------------
# Stodtmann et al. (2021) Clin Transl Sci. 14:1611-1619
# Peterson & Riggs (2010) + Riggs et al. (2012) estrogen extension
# Ref implementation: MetrumRG OpenBoneMin (faithful R port)
###############################################################################

# === Packages ================================================================
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# === Constants ===============================================================
HOURS_PER_MONTH <- 730.5

###############################################################################
#
#  1.  Module A --- Dose -> E2  (Scaled Logistic, Eq. 3)
#
###############################################################################

dose_e2_params <- list(slope = 0.00894, logE2max = 5.20, logE2min = 2.14)

predict_E2 <- function(daily_dose, p = dose_e2_params) {
  E2min <- exp(p$logE2min)   # 8.50
  E2max <- exp(p$logE2max)   # 181.3
  E2min + (E2max - E2min) / (1 + exp(p$slope * daily_dose))
}

E2_baseline <- predict_E2(0)   # ~94.9 pg/mL

###############################################################################
#
#  2.  QSP Parameters (faithful to OpenBoneMin.cpp)
#
###############################################################################

make_params <- function() {
  list(
    # -- Bone cell dynamics --
    OBtot0      = 0.00501324,
    k1          = 6.24e-6,
    k2          = 0.112013,
    k3          = 6.24e-6,
    k4          = 0.112013,
    kO          = 15.8885,
    kb          = 0.000605516,
    Pic0        = 0.228142,
    FracOBfast  = 0.797629,
    Frackb      = 0.313186,
    Da          = 0.7 / 24,

    # -- TGF-beta --
    OBtgfGAM    = 0.0111319,
    koutTGF0    = 2.98449e-5,
    koutTGFGam  = 0.919131,
    OCtgfGAM    = 0.593891,

    # -- ROB (responding osteoblasts) differentiation --
    EmaxPicROB   = 3.9745,
    PicROBgam    = 1.80968,
    FracPicROB   = 0.883824,

    # -- OB differentiation (PicOB) --
    PicOBgam     = 0.122313,
    FracPicOB    = 0.000244818,
    EmaxPicOB    = 0.251636,

    # -- OB apoptosis (TGF-beta -> kbprime) --
    PicOBgamkb      = 2.92375,
    MultPicOBkb     = 3.11842,
    FracPic0kb      = 0.764028,
    E0RUNX2kbEffFACT = 1.01,
    RUNkbGAM        = 3.81644,
    RUNkbMaxFact    = 0.638114,

    # -- OC (osteoclast) --
    E0Meff       = 0.388267,
    EmaxMeffOC   = 3.15667,
    kinOCgam     = 8.53065,
    EmaxPicOC    = 1.9746,
    FracPicOC    = 0.878215,
    PicOCgam     = 1.0168,
    MOCratioGam  = 0.603754,

    # -- OC survival --
    E0RANKL      = 3.80338,
    EmaxL        = 0.469779,
    LsurvOCgam   = 3.09023,
    LsurvOCCgam  = 3.0923,

    # -- Ca/Phos homeostasis --
    CaDay        = 88.0,
    V1           = 14.0,
    FracJ14      = 0.107763,
    J14OCmax     = 0.543488,
    J14OCgam     = 1.6971,
    FracJ15      = 0.114376,
    k14a         = 2.44437e-5,
    HApMRT       = 3.60609,

    # -- Oral Ca/Phos --
    OralCa       = 24.055 / 24,
    OralPhos     = 10.5 / 24,
    F12          = 0.7,

    # -- Gut Ca absorption --
    T28          = 0.9,
    T310         = 0.105929,
    T81          = 0.75,
    T87          = 0.0495,
    T77          = 0.909359,
    T80          = 4.0,
    T43          = 1.03856,
    T45          = 3.85,
    T84          = 0.25,

    # -- Renal Ca --
    Reabs50      = 1.57322,
    T16          = 1.06147,
    maxTmESTkid  = 0.923737,
    T7           = 2.0,
    T9           = 90.0,
    T70          = 0.01,
    T71          = 0.03,

    # -- Calcitriol / 1,25(OH)2D --
    T33          = 0.003,
    T34          = 0.037,
    T35          = 90.0,
    CaPOgam      = 1.0,
    AlphOHgam    = 0.111241,
    T64          = 0.05,
    T65          = 6.3,
    T67          = 1.54865,
    T69          = 0.10,
    CtriolPTgam  = 12.5033,
    CtriolMax    = 4.1029,
    CtriolMin    = 0.9,

    # -- Phosphate effects --
    ScaEffGam    = 0.9,
    PhosEff0     = 1.52493,
    PhosEff50    = 1.3021,
    PhosEffGam   = 8.25229,
    PO4inhPTHgam = 0.0,

    # -- Phosphate fluxes --
    T46          = 1.142,
    T49          = 51.8,
    T52          = 0.365,
    T55          = 0.019268,

    # -- PTH --
    kout         = 100 / 14,
    PTout        = 1.604e-4,
    T58          = 6249.09,
    T59          = 11.7387,
    T61          = 96.25,

    # -- RANKL / OPG --
    koutL        = 0.00293273,
    kinRNKgam    = 0.151825,
    koutRNK      = 0.00323667,
    opgPTH50     = 3.85,
    EmaxLpth     = 1.30721,
    OsteoEffectGam = 0.173833,

    # -- Intracellular signaling (RUNX2, CREB, BCL2) --
    RUNX20       = 10.0,
    RX2Kout0     = 0.693,
    E0rx2Kout    = 0.125,
    EmaxPTHRX2x  = 5.0,
    E0crebKin    = 0.5,
    EmaxPTHcreb  = 3.39745,
    crebKout     = 0.00279513,
    bcl2Kout     = 0.693,

    # -- BMD --
    koutBMDls    = 0.000397,
    gamOB        = 0.0793,
    gamOCls      = 0.14,

    # -- Estrogen extension (Riggs 2012) --
    ESTON        = 0.0,      # OFF for Elagolix (direct EST forcing)
    koutEST      = 0.05776227,
    tgfbGAM      = 0.0374,
    tgfbactGAM   = 0.045273,
    robGAM       = 0.16,
    obGAM        = 0.000012,
    E2scalePicB1 = 1.16832e-5,

    # -- GFR --
    GFR0         = 100 / 16.667,

    # -- Denosumab (unused for elagolix, set to zero) --
    kdenosl      = 2.0e-06
  )
}

###############################################################################
#
#  3.  Initial Conditions (31 state variables)
#      Matches OpenBoneMin.cpp $INIT + [MAIN] block
#
###############################################################################

make_init <- function(p) {
  L0   <- 0.4
  RNK0 <- 10.0
  O0   <- 4.0
  M0   <- p$k3 * RNK0 * L0 / p$k4
  N0   <- p$k1 * O0   * L0 / p$k2

  c(
    PTH     = 53.90,
    S       = 0.5,
    PTmax   = 1.0,
    B       = 1260.0,
    AOH     = 126.0,      # B_0 / 10
    P       = 32.90,
    ECCPhos = 16.8,
    Tgut    = 1.58471,
    Rkid    = 0.50,
    HAp     = 1.0,
    PhosGut = 0.839,
    IntraPO = 3226.0,
    Q       = 100.0,
    Qbone   = 24900.0,
    UCA     = 0.0,
    ROB1    = 0.00104122,
    OBfast  = p$OBtot0 * p$FracOBfast,
    OBslow  = p$OBtot0 * (1 - p$FracOBfast),
    OC      = 0.00115398,
    L       = L0,
    RNK     = RNK0,
    Oopg    = O0,
    Mcmplx  = M0,
    Ncmplx  = N0,
    TGFB    = p$Pic0 * 1000,
    TGFBact = p$Pic0,
    RX2     = 10.0,
    CREB    = 10.0,
    BCL2    = 100.0,
    EST     = 1.0,
    BMDls   = 1.0
  )
}

###############################################################################
#
#  4.  ODE System  (faithful port of OpenBoneMin.cpp [ODE] block)
#
###############################################################################

qsp_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {

    # --- Baseline reference values ---
    OB0   <- OBtot0
    OC0   <- 0.00115398
    ROB0  <- 0.00104122
    L0    <- 0.4
    RNK0  <- 10.0
    O0    <- 4.0
    PTH0  <- 53.90
    B0    <- 1260.0
    P0    <- 32.90
    M0    <- k3 * RNK0 * L0 / k4
    OBfast0 <- OBtot0 * FracOBfast
    OBslow0 <- OBtot0 * (1 - FracOBfast)
    TGFBact0 <- Pic0
    TGFB0    <- Pic0 * 1000
    RX20v    <- RUNX20
    Q0       <- 100.0
    Qbone0   <- 24900.0
    S0       <- 0.5
    T0       <- 1.58471   # Tgut initial
    R0       <- 0.5       # Rkid initial
    Calcitriol0 <- B0 / V1

    OB    <- OBfast + OBslow
    EST_s <- max(EST, 1e-6)

    # Concentrations
    CaConc0  <- P0 / V1
    PTHconc0 <- PTH0 / V1
    PTHconc  <- PTH / V1
    CaConc   <- P / V1
    C1       <- P / V1
    C2       <- ECCPhos / V1   # extracellular phosphate concentration
    C8       <- B / V1         # calcitriol concentration

    #=================================================================
    # Derived parameters from SS initial conditions
    #=================================================================
    T13 <- (CaDay / 24) / Q0

    T15 <- CaDay / (CaConc0 * V1 * 24)

    J14OC50 <- exp(log((J14OCmax * OC0^J14OCgam / T13) -
                        OC0^J14OCgam) / J14OCgam)

    OCeqn <- (J14OCmax * OC^J14OCgam) /
              (OC^J14OCgam + J14OC50^J14OCgam)

    kinRNKval <- (koutRNK * RNK0 + k3 * RNK0 * L0 - k4 * M0) /
                  TGFBact0^kinRNKgam

    MOCratio  <- Mcmplx / max(OC, 1e-15)
    MOCratio0 <- M0 / OC0
    MOCratioEff <- (MOCratio / MOCratio0)^MOCratioGam

    # J14: calcium flux from bone to plasma
    J14OCdepend <- OCeqn * Q0 * FracJ14 * MOCratioEff
    J14val <- T13 * Q0 * (1 - FracJ14) + J14OCdepend

    # J41: phosphate from bone (stoichiometric 0.464)
    J41 <- 0.464 * J14val

    bigDb <- kb * OB0 * Pic0 / ROB0

    #=================================================================
    # TGF-beta
    #=================================================================
    kinTGF   <- koutTGF0 * TGFB0
    koutTGF  <- koutTGF0 * (TGFB / TGFB0)^koutTGFGam
    koutTGFact_val <- koutTGF0 * 1000
    koutTGFeqn <- koutTGF * TGFB * (OC / OC0)^OCtgfGAM

    #=================================================================
    # ROB differentiation (PicROB)
    #=================================================================
    E0PicROB <- FracPicROB * Pic0
    EC50PicROBparen <- (EmaxPicROB * TGFBact0^PicROBgam /
                        (Pic0 - E0PicROB)) - TGFBact0^PicROBgam
    EC50PicROB <- exp(log(EC50PicROBparen) / PicROBgam)

    Dr <- kb * OB0 / Pic0

    PicROB_val <- E0PicROB + EmaxPicROB * TGFBact^PicROBgam /
                   (TGFBact^PicROBgam + EC50PicROB^PicROBgam)
    ROBin <- Dr * PicROB_val

    #=================================================================
    # OB differentiation (PicOB)
    #=================================================================
    E0PicOB <- FracPicOB * Pic0
    EC50PicOBparen <- (EmaxPicOB * TGFBact0^PicOBgam /
                       (Pic0 - E0PicOB)) - TGFBact0^PicOBgam
    EC50PicOB <- exp(log(EC50PicOBparen) / PicOBgam)

    PicOB_val <- E0PicOB + EmaxPicOB * TGFBact^PicOBgam /
                  (TGFBact^PicOBgam + EC50PicOB^PicOBgam)

    KPT <- 1 * (bigDb / PicOB_val)

    #=================================================================
    # OC differentiation (MeffOC)
    #=================================================================
    EC50MeffOC <- exp(log(M0^kinOCgam * EmaxMeffOC / (1 - E0Meff) -
                          M0^kinOCgam) / kinOCgam)
    MeffOC <- E0Meff + EmaxMeffOC * Mcmplx^kinOCgam /
               (Mcmplx^kinOCgam + EC50MeffOC^kinOCgam)
    kinOC2 <- Da * Pic0 * MeffOC * OC0

    #=================================================================
    # OC (PicOC for TGFBact effect)
    #=================================================================
    E0PicOC <- FracPicOC * Pic0
    EC50PicOCparen <- (EmaxPicOC * TGFBact0^PicOCgam /
                       (Pic0 - E0PicOC)) - TGFBact0^PicOCgam
    EC50PicOC <- exp(log(EC50PicOCparen) / PicOCgam)
    PicOC_val <- E0PicOC + EmaxPicOC * TGFBact^PicOCgam /
                  (TGFBact^PicOCgam + EC50PicOC^PicOCgam)

    #=================================================================
    # OC survival (LsurvOC)
    #=================================================================
    PiL0     <- (k3 / k4) * L0
    PiL      <- Mcmplx / 10
    EC50survInPar <- (E0RANKL - EmaxL) * PiL0^LsurvOCgam /
                      (E0RANKL - 1) - PiL0^LsurvOCgam
    EC50surv <- exp(log(EC50survInPar) / LsurvOCgam)
    LsurvOC  <- E0RANKL - (E0RANKL - EmaxL) * PiL^LsurvOCgam /
                 (PiL^LsurvOCgam + EC50surv^LsurvOCgam)
    KLSoc <- Da * PicOC_val * LsurvOC

    #=================================================================
    # 1-alpha-hydroxylase (AOH)
    #=================================================================
    T66 <- (T67^AlphOHgam + PTHconc0^AlphOHgam) / PTHconc0^AlphOHgam

    # Hydroxyapatite
    kLShap <- 1 / HApMRT
    kHApIn <- kLShap / OB0

    # Bone Ca exchange (slow)
    k15a <- k14a * Qbone0 / Q0
    J14a <- k14a * Qbone
    J15a <- k15a * Q

    # J15: calcium flux from plasma into bone
    J15val <- T15 * P * (1 - FracJ15) + T15 * P * FracJ15 * HAp

    # J42: phosphate into bone (stoichiometric)
    J42 <- 0.464 * J15val

    #=================================================================
    # RANKL production (OsteoEffect + LpthEff)
    #=================================================================
    kinLbase <- koutL * L0
    OsteoEffect <- (OB / OB0)^OsteoEffectGam
    PTH50 <- EmaxLpth * PTHconc0 - PTHconc0
    LpthEff <- EmaxLpth * PTHconc / (PTH50 * OsteoEffect + PTHconc)
    kinL_val <- kinLbase * OsteoEffect * LpthEff

    #=================================================================
    # OPG production
    #=================================================================
    pObase <- kO * O0
    pO_val <- pObase * (ROB1 / ROB0) *
              ((PTHconc + opgPTH50 * (ROB1 / ROB0)) /
               (2 * PTHconc))

    #=================================================================
    # RUNX2, CREB, BCL2
    #=================================================================
    RX2Kin <- RX2Kout0 * RX20v
    EC50PTHRX2x <- (EmaxPTHRX2x * PTHconc0 / (RX2Kout0 - E0rx2Kout)) -
                    PTHconc0
    RX2Kout_val <- E0rx2Kout + EmaxPTHRX2x * PTHconc /
                    (PTHconc + EC50PTHRX2x)

    EC50PTHcreb <- (EmaxPTHcreb * PTHconc0 / (1 - E0crebKin)) -
                    PTHconc0
    crebKin0 <- crebKout * 10.0  # CREB_0 = 10
    crebKin_val <- crebKin0 * (E0crebKin + EmaxPTHcreb * PTHconc /
                                (PTHconc + EC50PTHcreb))
    bcl2Kin <- RX2 * CREB * bcl2Kout

    #=================================================================
    # Phosphate effects
    #=================================================================
    PO4inhPTH <- (C2 / 1.2)^PO4inhPTHgam

    PhosEffTop <- (PhosEff0 - 1) * (1.2^PhosEffGam + PhosEff50^PhosEffGam)
    PhosEffBot <- PhosEff0 * 1.2^PhosEffGam
    PhosEffMax <- PhosEffTop / PhosEffBot
    PhosEff <- PhosEff0 - PhosEffMax * PhosEff0 * C2^PhosEffGam /
                (C2^PhosEffGam + PhosEff50^PhosEffGam)
    PhosEffect <- ifelse(C2 > 1.2, PhosEff, 1.0)

    T68 <- T66 * PTHconc^AlphOHgam /
            (T67^AlphOHgam * PO4inhPTH + PTHconc^AlphOHgam)

    SE <- T65 * T68 * PhosEffect

    #=================================================================
    # Calcitriol-dependent calcium absorption
    #=================================================================
    T36 <- T33 + (T34 - T33) * C8^CaPOgam / (T35^CaPOgam + C8^CaPOgam)
    T37 <- T34 - (T34 - T33) * C8^CaPOgam / (T35^CaPOgam + C8^CaPOgam)

    #=================================================================
    # RENAL CALCIUM HANDLING
    #=================================================================
    CaFilt <- 0.6 * 0.5 * GFR0 * CaConc

    # Estrogen effect on renal Tm
    mtmEST <- (1 - maxTmESTkid) / (1 - 0.1)
    tmEST  <- 1 - mtmEST + mtmEST * EST_s

    ReabsMax <- tmEST * (0.3 * GFR0 * CaConc0 - 0.149997) *
                (Reabs50 + CaConc0) / CaConc0

    # PTH effect on reabsorption
    T17 <- PTHconc0 * T16 - PTHconc0
    ReabsPTHeff <- T16 * PTHconc / (PTHconc + T17)

    # PTH-sensitive reabsorption
    CaReabsActive <- ReabsMax * C1 / (Reabs50 + C1) * ReabsPTHeff

    T20 <- CaFilt - CaReabsActive

    # Calcitriol effect
    T10 <- T7 * C8 / (C8 + T9)

    # Urinary calcium excretion
    J27a <- (2 - T10) * T20
    J27 <- max(J27a, 0)

    ScaEff <- (CaConc0 / CaConc)^ScaEffGam
    T72 <- 90 * ScaEff
    T73 <- T71 * (C8 - T72)
    T74 <- (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73))  # tanh
    T75 <- T70 * (0.85 * (1 + T74) + 0.15)
    T76 <- T70 * (0.85 * (1 - T74) + 0.15)

    #=================================================================
    # Phosphate renal excretion
    #=================================================================
    T47 <- T46 * 0.88 * GFR0
    J48a <- 0.88 * GFR0 * C2 - T47
    J48 <- max(J48a, 0)

    # Phosphate oral absorption and fluxes
    J53 <- T52 * PhosGut
    J54 <- T49 * C2
    J56 <- T55 * IntraPO

    #=================================================================
    # OB apoptosis: PicOBkb -> kbprime -> kbfast, kbslow
    #=================================================================
    E0PicOBkb <- MultPicOBkb * Pic0
    EmaxPicOBkb <- FracPic0kb * Pic0
    EC50PicOBparenKb <- ((E0PicOBkb - EmaxPicOBkb) *
                          TGFBact0^PicOBgamkb) /
                         (E0PicOBkb - Pic0) - TGFBact0^PicOBgamkb
    EC50PicOBkb <- exp(log(EC50PicOBparenKb) / PicOBgamkb)

    PicOBkb <- E0PicOBkb - (E0PicOBkb - EmaxPicOBkb) *
                TGFBact^PicOBgamkb /
                (TGFBact^PicOBgamkb + EC50PicOBkb^PicOBgamkb)

    # Estrogen effect on OB apoptosis
    PicOBkbEff <- (PicOBkb / Pic0) * (1 / EST_s^E2scalePicB1)

    # RUNX2/BCL2 effect on OB apoptosis
    E0RUNX2kbEff <- E0RUNX2kbEffFACT * kb
    RUNX2_val <- ifelse(BCL2 > 105, BCL2 - 90, 10)
    RUNkbMax <- E0RUNX2kbEff * RUNkbMaxFact

    INparen_kb <- (RUNkbMax * RUNX20^RUNkbGAM) /
                   (E0RUNX2kbEff - kb) - RUNX20^RUNkbGAM
    RUNkb50 <- exp(log(INparen_kb) / RUNkbGAM)
    RUNX2kbPrimeEff <- RUNkbMax * RUNX2_val^RUNkbGAM /
                        (RUNX2_val^RUNkbGAM + RUNkb50^RUNkbGAM)

    kbprime <- E0RUNX2kbEff * PicOBkbEff - RUNX2kbPrimeEff
    kbslow  <- kbprime * Frackb
    kbfast  <- (kb * OB0 + kbslow * OBfast0 - kbslow * OB0) / OBfast0
    Frackb2 <- kbfast / kbprime

    #=================================================================
    # Gut calcium equations
    #=================================================================
    T29 <- (T28 * T0 - T310 * T0) / T310
    T31 <- T28 * Tgut / (Tgut + T29)

    # R: calcitriol-dependent gut Ca2+ absorption
    T83 <- Rkid / 0.5

    # J40: calcium flux from gut to plasma
    J40 <- T31 * Tgut * T83 / (Tgut + T81) + T87 * Tgut

    # T85: extent of absorption of orally-administered dose
    T85Rpart <- Rkid^T80 / (Rkid^T80 + T81^T80)
    T85_val  <- T77 * T85Rpart

    #=================================================================
    # Calcitriol equations
    #=================================================================
    INparenCtriol <- ((CtriolMax - CtriolMin) * Calcitriol0^CtriolPTgam) /
                      (CtriolMax - 1) - Calcitriol0^CtriolPTgam
    Ctriol50 <- exp(log(INparenCtriol) / CtriolPTgam)
    CtriolPTeff <- CtriolMax - (CtriolMax - CtriolMin) *
                    C8^CtriolPTgam / (C8^CtriolPTgam + Ctriol50^CtriolPTgam)
    PTin <- PTout * CtriolPTeff

    #=================================================================
    # PTH secretion (complex sigmoid)
    #=================================================================
    FCTD <- (S / 0.5) * PTmax

    INparenCa <- (T58 - T61) * CaConc0^T59 / (T58 - 385) -
                  CaConc0^T59
    T60 <- exp(log(INparenCa) / T59)
    T63 <- T58 - (T58 - T61) * CaConc^T59 /
            (CaConc^T59 + T60^T59)
    SPTH <- T63 * FCTD

    #=================================================================
    # RANK kinetics
    #=================================================================
    kinRNK_val <- kinRNKval * TGFBact^kinRNKgam

    #=================================================================
    # Estrogen (for menopause model, not used for Elagolix forcing)
    #=================================================================
    kinEST <- koutEST
    dEST   <- (kinEST - koutEST * EST) * ESTON

    #=================================================================
    #                  DIFFERENTIAL EQUATIONS
    #=================================================================

    # PTH (pmol)
    dPTH <- SPTH - kout * PTH

    # PTH gland pool
    dS <- (1 - S) * T76 - S * T75

    # PTH gland max capacity
    dPTmax <- PTin - PTout * PTmax

    # Calcitriol (1,25(OH)2D): B
    dB <- AOH - T69 * B

    # 1-alpha-hydroxylase (AOH)
    dAOH <- SE - T64 * AOH

    # Plasma calcium (P, mmol)
    dP <- J14val - J15val - J27 + J40

    # Extracellular phosphate (ECCPhos, mmol)
    dECCPhos <- J41 - J42 - J48 + J53 - J54 + J56

    # Oral calcium in gut (Tgut, mmol)
    dTgut <- OralCa * T85_val - J40

    # Calcitriol-dependent gut Ca absorption (Rkid / R)
    dRkid <- T36 * (1 - Rkid) - T37 * Rkid

    # Hydroxyapatite
    dHAp <- kHApIn * OB - kLShap * HAp

    # Dietary phosphate in gut
    dPhosGut <- OralPhos * F12 - J53

    # Intracellular phosphate
    dIntraPO <- J54 - J56

    # Exchangeable bone Ca
    dQ <- J15val - J14val + J14a - J15a

    # Non-exchangeable bone Ca
    dQbone <- J15a - J14a

    # Cumulative urinary Ca
    dUCA <- J27

    # ROB1: Responding Osteoblasts
    dROB1 <- ROBin * (1 / EST_s)^robGAM - KPT * ROB1

    # OBfast and OBslow
    dOBfast <- (bigDb / PicOB_val) * ROB1 * FracOBfast * Frackb2 -
                kbfast * OBfast
    dOBslow <- (bigDb / PicOB_val) * ROB1 * (1 - FracOBfast) * Frackb -
                kbslow * OBslow

    # OC: Active Osteoclasts
    dOC <- kinOC2 - KLSoc * OC

    # RANKL (L)
    dL <- kinL_val - koutL * L - k1 * Oopg * L + k2 * Ncmplx -
          k3 * RNK * L + k4 * Mcmplx

    # RANK (RNK)
    dRNK <- kinRNK_val - koutRNK * RNK - k3 * RNK * L + k4 * Mcmplx

    # OPG (Oopg / O)
    dOopg <- pO_val - k1 * Oopg * L + k2 * Ncmplx - kO * Oopg

    # RANK-RANKL complex (Mcmplx / M)
    dMcmplx <- k3 * RNK * L - k4 * Mcmplx

    # OPG-RANKL complex (Ncmplx / N)
    dNcmplx <- k1 * Oopg * L - k2 * Ncmplx

    # Latent TGF-beta
    dTGFB <- kinTGF * (OB / OB0)^OBtgfGAM *
              (1 / EST_s)^tgfbGAM -
             koutTGFeqn * EST_s^tgfbactGAM

    # Active TGF-beta
    dTGFBact <- koutTGFeqn * EST_s^tgfbactGAM -
                koutTGFact_val * TGFBact

    # RUNX2
    dRX2 <- RX2Kin - RX2Kout_val * RX2

    # CREB
    dCREB <- crebKin_val - crebKout * CREB

    # BCL2
    dBCL2 <- bcl2Kout * CREB * RX2 - bcl2Kout * BCL2

    # BMD (lumbar spine)
    OBr <- OB / OB0
    OCr <- OC / OC0
    kinBMDls <- koutBMDls   # SS: kinBMDls = koutBMDls * 1.0
    dBMDls <- kinBMDls * OBr^gamOB - koutBMDls * OCr^gamOCls * BMDls

    # ---- Return ---
    list(
      c(dPTH, dS, dPTmax, dB, dAOH, dP, dECCPhos, dTgut, dRkid,
        dHAp, dPhosGut, dIntraPO, dQ, dQbone, dUCA,
        dROB1, dOBfast, dOBslow, dOC,
        dL, dRNK, dOopg, dMcmplx, dNcmplx,
        dTGFB, dTGFBact, dRX2, dCREB, dBCL2,
        dEST, dBMDls),
      OB_total = OBfast + OBslow,
      CTX_pct  = (OC / OC0 - 1) * 100,
      P1NP_pct = (OB / OB0 - 1) * 100,
      BMD_pct  = (BMDls - 1) * 100
    )
  })
}

###############################################################################
#
#  5.  Simulation run function
#
###############################################################################

run_sim <- function(dose_mg, tx_months, fu_months = 0, ss_months = 120,
                    params = NULL) {

  p  <- if (is.null(params)) make_params() else params
  y0 <- make_init(p)

  ss_hr <- ss_months * HOURS_PER_MONTH
  tx_hr <- tx_months * HOURS_PER_MONTH
  fu_hr <- fu_months * HOURS_PER_MONTH

  EST_tx <- predict_E2(dose_mg) / E2_baseline

  # ---- Phase 1 : Steady-state (120 months) ---
  message(sprintf("[SS]  dose=%d mg  EST_tx=%.4f", dose_mg, EST_tx))
  sol_ss <- ode(y0, seq(0, ss_hr, by = 24), qsp_ode, p,
                method = "lsoda", atol = 1e-8, rtol = 1e-8,
                maxsteps = 1e5)
  n_state <- length(y0)
  y_ss <- sol_ss[nrow(sol_ss), 2:(n_state + 1)]
  names(y_ss) <- names(y0)

  OB_ss  <- as.numeric(y_ss["OBfast"] + y_ss["OBslow"])
  OC_ss  <- as.numeric(y_ss["OC"])
  BMD_ss <- as.numeric(y_ss["BMDls"])

  # ---- Phase 2 : Treatment ---
  y_tx <- y_ss
  y_tx["EST"] <- EST_tx

  ev_tx <- function(t, y, parms) { y["EST"] <- EST_tx; y }
  times_tx <- seq(0, tx_hr, by = 24)

  message(sprintf("[TX]  %d months ...", tx_months))
  sol_tx <- ode(y_tx, times_tx, qsp_ode, p,
                method = "lsoda", atol = 1e-8, rtol = 1e-8,
                maxsteps = 1e5,
                events = list(func = ev_tx, time = times_tx))
  res_tx <- as.data.frame(sol_tx)

  # ---- Phase 3 : Follow-up (EST -> 1.0) ---
  if (fu_months > 0) {
    y_fu <- sol_tx[nrow(sol_tx), 2:(n_state + 1)]
    names(y_fu) <- names(y0)
    y_fu["EST"] <- 1.0

    times_fu <- seq(0, fu_hr, by = 24)
    message(sprintf("[FU]  %d months ...", fu_months))
    sol_fu <- ode(y_fu, times_fu, qsp_ode, p,
                  method = "lsoda", atol = 1e-8, rtol = 1e-8,
                  maxsteps = 1e5)
    res_fu <- as.data.frame(sol_fu)
    res_fu$time <- res_fu$time + tx_hr
    res <- rbind(res_tx, res_fu[-1, ])
  } else {
    res <- res_tx
  }

  # ---- Derived columns ---
  res$months       <- res$time / HOURS_PER_MONTH
  res$OB_total     <- res$OBfast + res$OBslow
  res$BMD_chg_pct  <- (res$BMDls / BMD_ss - 1) * 100
  res$CTX_chg_pct  <- (res$OC / OC_ss - 1) * 100
  res$P1NP_chg_pct <- (res$OB_total / OB_ss - 1) * 100
  res$dose_mg      <- dose_mg
  res$dose_label   <- dplyr::case_when(
    dose_mg == 150 ~ "150 mg QD",
    dose_mg == 400 ~ "200 mg BID",
    dose_mg == 600 ~ "300 mg BID",
    TRUE           ~ paste0(dose_mg, " mg")
  )
  res
}

# Helper: extract value at a given month
at_month <- function(df, m, col) {
  idx <- which.min(abs(df$months - m))
  df[[col]][idx]
}

###############################################################################
#
#  6.  Scenario runs
#
###############################################################################

message("\n========== Scenario 1 : 6-month treatment ==========")
s150_6  <- run_sim(150, tx_months = 6)
s400_6  <- run_sim(400, tx_months = 6)

message("\n========== Scenario 2 : 12-month + 6-month FU ==========")
s150_12fu <- run_sim(150, tx_months = 12, fu_months = 6)
s400_12fu <- run_sim(400, tx_months = 12, fu_months = 6)

message("\n========== Scenario 3 : 24-month continuous ==========")
s150_24 <- run_sim(150, tx_months = 24)
s400_24 <- run_sim(400, tx_months = 24)

message("\n========== Scenario 4 : 12-month + 12-month FU ==========")
s150_12fu12 <- run_sim(150, tx_months = 12, fu_months = 12)
s400_12fu12 <- run_sim(400, tx_months = 12, fu_months = 12)

###############################################################################
#
#  7.  Table 3 : Predicted BMD Change (%)
#
###############################################################################

tbl3 <- data.frame(
  Dose  = c("150 mg QD", "200 mg BID"),
  `6mo`  = c(at_month(s150_24, 6,  "BMD_chg_pct"),
             at_month(s400_24, 6,  "BMD_chg_pct")),
  `12mo` = c(at_month(s150_24, 12, "BMD_chg_pct"),
             at_month(s400_24, 12, "BMD_chg_pct")),
  `18mo` = c(at_month(s150_24, 18, "BMD_chg_pct"),
             at_month(s400_24, 18, "BMD_chg_pct")),
  `24mo` = c(at_month(s150_24, 24, "BMD_chg_pct"),
             at_month(s400_24, 24, "BMD_chg_pct")),
  check.names = FALSE
)

cat("\n=== Predicted BMD Change (%) -- Table 3 ===\n")
print(tbl3, digits = 3)
cat("\nTarget (paper):\n")
cat("150 mg QD : -0.61  -0.91  -0.96  -0.91\n")
cat("200 mg BID: -3.47  -4.95  -5.15  -4.97\n")

###############################################################################
#
#  8.  Figure 1 -- Dose vs E2
#
###############################################################################

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

d_seq <- seq(0, 800, 1)
fig1  <- data.frame(dose = d_seq, E2 = sapply(d_seq, predict_E2))

key <- data.frame(
  dose  = c(0, 150, 400, 600),
  label = c("Baseline", "150 mg QD", "200 mg BID", "300 mg BID")
)
key$E2 <- sapply(key$dose, predict_E2)

p1 <- ggplot(fig1, aes(dose, E2)) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  geom_point(data = key, colour = "red", size = 3) +
  geom_text(data = key, aes(label = label), vjust = -1.3, size = 3.2) +
  geom_hline(yintercept = c(exp(2.14), exp(5.20)),
             linetype = "dashed", colour = "grey60") +
  annotate("text", x = 760, y = exp(5.20), label = "E2max", colour = "grey40") +
  annotate("text", x = 760, y = exp(2.14), label = "E2min", colour = "grey40") +
  scale_y_continuous(limits = c(0, 200)) +
  labs(title = "Figure 1 -- Dose-E2 Relationship (Scaled Logistic)",
       x = "Daily Elagolix Dose (mg)", y = "Predicted E2 (pg/mL)") +
  theme_bw(base_size = 13)

ggsave("output/figures/fig1_dose_E2.png", p1, width = 8, height = 5, dpi = 300)
message("Saved: output/figures/fig1_dose_E2.png")

###############################################################################
#
#  9.  Figure 2 -- E2 Suppression vs BMD / CTX / P1NP
#
###############################################################################

doses_f2 <- c(50, 100, 150, 200, 300, 400, 500, 600)
f2 <- do.call(rbind, lapply(doses_f2, function(d) {
  sim <- run_sim(d, tx_months = 12)
  supp <- (1 - predict_E2(d) / E2_baseline) * 100
  data.frame(
    dose           = d,
    E2_suppression = supp,
    BMD_6m   = at_month(sim, 6,  "BMD_chg_pct"),
    BMD_12m  = at_month(sim, 12, "BMD_chg_pct"),
    CTX_6m   = at_month(sim, 6,  "CTX_chg_pct"),
    CTX_12m  = at_month(sim, 12, "CTX_chg_pct"),
    P1NP_6m  = at_month(sim, 6,  "P1NP_chg_pct"),
    P1NP_12m = at_month(sim, 12, "P1NP_chg_pct")
  )
}))

f2_long <- f2 %>%
  select(E2_suppression, BMD_6m, BMD_12m) %>%
  pivot_longer(-E2_suppression, names_to = "tp", values_to = "BMD") %>%
  mutate(tp = ifelse(tp == "BMD_6m", "6 Months", "12 Months"))

p2 <- ggplot(f2_long, aes(E2_suppression, BMD, colour = tp, shape = tp)) +
  geom_point(size = 3) + geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Figure 2 -- E2 Suppression vs BMD Change",
       x = "E2 Suppression (%)", y = "BMD Change (%)",
       colour = "Timepoint", shape = "Timepoint") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig2_E2_vs_BMD.png", p2, width = 8, height = 5, dpi = 300)
message("Saved: output/figures/fig2_E2_vs_BMD.png")

# CTX & P1NP supplementary plots
f2_ctx <- f2 %>%
  select(E2_suppression, CTX_6m, CTX_12m) %>%
  pivot_longer(-E2_suppression, names_to = "tp", values_to = "CTX") %>%
  mutate(tp = ifelse(tp == "CTX_6m", "6 Months", "12 Months"))

p2b <- ggplot(f2_ctx, aes(E2_suppression, CTX, colour = tp)) +
  geom_point(size = 3) + geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Figure 2b -- E2 Suppression vs CTX Change",
       x = "E2 Suppression (%)", y = "CTX Change (%)",
       colour = "Timepoint") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

f2_p1np <- f2 %>%
  select(E2_suppression, P1NP_6m, P1NP_12m) %>%
  pivot_longer(-E2_suppression, names_to = "tp", values_to = "P1NP") %>%
  mutate(tp = ifelse(tp == "P1NP_6m", "6 Months", "12 Months"))

p2c <- ggplot(f2_p1np, aes(E2_suppression, P1NP, colour = tp)) +
  geom_point(size = 3) + geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Figure 2c -- E2 Suppression vs P1NP Change",
       x = "E2 Suppression (%)", y = "P1NP Change (%)",
       colour = "Timepoint") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig2b_E2_vs_CTX.png",  p2b, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig2c_E2_vs_P1NP.png", p2c, width = 8, height = 5, dpi = 300)

###############################################################################
#
# 10.  Figure 3 -- Time-course (12 mo Tx + 6 mo FU)
#
###############################################################################

cols_sel <- c("months", "BMD_chg_pct", "CTX_chg_pct", "P1NP_chg_pct",
              "dose_label")

comb3 <- rbind(
  s150_12fu[, cols_sel],
  s400_12fu[, cols_sel]
)

make_tc <- function(df, yvar, ylab, title) {
  ggplot(df, aes(months, .data[[yvar]], colour = dose_label)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 12, linetype = "dashed", colour = "grey50") +
    annotate("text", x = 12.3, y = max(df[[yvar]]) * 0.9,
             label = "Tx end", size = 3, hjust = 0) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(title = title, x = "Time (months)", y = ylab, colour = "Dose") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
}

p3a <- make_tc(comb3, "BMD_chg_pct",  "BMD Change (%)",
               "Figure 3A -- BMD (12-mo Tx + 6-mo FU)")
p3b <- make_tc(comb3, "CTX_chg_pct",  "CTX Change (%)",
               "Figure 3B -- CTX (12-mo Tx + 6-mo FU)")
p3c <- make_tc(comb3, "P1NP_chg_pct", "P1NP Change (%)",
               "Figure 3C -- P1NP (12-mo Tx + 6-mo FU)")

ggsave("output/figures/fig3a_BMD_tc.png",  p3a, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig3b_CTX_tc.png",  p3b, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig3c_P1NP_tc.png", p3c, width = 8, height = 5, dpi = 300)
message("Saved: fig3a/b/c")

###############################################################################
#
# 11.  Figure 4 -- 24-month continuous
#
###############################################################################

comb4 <- rbind(
  s150_24[, cols_sel],
  s400_24[, cols_sel]
)

p4a <- ggplot(comb4, aes(months, BMD_chg_pct, colour = dose_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  labs(title = "Figure 4A -- BMD: 24-Month Continuous Treatment",
       x = "Time (months)", y = "BMD Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

p4b <- ggplot(comb4, aes(months, CTX_chg_pct, colour = dose_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Figure 4B -- CTX: 24-Month",
       x = "Time (months)", y = "CTX Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

p4c <- ggplot(comb4, aes(months, P1NP_chg_pct, colour = dose_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Figure 4C -- P1NP: 24-Month",
       x = "Time (months)", y = "P1NP Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig4a_BMD_24m.png",  p4a, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig4b_CTX_24m.png",  p4b, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig4c_P1NP_24m.png", p4c, width = 8, height = 5, dpi = 300)
message("Saved: fig4a/b/c")

###############################################################################
#
# 12.  Recovery Figure -- 12 mo Tx + 12 mo FU
#
###############################################################################

comb5 <- rbind(
  s150_12fu12[, cols_sel],
  s400_12fu12[, cols_sel]
)

p5 <- ggplot(comb5, aes(months, BMD_chg_pct, colour = dose_label)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 12, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 12.3, y = 0.3, label = "Tx end", size = 3, hjust = 0) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  labs(title = "Figure 5 -- BMD Recovery (12-mo Tx + 12-mo FU)",
       x = "Time (months)", y = "BMD Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig5_BMD_recovery.png", p5, width = 8, height = 5, dpi = 300)
message("Saved: fig5_BMD_recovery.png")

###############################################################################
#
# 13.  Summary
#
###############################################################################

cat("\n")
cat(strrep("=", 62), "\n")
cat("  SUMMARY\n")
cat(strrep("=", 62), "\n\n")

cat(sprintf("State variables : %d\n", length(make_init(make_params()))))
cat(sprintf("Parameters      : %d\n", length(make_params())))
cat(sprintf("Baseline E2     : %.1f pg/mL\n", E2_baseline))
cat(sprintf("150 mg QD  E2   : %.1f pg/mL  (suppression %.1f%%)\n",
            predict_E2(150), (1 - predict_E2(150)/E2_baseline)*100))
cat(sprintf("200 mg BID E2   : %.1f pg/mL  (suppression %.1f%%)\n",
            predict_E2(400), (1 - predict_E2(400)/E2_baseline)*100))
cat(sprintf("300 mg BID E2   : %.1f pg/mL  (suppression %.1f%%)\n",
            predict_E2(600), (1 - predict_E2(600)/E2_baseline)*100))

cat("\nFigures saved to output/figures/:\n")
cat("  fig1_dose_E2.png        Dose-E2 relationship\n")
cat("  fig2_E2_vs_BMD.png      E2 suppression vs BMD\n")
cat("  fig2b_E2_vs_CTX.png     E2 suppression vs CTX\n")
cat("  fig2c_E2_vs_P1NP.png    E2 suppression vs P1NP\n")
cat("  fig3a_BMD_tc.png        BMD time-course (12+6)\n")
cat("  fig3b_CTX_tc.png        CTX time-course (12+6)\n")
cat("  fig3c_P1NP_tc.png       P1NP time-course (12+6)\n")
cat("  fig4a_BMD_24m.png       BMD 24-month\n")
cat("  fig4b_CTX_24m.png       CTX 24-month\n")
cat("  fig4c_P1NP_24m.png      P1NP 24-month\n")
cat("  fig5_BMD_recovery.png   BMD recovery (12+12)\n")
cat("\nDone.\n")
