###############################################################################
# app.R  —  Elagolix QSP Model  ·  Interactive Shiny Dashboard
# -----------------------------------------------------------------------
# Model core: faithful R port of MetrumRG OpenBoneMin.cpp
# Stodtmann et al. (2021) Clin Transl Sci. 14:1611-1619
###############################################################################

library(shiny)
library(bslib)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ── Constants ────────────────────────────────────────────────────────────────
HOURS_PER_MONTH <- 730.5

# ── Module A : Dose -> E2 (Scaled Logistic, Eq.3) ───────────────────────────
dose_e2_params <- list(slope = 0.00894, logE2max = 5.20, logE2min = 2.14)

predict_E2 <- function(daily_dose, p = dose_e2_params) {
  E2min <- exp(p$logE2min)
  E2max <- exp(p$logE2max)
  E2min + (E2max - E2min) / (1 + exp(p$slope * daily_dose))
}

E2_baseline <- predict_E2(0)

# ── Module B : Parameters (faithful to OpenBoneMin.cpp) ──────────────────────
make_params <- function() {
  list(
    OBtot0=0.00501324, k1=6.24e-6, k2=0.112013, k3=6.24e-6, k4=0.112013,
    kO=15.8885, kb=0.000605516, Pic0=0.228142,
    FracOBfast=0.797629, Frackb=0.313186, Da=0.7/24,
    # TGF-beta
    OBtgfGAM=0.0111319, koutTGF0=2.98449e-5, koutTGFGam=0.919131,
    OCtgfGAM=0.593891,
    # ROB differentiation
    EmaxPicROB=3.9745, PicROBgam=1.80968, FracPicROB=0.883824,
    # OB differentiation
    PicOBgam=0.122313, FracPicOB=0.000244818, EmaxPicOB=0.251636,
    # OB apoptosis (TGF-beta -> kbprime)
    PicOBgamkb=2.92375, MultPicOBkb=3.11842, FracPic0kb=0.764028,
    E0RUNX2kbEffFACT=1.01, RUNkbGAM=3.81644, RUNkbMaxFact=0.638114,
    # OC
    E0Meff=0.388267, EmaxMeffOC=3.15667, kinOCgam=8.53065,
    EmaxPicOC=1.9746, FracPicOC=0.878215, PicOCgam=1.0168,
    MOCratioGam=0.603754,
    # OC survival
    E0RANKL=3.80338, EmaxL=0.469779, LsurvOCgam=3.09023, LsurvOCCgam=3.0923,
    # Ca/Phos homeostasis
    CaDay=88.0, V1=14.0,
    FracJ14=0.107763, J14OCmax=0.543488, J14OCgam=1.6971,
    FracJ15=0.114376, k14a=2.44437e-5, HApMRT=3.60609,
    # Oral Ca/Phos
    OralCa=24.055/24, OralPhos=10.5/24, F12=0.7,
    # Gut Ca absorption
    T28=0.9, T310=0.105929, T81=0.75, T87=0.0495,
    T77=0.909359, T80=4.0, T43=1.03856, T45=3.85, T84=0.25,
    # Renal Ca
    Reabs50=1.57322, T16=1.06147, maxTmESTkid=0.923737,
    T7=2.0, T9=90.0, T70=0.01, T71=0.03,
    # Calcitriol / 1,25(OH)2D
    T33=0.003, T34=0.037, T35=90.0, CaPOgam=1.0,
    AlphOHgam=0.111241, T64=0.05, T65=6.3, T67=1.54865, T69=0.10,
    CtriolPTgam=12.5033, CtriolMax=4.1029, CtriolMin=0.9,
    # Phosphate effects
    ScaEffGam=0.9, PhosEff0=1.52493, PhosEff50=1.3021,
    PhosEffGam=8.25229, PO4inhPTHgam=0.0,
    # Phosphate fluxes
    T46=1.142, T49=51.8, T52=0.365, T55=0.019268,
    # PTH
    kout=100/14, PTout=1.604e-4, T58=6249.09, T59=11.7387, T61=96.25,
    # RANKL / OPG
    koutL=0.00293273, kinRNKgam=0.151825, koutRNK=0.00323667,
    opgPTH50=3.85, EmaxLpth=1.30721, OsteoEffectGam=0.173833,
    # Intracellular signaling
    RUNX20=10.0, RX2Kout0=0.693, E0rx2Kout=0.125, EmaxPTHRX2x=5.0,
    E0crebKin=0.5, EmaxPTHcreb=3.39745, crebKout=0.00279513, bcl2Kout=0.693,
    # BMD
    koutBMDls=0.000397, gamOB=0.0793, gamOCls=0.14,
    # Estrogen extension
    ESTON=0.0, koutEST=0.05776227,
    tgfbGAM=0.0374, tgfbactGAM=0.045273, robGAM=0.16, obGAM=0.000012,
    E2scalePicB1=1.16832e-5,
    # GFR
    GFR0=100/16.667,
    # Denosumab (unused)
    kdenosl=2.0e-06
  )
}

# ── Initial Conditions (31 states, matches OpenBoneMin.cpp) ──────────────────
make_init <- function(p) {
  L0 <- 0.4; RNK0 <- 10; O0 <- 4
  c(PTH=53.9, S=0.5, PTmax=1, B=1260, AOH=126,
    P=32.9, ECCPhos=16.8, Tgut=1.58471, Rkid=0.5, HAp=1,
    PhosGut=0.839, IntraPO=3226, Q=100, Qbone=24900, UCA=0,
    ROB1=0.00104122,
    OBfast=p$OBtot0*p$FracOBfast,
    OBslow=p$OBtot0*(1-p$FracOBfast),
    OC=0.00115398,
    L=L0, RNK=RNK0, Oopg=O0,
    Mcmplx=p$k3*RNK0*L0/p$k4,
    Ncmplx=p$k1*O0*L0/p$k2,
    TGFB=p$Pic0*1000, TGFBact=p$Pic0,
    RX2=10, CREB=10, BCL2=100,
    EST=1, BMDls=1)
}

# ── ODE System (faithful port of OpenBoneMin.cpp [ODE] block) ────────────────
qsp_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {

    # Baseline reference values
    OB0 <- OBtot0; OC0 <- 0.00115398; ROB0 <- 0.00104122
    L0 <- 0.4; RNK0 <- 10; O0 <- 4; PTH0 <- 53.9; B0 <- 1260; P0 <- 32.9
    M0 <- k3*RNK0*L0/k4
    OBfast0 <- OBtot0*FracOBfast; OBslow0 <- OBtot0*(1-FracOBfast)
    TGFBact0 <- Pic0; TGFB0 <- Pic0*1000; RX20v <- RUNX20
    Q0 <- 100; Qbone0 <- 24900; T0 <- 1.58471
    Calcitriol0 <- B0/V1

    OB <- OBfast + OBslow
    EST_s <- max(EST, 1e-6)

    # Concentrations
    CaConc0 <- P0/V1; PTHconc0 <- PTH0/V1
    PTHconc <- PTH/V1; CaConc <- P/V1
    C1 <- P/V1; C2 <- ECCPhos/V1; C8 <- B/V1

    # Derived parameters from SS
    T13 <- (CaDay/24)/Q0
    T15 <- CaDay/(CaConc0*V1*24)
    J14OC50 <- exp(log((J14OCmax*OC0^J14OCgam/T13) - OC0^J14OCgam)/J14OCgam)
    OCeqn <- (J14OCmax*OC^J14OCgam)/(OC^J14OCgam + J14OC50^J14OCgam)
    kinRNKval <- (koutRNK*RNK0 + k3*RNK0*L0 - k4*M0)/TGFBact0^kinRNKgam

    MOCratio <- Mcmplx/max(OC, 1e-15)
    MOCratio0 <- M0/OC0
    MOCratioEff <- (MOCratio/MOCratio0)^MOCratioGam

    # J14: Ca from bone to plasma
    J14OCdepend <- OCeqn*Q0*FracJ14*MOCratioEff
    J14val <- T13*Q0*(1-FracJ14) + J14OCdepend
    J41 <- 0.464*J14val

    bigDb <- kb*OB0*Pic0/ROB0

    # TGF-beta
    kinTGF <- koutTGF0*TGFB0
    koutTGF_v <- koutTGF0*(TGFB/TGFB0)^koutTGFGam
    koutTGFact_v <- koutTGF0*1000
    koutTGFeqn <- koutTGF_v*TGFB*(OC/OC0)^OCtgfGAM

    # ROB differentiation (PicROB)
    E0PicROB <- FracPicROB*Pic0
    EC50PicROBp <- (EmaxPicROB*TGFBact0^PicROBgam/(Pic0-E0PicROB)) - TGFBact0^PicROBgam
    EC50PicROB <- exp(log(EC50PicROBp)/PicROBgam)
    Dr <- kb*OB0/Pic0
    PicROB_v <- E0PicROB + EmaxPicROB*TGFBact^PicROBgam/(TGFBact^PicROBgam + EC50PicROB^PicROBgam)
    ROBin <- Dr*PicROB_v

    # OB differentiation (PicOB)
    E0PicOB <- FracPicOB*Pic0
    EC50PicOBp <- (EmaxPicOB*TGFBact0^PicOBgam/(Pic0-E0PicOB)) - TGFBact0^PicOBgam
    EC50PicOB <- exp(log(EC50PicOBp)/PicOBgam)
    PicOB_v <- E0PicOB + EmaxPicOB*TGFBact^PicOBgam/(TGFBact^PicOBgam + EC50PicOB^PicOBgam)
    KPT <- bigDb/PicOB_v

    # OC differentiation (MeffOC)
    EC50MeffOC <- exp(log(M0^kinOCgam*EmaxMeffOC/(1-E0Meff) - M0^kinOCgam)/kinOCgam)
    MeffOC <- E0Meff + EmaxMeffOC*Mcmplx^kinOCgam/(Mcmplx^kinOCgam + EC50MeffOC^kinOCgam)
    kinOC2 <- Da*Pic0*MeffOC*OC0

    # OC PicOC
    E0PicOC <- FracPicOC*Pic0
    EC50PicOCp <- (EmaxPicOC*TGFBact0^PicOCgam/(Pic0-E0PicOC)) - TGFBact0^PicOCgam
    EC50PicOC <- exp(log(EC50PicOCp)/PicOCgam)
    PicOC_v <- E0PicOC + EmaxPicOC*TGFBact^PicOCgam/(TGFBact^PicOCgam + EC50PicOC^PicOCgam)

    # OC survival (LsurvOC)
    PiL0 <- (k3/k4)*L0; PiL <- Mcmplx/10
    EC50survInPar <- (E0RANKL-EmaxL)*PiL0^LsurvOCgam/(E0RANKL-1) - PiL0^LsurvOCgam
    EC50surv <- exp(log(EC50survInPar)/LsurvOCgam)
    LsurvOC <- E0RANKL - (E0RANKL-EmaxL)*PiL^LsurvOCgam/(PiL^LsurvOCgam + EC50surv^LsurvOCgam)
    KLSoc <- Da*PicOC_v*LsurvOC

    # 1-alpha-hydroxylase
    T66 <- (T67^AlphOHgam + PTHconc0^AlphOHgam)/PTHconc0^AlphOHgam

    # Hydroxyapatite
    kLShap <- 1/HApMRT; kHApIn <- kLShap/OB0

    # Bone Ca exchange (slow)
    k15a <- k14a*Qbone0/Q0
    J14a <- k14a*Qbone; J15a <- k15a*Q

    # J15: Ca plasma to bone
    J15val <- T15*P*(1-FracJ15) + T15*P*FracJ15*HAp
    J42 <- 0.464*J15val

    # RANKL production
    kinLbase <- koutL*L0
    OsteoEffect <- (OB/OB0)^OsteoEffectGam
    PTH50 <- EmaxLpth*PTHconc0 - PTHconc0
    LpthEff <- EmaxLpth*PTHconc/(PTH50*OsteoEffect + PTHconc)
    kinL_v <- kinLbase*OsteoEffect*LpthEff

    # OPG production
    pObase <- kO*O0
    pO_v <- pObase*(ROB1/ROB0)*((PTHconc + opgPTH50*(ROB1/ROB0))/(2*PTHconc))

    # RUNX2, CREB, BCL2
    RX2Kin <- RX2Kout0*RX20v
    EC50PTHRX2x <- (EmaxPTHRX2x*PTHconc0/(RX2Kout0-E0rx2Kout)) - PTHconc0
    RX2Kout_v <- E0rx2Kout + EmaxPTHRX2x*PTHconc/(PTHconc + EC50PTHRX2x)
    EC50PTHcreb <- (EmaxPTHcreb*PTHconc0/(1-E0crebKin)) - PTHconc0
    crebKin0 <- crebKout*10
    crebKin_v <- crebKin0*(E0crebKin + EmaxPTHcreb*PTHconc/(PTHconc + EC50PTHcreb))

    # Phosphate effects
    PO4inhPTH <- (C2/1.2)^PO4inhPTHgam
    PhosEffTop <- (PhosEff0-1)*(1.2^PhosEffGam + PhosEff50^PhosEffGam)
    PhosEffBot <- PhosEff0*1.2^PhosEffGam
    PhosEffMax <- PhosEffTop/PhosEffBot
    PhosEff_v <- PhosEff0 - PhosEffMax*PhosEff0*C2^PhosEffGam/(C2^PhosEffGam + PhosEff50^PhosEffGam)
    PhosEffect <- ifelse(C2 > 1.2, PhosEff_v, 1.0)

    T68 <- T66*PTHconc^AlphOHgam/(T67^AlphOHgam*PO4inhPTH + PTHconc^AlphOHgam)
    SE <- T65*T68*PhosEffect

    # Calcitriol-dependent Ca absorption
    T36 <- T33 + (T34-T33)*C8^CaPOgam/(T35^CaPOgam + C8^CaPOgam)
    T37 <- T34 - (T34-T33)*C8^CaPOgam/(T35^CaPOgam + C8^CaPOgam)

    # Renal calcium handling
    CaFilt <- 0.6*0.5*GFR0*CaConc
    mtmEST <- (1-maxTmESTkid)/(1-0.1)
    tmEST <- 1 - mtmEST + mtmEST*EST_s
    ReabsMax <- tmEST*(0.3*GFR0*CaConc0 - 0.149997)*(Reabs50 + CaConc0)/CaConc0
    T17 <- PTHconc0*T16 - PTHconc0
    ReabsPTHeff <- T16*PTHconc/(PTHconc + T17)
    CaReabsActive <- ReabsMax*C1/(Reabs50 + C1)*ReabsPTHeff
    T20 <- CaFilt - CaReabsActive
    T10 <- T7*C8/(C8 + T9)
    J27a <- (2-T10)*T20
    J27 <- max(J27a, 0)

    ScaEff <- (CaConc0/CaConc)^ScaEffGam
    T72 <- 90*ScaEff; T73 <- T71*(C8-T72)
    T74 <- (exp(T73)-exp(-T73))/(exp(T73)+exp(-T73))
    T75 <- T70*(0.85*(1+T74)+0.15)
    T76 <- T70*(0.85*(1-T74)+0.15)

    # Phosphate renal
    T47 <- T46*0.88*GFR0
    J48a <- 0.88*GFR0*C2 - T47; J48 <- max(J48a, 0)
    J53 <- T52*PhosGut; J54 <- T49*C2; J56 <- T55*IntraPO

    # OB apoptosis: PicOBkb -> kbprime -> kbfast, kbslow
    E0PicOBkb <- MultPicOBkb*Pic0
    EmaxPicOBkb_v <- FracPic0kb*Pic0
    EC50PicOBparenKb <- ((E0PicOBkb-EmaxPicOBkb_v)*TGFBact0^PicOBgamkb)/(E0PicOBkb-Pic0) - TGFBact0^PicOBgamkb
    EC50PicOBkb <- exp(log(EC50PicOBparenKb)/PicOBgamkb)
    PicOBkb <- E0PicOBkb - (E0PicOBkb-EmaxPicOBkb_v)*TGFBact^PicOBgamkb/(TGFBact^PicOBgamkb + EC50PicOBkb^PicOBgamkb)
    PicOBkbEff <- (PicOBkb/Pic0)*(1/EST_s^E2scalePicB1)

    E0RUNX2kbEff <- E0RUNX2kbEffFACT*kb
    RUNX2_v <- ifelse(BCL2 > 105, BCL2 - 90, 10)
    RUNkbMax <- E0RUNX2kbEff*RUNkbMaxFact
    INparen_kb <- (RUNkbMax*RUNX20^RUNkbGAM)/(E0RUNX2kbEff-kb) - RUNX20^RUNkbGAM
    RUNkb50 <- exp(log(INparen_kb)/RUNkbGAM)
    RUNX2kbPrimeEff <- RUNkbMax*RUNX2_v^RUNkbGAM/(RUNX2_v^RUNkbGAM + RUNkb50^RUNkbGAM)

    kbprime <- E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff
    kbslow_v <- kbprime*Frackb
    kbfast_v <- (kb*OB0 + kbslow_v*OBfast0 - kbslow_v*OB0)/OBfast0
    Frackb2 <- kbfast_v/kbprime

    # Gut calcium
    T29 <- (T28*T0 - T310*T0)/T310
    T31 <- T28*Tgut/(Tgut + T29)
    T83 <- Rkid/0.5
    J40 <- T31*Tgut*T83/(Tgut + T81) + T87*Tgut
    T85Rpart <- Rkid^T80/(Rkid^T80 + T81^T80)
    T85_v <- T77*T85Rpart

    # Calcitriol equations
    INparenCtriol <- ((CtriolMax-CtriolMin)*Calcitriol0^CtriolPTgam)/(CtriolMax-1) - Calcitriol0^CtriolPTgam
    Ctriol50 <- exp(log(INparenCtriol)/CtriolPTgam)
    CtriolPTeff <- CtriolMax - (CtriolMax-CtriolMin)*C8^CtriolPTgam/(C8^CtriolPTgam + Ctriol50^CtriolPTgam)
    PTin <- PTout*CtriolPTeff

    # PTH secretion (complex sigmoid)
    FCTD <- (S/0.5)*PTmax
    INparenCa <- (T58-T61)*CaConc0^T59/(T58-385) - CaConc0^T59
    T60 <- exp(log(INparenCa)/T59)
    T63 <- T58 - (T58-T61)*CaConc^T59/(CaConc^T59 + T60^T59)
    SPTH <- T63*FCTD

    # RANK kinetics
    kinRNK_v <- kinRNKval*TGFBact^kinRNKgam

    # Estrogen
    kinEST <- koutEST
    dEST <- (kinEST - koutEST*EST)*ESTON

    #=============== DIFFERENTIAL EQUATIONS ===============

    dPTH <- SPTH - kout*PTH
    dS <- (1-S)*T76 - S*T75
    dPTmax <- PTin - PTout*PTmax
    dB <- AOH - T69*B
    dAOH <- SE - T64*AOH
    dP <- J14val - J15val - J27 + J40
    dECCPhos <- J41 - J42 - J48 + J53 - J54 + J56
    dTgut <- OralCa*T85_v - J40
    dRkid <- T36*(1-Rkid) - T37*Rkid
    dHAp <- kHApIn*OB - kLShap*HAp
    dPhosGut <- OralPhos*F12 - J53
    dIntraPO <- J54 - J56
    dQ <- J15val - J14val + J14a - J15a
    dQbone <- J15a - J14a
    dUCA <- J27
    dROB1 <- ROBin*(1/EST_s)^robGAM - KPT*ROB1
    dOBfast <- (bigDb/PicOB_v)*ROB1*FracOBfast*Frackb2 - kbfast_v*OBfast
    dOBslow <- (bigDb/PicOB_v)*ROB1*(1-FracOBfast)*Frackb - kbslow_v*OBslow
    dOC <- kinOC2 - KLSoc*OC
    dL <- kinL_v - koutL*L - k1*Oopg*L + k2*Ncmplx - k3*RNK*L + k4*Mcmplx
    dRNK <- kinRNK_v - koutRNK*RNK - k3*RNK*L + k4*Mcmplx
    dOopg <- pO_v - k1*Oopg*L + k2*Ncmplx - kO*Oopg
    dMcmplx <- k3*RNK*L - k4*Mcmplx
    dNcmplx <- k1*Oopg*L - k2*Ncmplx
    dTGFB <- kinTGF*(OB/OB0)^OBtgfGAM*(1/EST_s)^tgfbGAM - koutTGFeqn*EST_s^tgfbactGAM
    dTGFBact <- koutTGFeqn*EST_s^tgfbactGAM - koutTGFact_v*TGFBact
    dRX2 <- RX2Kin - RX2Kout_v*RX2
    dCREB <- crebKin_v - crebKout*CREB
    dBCL2 <- bcl2Kout*CREB*RX2 - bcl2Kout*BCL2
    dBMDls <- koutBMDls*(OB/OB0)^gamOB - koutBMDls*(OC/OC0)^gamOCls*BMDls

    list(c(dPTH,dS,dPTmax,dB,dAOH,dP,dECCPhos,dTgut,dRkid,
           dHAp,dPhosGut,dIntraPO,dQ,dQbone,dUCA,
           dROB1,dOBfast,dOBslow,dOC,
           dL,dRNK,dOopg,dMcmplx,dNcmplx,
           dTGFB,dTGFBact,dRX2,dCREB,dBCL2,dEST,dBMDls),
         OB_total=OBfast+OBslow,
         CTX_pct=(OC/OC0-1)*100,
         P1NP_pct=(OB/OB0-1)*100,
         BMD_pct=(BMDls-1)*100)
  })
}

# ── run_sim ──────────────────────────────────────────────────────────────────
run_sim <- function(dose_mg, tx_months, fu_months=0, ss_months=120) {
  p <- make_params(); y0 <- make_init(p)
  ss_hr <- ss_months*HOURS_PER_MONTH
  tx_hr <- tx_months*HOURS_PER_MONTH
  fu_hr <- fu_months*HOURS_PER_MONTH
  EST_tx <- predict_E2(dose_mg)/E2_baseline

  sol_ss <- ode(y0, seq(0,ss_hr,by=24), qsp_ode, p,
                method="lsoda", atol=1e-8, rtol=1e-8, maxsteps=1e5)
  n_state <- length(y0)
  y_ss <- sol_ss[nrow(sol_ss), 2:(n_state+1)]; names(y_ss) <- names(y0)
  OB_ss <- as.numeric(y_ss["OBfast"]+y_ss["OBslow"])
  OC_ss <- as.numeric(y_ss["OC"]); BMD_ss <- as.numeric(y_ss["BMDls"])

  y_tx <- y_ss; y_tx["EST"] <- EST_tx
  ev_tx <- function(t,y,parms){y["EST"] <- EST_tx; y}
  times_tx <- seq(0,tx_hr,by=24)
  sol_tx <- ode(y_tx, times_tx, qsp_ode, p, method="lsoda",
                atol=1e-8, rtol=1e-8, maxsteps=1e5,
                events=list(func=ev_tx, time=times_tx))
  res <- as.data.frame(sol_tx)

  if(fu_months > 0){
    y_fu <- sol_tx[nrow(sol_tx), 2:(n_state+1)]; names(y_fu) <- names(y0)
    y_fu["EST"] <- 1
    sol_fu <- ode(y_fu, seq(0,fu_hr,by=24), qsp_ode, p,
                  method="lsoda", atol=1e-8, rtol=1e-8, maxsteps=1e5)
    res_fu <- as.data.frame(sol_fu); res_fu$time <- res_fu$time + tx_hr
    res <- rbind(res, res_fu[-1,])
  }
  res$months       <- res$time/HOURS_PER_MONTH
  res$OB_total     <- res$OBfast + res$OBslow
  res$BMD_chg_pct  <- (res$BMDls/BMD_ss-1)*100
  res$CTX_chg_pct  <- (res$OC/OC_ss-1)*100
  res$P1NP_chg_pct <- (res$OB_total/OB_ss-1)*100
  res$dose_mg      <- dose_mg
  res$dose_label   <- dplyr::case_when(
    dose_mg==150~"150 mg QD", dose_mg==400~"200 mg BID",
    dose_mg==600~"300 mg BID", TRUE~paste0(dose_mg," mg"))
  res
}

at_month <- function(df,m,col) df[[col]][which.min(abs(df$months-m))]

# ── Paper target (Table 3) ──────────────────────────────────────────────────
paper_target <- data.frame(
  Dose = c("150 mg QD","200 mg BID"),
  `6mo`  = c(-0.61, -3.47),
  `12mo` = c(-0.91, -4.95),
  `18mo` = c(-0.96, -5.15),
  `24mo` = c(-0.91, -4.97),
  check.names = FALSE
)

###############################################################################
#                         Shiny UI
###############################################################################

ui <- page_navbar(
  title = "Elagolix QSP Model",
  theme = bs_theme(bootswatch = "flatly", base_font = font_google("Noto Sans KR")),
  window_title = "Elagolix QSP | Stodtmann 2021",

  # Tab 1 : Dose-E2
  nav_panel("Dose-E2 Model",
    layout_sidebar(
      sidebar = sidebar(
        title = "Module A Parameters",
        numericInput("slope", "slope", 0.00894, step = 0.0001),
        numericInput("logE2max", "log(E2max)", 5.20, step = 0.01),
        numericInput("logE2min", "log(E2min)", 2.14, step = 0.01),
        hr(),
        sliderInput("dose_highlight", "Highlight dose (mg)",
                     min = 0, max = 800, value = 150, step = 10),
        hr(),
        helpText("Scaled Logistic (Eq.3, Stodtmann 2021)")
      ),
      card(card_header("Figure 1 -- Dose vs E2"), plotOutput("fig1_plot", height = "450px")),
      layout_column_wrap(
        width = 1/3,
        value_box("Baseline E2", textOutput("vb_baseline"), theme = "primary"),
        value_box("Predicted E2", textOutput("vb_pred_e2"), theme = "info"),
        value_box("Suppression", textOutput("vb_supp"), theme = "warning")
      )
    )
  ),

  # Tab 2 : QSP Simulation
  nav_panel("QSP Simulation",
    layout_sidebar(
      sidebar = sidebar(
        title = "Simulation Settings", width = 300,
        selectInput("dose_sel", "Elagolix Dose",
                     choices = c("150 mg QD"=150, "200 mg BID"=400, "300 mg BID"=600),
                     selected = 150),
        sliderInput("tx_mo", "Treatment (months)", 1, 36, 12, step = 1),
        sliderInput("fu_mo", "Follow-up (months)", 0, 24, 6, step = 1),
        hr(),
        actionButton("run_btn", "Run Simulation",
                     class = "btn-primary btn-lg w-100", icon = icon("play")),
        hr(),
        helpText("Steady-state warm-up: 120 months (auto).")
      ),
      navset_card_tab(
        nav_panel("BMD",  plotOutput("plot_bmd",  height = "420px")),
        nav_panel("CTX",  plotOutput("plot_ctx",  height = "420px")),
        nav_panel("P1NP", plotOutput("plot_p1np", height = "420px")),
        nav_panel("All Markers", plotOutput("plot_all", height = "480px"))
      ),
      card(card_header("Biomarker Values at Selected Timepoints"),
           DTOutput("tbl_biomarker"))
    )
  ),

  # Tab 3 : Multi-dose Comparison
  nav_panel("Multi-Dose Comparison",
    layout_sidebar(
      sidebar = sidebar(
        title = "Comparison Settings", width = 300,
        checkboxGroupInput("doses_cmp", "Select Doses",
                            choices = c("150 mg QD"=150, "200 mg BID"=400, "300 mg BID"=600),
                            selected = c(150, 400)),
        sliderInput("tx_mo_cmp", "Treatment (months)", 1, 36, 24, step = 1),
        sliderInput("fu_mo_cmp", "Follow-up (months)", 0, 24, 0, step = 1),
        hr(),
        actionButton("run_cmp", "Run Comparison",
                     class = "btn-success btn-lg w-100", icon = icon("layer-group")),
        hr(),
        helpText("Compare multiple regimens side-by-side.")
      ),
      navset_card_tab(
        nav_panel("BMD",  plotOutput("cmp_bmd",  height = "420px")),
        nav_panel("CTX",  plotOutput("cmp_ctx",  height = "420px")),
        nav_panel("P1NP", plotOutput("cmp_p1np", height = "420px"))
      ),
      card(card_header("Table 3 Comparison -- BMD Change (%)"),
           DTOutput("tbl3_cmp"))
    )
  ),

  # Tab 4 : Validation Table
  nav_panel("Validation (Table 3)",
    card(card_header("Paper Target: Predicted BMD Change (%)"),
         DTOutput("tbl_paper_target")),
    card(card_header("Model Prediction (run Comparison first)"),
         DTOutput("tbl_model_pred"))
  ),

  # Tab 5 : About
  nav_panel("About",
    card(
      card_header("Model Information"),
      tags$ul(
        tags$li(tags$b("Paper:"), " Stodtmann et al. (2021) Clin Transl Sci. 14:1611-1619"),
        tags$li(tags$b("Foundation:"), " Peterson & Riggs (2010) + Riggs et al. (2012)"),
        tags$li(tags$b("Reference impl.:"), " MetrumRG OpenBoneMin (faithful R port)"),
        tags$li(tags$b("State variables:"), " 31"),
        tags$li(tags$b("Parameters:"), sprintf(" %d", length(make_params()))),
        tags$li(tags$b("ODE solver:"), " deSolve::lsoda (atol/rtol = 1e-8)")
      ),
      hr(),
      tags$h5("Module A -- Dose-E2 (Scaled Logistic)"),
      tags$pre("E2 = exp(logE2min) + (exp(logE2max) - exp(logE2min)) / (1 + exp(slope * Dose))"),
      tags$h5("Module B -- Calcium Homeostasis & Bone Remodeling QSP"),
      tags$p("31 ODEs covering: PTH, calcitriol, Ca/Phos homeostasis,
              ROB, OBfast, OBslow, OC, RANKL/RANK/OPG, TGF-beta,
              RUNX2, CREB, BCL2, estrogen, BMD (lumbar spine)."),
      tags$p("All equations faithfully ported from OpenBoneMin.cpp source.")
    )
  )
)

###############################################################################
#                         Shiny Server
###############################################################################

server <- function(input, output, session) {

  # ========== Tab 1 : Dose-E2 =============================================

  output$fig1_plot <- renderPlot({
    p <- list(slope = input$slope, logE2max = input$logE2max,
              logE2min = input$logE2min)
    d_seq <- seq(0, 800, 1)
    df <- data.frame(dose = d_seq, E2 = sapply(d_seq, predict_E2, p = p))
    key <- data.frame(
      dose = c(0, 150, 400, 600),
      label = c("Baseline", "150 mg QD", "200 mg BID", "300 mg BID"))
    key$E2 <- sapply(key$dose, predict_E2, p = p)
    hl <- data.frame(dose = input$dose_highlight,
                     E2 = predict_E2(input$dose_highlight, p))
    ggplot(df, aes(dose, E2)) +
      geom_line(colour = "steelblue", linewidth = 1.2) +
      geom_point(data = key, colour = "red", size = 3) +
      geom_text(data = key, aes(label = label), vjust = -1.3, size = 3.5) +
      geom_point(data = hl, colour = "orange", size = 5, shape = 18) +
      geom_segment(data = hl, aes(x = dose, xend = dose, y = 0, yend = E2),
                   linetype = "dotted", colour = "orange") +
      geom_hline(yintercept = c(exp(p$logE2min), exp(p$logE2max)),
                 linetype = "dashed", colour = "grey60") +
      annotate("text", x = 760, y = exp(p$logE2max), label = "E2max", colour = "grey40") +
      annotate("text", x = 760, y = exp(p$logE2min), label = "E2min", colour = "grey40") +
      scale_y_continuous(limits = c(0, 200)) +
      labs(title = "Figure 1 -- Dose-E2 Relationship (Scaled Logistic)",
           x = "Daily Elagolix Dose (mg)", y = "Predicted E2 (pg/mL)") +
      theme_bw(base_size = 14)
  })

  output$vb_baseline <- renderText(sprintf("%.1f pg/mL", predict_E2(0,
    list(slope=input$slope, logE2max=input$logE2max, logE2min=input$logE2min))))
  output$vb_pred_e2 <- renderText({
    e2 <- predict_E2(input$dose_highlight,
      list(slope=input$slope, logE2max=input$logE2max, logE2min=input$logE2min))
    sprintf("%.1f pg/mL", e2)
  })
  output$vb_supp <- renderText({
    p <- list(slope=input$slope, logE2max=input$logE2max, logE2min=input$logE2min)
    bl <- predict_E2(0, p); e2 <- predict_E2(input$dose_highlight, p)
    sprintf("%.1f %%", (1-e2/bl)*100)
  })

  # ========== Tab 2 : QSP Simulation ======================================

  sim_data <- reactiveVal(NULL)

  observeEvent(input$run_btn, {
    showNotification("Simulation running ... please wait.", type = "message",
                     duration = NULL, id = "sim_note")
    dose <- as.numeric(input$dose_sel)
    res <- tryCatch(
      run_sim(dose, input$tx_mo, input$fu_mo),
      error = function(e) { showNotification(paste("Error:", e$message), type="error"); NULL }
    )
    removeNotification("sim_note")
    if (!is.null(res)) {
      sim_data(res)
      showNotification("Simulation complete!", type = "message", duration = 3)
    }
  })

  make_tc_plot <- function(df, yvar, ylab, title_str) {
    if (is.null(df)) return(NULL)
    tx_end <- input$tx_mo
    ggplot(df, aes(months, .data[[yvar]])) +
      geom_line(colour = "steelblue", linewidth = 1) +
      {if (input$fu_mo > 0) geom_vline(xintercept = tx_end,
                                         linetype = "dashed", colour = "grey50")} +
      {if (input$fu_mo > 0)
        annotate("text", x = tx_end + 0.3, y = max(df[[yvar]], na.rm=T)*0.9,
                 label = "Tx end", size = 3.5, hjust = 0)} +
      geom_hline(yintercept = 0, linetype = "dotted") +
      labs(title = title_str, x = "Time (months)", y = ylab) +
      theme_bw(base_size = 14)
  }

  output$plot_bmd  <- renderPlot(make_tc_plot(sim_data(), "BMD_chg_pct",
                                              "BMD Change (%)", "BMD Change Over Time"))
  output$plot_ctx  <- renderPlot(make_tc_plot(sim_data(), "CTX_chg_pct",
                                              "CTX Change (%)", "CTX (Bone Resorption) Over Time"))
  output$plot_p1np <- renderPlot(make_tc_plot(sim_data(), "P1NP_chg_pct",
                                              "P1NP Change (%)", "P1NP (Bone Formation) Over Time"))

  output$plot_all <- renderPlot({
    df <- sim_data(); if (is.null(df)) return(NULL)
    tx_end <- input$tx_mo
    long <- df %>%
      select(months, BMD_chg_pct, CTX_chg_pct, P1NP_chg_pct) %>%
      pivot_longer(-months, names_to = "Marker", values_to = "Change_pct") %>%
      mutate(Marker = recode(Marker,
        BMD_chg_pct = "BMD", CTX_chg_pct = "CTX", P1NP_chg_pct = "P1NP"))
    ggplot(long, aes(months, Change_pct, colour = Marker)) +
      geom_line(linewidth = 1) +
      {if (input$fu_mo > 0) geom_vline(xintercept = tx_end,
                                         linetype = "dashed", colour = "grey50")} +
      geom_hline(yintercept = 0, linetype = "dotted") +
      scale_colour_manual(values = c(BMD = "#E41A1C", CTX = "#377EB8", P1NP = "#4DAF4A")) +
      labs(title = "All Biomarkers", x = "Time (months)", y = "Change from Baseline (%)") +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
  })

  output$tbl_biomarker <- renderDT({
    df <- sim_data(); if (is.null(df)) return(NULL)
    pts <- c(3, 6, 9, 12, 18, 24)
    pts <- pts[pts <= max(df$months)]
    tbl <- data.frame(
      `Month` = pts,
      `BMD (%)` = sapply(pts, function(m) round(at_month(df, m, "BMD_chg_pct"), 3)),
      `CTX (%)` = sapply(pts, function(m) round(at_month(df, m, "CTX_chg_pct"), 3)),
      `P1NP (%)` = sapply(pts, function(m) round(at_month(df, m, "P1NP_chg_pct"), 3)),
      check.names = FALSE
    )
    datatable(tbl, rownames = FALSE,
              options = list(dom = "t", pageLength = 20)) %>%
      formatRound(2:4, 2)
  })

  # ========== Tab 3 : Multi-dose Comparison ================================

  cmp_data <- reactiveVal(NULL)

  observeEvent(input$run_cmp, {
    doses <- as.numeric(input$doses_cmp)
    if (length(doses) == 0) {
      showNotification("Select at least one dose.", type = "warning"); return()
    }
    showNotification("Running comparison ... please wait.", type = "message",
                     duration = NULL, id = "cmp_note")
    all_res <- tryCatch({
      do.call(rbind, lapply(doses, function(d) {
        run_sim(d, input$tx_mo_cmp, input$fu_mo_cmp)
      }))
    }, error = function(e) { showNotification(paste("Error:",e$message), type="error"); NULL })
    removeNotification("cmp_note")
    if (!is.null(all_res)) {
      cmp_data(all_res)
      showNotification("Comparison complete!", type = "message", duration = 3)
    }
  })

  make_cmp_plot <- function(yvar, ylab, title_str) {
    df <- cmp_data(); if (is.null(df)) return(NULL)
    tx_end <- input$tx_mo_cmp
    ggplot(df, aes(months, .data[[yvar]], colour = dose_label)) +
      geom_line(linewidth = 1) +
      {if (input$fu_mo_cmp > 0) geom_vline(xintercept = tx_end,
                                             linetype = "dashed", colour = "grey50")} +
      geom_hline(yintercept = 0, linetype = "dotted") +
      labs(title = title_str, x = "Time (months)", y = ylab, colour = "Dose") +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
  }

  output$cmp_bmd  <- renderPlot(make_cmp_plot("BMD_chg_pct",  "BMD Change (%)",
                                               "BMD -- Multi-dose Comparison"))
  output$cmp_ctx  <- renderPlot(make_cmp_plot("CTX_chg_pct",  "CTX Change (%)",
                                               "CTX -- Multi-dose Comparison"))
  output$cmp_p1np <- renderPlot(make_cmp_plot("P1NP_chg_pct", "P1NP Change (%)",
                                               "P1NP -- Multi-dose Comparison"))

  output$tbl3_cmp <- renderDT({
    df <- cmp_data(); if (is.null(df)) return(NULL)
    pts <- c(6, 12, 18, 24)
    pts <- pts[pts <= max(df$months)]
    doses <- unique(df$dose_label)
    tbl <- do.call(rbind, lapply(doses, function(dl) {
      sub <- df[df$dose_label == dl, ]
      row <- data.frame(Dose = dl)
      for (m in pts) row[[paste0(m, "mo")]] <- round(at_month(sub, m, "BMD_chg_pct"), 2)
      row
    }))
    datatable(tbl, rownames = FALSE, options = list(dom = "t")) %>%
      formatRound(2:ncol(tbl), 2)
  })

  # ========== Tab 4 : Validation ===========================================

  output$tbl_paper_target <- renderDT({
    datatable(paper_target, rownames = FALSE, options = list(dom = "t")) %>%
      formatRound(2:5, 2)
  })

  output$tbl_model_pred <- renderDT({
    df <- cmp_data(); if (is.null(df)) return(NULL)
    pts <- c(6, 12, 18, 24)
    pts <- pts[pts <= max(df$months)]
    doses <- unique(df$dose_label)
    tbl <- do.call(rbind, lapply(doses, function(dl) {
      sub <- df[df$dose_label == dl, ]
      row <- data.frame(Dose = dl)
      for (m in pts) row[[paste0(m, "mo")]] <- round(at_month(sub, m, "BMD_chg_pct"), 2)
      row
    }))
    datatable(tbl, rownames = FALSE, options = list(dom = "t")) %>%
      formatRound(2:ncol(tbl), 2)
  })
}

###############################################################################
#                         Launch
###############################################################################

shinyApp(ui, server)
