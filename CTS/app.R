###############################################################################
# app.R  —  Elagolix QSP Model  ·  Interactive Shiny Dashboard
# -----------------------------------------------------------------------
# Model core: faithful R port of MetrumRG OpenBoneMin.cpp
# Stodtmann et al. (2021) Clin Transl Sci. 14:1611-1619
# 8-tab dashboard with full biomarker visualization
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

# ── Biomarker Registry ───────────────────────────────────────────────────────
BIOMARKERS <- list(
  # Clinical Endpoints (% change)
  list(var="BMD_chg_pct",  label="BMD (Lumbar Spine)", unit="% change", group="clinical",    color="#E41A1C"),
  list(var="CTX_chg_pct",  label="CTX (Resorption)",   unit="% change", group="clinical",    color="#377EB8"),
  list(var="P1NP_chg_pct", label="P1NP (Formation)",   unit="% change", group="clinical",    color="#4DAF4A"),
  # Ca & Hormones
  list(var="PTHconc",      label="PTH",                unit="pg/mL",    group="ca_hormones", color="#E6550D"),
  list(var="Calcitriol",   label="Calcitriol (1,25D)", unit="pg/mL",    group="ca_hormones", color="#FD8D3C"),
  list(var="CaConc",       label="Plasma Ca",          unit="mg/dL",    group="ca_hormones", color="#FDAE6B"),
  list(var="ECCPhos",      label="ECC Phosphate",      unit="mg",       group="ca_hormones", color="#A1D99B"),
  list(var="AOH",          label="1α-Hydroxylase",     unit="AU",       group="ca_hormones", color="#74C476"),
  list(var="S",            label="CaSR Sensor",        unit="AU",       group="ca_hormones", color="#31A354"),
  list(var="Tgut",         label="Gut Ca Pool",        unit="mg",       group="ca_hormones", color="#006D2C"),
  list(var="Rkid",         label="Renal VDR",          unit="AU",       group="ca_hormones", color="#9E9AC8"),
  list(var="HAp",          label="Hydroxyapatite",     unit="AU",       group="ca_hormones", color="#756BB1"),
  list(var="Q",            label="Exchangeable Ca",    unit="mg",       group="ca_hormones", color="#54278F"),
  list(var="Qbone",        label="Bone Ca Pool",       unit="mg",       group="ca_hormones", color="#636363"),
  # Bone Cells
  list(var="ROB1",         label="ROB (Responding OB)",unit="cells",    group="bone_cells",  color="#1B9E77"),
  list(var="OBfast",       label="OB fast (Active)",   unit="cells",    group="bone_cells",  color="#D95F02"),
  list(var="OBslow",       label="OB slow (Lining)",   unit="cells",    group="bone_cells",  color="#7570B3"),
  list(var="OB_total",     label="OB Total",           unit="cells",    group="bone_cells",  color="#E7298A"),
  list(var="OC",           label="Osteoclasts (OC)",   unit="cells",    group="bone_cells",  color="#66A61E"),
  # Signaling Pathways
  list(var="L",            label="RANKL",              unit="pM",       group="signaling",   color="#E41A1C"),
  list(var="RNK",          label="RANK",               unit="pM",       group="signaling",   color="#377EB8"),
  list(var="Oopg",         label="OPG",                unit="pM",       group="signaling",   color="#4DAF4A"),
  list(var="Mcmplx",       label="RANKL-RANK Complex", unit="pM",       group="signaling",   color="#984EA3"),
  list(var="Ncmplx",       label="RANKL-OPG Complex",  unit="pM",       group="signaling",   color="#FF7F00"),
  list(var="TGFB",         label="TGF-β (Latent)",     unit="AU",       group="signaling",   color="#A65628"),
  list(var="TGFBact",      label="TGF-β (Active)",     unit="AU",       group="signaling",   color="#F781BF"),
  list(var="RX2",          label="RUNX2",              unit="AU",       group="signaling",   color="#1F78B4"),
  list(var="CREB",         label="CREB",               unit="AU",       group="signaling",   color="#33A02C"),
  list(var="BCL2",         label="BCL2",               unit="AU",       group="signaling",   color="#FB9A99"),
  # Auxiliary
  list(var="EST",          label="Estrogen (EST)",     unit="ratio",    group="auxiliary",    color="#B15928"),
  list(var="BMDls",        label="BMD (raw)",          unit="ratio",    group="auxiliary",    color="#E41A1C"),
  list(var="UCA",          label="Urinary Ca",         unit="mg",       group="auxiliary",    color="#636363"),
  list(var="PhosGut",      label="Gut Phosphate",      unit="mg",       group="auxiliary",    color="#969696"),
  list(var="IntraPO",      label="Intracellular PO₄",  unit="mg",       group="auxiliary",    color="#BDBDBD"),
  list(var="PTmax",        label="PTH Max Capacity",   unit="AU",       group="auxiliary",    color="#D9D9D9")
)

GROUP_LABELS <- c(
  clinical    = "Clinical Endpoints",
  ca_hormones = "Calcium & Hormones",
  bone_cells  = "Bone Cells",
  signaling   = "Signaling Pathways",
  auxiliary   = "Auxiliary States"
)

GROUP_COLORS <- c(
  clinical    = "#E41A1C",
  ca_hormones = "#E6550D",
  bone_cells  = "#1B9E77",
  signaling   = "#984EA3",
  auxiliary   = "#636363"
)

# ── Custom ggplot2 theme ─────────────────────────────────────────────────────
theme_qsp <- function(base_size = 13) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.1), margin = margin(b = 8)),
      plot.subtitle    = element_text(color = "grey40", size = rel(0.85)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      legend.position  = "bottom",
      axis.title       = element_text(size = rel(0.9)),
      plot.margin      = margin(10, 12, 10, 12)
    )
}

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

# ── Pre-compute Steady State (performance optimization) ──────────────────────
message("Computing steady-state (120-month warm-up) ...")
SS_PARAMS  <- make_params()
SS_Y0      <- make_init(SS_PARAMS)
SS_SOL     <- ode(SS_Y0, seq(0, 120 * HOURS_PER_MONTH, by = 48),
                  qsp_ode, SS_PARAMS,
                  method = "lsoda", atol = 1e-8, rtol = 1e-8, maxsteps = 1e5)
n_state    <- length(SS_Y0)
SS_VEC     <- SS_SOL[nrow(SS_SOL), 2:(n_state + 1)]
names(SS_VEC) <- names(SS_Y0)
SS_OB      <- as.numeric(SS_VEC["OBfast"] + SS_VEC["OBslow"])
SS_OC      <- as.numeric(SS_VEC["OC"])
SS_BMD     <- as.numeric(SS_VEC["BMDls"])
message("Steady-state ready.")

# ── Enhanced run_sim (uses cached SS) ────────────────────────────────────────
run_sim <- function(dose_mg, tx_months, fu_months = 0) {
  p <- SS_PARAMS
  y_ss <- SS_VEC
  tx_hr <- tx_months * HOURS_PER_MONTH
  fu_hr <- fu_months * HOURS_PER_MONTH
  EST_tx <- predict_E2(dose_mg) / E2_baseline

  y_tx <- y_ss; y_tx["EST"] <- EST_tx
  ev_tx <- function(t, y, parms) { y["EST"] <- EST_tx; y }
  times_tx <- seq(0, tx_hr, by = 48)
  sol_tx <- ode(y_tx, times_tx, qsp_ode, p, method = "lsoda",
                atol = 1e-8, rtol = 1e-8, maxsteps = 1e5,
                events = list(func = ev_tx, time = times_tx))
  res <- as.data.frame(sol_tx)

  if (fu_months > 0) {
    y_fu <- sol_tx[nrow(sol_tx), 2:(n_state + 1)]
    names(y_fu) <- names(SS_Y0)
    y_fu["EST"] <- 1
    sol_fu <- ode(y_fu, seq(0, fu_hr, by = 48), qsp_ode, p,
                  method = "lsoda", atol = 1e-8, rtol = 1e-8, maxsteps = 1e5)
    res_fu <- as.data.frame(sol_fu)
    res_fu$time <- res_fu$time + tx_hr
    res <- rbind(res, res_fu[-1, ])
  }

  # Derived columns
  res$months       <- res$time / HOURS_PER_MONTH
  res$OB_total     <- res$OBfast + res$OBslow
  res$CaConc       <- res$P / 14
  res$Calcitriol   <- res$B / 14
  res$PTHconc      <- res$PTH / 14
  res$BMD_chg_pct  <- (res$BMDls / SS_BMD - 1) * 100
  res$CTX_chg_pct  <- (res$OC / SS_OC - 1) * 100
  res$P1NP_chg_pct <- (res$OB_total / SS_OB - 1) * 100

  # % change from SS for all state variables
  for (v in names(SS_VEC)) {
    ss_val <- as.numeric(SS_VEC[v])
    if (abs(ss_val) > 1e-15) {
      res[[paste0(v, "_pct")]] <- (res[[v]] / ss_val - 1) * 100
    }
  }
  res$OB_total_pct <- (res$OB_total / SS_OB - 1) * 100

  res$dose_mg    <- dose_mg
  res$dose_label <- dplyr::case_when(
    dose_mg == 150 ~ "150 mg QD",
    dose_mg == 400 ~ "200 mg BID",
    dose_mg == 600 ~ "300 mg BID",
    TRUE ~ paste0(dose_mg, " mg")
  )
  res
}

at_month <- function(df, m, col) df[[col]][which.min(abs(df$months - m))]

# ── Paper target (Table 3) ──────────────────────────────────────────────────
paper_target <- data.frame(
  Dose = c("150 mg QD", "200 mg BID"),
  `6mo`  = c(-0.61, -3.47),
  `12mo` = c(-0.91, -4.95),
  `18mo` = c(-0.96, -5.15),
  `24mo` = c(-0.91, -4.97),
  check.names = FALSE
)

# ── Plot helper functions ────────────────────────────────────────────────────
make_detail_plot <- function(df, var, label, unit, color, tx_end, fu_months) {
  if (is.null(df) || !var %in% names(df)) return(NULL)
  ggplot(df, aes(months, .data[[var]])) +
    geom_line(colour = color, linewidth = 1.1) +
    {if (fu_months > 0) geom_vline(xintercept = tx_end,
                                    linetype = "dashed", colour = "grey50")} +
    {if (fu_months > 0) annotate("text", x = tx_end, y = Inf, label = " Tx end",
                                  hjust = 0, vjust = 1.5, size = 3.2, colour = "grey50")} +
    labs(title = label, x = "Time (months)", y = unit) +
    theme_qsp(base_size = 12)
}

make_dashboard_plot <- function(df, bm_list, tx_end, fu_months) {
  if (is.null(df)) return(NULL)
  # Build long-format data for faceting
  plot_vars <- sapply(bm_list, `[[`, "var")
  plot_labels <- sapply(bm_list, `[[`, "label")
  avail <- plot_vars[plot_vars %in% names(df)]
  if (length(avail) == 0) return(NULL)

  long <- df %>%
    select(months, all_of(avail)) %>%
    pivot_longer(-months, names_to = "variable", values_to = "value")

  # Map variable names to readable labels
  label_map <- setNames(plot_labels, plot_vars)
  long$facet_label <- label_map[long$variable]
  long$facet_label <- factor(long$facet_label, levels = unique(label_map[avail]))

  ggplot(long, aes(months, value)) +
    geom_line(colour = "#2171B5", linewidth = 0.6) +
    {if (fu_months > 0) geom_vline(xintercept = tx_end,
                                    linetype = "dashed", colour = "grey60", linewidth = 0.4)} +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 5) +
    labs(x = "Time (months)", y = NULL) +
    theme_qsp(base_size = 10) +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      axis.text  = element_text(size = 7),
      panel.spacing = unit(0.4, "lines")
    )
}

make_comparison_plot <- function(df, var, label, unit, tx_end, fu_months) {
  if (is.null(df) || !var %in% names(df)) return(NULL)
  ggplot(df, aes(months, .data[[var]], colour = dose_label)) +
    geom_line(linewidth = 1) +
    {if (fu_months > 0) geom_vline(xintercept = tx_end,
                                    linetype = "dashed", colour = "grey50")} +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
    labs(title = label, x = "Time (months)", y = unit, colour = "Dose") +
    theme_qsp(base_size = 13) +
    theme(legend.position = "bottom")
}

###############################################################################
#                         Shiny UI
###############################################################################

# Biomarker choices for comparison selector
bm_choices <- setNames(

  sapply(BIOMARKERS, `[[`, "var"),
  sapply(BIOMARKERS, `[[`, "label")
)

ui <- page_navbar(
  title = "Elagolix QSP Model",
  theme = bs_theme(
    bootswatch = "flatly",
    base_font = font_google("Noto Sans KR"),
    "navbar-bg" = "#2C3E50"
  ),
  window_title = "Elagolix QSP | Stodtmann 2021",

  # ═══════════ Tab 1: Dose-E2 ═══════════════════════════════════════════════
  nav_panel("Dose \u2192 E2",
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
      card(card_header("Dose-E2 Relationship (Scaled Logistic)"),
           plotOutput("fig1_plot", height = "420px")),
      layout_column_wrap(
        width = 1/3,
        value_box("Baseline E2", textOutput("vb_baseline"), theme = "primary"),
        value_box("Predicted E2", textOutput("vb_pred_e2"), theme = "info"),
        value_box("Suppression", textOutput("vb_supp"), theme = "warning")
      )
    )
  ),

  # ═══════════ Tab 2: Dashboard ═════════════════════════════════════════════
  nav_panel("Dashboard",
    layout_sidebar(
      sidebar = sidebar(
        title = "Simulation Settings", width = 280,
        selectInput("dash_dose", "Elagolix Dose",
                    choices = c("150 mg QD" = 150, "200 mg BID" = 400, "300 mg BID" = 600),
                    selected = 400),
        sliderInput("dash_tx", "Treatment (months)", 1, 36, 24, step = 1),
        sliderInput("dash_fu", "Follow-up (months)", 0, 24, 6, step = 1),
        hr(),
        actionButton("dash_run", "Run Simulation",
                     class = "btn-primary btn-lg w-100", icon = icon("play")),
        hr(),
        helpText("All 31 state variables at a glance.")
      ),
      card(
        card_header(class = "bg-primary text-white",
                    "All Biomarkers Overview"),
        plotOutput("dash_plot", height = "900px")
      )
    )
  ),

  # ═══════════ Tab 3: Clinical Endpoints ════════════════════════════════════
  nav_panel("Clinical Endpoints",
    layout_sidebar(
      sidebar = sidebar(
        title = "Simulation Settings", width = 280,
        selectInput("clin_dose", "Elagolix Dose",
                    choices = c("150 mg QD" = 150, "200 mg BID" = 400, "300 mg BID" = 600),
                    selected = 150),
        sliderInput("clin_tx", "Treatment (months)", 1, 36, 12, step = 1),
        sliderInput("clin_fu", "Follow-up (months)", 0, 24, 6, step = 1),
        hr(),
        actionButton("clin_run", "Run Simulation",
                     class = "btn-primary btn-lg w-100", icon = icon("play")),
        hr(),
        helpText("BMD, CTX, P1NP — clinical bone markers.")
      ),
      layout_column_wrap(
        width = 1/2,
        card(card_header("BMD Change (%)"),   plotOutput("clin_bmd",  height = "320px")),
        card(card_header("CTX Change (%)"),   plotOutput("clin_ctx",  height = "320px")),
        card(card_header("P1NP Change (%)"),  plotOutput("clin_p1np", height = "320px")),
        card(card_header("All Clinical Markers"),
             plotOutput("clin_all", height = "320px"))
      ),
      card(card_header("Biomarker Values at Selected Timepoints"),
           DTOutput("clin_table"))
    )
  ),

  # ═══════════ Tab 4: Ca & Hormones ════════════════════════════════════════
  nav_panel("Ca & Hormones",
    layout_sidebar(
      sidebar = sidebar(
        title = "Simulation Settings", width = 280,
        selectInput("ca_dose", "Elagolix Dose",
                    choices = c("150 mg QD" = 150, "200 mg BID" = 400, "300 mg BID" = 600),
                    selected = 400),
        sliderInput("ca_tx", "Treatment (months)", 1, 36, 24, step = 1),
        sliderInput("ca_fu", "Follow-up (months)", 0, 24, 6, step = 1),
        hr(),
        actionButton("ca_run", "Run Simulation",
                     class = "btn-warning btn-lg w-100", icon = icon("play")),
        hr(),
        helpText("PTH, Calcitriol, Ca, Phosphate, and related states.")
      ),
      layout_column_wrap(
        width = 1/2,
        card(card_header("PTH Concentration"),        plotOutput("ca_pth",     height = "300px")),
        card(card_header("Calcitriol (1,25D)"),       plotOutput("ca_ctriol",  height = "300px")),
        card(card_header("Plasma Calcium"),            plotOutput("ca_plasma",  height = "300px")),
        card(card_header("ECC Phosphate"),             plotOutput("ca_phos",    height = "300px")),
        card(card_header("1\u03b1-Hydroxylase (AOH)"),plotOutput("ca_aoh",     height = "300px")),
        card(card_header("CaSR Sensor (S)"),           plotOutput("ca_sensor",  height = "300px")),
        card(card_header("Gut Ca Pool (Tgut)"),        plotOutput("ca_tgut",    height = "300px")),
        card(card_header("Renal VDR (Rkid)"),          plotOutput("ca_rkid",    height = "300px")),
        card(card_header("Hydroxyapatite"),            plotOutput("ca_hap",     height = "300px")),
        card(card_header("Exchangeable Ca (Q)"),       plotOutput("ca_q",       height = "300px")),
        card(card_header("Bone Ca Pool (Qbone)"),      plotOutput("ca_qbone",   height = "300px")),
        card(card_header("Urinary Ca (UCA)"),          plotOutput("ca_uca",     height = "300px"))
      )
    )
  ),

  # ═══════════ Tab 5: Bone Cells ═══════════════════════════════════════════
  nav_panel("Bone Cells",
    layout_sidebar(
      sidebar = sidebar(
        title = "Simulation Settings", width = 280,
        selectInput("bone_dose", "Elagolix Dose",
                    choices = c("150 mg QD" = 150, "200 mg BID" = 400, "300 mg BID" = 600),
                    selected = 400),
        sliderInput("bone_tx", "Treatment (months)", 1, 36, 24, step = 1),
        sliderInput("bone_fu", "Follow-up (months)", 0, 24, 6, step = 1),
        hr(),
        actionButton("bone_run", "Run Simulation",
                     class = "btn-success btn-lg w-100", icon = icon("play")),
        hr(),
        helpText("ROB, OBfast, OBslow, OB total, OC.")
      ),
      layout_column_wrap(
        width = 1/2,
        card(card_header("ROB (Responding OB)"), plotOutput("bone_rob",    height = "320px")),
        card(card_header("OB fast (Active)"),    plotOutput("bone_obfast", height = "320px")),
        card(card_header("OB slow (Lining)"),    plotOutput("bone_obslow", height = "320px")),
        card(card_header("OB Total"),            plotOutput("bone_obtot",  height = "320px")),
        card(card_header("Osteoclasts (OC)"),    plotOutput("bone_oc",     height = "320px"))
      )
    )
  ),

  # ═══════════ Tab 6: Signaling Pathways ═══════════════════════════════════
  nav_panel("Signaling",
    layout_sidebar(
      sidebar = sidebar(
        title = "Simulation Settings", width = 280,
        selectInput("sig_dose", "Elagolix Dose",
                    choices = c("150 mg QD" = 150, "200 mg BID" = 400, "300 mg BID" = 600),
                    selected = 400),
        sliderInput("sig_tx", "Treatment (months)", 1, 36, 24, step = 1),
        sliderInput("sig_fu", "Follow-up (months)", 0, 24, 6, step = 1),
        hr(),
        actionButton("sig_run", "Run Simulation",
                     class = "btn-info btn-lg w-100", icon = icon("play")),
        hr(),
        helpText("RANKL/RANK/OPG, TGF-\u03b2, RUNX2, CREB, BCL2.")
      ),
      layout_column_wrap(
        width = 1/2,
        card(card_header("RANKL"),              plotOutput("sig_rankl",   height = "300px")),
        card(card_header("RANK"),               plotOutput("sig_rank",    height = "300px")),
        card(card_header("OPG"),                plotOutput("sig_opg",     height = "300px")),
        card(card_header("RANKL-RANK Complex"), plotOutput("sig_mcmplx", height = "300px")),
        card(card_header("RANKL-OPG Complex"),  plotOutput("sig_ncmplx", height = "300px")),
        card(card_header("TGF-\u03b2 (Latent)"),plotOutput("sig_tgfb",   height = "300px")),
        card(card_header("TGF-\u03b2 (Active)"),plotOutput("sig_tgfba",  height = "300px")),
        card(card_header("RUNX2"),              plotOutput("sig_runx2",   height = "300px")),
        card(card_header("CREB"),               plotOutput("sig_creb",    height = "300px")),
        card(card_header("BCL2"),               plotOutput("sig_bcl2",    height = "300px"))
      )
    )
  ),

  # ═══════════ Tab 7: Multi-Dose Comparison ════════════════════════════════
  nav_panel("Multi-Dose",
    layout_sidebar(
      sidebar = sidebar(
        title = "Comparison Settings", width = 300,
        checkboxGroupInput("cmp_doses", "Select Doses",
                           choices = c("150 mg QD" = 150, "200 mg BID" = 400, "300 mg BID" = 600),
                           selected = c(150, 400)),
        sliderInput("cmp_tx", "Treatment (months)", 1, 36, 24, step = 1),
        sliderInput("cmp_fu", "Follow-up (months)", 0, 24, 0, step = 1),
        hr(),
        selectInput("cmp_var", "Biomarker to Compare", choices = bm_choices,
                    selected = "BMD_chg_pct"),
        hr(),
        actionButton("cmp_run", "Run Comparison",
                     class = "btn-success btn-lg w-100", icon = icon("layer-group")),
        hr(),
        helpText("Compare multiple regimens side-by-side for any biomarker.")
      ),
      card(card_header("Multi-Dose Comparison"),
           plotOutput("cmp_plot", height = "450px")),
      layout_column_wrap(
        width = 1/3,
        card(card_header("BMD Comparison"),  plotOutput("cmp_bmd",  height = "300px")),
        card(card_header("CTX Comparison"),  plotOutput("cmp_ctx",  height = "300px")),
        card(card_header("P1NP Comparison"), plotOutput("cmp_p1np", height = "300px"))
      ),
      card(card_header("BMD Change (%) at Key Timepoints"),
           DTOutput("cmp_table"))
    )
  ),

  # ═══════════ Tab 8: Validation ═══════════════════════════════════════════
  nav_panel("Validation",
    layout_column_wrap(
      width = 1/2,
      card(card_header(class = "bg-primary text-white", "Paper Target (Table 3)"),
           DTOutput("val_paper")),
      card(card_header(class = "bg-success text-white", "Model Prediction"),
           DTOutput("val_model"))
    ),
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
      tags$h5("Module A \u2014 Dose-E2 (Scaled Logistic)"),
      tags$pre("E2 = exp(logE2min) + (exp(logE2max) - exp(logE2min)) / (1 + exp(slope * Dose))"),
      tags$h5("Module B \u2014 Calcium Homeostasis & Bone Remodeling QSP"),
      tags$p("31 ODEs covering: PTH, calcitriol, Ca/Phos homeostasis,
              ROB, OBfast, OBslow, OC, RANKL/RANK/OPG, TGF-\u03b2,
              RUNX2, CREB, BCL2, estrogen, BMD (lumbar spine)."),
      tags$p("All equations faithfully ported from OpenBoneMin.cpp source.")
    )
  )
)

###############################################################################
#                         Shiny Server
###############################################################################

server <- function(input, output, session) {

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 1: Dose-E2
  # ══════════════════════════════════════════════════════════════════════════

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
      geom_line(colour = "#2C3E50", linewidth = 1.3) +
      geom_point(data = key, colour = "#E74C3C", size = 3.5) +
      geom_text(data = key, aes(label = label), vjust = -1.3, size = 3.5) +
      geom_point(data = hl, colour = "#F39C12", size = 5, shape = 18) +
      geom_segment(data = hl, aes(x = dose, xend = dose, y = 0, yend = E2),
                   linetype = "dotted", colour = "#F39C12") +
      geom_hline(yintercept = c(exp(p$logE2min), exp(p$logE2max)),
                 linetype = "dashed", colour = "grey60") +
      annotate("text", x = 760, y = exp(p$logE2max), label = "E2max", colour = "grey40") +
      annotate("text", x = 760, y = exp(p$logE2min), label = "E2min", colour = "grey40") +
      scale_y_continuous(limits = c(0, 200)) +
      labs(title = "Dose-E2 Relationship (Scaled Logistic, Eq. 3)",
           subtitle = "Stodtmann et al. (2021)",
           x = "Daily Elagolix Dose (mg)", y = "Predicted E2 (pg/mL)") +
      theme_qsp(base_size = 14)
  })

  output$vb_baseline <- renderText(sprintf("%.1f pg/mL", predict_E2(0,
    list(slope = input$slope, logE2max = input$logE2max, logE2min = input$logE2min))))
  output$vb_pred_e2 <- renderText({
    e2 <- predict_E2(input$dose_highlight,
      list(slope = input$slope, logE2max = input$logE2max, logE2min = input$logE2min))
    sprintf("%.1f pg/mL", e2)
  })
  output$vb_supp <- renderText({
    p <- list(slope = input$slope, logE2max = input$logE2max, logE2min = input$logE2min)
    bl <- predict_E2(0, p); e2 <- predict_E2(input$dose_highlight, p)
    sprintf("%.1f %%", (1 - e2 / bl) * 100)
  })

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 2: Dashboard
  # ══════════════════════════════════════════════════════════════════════════

  dash_data <- eventReactive(input$dash_run, {
    showNotification("Running simulation ...", type = "message", duration = NULL, id = "dash_note")
    res <- tryCatch(
      run_sim(as.numeric(input$dash_dose), input$dash_tx, input$dash_fu),
      error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL }
    )
    removeNotification("dash_note")
    if (!is.null(res)) showNotification("Done!", type = "message", duration = 2)
    res
  })

  output$dash_plot <- renderPlot({
    make_dashboard_plot(dash_data(), BIOMARKERS, input$dash_tx, input$dash_fu)
  })

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 3: Clinical Endpoints
  # ══════════════════════════════════════════════════════════════════════════

  clin_data <- eventReactive(input$clin_run, {
    showNotification("Running simulation ...", type = "message", duration = NULL, id = "clin_note")
    res <- tryCatch(
      run_sim(as.numeric(input$clin_dose), input$clin_tx, input$clin_fu),
      error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL }
    )
    removeNotification("clin_note")
    if (!is.null(res)) showNotification("Done!", type = "message", duration = 2)
    res
  })

  output$clin_bmd <- renderPlot({
    make_detail_plot(clin_data(), "BMD_chg_pct", "BMD (Lumbar Spine)", "% change",
                     "#E41A1C", input$clin_tx, input$clin_fu)
  })
  output$clin_ctx <- renderPlot({
    make_detail_plot(clin_data(), "CTX_chg_pct", "CTX (Bone Resorption)", "% change",
                     "#377EB8", input$clin_tx, input$clin_fu)
  })
  output$clin_p1np <- renderPlot({
    make_detail_plot(clin_data(), "P1NP_chg_pct", "P1NP (Bone Formation)", "% change",
                     "#4DAF4A", input$clin_tx, input$clin_fu)
  })
  output$clin_all <- renderPlot({
    df <- clin_data(); if (is.null(df)) return(NULL)
    long <- df %>%
      select(months, BMD_chg_pct, CTX_chg_pct, P1NP_chg_pct) %>%
      pivot_longer(-months, names_to = "Marker", values_to = "Change") %>%
      mutate(Marker = recode(Marker,
        BMD_chg_pct = "BMD", CTX_chg_pct = "CTX", P1NP_chg_pct = "P1NP"))
    ggplot(long, aes(months, Change, colour = Marker)) +
      geom_line(linewidth = 1) +
      {if (input$clin_fu > 0) geom_vline(xintercept = input$clin_tx,
                                           linetype = "dashed", colour = "grey50")} +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
      scale_colour_manual(values = c(BMD = "#E41A1C", CTX = "#377EB8", P1NP = "#4DAF4A")) +
      labs(title = "All Clinical Markers", x = "Time (months)",
           y = "Change from Baseline (%)", colour = NULL) +
      theme_qsp(base_size = 12)
  })

  output$clin_table <- renderDT({
    df <- clin_data(); if (is.null(df)) return(NULL)
    pts <- c(3, 6, 9, 12, 18, 24)
    pts <- pts[pts <= max(df$months)]
    tbl <- data.frame(
      Month    = pts,
      `BMD (%)`  = sapply(pts, function(m) round(at_month(df, m, "BMD_chg_pct"), 3)),
      `CTX (%)`  = sapply(pts, function(m) round(at_month(df, m, "CTX_chg_pct"), 3)),
      `P1NP (%)` = sapply(pts, function(m) round(at_month(df, m, "P1NP_chg_pct"), 3)),
      check.names = FALSE
    )
    datatable(tbl, rownames = FALSE, options = list(dom = "t", pageLength = 20)) %>%
      formatRound(2:4, 2)
  })

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 4: Ca & Hormones
  # ══════════════════════════════════════════════════════════════════════════

  ca_data <- eventReactive(input$ca_run, {
    showNotification("Running simulation ...", type = "message", duration = NULL, id = "ca_note")
    res <- tryCatch(
      run_sim(as.numeric(input$ca_dose), input$ca_tx, input$ca_fu),
      error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL }
    )
    removeNotification("ca_note")
    if (!is.null(res)) showNotification("Done!", type = "message", duration = 2)
    res
  })

  output$ca_pth    <- renderPlot(make_detail_plot(ca_data(), "PTHconc",    "PTH",              "pg/mL", "#E6550D", input$ca_tx, input$ca_fu))
  output$ca_ctriol <- renderPlot(make_detail_plot(ca_data(), "Calcitriol", "Calcitriol (1,25D)","pg/mL","#FD8D3C", input$ca_tx, input$ca_fu))
  output$ca_plasma <- renderPlot(make_detail_plot(ca_data(), "CaConc",     "Plasma Ca",         "mg/dL","#FDAE6B", input$ca_tx, input$ca_fu))
  output$ca_phos   <- renderPlot(make_detail_plot(ca_data(), "ECCPhos",    "ECC Phosphate",     "mg",   "#A1D99B", input$ca_tx, input$ca_fu))
  output$ca_aoh    <- renderPlot(make_detail_plot(ca_data(), "AOH",        "1\u03b1-Hydroxylase","AU",  "#74C476", input$ca_tx, input$ca_fu))
  output$ca_sensor <- renderPlot(make_detail_plot(ca_data(), "S",          "CaSR Sensor",       "AU",   "#31A354", input$ca_tx, input$ca_fu))
  output$ca_tgut   <- renderPlot(make_detail_plot(ca_data(), "Tgut",       "Gut Ca Pool",       "mg",   "#006D2C", input$ca_tx, input$ca_fu))
  output$ca_rkid   <- renderPlot(make_detail_plot(ca_data(), "Rkid",       "Renal VDR",         "AU",   "#9E9AC8", input$ca_tx, input$ca_fu))
  output$ca_hap    <- renderPlot(make_detail_plot(ca_data(), "HAp",        "Hydroxyapatite",    "AU",   "#756BB1", input$ca_tx, input$ca_fu))
  output$ca_q      <- renderPlot(make_detail_plot(ca_data(), "Q",          "Exchangeable Ca",   "mg",   "#54278F", input$ca_tx, input$ca_fu))
  output$ca_qbone  <- renderPlot(make_detail_plot(ca_data(), "Qbone",      "Bone Ca Pool",      "mg",   "#636363", input$ca_tx, input$ca_fu))
  output$ca_uca    <- renderPlot(make_detail_plot(ca_data(), "UCA",        "Urinary Ca",        "mg",   "#969696", input$ca_tx, input$ca_fu))

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 5: Bone Cells
  # ══════════════════════════════════════════════════════════════════════════

  bone_data <- eventReactive(input$bone_run, {
    showNotification("Running simulation ...", type = "message", duration = NULL, id = "bone_note")
    res <- tryCatch(
      run_sim(as.numeric(input$bone_dose), input$bone_tx, input$bone_fu),
      error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL }
    )
    removeNotification("bone_note")
    if (!is.null(res)) showNotification("Done!", type = "message", duration = 2)
    res
  })

  output$bone_rob    <- renderPlot(make_detail_plot(bone_data(), "ROB1",     "ROB (Responding OB)","cells","#1B9E77", input$bone_tx, input$bone_fu))
  output$bone_obfast <- renderPlot(make_detail_plot(bone_data(), "OBfast",   "OB fast (Active)",   "cells","#D95F02", input$bone_tx, input$bone_fu))
  output$bone_obslow <- renderPlot(make_detail_plot(bone_data(), "OBslow",   "OB slow (Lining)",   "cells","#7570B3", input$bone_tx, input$bone_fu))
  output$bone_obtot  <- renderPlot(make_detail_plot(bone_data(), "OB_total", "OB Total",           "cells","#E7298A", input$bone_tx, input$bone_fu))
  output$bone_oc     <- renderPlot(make_detail_plot(bone_data(), "OC",       "Osteoclasts (OC)",   "cells","#66A61E", input$bone_tx, input$bone_fu))

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 6: Signaling Pathways
  # ══════════════════════════════════════════════════════════════════════════

  sig_data <- eventReactive(input$sig_run, {
    showNotification("Running simulation ...", type = "message", duration = NULL, id = "sig_note")
    res <- tryCatch(
      run_sim(as.numeric(input$sig_dose), input$sig_tx, input$sig_fu),
      error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL }
    )
    removeNotification("sig_note")
    if (!is.null(res)) showNotification("Done!", type = "message", duration = 2)
    res
  })

  output$sig_rankl  <- renderPlot(make_detail_plot(sig_data(), "L",       "RANKL",              "pM", "#E41A1C", input$sig_tx, input$sig_fu))
  output$sig_rank   <- renderPlot(make_detail_plot(sig_data(), "RNK",     "RANK",               "pM", "#377EB8", input$sig_tx, input$sig_fu))
  output$sig_opg    <- renderPlot(make_detail_plot(sig_data(), "Oopg",    "OPG",                "pM", "#4DAF4A", input$sig_tx, input$sig_fu))
  output$sig_mcmplx <- renderPlot(make_detail_plot(sig_data(), "Mcmplx",  "RANKL-RANK Complex", "pM", "#984EA3", input$sig_tx, input$sig_fu))
  output$sig_ncmplx <- renderPlot(make_detail_plot(sig_data(), "Ncmplx",  "RANKL-OPG Complex",  "pM", "#FF7F00", input$sig_tx, input$sig_fu))
  output$sig_tgfb   <- renderPlot(make_detail_plot(sig_data(), "TGFB",    "TGF-\u03b2 (Latent)","AU", "#A65628", input$sig_tx, input$sig_fu))
  output$sig_tgfba  <- renderPlot(make_detail_plot(sig_data(), "TGFBact", "TGF-\u03b2 (Active)","AU", "#F781BF", input$sig_tx, input$sig_fu))
  output$sig_runx2  <- renderPlot(make_detail_plot(sig_data(), "RX2",     "RUNX2",              "AU", "#1F78B4", input$sig_tx, input$sig_fu))
  output$sig_creb   <- renderPlot(make_detail_plot(sig_data(), "CREB",    "CREB",               "AU", "#33A02C", input$sig_tx, input$sig_fu))
  output$sig_bcl2   <- renderPlot(make_detail_plot(sig_data(), "BCL2",    "BCL2",               "AU", "#FB9A99", input$sig_tx, input$sig_fu))

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 7: Multi-Dose Comparison
  # ══════════════════════════════════════════════════════════════════════════

  cmp_data <- eventReactive(input$cmp_run, {
    doses <- as.numeric(input$cmp_doses)
    if (length(doses) == 0) {
      showNotification("Select at least one dose.", type = "warning"); return(NULL)
    }
    showNotification("Running comparison ...", type = "message", duration = NULL, id = "cmp_note")
    all_res <- tryCatch({
      do.call(rbind, lapply(doses, function(d) run_sim(d, input$cmp_tx, input$cmp_fu)))
    }, error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL })
    removeNotification("cmp_note")
    if (!is.null(all_res)) showNotification("Done!", type = "message", duration = 2)
    all_res
  })

  # Selected biomarker comparison
  output$cmp_plot <- renderPlot({
    df <- cmp_data(); if (is.null(df)) return(NULL)
    var <- input$cmp_var
    bm <- Filter(function(b) b$var == var, BIOMARKERS)
    label <- if (length(bm) > 0) bm[[1]]$label else var
    unit  <- if (length(bm) > 0) bm[[1]]$unit  else ""
    make_comparison_plot(df, var, label, unit, input$cmp_tx, input$cmp_fu)
  })

  # Fixed clinical comparisons
  output$cmp_bmd  <- renderPlot(make_comparison_plot(cmp_data(), "BMD_chg_pct",  "BMD",  "% change", input$cmp_tx, input$cmp_fu))
  output$cmp_ctx  <- renderPlot(make_comparison_plot(cmp_data(), "CTX_chg_pct",  "CTX",  "% change", input$cmp_tx, input$cmp_fu))
  output$cmp_p1np <- renderPlot(make_comparison_plot(cmp_data(), "P1NP_chg_pct", "P1NP", "% change", input$cmp_tx, input$cmp_fu))

  output$cmp_table <- renderDT({
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

  # ══════════════════════════════════════════════════════════════════════════
  #  Tab 8: Validation
  # ══════════════════════════════════════════════════════════════════════════

  # Auto-run validation on startup
  val_data <- reactive({
    tryCatch({
      do.call(rbind, lapply(c(150, 400), function(d) run_sim(d, 24, 0)))
    }, error = function(e) NULL)
  })

  output$val_paper <- renderDT({
    datatable(paper_target, rownames = FALSE, options = list(dom = "t")) %>%
      formatRound(2:5, 2)
  })

  output$val_model <- renderDT({
    df <- val_data(); if (is.null(df)) return(NULL)
    pts <- c(6, 12, 18, 24)
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
