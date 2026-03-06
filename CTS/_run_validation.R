###############################################################################
# _run_validation.R
# Elagolix Calcium Homeostasis & Bone Remodeling QSP Model
# -----------------------------------------------------------------------
# Stodtmann et al. (2021) Clin Transl Sci. 14:1611-1619
# Peterson & Riggs (2010) + Riggs et al. (2012) estrogen extension
# Ref implementation: MetrumRG OpenBoneMin
###############################################################################

# === 패키지 ================================================================
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# === 상수 ==================================================================
HOURS_PER_MONTH <- 730.5

###############################################################################
#
#  1.  모듈 A — Dose → E2  (Scaled Logistic, Eq. 3)
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
#  2.  QSP 파라미터 (~90 개)
#
###############################################################################

make_params <- function() {
  list(
    # -- 골세포 역학 --
    OBtot0      = 0.00501324,
    k1          = 6.24e-6,       k2 = 0.112013,
    k3          = 6.24e-6,       k4 = 0.112013,
    kO          = 15.8885,
    kb          = 0.000605516,
    Pic0        = 0.228142,
    FracOBfast  = 0.797629,
    Frackb      = 0.313186,

    # -- TGF-β --
    OBtgfGAM   = 0.0111319,
    koutTGF0    = 2.98449e-5,
    koutTGFGam  = 0.919131,
    OCtgfGAM   = 0.593891,

    # -- 파골세포 --
    kinOCgam    = 8.53065,
    EmaxPicOC   = 1.9746,
    FracPicOC   = 0.878215,
    PicOCgam    = 1.0168,
    MOCratioGam = 0.603754,
    LsurvOCgam  = 3.0923,
    LsurvOCCgam = 3.09023,

    # -- 조골세포 --
    EmaxPicROB  = 3.9745,
    PicROBgam   = 1.80968,
    FracPicROB  = 0.883824,
    PicOBgam    = 0.122313,
    FracPicOB   = 0.000244818,
    EmaxPicOB   = 0.251636,

    # -- Ca/Phos 항상성 --
    CaDay       = 88.0,
    OralCa      = 24.055 / 24,
    OralPhos    = 10.5 / 24,
    F12         = 0.7,
    V1          = 14.0,
    Da          = 0.7 / 24,
    Reabs50     = 1.57322,

    # -- PTH --
    kout_PTH    = 100 / 14,
    PTout       = 1.604e-4,
    opgPTH50    = 3.85,

    # -- 칼시트리올 / 비타민 D --
    AlphOHgam   = 0.111241,
    CtriolPTgam = 12.5033,
    CtriolMax   = 4.1029,
    CtriolMin   = 0.9,

    # -- RANKL / OPG --
    E0RANKL     = 3.80338,
    EmaxL       = 0.469779,
    koutL       = 0.00293273,
    kinRNKgam   = 0.151825,
    koutRNK     = 0.00323667,

    # -- 세포내 신호전달 (RUNX2, CREB, BCL2) --
    RUNX20      = 10.0,
    RX2Kout0    = 0.693,
    E0rx2Kout   = 0.125,
    EmaxPTHRX2x = 5.0,
    E0crebKin   = 0.5,
    EmaxPTHcreb = 3.39745,
    crebKout    = 0.00279513,
    bcl2Kout    = 0.693,

    # -- BMD --
    koutBMDls   = 0.000397,
    gamOB       = 0.0793,
    gamOCls     = 0.14,

    # -- 에스트로겐 확장 (Riggs 2012) --
    ESTON       = 1.0,
    koutEST     = 0.05776227,
    tgfbGAM     = 0.0374,
    tgfbactGAM  = 0.045273,
    robGAM      = 0.16,
    obGAM       = 0.000012,
    E2scalePicB1 = 1.16832e-5,
    maxTmESTkid = 0.923737,

    # -- 보조 --
    Smax = 8.0,  Smin = 0.2,  gamS = 3.0,
    GFR  = 6.0,
    TPmaxCA  = 2.15,
    CaPset   = 32.90,
    CaFrac   = 0.5,
    T85_base = 0.25,
    kinB0    = 126.0,
    T69_base = 0.1,
    T64_base = 1.0,
    PhosKid     = 0.8,
    IntraPO_kel = 0.01,
    IntraPO_kin = 32.26,
    HApForm     = 0.001,
    HApResorb   = 0.001
  )
}

###############################################################################
#
#  3.  초기 조건 (31 state variables)
#
###############################################################################

make_init <- function(p) {
  L0   <- 0.4;  RNK0 <- 10.0;  O0 <- 4.0
  M0   <- p$k3 * RNK0 * L0 / p$k4
  N0   <- p$k1 * O0   * L0 / p$k2

  c(
    PTH     = 53.90,
    S       = 0.5,
    PTmax   = 1.0,
    B       = 1260.0,
    AOH     = 126.0,
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
#  4.  ODE 시스템  (deSolve 호환)
#
###############################################################################

qsp_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {

    # --- 기저 참조값 ---
    OB0   <- OBtot0
    OC0   <- 0.00115398
    ROB0  <- 0.00104122
    L0    <- 0.4;  RNK0 <- 10.0;  O0 <- 4.0
    PTH0  <- 53.90;  B0 <- 1260.0;  P0 <- 32.90
    M0    <- k3 * RNK0 * L0 / k4
    RX20v <- RUNX20

    OB <- OBfast + OBslow
    EST_s <- max(EST, 1e-6)

    # =================================================================
    # (A) 에스트로겐
    # =================================================================
    kinEST <- koutEST          # SS: EST = 1
    dEST   <- (kinEST - koutEST * EST) * ESTON

    # =================================================================
    # (B) PTH
    # =================================================================
    CaRatio <- P / P0
    Sfun    <- Smin + (Smax - Smin) / (1 + CaRatio^gamS)
    SPTH    <- Sfun * kout_PTH * PTH0
    dPTH    <- SPTH - kout_PTH * PTH

    dS     <- 0
    dPTmax <- 0

    PTHr <- PTH / PTH0

    # =================================================================
    # (C) 칼시트리올 / 비타민 D
    # =================================================================
    SE   <- kinB0 * PTHr^AlphOHgam
    dAOH <- SE - T64_base * AOH

    dB   <- AOH - T69_base * B

    Br <- B / B0

    # =================================================================
    # (D) 칼슘 항상성
    # =================================================================
    T85  <- T85_base * (CtriolMin + (CtriolMax - CtriolMin) *
              Br^CtriolPTgam / (1 + Br^CtriolPTgam))
    J40  <- Da * Tgut
    dTgut <- OralCa * T85 - J40

    # 신장 재흡수 — 에스트로겐 효과
    ESTkid <- max(1.0 + maxTmESTkid * (EST_s - 1.0) * ESTON, 0.1)
    TPeff  <- TPmaxCA * ESTkid
    CaFilt <- GFR * CaFrac * P / V1
    T36    <- CaFilt * TPeff / (Reabs50 + CaFrac * P / V1)
    dRkid  <- T36 * (1 - Rkid) - CaFilt * Rkid * 0.01

    Jurine <- CaFilt * (1 - Rkid)

    J14  <- 0.02 * Q
    J15  <- 0.02 * P * 3
    J14a <- HApResorb * (OC / OC0) * Q
    J15a <- HApForm   * (OB / OB0) * P

    dP     <- (J14 - J15 - Jurine + J40 + J14a - J15a) / V1
    dQ     <- J15 - J14 + J14a - J15a
    dQbone <- J15a - J14a
    dUCA   <- Jurine
    dHAp   <- HApForm * (OB / OB0) - HApResorb * (OC / OC0)

    # =================================================================
    # (E) 인산 항상성
    # =================================================================
    J53 <- Da * PhosGut
    J41 <- OralPhos * F12
    J42 <- PhosKid * ECCPhos
    J48 <- 0.001 * ECCPhos
    J54 <- IntraPO_kin * ECCPhos / (ECCPhos + 10)
    J56 <- IntraPO_kel * IntraPO

    dECCPhos <- J41 + J53 - J42 - J48 - J54 + J56
    dPhosGut <- OralPhos * F12 - J53
    dIntraPO <- J54 - J56

    # =================================================================
    # (F) RUNX2, CREB, BCL2
    # =================================================================
    ptf  <- PTHr / (1 + PTHr)
    RX2Kin  <- RX2Kout0 * RX20v * (1 + E0rx2Kout + EmaxPTHRX2x * ptf)
    RX2Kout <- RX2Kout0 * (1 + E0rx2Kout + EmaxPTHRX2x * ptf)
    dRX2    <- RX2Kin - RX2Kout * RX2

    crebKin <- crebKout * 10.0 * (E0crebKin + EmaxPTHcreb * ptf)
    dCREB   <- crebKin - crebKout * CREB

    dBCL2 <- bcl2Kout * (CREB / 10.0) * (RX2 / RX20v) * 100.0 -
             bcl2Kout * BCL2

    # =================================================================
    # (G) TGF-β
    # =================================================================
    kinTGF  <- koutTGF0 * Pic0 * 1000     # SS production
    koutTGFeff <- koutTGF0 * (OC / OC0)^koutTGFGam

    dTGFB <- kinTGF * (OB / OB0)^OBtgfGAM *
               ((1 / EST_s)^tgfbGAM * ESTON + (1 - ESTON)) -
             koutTGFeff * TGFB *
               (EST_s^tgfbactGAM * ESTON + (1 - ESTON))

    koutTGFact <- koutTGF0 * 1000
    dTGFBact <- koutTGFeff * TGFB *
                  (EST_s^tgfbactGAM * ESTON + (1 - ESTON)) -
                koutTGFact * TGFBact

    PicR <- TGFBact / Pic0

    # =================================================================
    # (H) RANKL / RANK / OPG
    # =================================================================
    kinL <- koutL * L0 * (E0RANKL + EmaxL * ptf) / (E0RANKL + EmaxL * 0.5)
    dL   <- kinL - koutL * L - k1 * Oopg * L + k2 * Ncmplx -
            k3 * RNK * L + k4 * Mcmplx

    kinRNK <- koutRNK * RNK0 * (TGFBact / Pic0)^kinRNKgam
    dRNK   <- kinRNK - koutRNK * RNK - k3 * RNK * L + k4 * Mcmplx

    PTHinhOPG  <- opgPTH50 / (opgPTH50 + PTH)
    PTHinhOPG0 <- opgPTH50 / (opgPTH50 + PTH0)
    pO <- kO * O0 * PTHinhOPG / PTHinhOPG0 *
          (EST_s^0.05 * ESTON + (1 - ESTON))
    dOopg <- pO - k1 * Oopg * L + k2 * Ncmplx - kO * Oopg

    dMcmplx <- k3 * RNK * L - k4 * Mcmplx
    dNcmplx <- k1 * Oopg * L - k2 * Ncmplx

    # =================================================================
    # (I) ROB — Responding Osteoblasts
    # =================================================================
    PicROB  <- 1 + EmaxPicROB * FracPicROB * PicR^PicROBgam /
                   (1 + FracPicROB * PicR^PicROBgam)
    PicROB0 <- 1 + EmaxPicROB * FracPicROB /
                   (1 + FracPicROB)
    ROBin   <- kb * ROB0 * PicROB / PicROB0
    ESTrob  <- (1 / EST_s)^robGAM * ESTON + (1 - ESTON)
    dROB1   <- ROBin * ESTrob - kb * ROB1

    # =================================================================
    # (J) OBfast / OBslow
    # =================================================================
    PicOB  <- 1 - EmaxPicOB * FracPicOB * PicR^PicOBgam /
                  (1 + FracPicOB * PicR^PicOBgam)
    PicOB0 <- 1 - EmaxPicOB * FracPicOB /
                  (1 + FracPicOB)
    PicOB  <- max(PicOB, 1e-4)
    PicOB0 <- max(PicOB0, 1e-4)

    Drate  <- kb * ROB1
    bigDb  <- Drate * PicOB / PicOB0
    ESTob  <- EST_s^obGAM * ESTON + (1 - ESTON)

    kbfast <- kb * 10
    kbslow <- kb * 0.5
    dOBfast <- bigDb * FracOBfast * Frackb * ESTob - kbfast * OBfast
    dOBslow <- bigDb * (1 - FracOBfast) * Frackb - kbslow * OBslow

    # =================================================================
    # (K) 파골세포
    # =================================================================
    Mr <- max(Mcmplx / M0, 1e-6)
    Lr <- max(L / L0, 1e-6)

    PicOC  <- 1 + EmaxPicOC * FracPicOC * PicR^PicOCgam /
                  (1 + FracPicOC * PicR^PicOCgam)
    PicOC0 <- 1 + EmaxPicOC * FracPicOC /
                  (1 + FracPicOC)

    kinOC2 <- kb * OC0 * Mr^MOCratioGam * PicOC / PicOC0
    BCL2r  <- max(BCL2 / 100.0, 1e-6)
    KLSoc  <- kb * (1 / Lr)^LsurvOCgam / BCL2r^0.1

    dOC <- kinOC2 - KLSoc * OC

    # =================================================================
    # (L) BMD (lumbar spine)
    # =================================================================
    OBr <- max(OB / OB0, 1e-6)
    OCr <- max(OC / OC0, 1e-6)

    kinBMDls <- koutBMDls   # SS: kinBMDls = koutBMDls
    dBMDls   <- kinBMDls * OBr^gamOB - koutBMDls * OCr^gamOCls * BMDls

    # ---- 반환 ---
    list(
      c(dPTH, dS, dPTmax, dB, dAOH, dP, dECCPhos, dTgut, dRkid,
        dHAp, dPhosGut, dIntraPO, dQ, dQbone, dUCA,
        dROB1, dOBfast, dOBslow, dOC,
        dL, dRNK, dOopg, dMcmplx, dNcmplx,
        dTGFB, dTGFBact, dRX2, dCREB, dBCL2,
        dEST, dBMDls),
      OB_total = OBfast + OBslow,
      CTX_pct  = (OC / OC0 - 1) * 100,
      P1NP_pct = ((OBfast + OBslow) / OB0 - 1) * 100,
      BMD_pct  = (BMDls - 1) * 100
    )
  })
}

###############################################################################
#
#  5.  시뮬레이션 실행 함수
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

  # ---- Phase 1 : 정상상태 (120 개월) ---
  message(sprintf("[SS]  dose=%d mg  EST_tx=%.4f", dose_mg, EST_tx))
  sol_ss <- ode(y0, seq(0, ss_hr, by = 24), qsp_ode, p,
                method = "lsoda", atol = 1e-8, rtol = 1e-8,
                maxsteps = 1e5)
  n_state <- length(y0)
  y_ss <- sol_ss[nrow(sol_ss), 2:(n_state + 1)]
  names(y_ss) <- names(y0)

  OB_ss <- as.numeric(y_ss["OBfast"] + y_ss["OBslow"])
  OC_ss <- as.numeric(y_ss["OC"])
  BMD_ss <- as.numeric(y_ss["BMDls"])

  # ---- Phase 2 : 치료 ---
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

  # ---- Phase 3 : 후속관찰 (EST → 1.0) ---
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

  # ---- 파생 열 ---
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

# 월 시점에서 값 추출 헬퍼
at_month <- function(df, m, col) {
  idx <- which.min(abs(df$months - m))
  df[[col]][idx]
}

###############################################################################
#
#  6.  시나리오 실행
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
#  7.  Table 3 재현 : 예측 BMD 변화 (%)
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

cat("\n=== Predicted BMD Change (%) — Table 3 ===\n")
print(tbl3, digits = 3)
cat("\nTarget (paper):\n")
cat("150 mg QD : -0.61  -0.91  -0.96  -0.91\n")
cat("200 mg BID: -3.47  -4.95  -5.15  -4.97\n")

###############################################################################
#
#  8.  Figure 1 — Dose vs E2
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
  labs(title = "Figure 1 — Dose-E2 Relationship (Scaled Logistic)",
       x = "Daily Elagolix Dose (mg)", y = "Predicted E2 (pg/mL)") +
  theme_bw(base_size = 13)

ggsave("output/figures/fig1_dose_E2.png", p1, width = 8, height = 5, dpi = 300)
message("Saved: output/figures/fig1_dose_E2.png")

###############################################################################
#
#  9.  Figure 2 — E2 Suppression vs BMD / CTX / P1NP
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
  labs(title = "Figure 2 — E2 Suppression vs BMD Change",
       x = "E2 Suppression (%)", y = "BMD Change (%)",
       colour = "Timepoint", shape = "Timepoint") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig2_E2_vs_BMD.png", p2, width = 8, height = 5, dpi = 300)
message("Saved: output/figures/fig2_E2_vs_BMD.png")

# CTX & P1NP 보조 플롯
f2_ctx <- f2 %>%
  select(E2_suppression, CTX_6m, CTX_12m) %>%
  pivot_longer(-E2_suppression, names_to = "tp", values_to = "CTX") %>%
  mutate(tp = ifelse(tp == "CTX_6m", "6 Months", "12 Months"))

p2b <- ggplot(f2_ctx, aes(E2_suppression, CTX, colour = tp)) +
  geom_point(size = 3) + geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Figure 2b — E2 Suppression vs CTX Change",
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
  labs(title = "Figure 2c — E2 Suppression vs P1NP Change",
       x = "E2 Suppression (%)", y = "P1NP Change (%)",
       colour = "Timepoint") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig2b_E2_vs_CTX.png",  p2b, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig2c_E2_vs_P1NP.png", p2c, width = 8, height = 5, dpi = 300)

###############################################################################
#
# 10.  Figure 3 — Time-course (12 mo Tx + 6 mo FU)
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
               "Figure 3A — BMD (12-mo Tx + 6-mo FU)")
p3b <- make_tc(comb3, "CTX_chg_pct",  "CTX Change (%)",
               "Figure 3B — CTX (12-mo Tx + 6-mo FU)")
p3c <- make_tc(comb3, "P1NP_chg_pct", "P1NP Change (%)",
               "Figure 3C — P1NP (12-mo Tx + 6-mo FU)")

ggsave("output/figures/fig3a_BMD_tc.png",  p3a, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig3b_CTX_tc.png",  p3b, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig3c_P1NP_tc.png", p3c, width = 8, height = 5, dpi = 300)
message("Saved: fig3a/b/c")

###############################################################################
#
# 11.  Figure 4 — 24-month continuous
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
  labs(title = "Figure 4A — BMD: 24-Month Continuous Treatment",
       x = "Time (months)", y = "BMD Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

p4b <- ggplot(comb4, aes(months, CTX_chg_pct, colour = dose_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Figure 4B — CTX: 24-Month",
       x = "Time (months)", y = "CTX Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

p4c <- ggplot(comb4, aes(months, P1NP_chg_pct, colour = dose_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "Figure 4C — P1NP: 24-Month",
       x = "Time (months)", y = "P1NP Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig4a_BMD_24m.png",  p4a, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig4b_CTX_24m.png",  p4b, width = 8, height = 5, dpi = 300)
ggsave("output/figures/fig4c_P1NP_24m.png", p4c, width = 8, height = 5, dpi = 300)
message("Saved: fig4a/b/c")

###############################################################################
#
# 12.  Recovery Figure — 12 mo Tx + 12 mo FU
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
  labs(title = "Figure 5 — BMD Recovery (12-mo Tx + 12-mo FU)",
       x = "Time (months)", y = "BMD Change (%)", colour = "Dose") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

ggsave("output/figures/fig5_BMD_recovery.png", p5, width = 8, height = 5, dpi = 300)
message("Saved: fig5_BMD_recovery.png")

###############################################################################
#
# 13.  요약
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
