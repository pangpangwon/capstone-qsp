###############################################################################
# app.R  —  Elagolix QSP Model  ·  Interactive Shiny Dashboard
# -----------------------------------------------------------------------
# 모델 코어는 _run_validation.R 에서 source() 로 로드
# (predict_E2, make_params, make_init, qsp_ode, run_sim, at_month)
###############################################################################

library(shiny)
library(bslib)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ── 모델 코어 로드 (실행 코드 제외, 함수만 사용) ──────────────────────────
# _run_validation.R 의 6~13번 섹션(시나리오 자동 실행)은 source 시 실행되므로
# 함수만 별도로 가져오기 위해 여기서 직접 정의하지 않고 source 후
# 필요한 전역 객체만 사용한다.
# 아래는 _run_validation.R 섹션 1~5 (함수 정의부)만 재사용하는 방식이다.

# ---- 상수 -----------------------------------------------------------------
HOURS_PER_MONTH <- 730.5

# ---- 모듈 A : Dose → E2 ---------------------------------------------------
dose_e2_params <- list(slope = 0.00894, logE2max = 5.20, logE2min = 2.14)

predict_E2 <- function(daily_dose, p = dose_e2_params) {
  E2min <- exp(p$logE2min)
  E2max <- exp(p$logE2max)
  E2min + (E2max - E2min) / (1 + exp(p$slope * daily_dose))
}

E2_baseline <- predict_E2(0)

# ---- 모듈 B : 파라미터 / 초기조건 / ODE ------------------------------------
make_params <- function() {
  list(
    OBtot0=0.00501324, k1=6.24e-6, k2=0.112013, k3=6.24e-6, k4=0.112013,
    kO=15.8885, kb=0.000605516, Pic0=0.228142,
    FracOBfast=0.797629, Frackb=0.313186,
    OBtgfGAM=0.0111319, koutTGF0=2.98449e-5, koutTGFGam=0.919131,
    OCtgfGAM=0.593891,
    kinOCgam=8.53065, EmaxPicOC=1.9746, FracPicOC=0.878215,
    PicOCgam=1.0168, MOCratioGam=0.603754, LsurvOCgam=3.0923,
    LsurvOCCgam=3.09023,
    EmaxPicROB=3.9745, PicROBgam=1.80968, FracPicROB=0.883824,
    PicOBgam=0.122313, FracPicOB=0.000244818, EmaxPicOB=0.251636,
    CaDay=88, OralCa=24.055/24, OralPhos=10.5/24, F12=0.7,
    V1=14, Da=0.7/24, Reabs50=1.57322,
    kout_PTH=100/14, PTout=1.604e-4, opgPTH50=3.85,
    AlphOHgam=0.111241, CtriolPTgam=12.5033, CtriolMax=4.1029, CtriolMin=0.9,
    E0RANKL=3.80338, EmaxL=0.469779, koutL=0.00293273,
    kinRNKgam=0.151825, koutRNK=0.00323667,
    RUNX20=10, RX2Kout0=0.693, E0rx2Kout=0.125, EmaxPTHRX2x=5,
    E0crebKin=0.5, EmaxPTHcreb=3.39745, crebKout=0.00279513, bcl2Kout=0.693,
    koutBMDls=0.000397, gamOB=0.0793, gamOCls=0.14,
    ESTON=1, koutEST=0.05776227,
    tgfbGAM=0.0374, tgfbactGAM=0.045273, robGAM=0.16, obGAM=0.000012,
    E2scalePicB1=1.16832e-5, maxTmESTkid=0.923737,
    Smax=8, Smin=0.2, gamS=3, GFR=6,
    TPmaxCA=2.15, CaPset=32.9, CaFrac=0.5, T85_base=0.25,
    kinB0=126, T69_base=0.1, T64_base=1,
    PhosKid=0.8, IntraPO_kel=0.01, IntraPO_kin=32.26,
    HApForm=0.001, HApResorb=0.001
  )
}

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

qsp_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    OB0<-OBtot0; OC0<-0.00115398; ROB0<-0.00104122
    L0<-0.4; RNK0<-10; O0<-4; PTH0<-53.9; B0<-1260; P0<-32.9
    M0<-k3*RNK0*L0/k4; RX20v<-RUNX20
    OB<-OBfast+OBslow; EST_s<-max(EST,1e-6)

    kinEST<-koutEST; dEST<-(kinEST-koutEST*EST)*ESTON
    CaRatio<-P/P0
    Sfun<-Smin+(Smax-Smin)/(1+CaRatio^gamS)
    SPTH<-Sfun*kout_PTH*PTH0; dPTH<-SPTH-kout_PTH*PTH
    dS<-0; dPTmax<-0; PTHr<-PTH/PTH0
    SE<-kinB0*PTHr^AlphOHgam; dAOH<-SE-T64_base*AOH
    dB<-AOH-T69_base*B; Br<-B/B0
    T85<-T85_base*(CtriolMin+(CtriolMax-CtriolMin)*Br^CtriolPTgam/(1+Br^CtriolPTgam))
    J40<-Da*Tgut; dTgut<-OralCa*T85-J40
    ESTkid<-max(1+maxTmESTkid*(EST_s-1)*ESTON,0.1)
    TPeff<-TPmaxCA*ESTkid; CaFilt<-GFR*CaFrac*P/V1
    T36<-CaFilt*TPeff/(Reabs50+CaFrac*P/V1)
    dRkid<-T36*(1-Rkid)-CaFilt*Rkid*0.01
    Jurine<-CaFilt*(1-Rkid)
    J14<-0.02*Q; J15<-0.02*P*3
    J14a<-HApResorb*(OC/OC0)*Q; J15a<-HApForm*(OB/OB0)*P
    dP<-(J14-J15-Jurine+J40+J14a-J15a)/V1
    dQ<-J15-J14+J14a-J15a; dQbone<-J15a-J14a
    dUCA<-Jurine; dHAp<-HApForm*(OB/OB0)-HApResorb*(OC/OC0)
    J53<-Da*PhosGut; J41<-OralPhos*F12; J42<-PhosKid*ECCPhos
    J48<-0.001*ECCPhos; J54<-IntraPO_kin*ECCPhos/(ECCPhos+10)
    J56<-IntraPO_kel*IntraPO
    dECCPhos<-J41+J53-J42-J48-J54+J56
    dPhosGut<-OralPhos*F12-J53; dIntraPO<-J54-J56
    ptf<-PTHr/(1+PTHr)
    RX2Kin<-RX2Kout0*RX20v*(1+E0rx2Kout+EmaxPTHRX2x*ptf)
    RX2Kout<-RX2Kout0*(1+E0rx2Kout+EmaxPTHRX2x*ptf)
    dRX2<-RX2Kin-RX2Kout*RX2
    crebKin<-crebKout*10*(E0crebKin+EmaxPTHcreb*ptf)
    dCREB<-crebKin-crebKout*CREB
    dBCL2<-bcl2Kout*(CREB/10)*(RX2/RX20v)*100-bcl2Kout*BCL2
    kinTGF<-koutTGF0*Pic0*1000
    koutTGFeff<-koutTGF0*(OC/OC0)^koutTGFGam
    dTGFB<-kinTGF*(OB/OB0)^OBtgfGAM*((1/EST_s)^tgfbGAM*ESTON+(1-ESTON))-
      koutTGFeff*TGFB*(EST_s^tgfbactGAM*ESTON+(1-ESTON))
    koutTGFact<-koutTGF0*1000
    dTGFBact<-koutTGFeff*TGFB*(EST_s^tgfbactGAM*ESTON+(1-ESTON))-koutTGFact*TGFBact
    PicR<-TGFBact/Pic0
    kinL<-koutL*L0*(E0RANKL+EmaxL*ptf)/(E0RANKL+EmaxL*0.5)
    dL<-kinL-koutL*L-k1*Oopg*L+k2*Ncmplx-k3*RNK*L+k4*Mcmplx
    kinRNK<-koutRNK*RNK0*(TGFBact/Pic0)^kinRNKgam
    dRNK<-kinRNK-koutRNK*RNK-k3*RNK*L+k4*Mcmplx
    PTHinhOPG<-opgPTH50/(opgPTH50+PTH)
    PTHinhOPG0<-opgPTH50/(opgPTH50+PTH0)
    pO<-kO*O0*PTHinhOPG/PTHinhOPG0*(EST_s^0.05*ESTON+(1-ESTON))
    dOopg<-pO-k1*Oopg*L+k2*Ncmplx-kO*Oopg
    dMcmplx<-k3*RNK*L-k4*Mcmplx; dNcmplx<-k1*Oopg*L-k2*Ncmplx
    PicROB<-1+EmaxPicROB*FracPicROB*PicR^PicROBgam/(1+FracPicROB*PicR^PicROBgam)
    PicROB0<-1+EmaxPicROB*FracPicROB/(1+FracPicROB)
    ROBin<-kb*ROB0*PicROB/PicROB0
    ESTrob<-(1/EST_s)^robGAM*ESTON+(1-ESTON)
    dROB1<-ROBin*ESTrob-kb*ROB1
    PicOB<-1-EmaxPicOB*FracPicOB*PicR^PicOBgam/(1+FracPicOB*PicR^PicOBgam)
    PicOB0<-1-EmaxPicOB*FracPicOB/(1+FracPicOB)
    PicOB<-max(PicOB,1e-4); PicOB0<-max(PicOB0,1e-4)
    Drate<-kb*ROB1; bigDb<-Drate*PicOB/PicOB0
    ESTob<-EST_s^obGAM*ESTON+(1-ESTON)
    kbfast<-kb*10; kbslow<-kb*0.5
    dOBfast<-bigDb*FracOBfast*Frackb*ESTob-kbfast*OBfast
    dOBslow<-bigDb*(1-FracOBfast)*Frackb-kbslow*OBslow
    Mr<-max(Mcmplx/M0,1e-6); Lr<-max(L/L0,1e-6)
    PicOC<-1+EmaxPicOC*FracPicOC*PicR^PicOCgam/(1+FracPicOC*PicR^PicOCgam)
    PicOC0<-1+EmaxPicOC*FracPicOC/(1+FracPicOC)
    kinOC2<-kb*OC0*Mr^MOCratioGam*PicOC/PicOC0
    BCL2r<-max(BCL2/100,1e-6)
    KLSoc<-kb*(1/Lr)^LsurvOCgam/BCL2r^0.1
    dOC<-kinOC2-KLSoc*OC
    OBr<-max(OB/OB0,1e-6); OCr<-max(OC/OC0,1e-6)
    kinBMDls<-koutBMDls
    dBMDls<-kinBMDls*OBr^gamOB-koutBMDls*OCr^gamOCls*BMDls

    list(c(dPTH,dS,dPTmax,dB,dAOH,dP,dECCPhos,dTgut,dRkid,
           dHAp,dPhosGut,dIntraPO,dQ,dQbone,dUCA,
           dROB1,dOBfast,dOBslow,dOC,
           dL,dRNK,dOopg,dMcmplx,dNcmplx,
           dTGFB,dTGFBact,dRX2,dCREB,dBCL2,dEST,dBMDls),
         OB_total=OBfast+OBslow,
         CTX_pct=(OC/OC0-1)*100,
         P1NP_pct=((OBfast+OBslow)/OB0-1)*100,
         BMD_pct=(BMDls-1)*100)
  })
}

# ---- run_sim (from _run_validation.R) -------------------------------------
run_sim <- function(dose_mg, tx_months, fu_months=0, ss_months=120) {
  p  <- make_params(); y0 <- make_init(p)
  ss_hr <- ss_months*HOURS_PER_MONTH
  tx_hr <- tx_months*HOURS_PER_MONTH
  fu_hr <- fu_months*HOURS_PER_MONTH
  EST_tx <- predict_E2(dose_mg)/E2_baseline

  sol_ss <- ode(y0, seq(0,ss_hr,by=24), qsp_ode, p,
                method="lsoda", atol=1e-8, rtol=1e-8, maxsteps=1e5)
  n_state <- length(y0)
  y_ss <- sol_ss[nrow(sol_ss), 2:(n_state+1)]; names(y_ss)<-names(y0)
  OB_ss<-as.numeric(y_ss["OBfast"]+y_ss["OBslow"])
  OC_ss<-as.numeric(y_ss["OC"]); BMD_ss<-as.numeric(y_ss["BMDls"])

  y_tx<-y_ss; y_tx["EST"]<-EST_tx
  ev_tx<-function(t,y,parms){y["EST"]<-EST_tx;y}
  times_tx<-seq(0,tx_hr,by=24)
  sol_tx<-ode(y_tx,times_tx,qsp_ode,p,method="lsoda",
              atol=1e-8,rtol=1e-8,maxsteps=1e5,
              events=list(func=ev_tx,time=times_tx))
  res<-as.data.frame(sol_tx)

  if(fu_months>0){
    y_fu<-sol_tx[nrow(sol_tx), 2:(n_state+1)]; names(y_fu)<-names(y0); y_fu["EST"]<-1
    sol_fu<-ode(y_fu,seq(0,fu_hr,by=24),qsp_ode,p,
                method="lsoda",atol=1e-8,rtol=1e-8,maxsteps=1e5)
    res_fu<-as.data.frame(sol_fu); res_fu$time<-res_fu$time+tx_hr
    res<-rbind(res,res_fu[-1,])
  }
  res$months      <- res$time/HOURS_PER_MONTH
  res$OB_total    <- res$OBfast+res$OBslow
  res$BMD_chg_pct <- (res$BMDls/BMD_ss-1)*100
  res$CTX_chg_pct <- (res$OC/OC_ss-1)*100
  res$P1NP_chg_pct<- (res$OB_total/OB_ss-1)*100
  res$dose_mg     <- dose_mg
  res$dose_label  <- dplyr::case_when(
    dose_mg==150~"150 mg QD", dose_mg==400~"200 mg BID",
    dose_mg==600~"300 mg BID", TRUE~paste0(dose_mg," mg"))
  res
}

at_month <- function(df,m,col) df[[col]][which.min(abs(df$months-m))]

# ── 논문 목표값 (Table 3) ──────────────────────────────────────────────────
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
  window_title = "Elagolix QSP · Stodtmann 2021",

  # ── Tab 1 : Dose-E2 ──────────────────────────────────────────────────────

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
      card(card_header("Figure 1 — Dose vs E2"), plotOutput("fig1_plot", height = "450px")),
      layout_column_wrap(
        width = 1/3,
        value_box("Baseline E2", textOutput("vb_baseline"), theme = "primary"),
        value_box("Predicted E2", textOutput("vb_pred_e2"), theme = "info"),
        value_box("Suppression", textOutput("vb_supp"), theme = "warning")
      )
    )
  ),

  # ── Tab 2 : QSP Simulation ──────────────────────────────────────────────
  nav_panel("QSP Simulation",
    layout_sidebar(
      sidebar = sidebar(
        title = "Simulation Settings",
        width = 300,
        selectInput("dose_sel", "Elagolix Dose",
                     choices = c("150 mg QD" = 150, "200 mg BID" = 400,
                                 "300 mg BID" = 600),
                     selected = 150),
        sliderInput("tx_mo", "Treatment (months)", 1, 36, 12, step = 1),
        sliderInput("fu_mo", "Follow-up (months)", 0, 24, 6, step = 1),
        hr(),
        actionButton("run_btn", "Run Simulation",
                     class = "btn-primary btn-lg w-100",
                     icon = icon("play")),
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

  # ── Tab 3 : Multi-dose Comparison ───────────────────────────────────────
  nav_panel("Multi-Dose Comparison",
    layout_sidebar(
      sidebar = sidebar(
        title = "Comparison Settings",
        width = 300,
        checkboxGroupInput("doses_cmp", "Select Doses",
                            choices = c("150 mg QD" = 150,
                                        "200 mg BID" = 400,
                                        "300 mg BID" = 600),
                            selected = c(150, 400)),
        sliderInput("tx_mo_cmp", "Treatment (months)", 1, 36, 24, step = 1),
        sliderInput("fu_mo_cmp", "Follow-up (months)", 0, 24, 0, step = 1),
        hr(),
        actionButton("run_cmp", "Run Comparison",
                     class = "btn-success btn-lg w-100",
                     icon = icon("layer-group")),
        hr(),
        helpText("Compare multiple regimens side-by-side.")
      ),
      navset_card_tab(
        nav_panel("BMD",  plotOutput("cmp_bmd",  height = "420px")),
        nav_panel("CTX",  plotOutput("cmp_ctx",  height = "420px")),
        nav_panel("P1NP", plotOutput("cmp_p1np", height = "420px"))
      ),
      card(card_header("Table 3 Comparison — BMD Change (%)"),
           DTOutput("tbl3_cmp"))
    )
  ),

  # ── Tab 4 : Validation Table ────────────────────────────────────────────
  nav_panel("Validation (Table 3)",
    card(
      card_header("Paper Target: Predicted BMD Change (%)"),
      DTOutput("tbl_paper_target")
    ),
    card(
      card_header("Model Prediction (run Comparison first)"),
      DTOutput("tbl_model_pred")
    )
  ),

  # ── Tab 5 : About ───────────────────────────────────────────────────────
  nav_panel("About",
    card(
      card_header("Model Information"),
      tags$ul(
        tags$li(tags$b("Paper:"), " Stodtmann et al. (2021) Clin Transl Sci. 14:1611-1619"),
        tags$li(tags$b("Foundation:"), " Peterson & Riggs (2010) + Riggs et al. (2012)"),
        tags$li(tags$b("Reference impl.:"), " MetrumRG OpenBoneMin"),
        tags$li(tags$b("State variables:"), " 31"),
        tags$li(tags$b("Parameters:"), sprintf(" %d", length(make_params()))),
        tags$li(tags$b("ODE solver:"), " deSolve::lsoda (atol/rtol = 1e-8)")
      ),
      hr(),
      tags$h5("Module A — Dose-E2 (Scaled Logistic)"),
      tags$pre("E2 = exp(logE2min) + (exp(logE2max) - exp(logE2min)) / (1 + exp(slope * Dose))"),
      tags$h5("Module B — Calcium Homeostasis & Bone Remodeling QSP"),
      tags$p("31 ODEs covering: PTH, calcitriol, Ca/Phos homeostasis,
              ROB, OBfast, OBslow, OC, RANKL/RANK/OPG, TGF-β,
              RUNX2, CREB, BCL2, estrogen, BMD (lumbar spine).")
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
      dose  = c(0, 150, 400, 600),
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
      labs(title = "Figure 1 — Dose-E2 Relationship (Scaled Logistic)",
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
                                               "BMD — Multi-dose Comparison"))
  output$cmp_ctx  <- renderPlot(make_cmp_plot("CTX_chg_pct",  "CTX Change (%)",
                                               "CTX — Multi-dose Comparison"))
  output$cmp_p1np <- renderPlot(make_cmp_plot("P1NP_chg_pct", "P1NP Change (%)",
                                               "P1NP — Multi-dose Comparison"))

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
