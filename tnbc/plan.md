# mTNBC QSP 모델 R 구현 계획

## 프로젝트 개요

**연구자**: 박정민 (wjdals011106)
**참고 논문**: Arulraj et al. (2023) "A transcriptome-informed QSP model of metastatic triple-negative breast cancer identifies predictive biomarkers for PD-1 inhibition" *Science Advances*
**목표**: 논문의 MATLAB SimBiology 모델을 R로 재현하고 Shiny 웹앱 개발

---

## 1. 모델 구조 개요

### 1.1 구획 (Compartments) - 9개

```
┌─────────────────────────────────────────────────────────────┐
│                      Central (Blood)                        │
│                            ↕                                │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐              │
│  │Peripheral│    │ Primary  │    │  Other   │              │
│  │          │    │  Tumor   │    │Met Tumor │              │
│  └──────────┘    └────┬─────┘    └────┬─────┘              │
│                       ↓               ↓                     │
│                 ┌──────────┐    ┌──────────┐               │
│                 │Primary LN│    │Other LN  │               │
│                 └──────────┘    └──────────┘               │
│                                                             │
│  ┌──────────┐    ┌──────────┐                              │
│  │Lung Met 1│    │Lung Met 2│                              │
│  └────┬─────┘    └────┬─────┘                              │
│       └───────┬───────┘                                     │
│               ↓                                             │
│        ┌──────────┐                                         │
│        │ Lung LN  │                                         │
│        └──────────┘                                         │
└─────────────────────────────────────────────────────────────┘
```

### 1.2 세포 유형 (Cell Types)

| 세포 유형 | 약어 | 위치 |
|----------|------|------|
| Naive CD4+ T cells | nTCD4 | Central, Peripheral, LN |
| Naive CD8+ T cells | nTCD8 | Central, Peripheral, LN |
| Activated CD4+ T cells | aTCD4 | Central, LN |
| Activated CD8+ T cells | aTCD8 | Central, LN |
| Helper T cells | Th | Central, Tumor, LN |
| Regulatory T cells | Treg | Central, Tumor, LN |
| Cytotoxic T cells | Tcyt | Central, Tumor, LN |
| Exhausted T cells | Texh | Tumor |
| Antigen Presenting Cells | APC | Tumor |
| Mature APC | mAPC | Tumor, LN |
| M1 Macrophages | M1_Mac | Tumor |
| M2 Macrophages | M2_Mac | Tumor |
| MDSCs | MDSC | Tumor |
| Cancer Cells (5 clones) | C1-C5 | Tumor |

### 1.3 사이토카인 및 분자 (Cytokines & Molecules)

- **사이토카인**: IL-2, IL-10, IL-12, IFN-γ, TGF-β, CCL2
- **Checkpoint 분자**: PD-1, PD-L1, PD-L2, CTLA-4, CD28, CD80, CD86, SIRPα, CD47
- **기타**: Angiogenic factors, Arginase-I (Arg-I), Nitric Oxide (NO)

---

## 2. 핵심 수식 (Key Equations)

### 2.1 암세포 성장 (Cancer Cell Growth) - Modified Gompertzian

```
dC_i/dt = λ_i × C_i × log(K/C_total) - k_death × C_i - k_kill × f(Tcyt, C_i)
```

- `C_i`: i번째 암 클론의 세포 수
- `λ_i`: 성장률 상수
- `K`: Carrying capacity (혈관 신생에 의해 조절)
- `k_death`: 자연 사망률
- `k_kill`: T세포에 의한 사멸률

### 2.2 Carrying Capacity (Tumor Vasculature)

```
dK/dt = k_vas_Csec × C_total × (1 - K/K_max) - k_K_d × K
```

- Angiogenic factor에 의해 증가
- 최대 용량 `K_max`에 수렴

### 2.3 T세포 활성화 (T Cell Activation)

```
Activation = H_TCR × (pMHC / (pMHC + K_TCR))
```

- Hill function 기반
- pMHC: peptide-MHC 복합체 농도

### 2.4 T세포 증식 횟수

```
N_div = N_TCR + N_CD28 + N_IL2
```

- TCR ligation 기여분
- CD28 engagement (CTLA-4에 의해 경쟁적 억제)
- IL-2 기여분

### 2.5 PD-1/PD-L1 상호작용 (Checkpoint Inhibition)

```
PD1_inhibition = (PD1_PDL1 + PD1_PDL2) / (K_PD1 + PD1_PDL1 + PD1_PDL2)
```

- Pembrolizumab이 PD-1에 결합하여 이 상호작용 차단

### 2.6 Pembrolizumab PK (2-compartment model)

```
dA_central/dt = -CL × C_central - Q × (C_central - C_peripheral) + Rate_infusion
dA_peripheral/dt = Q × (C_central - C_peripheral)
```

### 2.7 Macrophage Polarization

```
M1 ⇌ M2

dM1/dt = k_recruit × f(CCL2) - k_M1_to_M2 × f(IL10, TGFβ) × M1
         + k_M2_to_M1 × f(IL12, IFNγ) × M2 - k_death × M1

dM2/dt = k_M1_to_M2 × f(IL10, TGFβ) × M1
         - k_M2_to_M1 × f(IL12, IFNγ) × M2 - k_death × M2
```

### 2.8 MDSC 모집

```
dMDSC/dt = k_MDSC_recruit × CCL2/(CCL2 + CCL2_50) - k_MDSC_death × MDSC
```

### 2.9 반응 평가 지표

**Responder Inclusion Score (RIS)**:
```
RIS = (Responders in subgroup / Total responders)
    - (Non-responders in subgroup / Total non-responders)
```

**Diversity Index (Shannon)**:
```
H = -Σ p_k × ln(p_k)
```

**Evenness**:
```
J = H / ln(S)
```

---

## 3. 파라미터 분류 (Data S1 기준)

### 3.1 고정 파라미터 (Fixed Parameters)

| 카테고리 | 개수 | 예시 |
|---------|------|------|
| PK parameters | ~15 | CL, V_central, V_peripheral, Q |
| Cell kinetics | ~50 | k_death, k_prolif, k_migrate |
| Cytokine dynamics | ~30 | k_secretion, k_degradation, IC50 |
| Checkpoint molecules | ~40 | Expression levels, Kd values |

### 3.2 가변 파라미터 (Virtual Patient Generation용, 127개)

| 파라미터 | 분포 | 의미 |
|---------|------|------|
| Initial tumor diameter | Log-normal | 초기 종양 크기 |
| PD-L1 expression | Log-normal | PD-L1 발현 수준 |
| Cancer clone growth rates | Uniform | 클론별 성장률 |
| T cell diversity | Normal | T세포 다양성 |
| Metastasis seeding time | Uniform | 전이 발생 시점 |

---

## 4. R 구현 전략

### 4.1 단계별 접근 (Bottom-up)

```
Phase 1: Core Modules (단일 종양)
├── 1.1 Cancer growth module
├── 1.2 T cell module
├── 1.3 APC/Antigen presentation module
├── 1.4 Macrophage module
└── 1.5 Cytokine module

Phase 2: Integration
├── 2.1 Checkpoint interactions
├── 2.2 Pembrolizumab PK/PD
└── 2.3 Single tumor + LN integration

Phase 3: Extension to Metastatic Setting
├── 3.1 Multiple tumor compartments
├── 3.2 Differential immune composition
└── 3.3 Virtual patient generation

Phase 4: Analysis & Visualization
├── 4.1 RECIST response evaluation
├── 4.2 Biomarker analysis
└── 4.3 Shiny app development
```

### 4.2 R 패키지 선택

```r
# ODE Solver
library(deSolve)      # 기본 ODE solver
library(mrgsolve)     # 대안 (더 빠름)

# Data manipulation
library(tidyverse)
library(data.table)

# Visualization
library(ggplot2)
library(plotly)

# Shiny app
library(shiny)
library(shinydashboard)
library(bslib)

# Parameter sampling
library(lhs)          # Latin Hypercube Sampling
```

---

## 5. 단순화 전략 (Simplification for Student Project)

원 모델 (641 equations, 737 parameters)을 학생 프로젝트로 단순화:

### 5.1 구획 축소
- **원본**: 9 compartments
- **단순화**: 4 compartments (Central, Tumor, LN, Peripheral)
- 폐 전이 2개를 하나로 통합

### 5.2 세포 유형 축소
- **원본**: 14+ cell types
- **단순화**: 8 essential cell types
  - Cancer (단일 클론으로 시작)
  - Tcyt, Treg, Th
  - APC, mAPC
  - M1_Mac, M2_Mac, MDSC

### 5.3 사이토카인 축소
- **원본**: 6+ cytokines
- **단순화**: 4 key cytokines (IL-2, IL-10, IL-12, IFN-γ)

### 5.4 예상 규모
- **방정식**: ~50-80개
- **파라미터**: ~100개

---

## 6. 구현 일정 (12주)

### Phase 1: 기초 (1-3주)

#### Week 1: 환경 설정 및 문헌 심층 분석
- [ ] R 환경 및 패키지 설치
- [ ] 논문 및 Supplement 정독
- [ ] Data S1 (파라미터), Data S2 (바이오마커) 엑셀 파일 분석
- [ ] 모델 구조 다이어그램 작성

#### Week 2: 단일 모듈 구현 시작
- [ ] Cancer growth module 구현
- [ ] Carrying capacity dynamics 구현
- [ ] 단위 테스트

#### Week 3: T세포 모듈
- [ ] Naive T cell dynamics
- [ ] T cell activation (simplified)
- [ ] Effector T cell (Tcyt, Th, Treg) 분화

### Phase 2: 핵심 모듈 (4-6주)

#### Week 4: 면역 세포 모듈
- [ ] APC/mAPC dynamics
- [ ] Macrophage (M1/M2) polarization
- [ ] MDSC recruitment

#### Week 5: 사이토카인 및 Checkpoint
- [ ] Cytokine network (IL-2, IL-10, IL-12, IFN-γ)
- [ ] PD-1/PD-L1 interactions
- [ ] T cell exhaustion

#### Week 6: Pembrolizumab PK/PD
- [ ] 2-compartment PK model
- [ ] Target-mediated drug disposition (TMDD) 고려
- [ ] PD-1 occupancy 계산

### Phase 3: 통합 및 검증 (7-9주)

#### Week 7: 시스템 통합
- [ ] 모든 모듈 연결
- [ ] 단일 종양 + LN 시뮬레이션
- [ ] Steady-state 확인

#### Week 8: 전이 설정 확장
- [ ] Multiple tumor compartments 추가
- [ ] Immune cell composition 차이 반영
- [ ] Virtual patient 생성 (Latin Hypercube Sampling)

#### Week 9: 모델 검증
- [ ] 논문 Figure 3 (response rates) 재현
- [ ] KEYNOTE-119 데이터와 비교
- [ ] 민감도 분석

### Phase 4: Shiny 앱 개발 (10-12주)

#### Week 10: 앱 UI 설계
- [ ] 입력 패널 (용량, 투여 간격, 환자 특성)
- [ ] 출력 패널 레이아웃
- [ ] 반응형 디자인

#### Week 11: 핵심 기능 구현
- [ ] 시뮬레이션 실행 및 결과 표시
- [ ] 종양 크기 변화 그래프
- [ ] 면역 세포 조성 시각화
- [ ] Waterfall plot

#### Week 12: 마무리
- [ ] 바이오마커 분석 기능
- [ ] 최적화 및 디버깅
- [ ] 배포 (shinyapps.io)
- [ ] 문서화

---

## 7. Shiny 앱 설계

### 7.1 입력 패널

```
┌─────────────────────────────────────┐
│ Simulation Settings                 │
├─────────────────────────────────────┤
│ Pembrolizumab Dose: [200] mg        │
│ Dosing Interval:   [3 weeks] ▼      │
│ Treatment Duration:[52] weeks       │
├─────────────────────────────────────┤
│ Patient Characteristics             │
├─────────────────────────────────────┤
│ Initial Tumor Size: [1.5] cm        │
│ PD-L1 Expression:   [300] mol/μm²   │
│ TIL Fraction:       [0.4]           │
│ Number of Metastases: [2]           │
├─────────────────────────────────────┤
│ [Run Single Patient]                │
│ [Run Virtual Trial (N=100)]         │
└─────────────────────────────────────┘
```

### 7.2 출력 패널

```
┌─────────────────────────────────────────────────────────┐
│ Tab 1: Tumor Response                                   │
│ ┌─────────────────────┐ ┌─────────────────────┐        │
│ │ Tumor Size vs Time  │ │ Waterfall Plot      │        │
│ │ (Spider Plot)       │ │                     │        │
│ └─────────────────────┘ └─────────────────────┘        │
├─────────────────────────────────────────────────────────┤
│ Tab 2: Immune Microenvironment                         │
│ ┌─────────────────────┐ ┌─────────────────────┐        │
│ │ Cell Composition    │ │ Cytokine Levels     │        │
│ │ (Stacked Area)      │ │ (Line Plot)         │        │
│ └─────────────────────┘ └─────────────────────┘        │
├─────────────────────────────────────────────────────────┤
│ Tab 3: PK/PD Profile                                   │
│ ┌─────────────────────┐ ┌─────────────────────┐        │
│ │ Drug Concentration  │ │ PD-1 Occupancy      │        │
│ └─────────────────────┘ └─────────────────────┘        │
├─────────────────────────────────────────────────────────┤
│ Tab 4: Biomarker Analysis                              │
│ ┌─────────────────────┐ ┌─────────────────────┐        │
│ │ Response by         │ │ RIS Score           │        │
│ │ Biomarker Level     │ │ Ranking             │        │
│ └─────────────────────┘ └─────────────────────┘        │
└─────────────────────────────────────────────────────────┘
```

### 7.3 주요 출력 지표

| 지표 | 설명 |
|------|------|
| Tumor diameter change (%) | RECIST 기준 |
| Response status | CR/PR/SD/PD |
| TIL fraction | 종양 내 림프구 비율 |
| M2/M1 ratio | 대식세포 극성화 |
| PD-1 occupancy (%) | 약물 점유율 |
| Cancer clone richness | 클론 다양성 |

---

## 8. 파일 구조

```
tnbc/
├── plan.md                      # 이 문서
├── model/
│   ├── parameters.R             # 파라미터 정의
│   ├── modules/
│   │   ├── cancer_growth.R      # 암세포 성장
│   │   ├── t_cell.R             # T세포 dynamics
│   │   ├── apc.R                # APC/항원 제시
│   │   ├── macrophage.R         # 대식세포
│   │   ├── mdsc.R               # MDSC
│   │   ├── cytokines.R          # 사이토카인
│   │   ├── checkpoint.R         # Checkpoint 상호작용
│   │   └── pembrolizumab_pk.R   # 약물 PK
│   ├── tnbc_qsp_model.R         # 통합 모델
│   └── virtual_patient.R        # 가상 환자 생성
├── shiny/
│   ├── app.R
│   ├── ui.R
│   ├── server.R
│   └── modules/
│       ├── simulation_module.R
│       ├── visualization_module.R
│       └── biomarker_module.R
├── data/
│   ├── parameters_s1.xlsx       # Data S1
│   └── biomarkers_s2.xlsx       # Data S2
├── tests/
│   └── test_modules.R
└── docs/
    ├── model_equations.md
    └── validation_report.md
```

---

## 9. 검증 목표

### 9.1 정량적 검증 (논문 Figure 재현)

| Figure | 내용 | 목표 |
|--------|------|------|
| Fig 3A | Response rates | ORR ~10%, CR/PR 재현 |
| Fig 3B | Waterfall plot | 분포 패턴 유사 |
| Fig 3C | Spider plot | 개별 환자 trajectory |
| Fig 5D | Biomarker ranking | Top biomarkers 일치 |

### 9.2 정성적 검증

- [ ] Pembrolizumab 투여 시 TIL 증가
- [ ] M2/M1 ratio 변화 패턴
- [ ] Cancer clone richness 감소 (responders)
- [ ] 폐 전이 vs 기타 전이 반응 차이

---

## 10. 참고 자료

### 10.1 원본 코드
- MATLAB code: http://dx.doi.org/10.17632/r46rk4vwdv.1

### 10.2 관련 논문
- Wang et al. (2022) iScience - Base TNBC QSP model
- Milberg et al. (2019) Scientific Reports - Checkpoint inhibitor QSP

### 10.3 R 패키지 문서
- deSolve: https://cran.r-project.org/web/packages/deSolve/
- mrgsolve: https://mrgsolve.org/

---

## 11. 위험 요소 및 대응

| 위험 | 가능성 | 대응 |
|------|--------|------|
| 모델 복잡성 | 높음 | 단순화 버전으로 시작 |
| 파라미터 불확실성 | 중간 | 논문 Data S1 활용 |
| 계산 시간 | 중간 | mrgsolve 사용, 병렬화 |
| Shiny 성능 | 낮음 | Pre-computed 결과 캐싱 |

---

*Last updated: 2026-01-23*
