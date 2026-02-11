# QSP 모델 구현 계획서
## Elagolix 칼슘 항상성 및 골 리모델링 모델 (R 코드 구현)

**논문**: Stodtmann et al. (2021) "Validation of a quantitative systems pharmacology model of calcium homeostasis using elagolix Phase 3 clinical trial data in women with endometriosis" *Clin Transl Sci.* 14:1611-1619.

**기반 모델**: Peterson & Riggs (2010) + Riggs et al. (2012) 칼슘 항상성 및 골 리모델링 QSP 모델

**참고 구현체**: MetrumRG OpenBoneMin R 패키지 (https://github.com/metrumresearchgroup/OpenBoneMin)

---

## 1. 모델 개요

이 QSP 모델은 두 개의 핵심 모듈로 구성된다:

### 모듈 A: Elagolix 용량-에스트라디올(E2) 모델
- Elagolix 일일 용량에서 평균 E2 수준을 예측하는 비선형 회귀 모델

### 모듈 B: 칼슘 항상성 및 골 리모델링 QSP 모델
- E2 수준 변화가 골 바이오마커(CTX, P1NP) 및 BMD에 미치는 영향을 예측
- Peterson & Riggs (2010)의 29개 ODE 기반 모델 + Riggs et al. (2012)의 에스트로겐 확장

---

## 2. 모듈 A: 용량-E2 모델 (Dose-Response)

### 2.1 수식

최종 선택 모델: Scaled Logistic Function (Equation 3)

```
E2 = exp(logE2min) + (exp(logE2max) - exp(logE2min)) / (1 + exp(slope * DailyDose))
```

비교 모델들:
- Linear: `E2 = intercept - slope * DailyDose`
- Exponential: `E2 = exp(-(log(intercept) - slope * DailyDose))`

### 2.2 최종 파라미터 (Table 2)

| 파라미터 | 추정값 | 표준오차 |
|----------|--------|----------|
| slope | 0.00894 | 0.000427 |
| log(E2max) | 5.20 | 0.0314 |
| log(E2min) | 2.14 | 0.104 |

- E2max = exp(5.20) = 181.3 pg/ml
- E2min = exp(2.14) = 8.50 pg/ml
- 기저 E2 (0 mg 용량 시) = (exp(logE2max) + exp(logE2min)) / 2 = 약 94.9 pg/ml

### 2.3 모델 적합도 비교 (Supplement Table S1)

| 모델 | Deviance | AIC | BIC |
|------|----------|-----|-----|
| logit | 1452 | 9326 | 9337 |
| exponential | 1911 | 9374 | 9382 |
| linear | 6639 | 9514 | 9522 |

### 2.4 Elagolix 용량별 예측 E2 수준

| 용량 | DailyDose (mg) | 예측 E2 (pg/ml) | E2 억제율 (%) |
|------|----------------|-----------------|---------------|
| 150 mg QD | 150 | ~47 pg/ml | ~50% |
| 200 mg BID | 400 | ~13 pg/ml | ~86% |
| 300 mg BID | 600 | ~10 pg/ml | ~89% |

---

## 3. 모듈 B: 칼슘 항상성 & 골 리모델링 QSP 모델

### 3.1 모델 구조 개요

모델은 다음 생리학적 구획(compartment)을 포함:
- **소화관 (Gut)**: 칼슘 및 인산 경구 흡수
- **혈관/세포외액 (Vasculature/ECC)**: 혈중 칼슘, 인산 농도
- **신장 (Kidney)**: 칼슘/인산 재흡수 및 배설
- **부갑상선 (Parathyroid Gland)**: PTH 분비 조절
- **뼈 (Bone)**: 골 리모델링 (형성 및 흡수)
- **조골세포/파골세포 (Osteoblasts/Osteoclasts)**: 골 세포 역학
- **세포 내 신호전달**: RANK/RANKL/OPG, TGF-β, RUNX2, CREB, BCL2

### 3.2 주요 상태 변수 (State Variables) 및 초기값

#### 칼슘/인산 항상성 관련
| 변수 | 설명 | 초기값 | 단위 |
|------|------|--------|------|
| PTH | 부갑상선 호르몬 | 53.90 | pM |
| S | PTH 분비 신호 | 0.5 | - |
| PTmax | 부갑상선 최대 용량 | 1.00 | - |
| B | 칼시트리올 (활성 비타민D) | 1260.0 | fM |
| AOH | 1-α-히드록실라제 | 126.0 | - |
| P | 세포외 칼슘 | 32.90 | (≈2.35 mM 환산) |
| ECCPhos | 세포외 인산 | 16.8 | - |
| T | 장관 칼슘 | 1.58471 | - |
| R | 신장 칼슘 재흡수율 | 0.50 | - |
| HAp | 골 히드록시아파타이트 | 1.00 | - |
| PhosGut | 장내 인산 | 0.839 | - |
| IntraPO | 세포내 인산 | 3226.0 | - |
| Q | 골 칼슘 교환가능 풀 | 100.0 | - |
| Qbone | 골 칼슘 저장 풀 | 24900.0 | - |
| UCA | 누적 뇨 칼슘 | 0.0 | - |

#### 골 세포 역학 관련
| 변수 | 설명 | 초기값 | 단위 |
|------|------|--------|------|
| ROB1 | 반응 조골세포 (Responding OB) | 0.00104122 | - |
| OBfast | 활성 조골세포 (빠른 전환) | OBtot0 * FracOBfast | - |
| OBslow | 활성 조골세포 (느린 전환) | OBtot0 * (1-FracOBfast) | - |
| OC | 파골세포 | 0.00115398 | - |

#### 신호전달 분자
| 변수 | 설명 | 초기값 | 단위 |
|------|------|--------|------|
| L | RANKL | 0.4 | pM |
| RNK | RANK | 10.0 | pM |
| O | OPG | 4.0 | pM |
| M | RANK-RANKL 복합체 | k3*RNK_0*L_0/k4 | - |
| N | OPG-RANKL 복합체 | k1*O_0*L_0/k2 | - |
| TGFB | 잠재 TGF-β | Pic0*1000 | - |
| TGFBact | 활성 TGF-β | Pic0 | - |
| RX2 | RUNX2 | 10.0 | - |
| CREB | CREB | 10.0 | - |
| BCL2 | BCL2 | 100.0 | - |

#### 에스트로겐 & BMD
| 변수 | 설명 | 초기값 | 단위 |
|------|------|--------|------|
| EST | 에스트로겐 (상대값) | 1.0 | - |
| BMDls | 요추 BMD | 1.0 | 상대값 |

### 3.3 핵심 파라미터 (~90개 이상)

#### 골세포 역학 파라미터
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| OBtot0 | 0.00501324 | 기저 총 조골세포 수 |
| k1 | 0.00000624 | OPG-RANKL 결합 속도 |
| k2 | 0.112013 | OPG-RANKL 해리 속도 |
| k3 | 0.00000624 | RANK-RANKL 결합 속도 |
| k4 | 0.112013 | RANK-RANKL 해리 속도 |
| kO | 15.8885 | OPG 제거 속도 |
| kb | 0.000605516 | 조골세포 기저 전환 |
| Pic0 | 0.228142 | 기저 TGF-β 활성값 |
| FracOBfast | 0.797629 | 빠른 전환 조골세포 분율 |
| Frackb | 0.313186 | 골 형성 분율 |

#### TGF-β 관련
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| OBtgfGAM | 0.0111319 | OB→TGF-β 생산 감마 |
| koutTGF0 | 0.0000298449 | TGF-β 제거 기저값 |
| koutTGFGam | 0.919131 | TGF-β 제거 감마 |
| OCtgfGAM | 0.593891 | OC의 TGF-β 영향 감마 |

#### 파골세포 관련
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| kinOCgam | 8.53065 | OC 유입 감마 |
| EmaxPicOC | 1.9746 | TGF-β→OC Emax |
| FracPicOC | 0.878215 | TGF-β→OC 분율 |
| PicOCgam | 1.0168 | TGF-β→OC 감마 |
| MOCratioGam | 0.603754 | RANK-RANKL→OC 비율 감마 |
| LsurvOCgam | 3.0923 | RANKL→OC 생존 감마 |
| LsurvOCCgam | 3.09023 | RANKL→OC 세포 감마 |

#### 조골세포 관련
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| EmaxPicROB | 3.9745 | TGF-β→ROB Emax |
| PicROBgam | 1.80968 | TGF-β→ROB 감마 |
| FracPicROB | 0.883824 | TGF-β→ROB 분율 |
| PicOBgam | 0.122313 | TGF-β→OB 감마 |
| FracPicOB | 0.000244818 | TGF-β→OB 분율 |
| EmaxPicOB | 0.251636 | TGF-β→OB Emax |

#### 칼슘/인산 항상성
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| CaDay | 88.0 | 일일 칼슘 총량 |
| OralCa | 24.055/24 | 시간당 경구 칼슘 |
| OralPhos | 10.5/24 | 시간당 경구 인산 |
| F12 | 0.7 | 인산 흡수 분율 |
| V1 | 14.0 | 세포외액 부피 (L) |
| Da | 0.7/24 | 장관 흡수율 |
| Reabs50 | 1.57322 | 신장 칼슘 재흡수 EC50 |

#### PTH 관련
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| kout | 100/14 | PTH 제거 속도 (/hr) |
| PTout | 0.0001604 | 부갑상선 제거 속도 |
| opgPTH50 | 3.85 | PTH→OPG 억제 EC50 |

#### 칼시트리올/비타민D
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| AlphOHgam | 0.111241 | 1-α-히드록실라제 감마 |
| CtriolPTgam | 12.5033 | PTH→칼시트리올 감마 |
| CtriolMax | 4.1029 | 칼시트리올 Emax |
| CtriolMin | 0.9 | 칼시트리올 Emin |

#### RANKL/OPG
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| E0RANKL | 3.80338 | RANKL 기저 생산율 |
| EmaxL | 0.469779 | PTH→RANKL Emax |
| koutL | 0.00293273 | RANKL 제거 속도 |
| kinRNKgam | 0.151825 | RANK 유입 감마 |
| koutRNK | 0.00323667 | RANK 제거 속도 |

#### 세포내 신호전달 (RUNX2, CREB, BCL2)
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| RUNX20 | 10.0 | RUNX2 초기값 |
| RX2Kout0 | 0.693 | RUNX2 제거 기저 |
| E0rx2Kout | 0.125 | RUNX2 제거 E0 |
| EmaxPTHRX2x | 5.0 | PTH→RUNX2 Emax |
| E0crebKin | 0.5 | CREB 유입 E0 |
| EmaxPTHcreb | 3.39745 | PTH→CREB Emax |
| crebKout | 0.00279513 | CREB 제거 속도 |
| bcl2Kout | 0.693 | BCL2 전환 속도 |

#### BMD 관련
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| koutBMDls | 0.000397 | 요추 BMD 제거 속도 (/hr) |
| gamOB | 0.0793 | OB→BMD 형성 감마 |
| gamOCls | 0.14 | OC→BMD 흡수 감마 (요추) |

#### 에스트로겐 관련 확장 (Riggs 2012)
| 파라미터 | 값 | 설명 |
|----------|-----|------|
| ESTON | 0.0 (1.0 활성화 시) | 에스트로겐 모듈 스위치 |
| koutEST | 0.05776227 | 에스트로겐 제거 속도 |
| tgfbGAM | 0.0374 | E→잠재 TGF-β 감마 (논문: 0.075) |
| tgfbactGAM | 0.045273 | E→활성 TGF-β 감마 (논문: 0.045) |
| robGAM | 0.16 | E→ROB 감마 |
| obGAM | 0.000012 | E→OB 감마 |
| E2scalePicB1 | 0.0000116832 | E2→OB 스케일링 |
| maxTmESTkid | 0.923737 | E→신장 칼슘 재흡수 최대 |

### 3.4 핵심 ODE 방정식

#### 에스트로겐 역학
```
dEST/dt = (kinEST - koutEST * EST) * ESTON
```
- `ESTON = 1`로 설정하여 에스트로겐 모듈 활성화
- Elagolix 치료 시, 용량-E2 모델로 계산한 E2 수준을 EST에 반영

#### TGF-β 역학
```
dTGFB/dt = kinTGF * (OB/OB0)^OBtgfGAM * (1/EST)^tgfbGAM - koutTGF * EST^tgfbactGAM
dTGFBact/dt = koutTGF * EST^tgfbactGAM - koutTGFact * TGFBact
```
- 에스트로겐 감소 → 잠재 TGF-β 생산 증가 (`(1/EST)^tgfbGAM`)
- 에스트로겐 감소 → 활성 TGF-β 전환 감소 (`EST^tgfbactGAM`)

#### 조골세포 (ROB) 역학
```
dROB1/dt = ROBin * (1/EST)^robGAM - KPT * ROB1
```
- 에스트로겐 감소 → ROB 생산 증가

#### 활성 조골세포 역학
```
dOBfast/dt = (bigDb/PicOB) * D * FracOBfast * Frackb2 - kbfast * OBfast
dOBslow/dt = (bigDb/PicOB) * D * (1-FracOBfast) * Frackb - kbslow * OBslow
```
- Osteoblast = OBfast + OBslow

#### 파골세포 역학
```
dOC/dt = kinOC2 - KLSoc * OC
```
- kinOC2는 RANK-RANKL 신호, TGF-β에 의해 조절

#### RANK/RANKL/OPG 시스템
```
dL/dt = kinL - koutL*L - k1*O*L + k2*N - k3*RNK*L + k4*M
dRNK/dt = kinRNK * TGFBact^kinRNKgam - koutRNK*RNK - k3*RNK*L + k4*M
dO/dt = pO - k1*O*L + k2*N - kO*O
dM/dt = k3*RNK*L - k4*M    (RANK-RANKL 복합체)
dN/dt = k1*O*L - k2*N       (OPG-RANKL 복합체)
```

#### PTH 역학
```
dPTH/dt = SPTH - kout * PTH
```
- SPTH는 칼슘 농도에 의해 시그모이드 함수로 조절

#### 칼시트리올/비타민D
```
dB/dt = AOH - T69 * B
dAOH/dt = SE - T64 * AOH
```

#### 칼슘 항상성
```
dP/dt = J14 - J15 - J27 + J40       (세포외 칼슘)
dT/dt = OralCa*T85 - J40              (장관 칼슘)
dR/dt = T36*(1-R) - T37*R             (신장 재흡수)
dQ/dt = J15 - J14 + J14a - J15a       (골 교환가능 칼슘)
dQbone/dt = J15a - J14a               (골 저장 칼슘)
```

#### 인산 항상성
```
dECCPhos/dt = J41 - J42 - J48 + J53 - J54 + J56
dPhosGut/dt = OralPhos*F12 - J53
dIntraPO/dt = J54 - J56
```

#### BMD (요추, Lumbar Spine)
```
dBMDls/dt = kinBMDls * (OB/OB0)^gamOB - koutBMDls * (OC/OC0)^gamOCls * BMDls
```
- kinBMDls = koutBMDls * (OC0/OC0)^gamOCls * BMDls0 / (OB0/OB0)^gamOB = koutBMDls (정상상태)
- 조골세포 증가 → BMD 증가, 파골세포 증가 → BMD 감소

#### 세포내 신호전달
```
dRX2/dt = RX2Kin - RX2Kout * RX2
dCREB/dt = crebKin - crebKout * CREB
dBCL2/dt = bcl2Kout * CREB * RX2 - bcl2Kout * BCL2
```

### 3.5 바이오마커 출력 (P1NP, CTX)

논문에서:
- **CTX** (골 흡수 마커): 파골세포(OC) 활성과 비례 → `CTX_change = (OC/OC0 - 1) * 100%`
- **P1NP** (골 형성 마커): 조골세포(OB) 활성과 비례 → `P1NP_change = (OB/OB0 - 1) * 100%`
- **BMDls** 변화율: `BMD_change = (BMDls/BMDls0 - 1) * 100%`

### 3.6 RMSE 검증 결과 (Supplement Table S2)

| 엔드포인트 | 6개월 RMSE | 6개월 최대편차 | 12개월 RMSE | 12개월 최대편차 |
|-----------|-----------|-------------|-----------|-------------|
| BMD | 0.780% | -1.46% | 0.684% | -1.25% |
| CTX | 13.10% | 33.69% | 12.34% | 22.21% |
| P1NP | 4.39% | 7.51% | 7.09% | -13.30% |

---

## 4. R 코드 구현 계획

### 4.1 기술 스택

| 도구 | 용도 |
|------|------|
| `R` (>= 4.0) | 주 프로그래밍 언어 |
| `deSolve` | ODE 수치 적분 (lsoda 솔버) |
| `ggplot2` | 시각화 |
| `tidyverse` | 데이터 처리 |

**대안**: `mrgsolve` 패키지 사용 시 C++ 기반 고속 ODE 솔버 활용 가능 (OpenBoneMin 참조)

### 4.2 파일 구조

```
CTS/
├── plan.md                      # 이 계획서
├── R/
│   ├── 01_dose_E2_model.R       # 모듈 A: 용량-E2 모델
│   ├── 02_qsp_parameters.R      # QSP 모델 파라미터 정의
│   ├── 03_qsp_ode_system.R      # QSP ODE 시스템 (29+ 미분방정식)
│   ├── 04_qsp_simulation.R      # 시뮬레이션 실행 스크립트
│   ├── 05_validation_plots.R    # 검증용 시각화
│   └── utils.R                  # 유틸리티 함수
├── output/
│   ├── figures/                  # 생성된 그래프
│   └── results/                 # 시뮬레이션 결과
└── data/
    └── observed_data.R          # 논문 Figure에서 추출한 관측 데이터
```

### 4.3 구현 단계

#### Step 1: 용량-E2 모델 (`01_dose_E2_model.R`)
- Scaled logistic function 구현
- 파라미터: slope=0.00894, logE2max=5.20, logE2min=2.14
- 용량별 E2 예측 함수 작성
- Figure 1 재현 (용량 vs. E2 곡선)

#### Step 2: QSP 파라미터 정의 (`02_qsp_parameters.R`)
- 약 90개 이상의 모델 파라미터를 named list로 정의
- 초기 조건 (약 35개 상태변수) 정의
- 에스트로겐 모듈 파라미터 포함
- 정상상태(baseline) 조건에서의 파생 변수 계산

#### Step 3: ODE 시스템 구현 (`03_qsp_ode_system.R`)
- `deSolve::ode()` 호환 형식의 ODE 함수 작성
- 약 35개 미분방정식 구현:
  - 칼슘/인산 항상성 (PTH, 칼시트리올, 칼슘, 인산)
  - 골세포 역학 (ROB, OBfast, OBslow, OC)
  - RANK/RANKL/OPG 시스템
  - TGF-β 역학
  - 세포내 신호전달 (RUNX2, CREB, BCL2)
  - BMD 역학
  - 에스트로겐 역학
- 보조 대수 방정식 (flux 계산, Hill function 등)
- 치료 ON/OFF 스위치 구현 (치료기간 & 후속관찰기간)

#### Step 4: 시뮬레이션 실행 (`04_qsp_simulation.R`)
- **시나리오 1**: Elagolix 150 mg QD, 6개월 → BMD, CTX, P1NP 변화 예측
- **시나리오 2**: Elagolix 200 mg BID, 6개월 → BMD, CTX, P1NP 변화 예측
- **시나리오 3**: 12개월 치료 + 6개월 후속관찰 (Figure 3 재현)
- **시나리오 4**: 24개월 연속 치료 (Table 3, Figure 4 재현)
- **시나리오 5**: 12개월 치료 + 12개월 후속관찰

#### Step 5: 검증 및 시각화 (`05_validation_plots.R`)
- Figure 1 재현: 용량-E2 관계
- Figure 2 재현: E2 억제율 vs. BMD/CTX/P1NP (6개월, 12개월)
- Figure 3 재현: 시간 경과에 따른 BMD/CTX/P1NP (치료+후속관찰)
- Figure 4 재현: 24개월 연속 치료 시뮬레이션
- Table 3 재현: 예측 BMD 변화율

### 4.4 주요 구현 고려사항

1. **정상상태 초기화**:
   - 약물 투여 전 모델을 충분히 오래 실행하여 정상상태 도달 확인
   - 초기값은 OpenBoneMin 참조값 사용

2. **에스트로겐 모듈 통합**:
   - Riggs (2012) 에스트로겐 확장을 기본 모델에 통합
   - `ESTON = 1`로 활성화
   - 용량-E2 모델의 출력을 EST 상태변수로 연결

3. **치료 시뮬레이션 방법**:
   - 기저선에서 EST = 1.0 (정상 에스트로겐)
   - 치료 시작: EST를 E2_predicted/E2_baseline으로 설정 (단계적 변화)
   - 치료 종료: EST를 다시 1.0으로 복귀

4. **수치 안정성**:
   - `deSolve`의 `lsoda` 적응 솔버 사용
   - 적절한 atol/rtol 설정 (1e-8 ~ 1e-10)
   - 상태변수 음수 방지 (`events` 또는 `rootfunc` 활용)

5. **시간 단위**:
   - 모델 내부: 시간(hours) 단위
   - 출력/시각화: 월(months) 단위로 변환 (1개월 = 730.5시간)

---

## 5. 예상 결과 (Table 3 기준)

| 용량 | 6개월 | 12개월 | 18개월 | 24개월 |
|------|-------|--------|--------|--------|
| 150 mg QD | -0.61% | -0.91% | -0.96% | -0.91% |
| 200 mg BID | -3.47% | -4.95% | -5.15% | -4.97% |

- 150 mg QD: BMD 감소가 약 -1% 수준에서 안정화
- 200 mg BID: BMD 감소가 약 -5% 수준에서 안정화
- 치료 중단 후 12개월 이내 기저선 근처로 회복 예측

---

## 6. 참고문헌

1. Stodtmann S, et al. *Clin Transl Sci.* 2021;14:1611-1619. (주 논문)
2. Peterson MC, Riggs MM. *Bone.* 2010;46:49-63. (기반 QSP 모델)
3. Riggs MM, et al. *CPT Pharmacometrics Syst Pharmacol.* 2012;1:e11. (에스트로겐 확장)
4. MetrumRG OpenBoneMin: https://github.com/metrumresearchgroup/OpenBoneMin (R 구현 참조)
5. BioModels BIOMD0000000613: https://www.ebi.ac.uk/biomodels/BIOMD0000000613 (SBML 모델)

---

## 7. 작업 우선순위 및 일정

| 순서 | 작업 | 의존성 |
|------|------|--------|
| 1 | 용량-E2 모델 구현 및 검증 | 없음 |
| 2 | QSP 파라미터 정의 | 없음 |
| 3 | ODE 시스템 구현 | 2 |
| 4 | 정상상태 초기화 검증 | 2, 3 |
| 5 | 에스트로겐 모듈 통합 | 1, 3, 4 |
| 6 | 시뮬레이션 실행 | 5 |
| 7 | 검증 시각화 | 6 |
