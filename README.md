# 의과대학 4학년 QSP 캡스톤 프로젝트 (2025)

## 프로젝트 개요

본 프로젝트는 가톨릭대학교 의과대학 4학년 학생들이 **Quantitative Systems Pharmacology (QSP)** 모델링을 학습하고, 각자의 연구 주제에 맞는 웹 기반 시뮬레이션 앱을 개발하는 것을 목표로 합니다.

**최종 산출물**: [Merigolix QSP 대시보드](https://pipetqsp.shinyapps.io/merigolix/)와 같은 R Shiny 기반 웹 애플리케이션

---

## 연구 주제 및 참여자

| 참여자 | GitHub ID | 연구 주제 | 참고 논문 |
|--------|-----------|-----------|-----------|
| 박광원 | pangpangwon | Elagolix 임상 3상 데이터 기반 칼슘 항상성/골대사 QSP 모델 검증 | [Validation of a QSP model of calcium homeostasis](https://pubmed.ncbi.nlm.nih.gov/33963686/) |
| 박정민 | wjdals011106 | 전이성 삼중음성유방암(mTNBC) 면역항암제 반응 예측 QSP 모델 | [A transcriptome-informed QSP model of mTNBC](https://www.science.org/doi/10.1126/sciadv.adg0289) |
| 배정현 | john83072 | 알츠하이머병 Aβ 플라크 표적 치료제 임상실패 원인 분석 QSP 모델 | [In silico analysis for reducing Aβ plaque](https://pubmed.ncbi.nlm.nih.gov/33938131/) |
| 한상하 | magen2001 | 조현병 항정신병약 반응 예측 QSP 모델 | [Computer-based mechanistic schizophrenia disease model](https://doi.org/10.1371/journal.pone.0049732) |

---

## 3개월 연구 일정 (12주)

### Phase 1: 기초 학습 및 모델 이해 (1-4주)

#### 1주차: 프로젝트 킥오프 및 환경 설정
- [ ] 개인 GitHub 저장소 생성 및 프로젝트 구조 설정
- [ ] R/RStudio 및 필수 패키지 설치 (deSolve, shiny, ggplot2, plotly)
- [ ] 참고 논문 정독 및 핵심 내용 정리
- [ ] 연구 배경 및 목적 문서화

#### 2주차: QSP 모델링 기초 학습
- [ ] ODE (Ordinary Differential Equations) 기초 복습
- [ ] deSolve 패키지를 이용한 간단한 PK/PD 모델 구현 실습
- [ ] 참고 논문의 수학적 모델 구조 분석
- [ ] 모델 변수 및 파라미터 목록 작성

#### 3주차: 참고 모델 심층 분석
| 참여자 | 분석 내용 |
|--------|-----------|
| 박광원 | HPG축-칼슘 항상성-골대사 네트워크 구조, E2/PTH/비타민D 상호작용 |
| 박정민 | 종양 면역 미세환경 구성요소, PD-1/PD-L1 경로, 바이오마커 |
| 배정현 | Aβ 생성-응집-제거 경로, BACE/γ-secretase 억제 기전, ADCP |
| 한상하 | 도파민-세로토닌 수용체 경쟁 모델, MSN 모델, PANSS 예측 |

#### 4주차: 모델 구현 준비
- [ ] ODE 시스템 수식 정리 (LaTeX 형식)
- [ ] 초기 조건 및 파라미터 값 문헌 조사
- [ ] 모델 검증용 임상 데이터 수집
- [ ] Phase 1 중간 발표 및 피드백

---

### Phase 2: 모델 구현 및 검증 (5-8주)

#### 5주차: 핵심 ODE 시스템 구현
- [ ] R에서 ODE 시스템 코딩
- [ ] 기본 시뮬레이션 실행 및 디버깅
- [ ] 시간에 따른 주요 변수 변화 그래프 생성

#### 6주차: 약물 효과 모듈 추가
| 참여자 | 구현 내용 |
|--------|-----------|
| 박광원 | Elagolix PK 모델 + GnRH 수용체 길항 효과 |
| 박정민 | Pembrolizumab PK + PD-1 수용체 점유율 모델 |
| 배정현 | 항-Aβ 항체 (Aducanumab 등) PK + 플라크 결합 모델 |
| 한상하 | 항정신병약 PK + 다중 수용체 점유율 모델 |

#### 7주차: 모델 캘리브레이션
- [ ] 문헌 데이터와 시뮬레이션 결과 비교
- [ ] 파라미터 민감도 분석 수행
- [ ] 필요시 파라미터 조정 및 최적화
- [ ] 모델 예측값 vs 임상 관측값 그래프 생성

#### 8주차: 모델 검증 및 문서화
- [ ] 검증 결과 정리 (정확도, 한계점)
- [ ] 모델 구조도 및 플로우차트 작성
- [ ] Phase 2 중간 발표 및 피드백
- [ ] 코드 리뷰 및 정리

---

### Phase 3: Shiny 웹앱 개발 (9-12주)

#### 9주차: Shiny 기초 및 UI 설계
- [ ] Shiny 기초 학습 (reactive programming)
- [ ] [Merigolix 앱](https://pipetqsp.shinyapps.io/merigolix/) 구조 분석
- [ ] UI 와이어프레임 설계
- [ ] 사용자 입력 패널 구현 (용량, 투여 간격, 시뮬레이션 기간)

#### 10주차: 핵심 기능 구현
| 구성요소 | 기능 |
|----------|------|
| 입력 패널 | 약물 용량, 투여 빈도, 환자 특성 |
| 출력 패널 | PK/PD 프로파일, 바이오마커 변화, 임상 지표 |
| 시각화 | 인터랙티브 그래프 (plotly), 데이터 테이블 |

#### 11주차: 고급 기능 및 최적화
- [ ] 가상 환자 시뮬레이션 기능 추가
- [ ] 용량 최적화 도구 구현
- [ ] 앱 성능 최적화
- [ ] 사용자 가이드 작성

#### 12주차: 배포 및 최종 발표
- [ ] shinyapps.io 배포
- [ ] 최종 버그 수정 및 테스트
- [ ] 최종 발표 자료 준비
- [ ] 프로젝트 보고서 작성

---

## 웹앱 주요 기능 명세

### 공통 기능
```
1. 시뮬레이션 설정
   - 약물 용량 (mg)
   - 투여 빈도 (QD/BID/QW 등)
   - 시뮬레이션 기간 (일/주)

2. PK/PD 프로파일
   - 혈중 약물 농도 시간 곡선
   - 타겟 점유율/억제율

3. 바이오마커 변화
   - 질환별 핵심 바이오마커 시계열 그래프

4. 임상 지표 예측
   - 질환별 임상 평가 지표
```

### 개인별 특화 기능

| 참여자 | 바이오마커 | 임상 지표 |
|--------|-----------|-----------|
| 박광원 | E2, PTH, CTX, P1NP | BMD 변화율 (%) |
| 박정민 | TIL, APC 밀도, 암클론 다양성 | 반응률, 반응 지속 기간 |
| 배정현 | CSF Aβ, 플라크 부하 (PET-SUVR) | 플라크 감소율 (%) |
| 한상하 | D2 점유율, 5-HT2A 점유율 | PANSS 점수 변화, EPS 위험도 |

---

## 기술 스택

- **모델링**: R, deSolve, mrgsolve
- **웹앱**: Shiny, shinydashboard, bslib
- **시각화**: ggplot2, plotly
- **배포**: shinyapps.io
- **버전관리**: Git, GitHub

---

## 디렉토리 구조

```
capstone-qsp/
├── README.md
├── parkgwangwon/          # 박광원: 칼슘 항상성/골대사
│   ├── model/
│   │   └── elagolix_bone_qsp.R
│   ├── shiny/
│   │   ├── app.R
│   │   ├── ui.R
│   │   └── server.R
│   └── docs/
│       └── model_description.md
├── parkjungmin/           # 박정민: mTNBC 면역항암제
│   ├── model/
│   │   └── mtnbc_io_qsp.R
│   ├── shiny/
│   └── docs/
├── baejunghyun/           # 배정현: 알츠하이머 Aβ
│   ├── model/
│   │   └── alzheimer_abeta_qsp.R
│   ├── shiny/
│   └── docs/
└── hansangha/             # 한상하: 조현병 항정신병약
    ├── model/
    │   └── schizophrenia_qsp.R
    ├── shiny/
    └── docs/
```

---

## 주간 미팅 일정

- **매주 금요일 오후 4시**: 진행 상황 공유 및 피드백
- **월 1회**: 전체 발표 및 코드 리뷰

---

## 평가 기준

| 항목 | 비중 | 내용 |
|------|------|------|
| 모델 정확성 | 30% | 문헌 데이터와의 일치도 |
| 코드 품질 | 20% | 가독성, 문서화, 재현성 |
| 앱 완성도 | 30% | UI/UX, 기능 구현 |
| 발표 및 문서 | 20% | 이해도, 설명력 |

---

## 참고 자료

### QSP 모델링 기초
- [QSP Modeling Workflow (CPT:PSP)](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12071)
- [deSolve Package Documentation](https://cran.r-project.org/web/packages/deSolve/)

### Shiny 개발
- [Mastering Shiny (Hadley Wickham)](https://mastering-shiny.org/)
- [Shiny Dashboard](https://rstudio.github.io/shinydashboard/)

### 예시 앱
- [Merigolix QSP Dashboard](https://pipetqsp.shinyapps.io/merigolix/)

---

## 문의

- **지도교수**: 한성필 (PIPET 가톨릭대학교 계량약리학연구소)
- **연구실**: 가톨릭대학교 성의교정 약리학교실

---

*Last updated: 2026-01*
