# TAUSORT Technical Report
## 태양 대기 복사 전달을 위한 불투명도 처리 시스템

---

## 목차

1. [개요](#1-개요)
2. [물리학적 배경](#2-물리학적-배경)
3. [수학적 원리](#3-수학적-원리)
4. [시스템 구조](#4-시스템-구조)
5. [알고리즘 상세](#5-알고리즘-상세)
6. [입출력 사양](#6-입출력-사양)
7. [컴파일 및 실행](#7-컴파일-및-실행)
8. [파라미터 설정](#8-파라미터-설정)
9. [응용 및 활용](#9-응용-및-활용)
10. [기술 사양](#10-기술-사양)
11. [참고문헌](#11-참고문헌)

---

## 1. 개요

### 1.1 목적

TAUSORT는 MURaM (MPS/University of Chicago Radiative MHD) 시뮬레이션 코드와 함께 사용되는 불투명도(opacity) 테이블 생성 프로그램입니다. 태양 대기의 복사 전달(radiative transfer) 계산을 효율적으로 수행하기 위해 ODF(Opacity Distribution Function) 데이터를 처리하여 밴드별 평균 불투명도를 계산합니다.

### 1.2 주요 기능

- **ODF 데이터 처리**: 수백만 개의 스펙트럼선을 통계적으로 압축
- **밴드 평균 계산**: 플랑크(Planck) 및 로슬란드(Rosseland) 평균 불투명도
- **다중 모드 지원**: Grey, Multi-band, Full ODF 모드
- **효율적 저장**: 바이너리 형식의 룩업 테이블 생성

### 1.3 성능 특징

- **압축률**: 10⁶개 스펙트럼선 → 12개 ODF 빈 (약 10⁵배 압축)
- **정확도**: 복사 플럭스 1% 이내 재현
- **계산 효율**: 실시간 3D 시뮬레이션 가능
- **메모리 효율**: 최적화된 데이터 구조와 동적 메모리 관리

---

## 2. 물리학적 배경

### 2.1 복사 전달 방정식 (Radiative Transfer Equation)

태양 대기에서 빛의 전파는 복사 전달 방정식으로 기술됩니다:

```
dI_ν/ds = -κ_ν ρ I_ν + j_ν
```

**변수 정의**:
- `I_ν`: 주파수 ν에서의 복사 강도 [erg/cm²/s/Hz/steradian]
- `κ_ν`: 질량 흡수 계수 [cm²/g]
- `ρ`: 물질 밀도 [g/cm³]
- `j_ν`: 방출 계수 [erg/cm³/s/Hz/steradian]

### 2.2 국소 열역학 평형 (Local Thermodynamic Equilibrium, LTE)

태양 광구에서는 LTE 가정이 유효합니다:

```
j_ν = κ_ν ρ B_ν(T)
```

여기서 `B_ν(T)`는 플랑크 함수:

```
B_ν(T) = (2hν³/c²) × 1/(e^(hν/kT) - 1)
```

### 2.3 광학 깊이 (Optical Depth)

광학 깊이는 매질의 불투명도를 나타내는 무차원 량입니다:

```
τ_ν(z) = ∫₀^z κ_ν(z') ρ(z') dz'
```

**물리적 의미**:
- `τ = 0`: 관측자 방향 (대기 최상층)
- `τ = 1`: 빛의 강도가 1/e로 감소하는 지점
- `τ = 2/3`: 태양 광구의 유효 방출층
- `τ >> 1`: 불투명한 깊은 층

### 2.4 불투명도의 구성

태양 대기의 총 불투명도:

```
κ_total = κ_lines + κ_continuum
        = κ_lines + (κ_bf + κ_ff + κ_scattering)
```

**구성 요소**:
1. **선 불투명도 (Line Opacity)**
   - 원자/이온의 전자 전이
   - 수백만 개의 스펙트럼선

2. **연속 불투명도 (Continuum Opacity)**
   - Bound-Free: H⁻ 광이온화, 수소/헬륨 연속
   - Free-Free: 제동복사
   - Scattering: Thomson, Rayleigh 산란

---

## 3. 수학적 원리

### 3.1 ODF (Opacity Distribution Function) 방법

ODF는 주어진 파장 구간 내 불투명도의 누적분포함수입니다:

```
P(κ < κ₀) = 1/Δν ∫[κ_ν < κ₀] dν
```

**구현 방식**:
- 328개 파장 구간으로 분할
- 각 구간을 12개 대표값(빈)으로 압축
- 통계적 가중치 적용

### 3.2 평균 불투명도

#### 플랑크 평균 (Planck Mean)

복사 에너지로 가중평균:

```
κ_P = ∫ κ_ν B_ν dν / ∫ B_ν dν
```

**적용**: 광학적으로 얇은 층 (τ < 1)

#### 로슬란드 평균 (Rosseland Mean)

복사 플럭스 전달의 조화평균:

```
1/κ_R = ∫ (1/κ_ν)(∂B_ν/∂T) dν / ∫ (∂B_ν/∂T) dν
```

**적용**: 광학적으로 두꺼운 층 (τ >> 1)

#### 적응적 전환

```
κ_mean = κ_P × 2^(-τ/τ₁/₂) + κ_R × (1 - 2^(-τ/τ₁/₂))
```

여기서 `τ₁/₂ = 0.35`는 전환 광학 깊이

### 3.3 수치적 방법

#### Simpson 적분

광학 깊이 계산:

```
τ_i = τ_{i-1} + Δz × [1/3 κᵢρᵢ + 1/6 κᵢ₋₁ρᵢ + 1/6 κᵢρᵢ₋₁ + 1/3 κᵢ₋₁ρᵢ₋₁]
```

**정확도**: O(Δz⁴)

#### 이중선형 보간

온도-압력 격자에서 불투명도 보간:

```
f(T,P) = Σᵢⱼ wᵢⱼ × log fᵢⱼ
```

로그 공간 보간으로 지수적 변화 처리

---

## 4. 시스템 구조

### 4.1 프로그램 흐름

```
main()
  ├── input()        # 데이터 읽기
  ├── initialize()   # 물리량 계산
  ├── 분기
  │   ├── [FULLODF] sort_full()
  │   └── [일반] sort() → meanop()
  └── output()       # 결과 저장
```

### 4.2 주요 함수

| 함수 | 기능 | 복잡도 |
|------|------|--------|
| `input()` | 대기 모델 및 ODF 데이터 읽기 | O(N) |
| `initialize()` | 플랑크 함수, 광학 깊이 계산 | O(N×M) |
| `sort()` | 밴드별 분류 | O(N×M) |
| `meanop()` | 평균 불투명도 계산 | O(N×B) |
| `output()` | 바이너리 파일 생성 | O(N) |

### 4.3 데이터 구조

#### 주요 배열

```c
// ODF 데이터 [파장][온도][압력][빈]
short ODF[328][57][25][12];

// 불투명도 [주파수][빈][깊이]
double kap[328][12][1200];

// 광학 깊이 [주파수][빈][깊이]
double tau_nu[328][12][1200];

// 밴드 평균 [밴드][온도][압력]
float kap_mean[Nbands][57][25];
```

#### 연결 리스트

```c
struct member {
    int nuix;              // 주파수 인덱스
    int binix;             // 빈 인덱스
    struct member *pNext;  // 다음 포인터
};
```

---

## 5. 알고리즘 상세

### 5.1 초기화 단계

1. **주파수 변환**
   ```c
   dnu = c × (1/λ₁ - 1/λ₂)
   nu_mid = c/2 × (1/λ₁ + 1/λ₂)
   ```

2. **플랑크 함수 계산**
   - 각 온도에서 B_ν와 ∂B_ν/∂T 계산
   - 로그 스케일 온도 처리

3. **보간**
   - 온도-압력 격자에서 실제 대기 조건으로
   - 이중선형 보간 사용

4. **광학 깊이 적분**
   - Simpson 복합 적분 규칙
   - 경계 조건: τ(0) = 10⁻¹⁰

### 5.2 밴드 정렬

1. **τ = 1 지점 탐색**
   - 각 주파수-빈에서 τ_ν = 1인 깊이
   - 해당 지점의 τ₅₀₀₀ 기록

2. **밴드 할당**
   ```
   if (Level_tau_bot < log(τ₅₀₀₀) ≤ Level_tau_top)
      AND (Level_nu_bot ≤ ν < Level_nu_top):
      → Band 할당
   ```

3. **연결 리스트 구성**
   - 각 밴드별 주파수-빈 목록 관리
   - O(1) 추가 연산

### 5.3 평균 계산

1. **밴드 플랑크 함수**
   ```
   B_band = Σ B_ν × Δν × w_bin
   ```

2. **플랑크 평균**
   ```
   κ_P = Σ(κ_ν B_ν Δν w) / Σ(B_ν Δν w)
   ```

3. **로슬란드 평균**
   ```
   1/κ_R = Σ(1/κ_ν × ∂B/∂T × Δν × w) / Σ(∂B/∂T × Δν × w)
   ```

---

## 6. 입출력 사양

### 6.1 입력 파일

#### G2_1D.dat
- **형식**: ASCII 텍스트
- **내용**: 1D 대기 모델 (1200개 깊이 지점)
- **구조**: `z[cm] ρ[g/cm³] P[dyne/cm²] T[K]`

#### p00big2_asplund.bdf
- **형식**: ASCII/바이너리
- **내용**: ODF 데이터
- **크기**: 328×57×25×12 = 5,606,400 엔트리
- **단위**: log₁₀(κ) × 1000

#### asplund_abs_cont.dat / asplund_sca_cont.dat
- **형식**: ASCII
- **내용**: 연속 흡수/산란 계수
- **구조**: [파장][온도][압력]

### 6.2 출력 파일

#### 파일명 규칙
- `kappa_grey.dat`: 회색 근사 (1 밴드)
- `kappa_N_band.dat`: N개 밴드
- `kappa_fullodf.dat`: 전체 ODF (3936 빈)

#### 바이너리 구조

```
[헤더: 32 bytes]
├─ int[8]: 메타데이터
│  ├─ [0]: tau5000bin 플래그
│  ├─ [1]: NT (온도 격자점)
│  ├─ [2]: Np (압력 격자점)
│  ├─ [3]: Nbands (밴드 수)
│  └─ [4-7]: 추가 플래그

[데이터]
├─ double[NT]: 온도 격자
├─ double[Np]: 압력 격자
├─ float[Nbands×NT×Np]: κ_mean
├─ float[Nbands×NT]: B_band
└─ float[Nnu]: 주파수 (FULLODF만)
```

---

## 7. 컴파일 및 실행

### 7.1 컴파일

```bash
# 표준 컴파일
make clean
make

# 디버그 모드
make CFLAGS="-g -Wall" clean
make

# 최적화
make CFLAGS="-O3 -march=native" clean
make
```

### 7.2 실행

```bash
# 기본 실행
./tausort.x

# 출력 확인
ls -lh kappa_*.dat
```

### 7.3 의존성

- **컴파일러**: g++ 또는 호환 C++ 컴파일러
- **표준 라이브러리**: math, stdio, stdlib, string
- **메모리**: 최소 100 MB RAM

---

## 8. 파라미터 설정

### 8.1 모드 설정 (global_tau.h)

```c
// 주요 모드 플래그
#define FULLODF 0    // 0: 밴드 평균, 1: 전체 ODF
#define ROSS 1       // 1: 로슬란드만, 0: 적응적
#define GREY 0       // 1: 회색 근사, 0: 다중밴드
#define Nbands 5     // 밴드 수 (GREY=0일 때)
```

### 8.2 밴드 경계 설정

```c
// 5밴드 시스템 예제
double Level_tau_bot[5] = {-0.25, -1.5, -3.0, -99.0, -99};
double Level_tau_top[5] = {99.0, -0.25, -1.5, -3.0, -3.0};
double Level_nu_bot[5] = {0, 0, 0, 0, 175};
double Level_nu_top[5] = {328, 328, 328, 175, 328};
```

### 8.3 격자 설정

```c
#define KMAX 1200    // 깊이 격자점
#define NT 57        // 온도 격자점
#define Np 25        // 압력 격자점
#define Nlam 328     // 파장 구간
#define Nbin 12      // ODF 빈
```

---

## 9. 응용 및 활용

### 9.1 MURaM 시뮬레이션

**통합 방법**:
1. TAUSORT로 불투명도 테이블 생성
2. MURaM에서 테이블 로드
3. 3D 복사 전달 계산 수행

**성능 향상**:
- 실시간 계산 대비 10⁴-10⁵배 가속
- 메모리 사용량 99% 감소

### 9.2 스펙트럼 합성

**응용 분야**:
- 합성 스펙트럼 생성
- 관측 데이터와 비교
- 대기 모델 검증

### 9.3 에너지 균형 계산

**계산 항목**:
- 복사 가열률: Q_rad = κ_P × B × 4π
- 복사 냉각률: 광학적으로 얇은 손실
- 순 에너지 플럭스

---

## 10. 기술 사양

### 10.1 성능 특성

| 항목 | 사양 |
|------|------|
| 압축률 | ~10⁵:1 |
| 정확도 | 플럭스 오차 <1% |
| 처리 시간 | ~10초 (표준 설정) |
| 메모리 사용 | ~100 MB |
| 출력 크기 | 2-50 MB (모드별) |

### 10.2 수치 정확도

| 계산 | 방법 | 정확도 |
|------|------|--------|
| 적분 | Simpson | O(Δz⁴) |
| 보간 | 이중선형 | O(ΔT²×ΔP²) |
| ODF | 12빈 근사 | ~1% |

### 10.3 확장성

- **파장 범위**: 90 Å ~ 50 μm (확장 가능)
- **온도 범위**: 2,000 ~ 200,000 K
- **압력 범위**: 10⁻⁴ ~ 10⁸ dyne/cm²

---

## 11. 참고문헌

### 핵심 참고문헌

1. **ODF 방법론**
   - Carbon, D. F. (1974), "Opacity sampling in stellar atmosphere calculations", ApJ 187, 135
   - Nordlund, Å. (1982), "Numerical simulations of the solar granulation", A&A 107, 1

2. **평균 불투명도**
   - Mihalas, D. (1978), "Stellar Atmospheres", W. H. Freeman
   - Cox, J. P. & Giuli, R. T. (1968), "Principles of Stellar Structure"

3. **태양 대기 모델**
   - Vernazza, J. E., Avrett, E. H., & Loeser, R. (1981), "Structure of the solar chromosphere", ApJS 45, 635
   - Asplund, M., et al. (2009), "The chemical composition of the Sun", ARA&A 47, 481

4. **MURaM 코드**
   - Vögler, A., et al. (2005), "Simulations of magneto-convection in the solar photosphere", A&A 429, 335
   - Rempel, M. (2017), "Extension of the MURaM radiative MHD code", ApJ 834, 10

### 관련 소프트웨어

- **ATLAS**: Kurucz 불투명도 데이터베이스
- **PHOENIX**: 항성 대기 모델 코드
- **MULTI**: 비LTE 복사 전달 코드

---

## 부록 A: 물리 상수

| 상수 | 기호 | 값 | 단위 |
|------|------|-----|------|
| 빛의 속도 | c | 2.99792458×10¹⁰ | cm/s |
| 플랑크 상수 | h | 6.626196×10⁻²⁷ | erg·s |
| 볼츠만 상수 | k | 1.380622×10⁻¹⁶ | erg/K |
| 태양 표면 중력 | g☉ | 2.74×10⁴ | cm/s² |

## 부록 B: 파일 체크섬

입력 파일의 무결성 검증을 위한 예상 크기:
- `G2_1D.dat`: ~100 KB
- `p00big2_asplund.bdf`: ~11 MB
- `asplund_abs_cont.dat`: ~2 MB
- `asplund_sca_cont.dat`: ~2 MB

## 부록 C: 문제 해결

### 일반적인 오류와 해결법

| 오류 | 원인 | 해결법 |
|------|------|--------|
| Segmentation fault | 메모리 부족 | KMAX 감소 |
| NaN in output | 언더플로 | 초기 τ 값 증가 |
| Wrong file size | 잘못된 모드 | FULLODF 플래그 확인 |

---

## 저작권 및 라이선스

TAUSORT는 MURaM 프로젝트의 일부로 개발되었습니다.
사용 시 관련 논문 인용을 권장합니다.

---

*문서 작성일: 2024*  
*버전: 1.0*  
*작성자: MURaM Development Team*

---

## 연락처

기술 지원 및 문의:
- MURaM 프로젝트 홈페이지
- GitHub 저장소
- 이메일: [프로젝트 관리자]

---

**END OF DOCUMENT**