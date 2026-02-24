# Photosphere `parameters.dat` 기반 경계조건 해석 보고서

## 0. 범위와 기준

이 문서는 다음 파일을 기준으로, Photosphere 설정에서 경계조건이 코드에서 어떻게 실제 실행되는지 정적 추적으로 해석한 기술 보고서이다.

- 입력 파라미터: `python_codes/SetupNew/Photosphere/parameters.dat`
- 해석 범위: MHD 물리 경계 + MPI 고스트/주기 교환 + TVD 경계층 + divB 경계 처리 + RT 경계
- 코드 수정: 없음

핵심 가정:

1. 빌드 플래그는 `FFT_MODE=FFTW` (`Make_defs:8`)이다.
2. 수직축은 `x` 방향(index 0)이다 (`clean_divB.C:198` 코멘트, `Physics.g = (-g,0,0)` in `physics.H:201-203`).
3. `backup.dat`는 `0 0 -1 -1`이며, 첫 경계 호출에서 `p_bc`, `s_bc`가 자동 계산된다.


## 1. 입력 파라미터의 런타임 고정

### 1.1 파싱과 저장 구조

`Initialize()`에서 `parameters.dat`를 읽어 `Run`, `Grid`, `Physics` 구조체로 주입한다 (`init.C:17-208`).

- 격자/주기성:
  - `periods = 0 1 1` (`parameters.dat:7`)
  - 파싱: `init.C:69`
  - 저장: `Grid.periods[3]` (`grid.H:32`)
- 경계 파라미터:
  - `bnd_top=1`, `bnd_pot=1`, `bnd_bcrit=1e10`, `bnd_eps_top=0` (`parameters.dat:57-61`)
  - 파싱: `init.C:113-116`
  - 인덱스: `i_bnd_top..i_bnd_fem` (`physics.H:25-29`)
- TVD 경계층:
  - `tvd_h_bnd`, `tvd_visc_bnd`, `tvd_eta_bnd` (`parameters.dat:83-85`)
  - 파싱: `init.C:132-134`
- RT:
  - `rt_*` (`parameters.dat:108-118`)
  - 파싱: `init.C:151-160`

`getvar_s()`로 읽은 항목은 누락 시 abort, `getvar()` 항목은 누락 시 기본값 유지이다 (`dfparser.C:59-89`).  
즉 `bnd_top`, `bnd_pot`, `bnd_bcrit`는 필수이고, `bnd_eps_top`는 선택 입력이다.

### 1.2 현재 설정에서 고정되는 수치

`parameters.dat` 값으로부터:

- 도메인: `Lx=3.072e8 cm`, `Ly=Lz=9.216e8 cm`
- 격자: `Nx=96`, `Ny=Nz=192`
- 격자 간격:
  - `dx = 3.2e6 cm`
  - `dy = dz = 4.8e6 cm`
  - `dx/dy = 0.6666667`

이 비율은 potential 경계 커널의 수직 감쇠 스케일(`scale=dx_vertical/dx_horizontal`)과 직접 연결된다 (`potential_sd_heffte_kernels.C:181`).

### 1.3 `backup.dat`와 `p_bc`, `s_bc`

`backup.dat` 포맷은 `iter time p_bc s_bc`다 (`boundary_pdmp_1_fftw3.C:44-59`).  
복원 시 `ReadBackupFile(..., bc=1)`로 `p_bc`, `s_bc`를 로드한다 (`io_xysl.C:330-331`).

현재 `backup.dat` 값:

```
0 0 -1 -1
```

따라서 첫 경계 호출에서 `set_p_s_bc()`가 실행되어 초기 `p_bc`, `s_bc`를 내부 상태로부터 자동 산출한다 (`boundary_pdmp_1_fftw3.C:117-121`, `496-564`).


## 2. 타임스텝 내 경계 호출 타이밍

### 2.1 실행 루프의 경계 호출 지점

`ComputeSolution()`에서 경계 호출 순서 (`solver.C`):

1. 시작 직후 1회:
   - `exchange_grid_acc(...,1)` 후 `SetBoundaryConditions(..., stage=0, pt_update=1)` (`solver.C:80-82`)
2. 각 RK stage마다:
   - `exchange_grid_acc(...,1)` 후 `SetBoundaryConditions(..., stage, pt_update=0)` (`solver.C:316-321`)
3. 스테이지 완료 후 최종 보정:
   - `exchange_grid_acc(...,0)` 후 `SetBoundaryConditions(..., stage=maxstage+1, pt_update)` (`solver.C:352-364`)

`maxstage=4` (`src_int_tck.H:5`)이므로 1 iteration당 경계 호출 수는:

- stage 호출 4회 + 최종 1회 = 5회
- 초기화 호출 포함하면 첫 iteration 시작 전 추가 1회

### 2.2 `bnd_pot=1`에서 potential 업데이트 빈도

potential 업데이트 조건:

```cpp
if( (POTENTIAL_BC == 1) or (POTENTIAL_BC > 1 and pt_update == 1) )
```

(`boundary_pdmp_1_fftw3.C:336`)

현재 `bnd_pot=1`이므로 `SetBoundaryConditions()`가 불릴 때마다 potential 재계산이 실행된다.


## 3. 수평 주기 경계 vs 수직 물리 경계

### 3.1 격자 관점

- 전역 하부/상부 여부: `is_gbeg[d]`, `is_gend[d]` (`grid.C:211-256`)
- 주기성: `periods[d]`는 MPI Cartesian topology에 직접 주입 (`grid.C:154`)

현재 `periods=0 1 1`이므로:

- `d=0 (x)` 비주기
- `d=1,2 (y,z)` 주기

### 3.2 MPI 교환 분기

`exchange_grid_acc()`에서 고스트 채우기 조건:

- 우측(상단) 고스트:
  - `if( !Grid.is_gend[d1] or Grid.periods[d1] )` (`exchange.C:711`)
- 좌측(하단) 고스트:
  - `if( !Grid.is_gbeg[d1] or Grid.periods[d1] )` (`exchange.C:744`)

즉 주기 방향(`periods[d1]=1`)은 물리 경계에서도 이웃 rank와 래핑 교환을 수행하고,  
비주기 방향(`periods[d1]=0`)의 전역 경계면은 MPI 교환으로 고스트를 채우지 않으며 이후 `SetBoundaryConditions()`가 물리 경계값을 넣는다.

결론:

- `x` 방향: 물리 경계 적용(하부/상부)
- `y,z` 방향: 주기 경계 적용


## 4. 하부 경계(`is_gbeg[0]`) 이산식

코드 본체: `boundary_pdmp_1_fftw3.C:160-331`

### 4.1 상태 변수 정의

코드의 보존 변수:

\[
U=(\rho,\rho v_x,\rho v_y,\rho v_z,e_t,B_x,B_y,B_z)
\]

(`mhd3d.H:24-27`)

내부에너지 밀도:

\[
e_i = e_t - \frac{1}{2}\frac{M^2}{\rho}
,\quad
\epsilon = \frac{e_i}{\rho}
\]

(`boundary_pdmp_1_fftw3.C:183-188`, `250-254`)

### 4.2 평균압/요동 분해와 `p_bc` 외삽

하부 첫 두 내부층에서

\[
p_1=\mathrm{EOS}_p(\epsilon_1,\rho_1),\quad p_2=\mathrm{EOS}_p(\epsilon_2,\rho_2)
\]

의 수평 평균 \( \bar p_1,\bar p_2 \)를 구한다 (`216-224`).

경계 목표 압력 `p_bc`를 쓰는 외삽 계수:

\[
c_p=\frac{p_{bc}}{\sqrt{\bar p_1\bar p_2}}
\]

(`273`)

평균장 외삽:

\[
p_{g1}^{(mean)}=\bar p_1 c_p,\quad
p_{g2}^{(mean)}=\bar p_1 c_p^2
\]

(`274-275`)

여기에 1층 압력요동 \(p'_1=p_1-\bar p_1\)을 감쇠 계수 \(p1\_dmp=0.8\)로 추가:

\[
p_{g1}\leftarrow p_{g1}^{(mean)} + p'_1 p1\_dmp
\]
\[
p_{g2}\leftarrow p_{g2}^{(mean)} + p'_1 p1\_dmp^2
\]

(`277-278`)

정규자기압 보정:

\[
p_{mag}=\frac{1}{2}B_x^2
\]
\[
p_{g1}\leftarrow p_{g1} - p_{mag}(1-p1\_dmp),\;
p_{g2}\leftarrow p_{g2} - p_{mag}(1-p1\_dmp^2)
\]

(`256`, `280-281`)

### 4.3 엔트로피 선택 규칙

\[
s_1,s_2=
\begin{cases}
s_{bc}, & M_x\ge 0 \\
s(\epsilon,\rho), & M_x<0
\end{cases}
\]

(`264-270`)

주의: 하부 경계에서 \(M_x\ge 0\)는 도메인 유입(upflow) 조건이므로, 유입 엔트로피를 `s_bc`로 고정하고 유출은 내부값을 통과시키는 전형적인 convection-driving 하부 경계다.

### 4.4 강자기장 하부 속도 억제

`bnd_bcrit`(Gauss)를 코드 임계값으로 변환:

\[
B_{crit,code}=\frac{bnd\_bcrit}{\sqrt{8\pi}}
\]

(`92`)

억제 계수:

\[
\xi=\frac{B_x}{B_{crit,code}},\quad
c_m=\frac{1-\xi^4}{1+\xi^4}
\]
\[
\mathbf{M}\leftarrow c_m\,\mathbf{M}
\]

(`284-289`)

현재 설정 수치:

- \(\sqrt{8\pi}=5.013256549\)
- \(B_{crit,code}=1.994711\times10^9\)

### 4.5 EOS 역변환으로 ghost 상태 구성

\[
\rho_{g} = d3\_interp(p_g,s_g),\quad
\epsilon_g = eps3\_interp(p_g,s_g)
\]
\[
e_{t,g}=\epsilon_g \rho_g + \frac{1}{2}\frac{M_g^2}{\rho_g}
\]

(`293-303`)

그리고 `lbeg-1`, `lbeg-2` ghost에 저장 (`304-305`).

### 4.6 `p_bc`,`s_bc` 자동 산출식

`set_p_s_bc()` (`496-564`)에서:

- \( \bar p_1, \bar p_2 \) 수평 평균 계산
- 초기값이 음수일 때

\[
p_{bc}=\bar p_1\sqrt{\frac{\bar p_1}{\bar p_2}}
\]

(`551-553`)

\[
s_{bc}=
\begin{cases}
\langle s\rangle_{M_x>0}, & \#(M_x>0)>0\\
\langle s\rangle_{all}, & \text{otherwise}
\end{cases}
\]

(`554-559`)

이 값은 `XCOL_COMM`으로 방송되어 수직 컬럼 전체가 동일 `p_bc,s_bc`를 사용한다 (`562-563`).


## 5. 상부 경계(`is_gend[0]`) 이산식

코드 본체: `boundary_pdmp_1_fftw3.C:385-487`

상부 내부 셀 \(U_0\), 1/2 ghost \(U_1,U_2\)를 구성한다.

### 5.1 속도(법선 모멘텀) 규칙: `bnd_top`

현재 `bnd_top=1` (half-open):

- \(M_x<0\) (상부에서 유입)일 때 부호 반전으로 유입 차단
- \(M_x\ge 0\) (유출)은 통과

(`415-419`)

`bnd_top=0`이면 항상 반전되어 닫힌 경계처럼 동작 (`420-423`).

### 5.2 내부에너지 상부 고정: `bnd_eps_top`

현재 `bnd_eps_top=0`이라 비활성 (`425` 조건 미충족).  
활성일 때는 ghost에

\[
e_t = \epsilon_{top}\rho + \frac{1}{2}\frac{M^2}{\rho}
\]

을 강제한다 (`425-430`).

### 5.3 자기장 규칙: `bnd_pot`

현재 `bnd_pot=1`이므로 potential 외삽 결과 `b_ext`를 사용:

- top 내부 셀 \(U_0\): \(B_y,B_z\)만 갱신
- ghost 1층 \(U_1\): \(B_x,B_y,B_z\) 모두 갱신
- ghost 2층 \(U_2\): \(B_x,B_y,B_z\) 모두 갱신

인덱스 매핑 (`432-442`):

- \(U_0.B_y \leftarrow b\_ext[2]\), \(U_0.B_z \leftarrow b\_ext[5]\)
- \(U_1.B_x \leftarrow b\_ext[0]\), \(U_1.B_y \leftarrow b\_ext[3]\), \(U_1.B_z \leftarrow b\_ext[6]\)
- \(U_2.B_x \leftarrow b\_ext[1]\), \(U_2.B_y \leftarrow b\_ext[4]\), \(U_2.B_z \leftarrow b\_ext[7]\)

`bnd_pot=0`이면 접선성분 반대칭(odd) 처리로 대체 (`443-448`).


## 6. `bnd_pot=1` potential-field 경계의 푸리에 해석

### 6.1 실행 경로

`SetBoundaryConditions()`는 상부 내부면 \(B_x\)를 추출해 `bz0`에 넣고 (`338-349`),  
`potential_ext_fftw3()` 또는 병렬판을 호출해 `b_ext`를 만든다 (`355-377`).

### 6.2 수평 주기 가정과 커널 컨볼루션

FFTW 경로는 `PSF-kernel-*.dat`를 읽어 FFT 공간에서 컨볼루션한다 (`potential_sd_fftw3.C:77-229`).

수학적 형태는 HeFFTe 분석식 구현과 같다 (`potential_sd_heffte_kernels.C:175-224`):

\[
\hat B(\mathbf{k},z) \propto e^{-|\mathbf{k}|z}
\]

\[
k_x = \frac{2\pi n_x}{N_x},\quad
k_y = \frac{2\pi n_y}{N_y},\quad
|\mathbf{k}|=\sqrt{k_x^2+k_y^2}
\]

감쇠 계수:

\[
e_0=\frac{e^{-s k^2}}{N_xN_y},\;
e_1=\frac{e^{-|\mathbf{k}|\,\mathrm{scale}-s k^2}}{N_xN_y},\;
e_2=\frac{e^{-2|\mathbf{k}|\,\mathrm{scale}-s k^2}}{N_xN_y}
\]

여기서 `scale = dx_vertical/dx_horizontal` (`181`), `s=0.25` (`178`).

현재 설정에서는:

- 수평 isotropic 조건 `dy==dz`를 만족해야 하며 (`potential_sd_fftw3.C:69-73`)
- `scale = dx/dy = 0.6666667`

즉 수직 감쇠 길이가 수평 격자보다 짧은 비등방 격자 보정이 들어간다.

### 6.3 경계 적용 물리 의미

이 처리는 상부에서 \(\nabla\times \mathbf{B}=0\) (무전류 potential) 가정으로,  
주어진 상부면 \(B_x\)를 경계조건으로 하는 외부 자기장을 재구성한다.  
결과적으로 상부에서 인공 전류 축적을 줄이고 개방 자기장 구조를 허용한다.


## 7. TVD 경계층 파라미터 영향 (`tvdlimit_SR.C`)

코드 본체: `tvdlimit_SR.C:50-695`

### 7.1 경계 마스크 두께 `tvd_h_bnd`

`tvd_h_bnd<1`이면 도메인 높이 정규화 두께로 해석 (`74-86`).

현재 `tvd_h_bnd = 0.01, 0.01`이므로:

- 하부/상부 각각 `2*h_bnd = 0.02` 구간에서 마스크가 작동
- 물리 길이: \(0.02L_x=6.144\times10^6\) cm = 61.44 km

마스크:

\[
h_{ft}, h_{fb} \in [0,1]
\]

를 포물선 프로파일로 만든다 (`137-149`).

### 7.2 `tvd_visc_bnd=tvd_eta_bnd=1`의 효과

현재 값:

- `visc_coeff = 1`, `visc_slope = 0`
- `eta_coeff = 1`, `eta_slope = 0`

(`50-58`)

`slope=0`이므로 경계 근방 slope 추가 감쇠 블록은 비활성 (`400-415`).  
또한 `tvd_coeff`가 전부 1이라 `need_tvd_coeff=0` (`121-130`), 계수 재스케일 블록도 실행되지 않는다 (`470-487`).

결론: 현재 설정에서는 경계 점성/저항 계수 강화가 사실상 중립이다.

### 7.3 `tvd_Qdiff_bnd=0`의 상부 가열 영향

`qft` 정의:

- 상부 경계층 내부: `qft = Qdiff_bnd`
- 그 외: `qft = 1`

(`145-149`)

현재 `Qdiff_bnd=0`이므로 상부 경계층에서:

- 저항 가열 항 \(qres\)이 `*qft`로 억제 (`618-623`)
- 점성 관련 에너지 항도 `qft` 경로로 억제/보정 (`549`, `572-576`, `665`)

즉 상부 경계 근처 수치적 가열 누적을 의도적으로 줄인다.


## 8. divB 클리닝의 수직 경계 처리

코드 본체: `clean_divB.C:35-280`

핵심:

1. `exchange_B_acc()`로 B 고스트 동기화 (`119`)
2. \(\nabla\cdot B\) 계산 후 보조장 \(\phi\) 반복 업데이트 (`125-185`)
3. 수직 경계에서 \(\phi\) 반대칭 조건:

\[
\phi_{-1}=-\phi_0,\;\phi_{-2}=-\phi_{+1}
\]
\[
\phi_{+1}=-\phi_0,\;\phi_{+2}=-\phi_{-1}
\]

(`200-224`)

4. \(\nabla\phi\)로 B 교정, \(\phi\,\nabla\cdot B\)를 에너지에 반영 (`252-269`)

이 경계 반대칭은 수직 경계에서 divB 클리닝 파동의 반사/흡수 성질을 제어해,  
경계 근방 비물리적 \(\nabla\cdot B\) 재증폭을 완화한다.


## 9. RT 경계조건

코드 본체: `rt.cc`

### 9.1 하부 확산근사 경계

하부 경계 rank(`isgbeg[0]==1`)에서 incoming intensity 버퍼를 Planck 함수 \(B\)로 설정 (`1007-1024`):

\[
I_{in,\;bottom} \leftarrow B(T)
\]

이는 광학적으로 두꺼운 하부 경계의 확산근사 경계와 같다.

### 9.2 상부 무입사 경계

상부 경계 rank(`isgend[0]==1`)에서 downwind incoming intensity를 0으로 강제 (`1027-1040`):

\[
I_{in,\;top}=0
\]

즉 외부 복사 유입이 없는 개방 상부 복사 경계다.

### 9.3 `rt_tau_min` 기반 \(Q_{rad}\) 스케일 제한

`calc_Qtot_and_Tau()`에서:

\[
\tau_{\min} = (\texttt{rt\_tau\_min})^2
\]
\[
\text{scale}=\frac{\tau^2}{\tau^2+\tau_{\min}}\cdot \text{TRSW}
\]
\[
Q_{tot}=Q_t\cdot \text{scale}
\]

(`1121`, `1157-1163`)

현재 `rt_tau_min=1e-8`이므로 \(\tau\ll 10^{-8}\) 극한에서 복사 가열/냉각 기여를 부드럽게 제한한다.


## 10. 현재 파라미터에서 활성/비활성 분기 요약

| 항목 | 값 | 코드 분기 | 상태 |
|---|---:|---|---|
| 수직 주기성 | `periods[0]=0` | `exchange.C:711,744` | 물리 경계 활성 |
| 수평 주기성 | `periods[1]=periods[2]=1` | `exchange.C:711,744` | 주기 경계 활성 |
| 상부 유체 경계 | `bnd_top=1` | `boundary_pdmp_1_fftw3.C:415-423` | half-open 활성 |
| 상부 자기 경계 | `bnd_pot=1` | `boundary_pdmp_1_fftw3.C:336` | potential 업데이트 매 호출 |
| 하부 강자기장 속도 억제 | `bnd_bcrit=1e10` | `boundary_pdmp_1_fftw3.C:284-289` | 활성 |
| 상부 hot plate | `bnd_eps_top=0` | `boundary_pdmp_1_fftw3.C:425` | 비활성 |
| 전도(Spitzer) 경계 | `param_spitzer=0` | `boundary_pdmp_1_fftw3.C:307-311,454-467` | 비활성 |
| ambipolar 경계 | `param_ambipolar=0` | `boundary_pdmp_1_fftw3.C:313-317,469-473` | 비활성 |
| `bnd_fem` | 기본 0, 사용처 없음 | 검색 결과 `init/physics` 출력만 | 비활성 |
| TVD 경계계수 강화 | `tvd_visc_bnd=tvd_eta_bnd=1` | `tvdlimit_SR.C:50-58,400-487` | 실질 중립 |
| 상부 TVD 가열 억제 | `tvd_Qdiff_bnd=0` | `tvdlimit_SR.C:145-149,549,618` | 활성 |
| RT 하부/상부 경계 | `rt_type=0` | `rt.cc:1007-1040` | 활성 |


## 11. 검증 시나리오 결과

### 11.1 분기 검증

`parameters.dat`의 각 경계 파라미터가 `if` 조건으로 연결되는 위치를 모두 추적했고, 현재 값 기준 활성/비활성 표와 일치함을 확인했다.

### 11.2 단위/스케일 검증

- 격자: `dx=3.2e6`, `dy=dz=4.8e6`로 potential 경계의 `dy==dz` 제약을 만족.
- `bnd_bcrit`는 코드에서 `sqrt(8*pi)`로 변환되어 비교됨 (`boundary_pdmp_1_fftw3.C:92`).

### 11.3 주기성-잠재장 일관성

MPI 교환에서 `y,z` 주기 래핑이 보장되고(`exchange.C`), potential 외삽도 수평 주기 FFT를 전제로 한다(`potential_sd_fftw3.C`). 두 가정은 일관된다.

### 11.4 타이밍 검증

`SetBoundaryConditions()` 호출 위치를 기준으로 per-iteration 5회 호출이 성립하고, `bnd_pot=1`에서는 매 호출 potential 업데이트가 실행되는 흐름과 일치한다.

### 11.5 RT 독립 계층 검증

RT 경계(`rt.cc`)는 복사강도 버퍼 레벨에서 적용되고, MHD 경계(`boundary_pdmp_1_fftw3.C`)와 별도 루틴/호출 단계로 분리되어 동작한다.


## 12. 구현상 주의점(운영 관점)

1. `bnd_pot=1` 실행에는 `PSF-kernel-*.dat` 커널 파일이 런타임 경로에 필요하다 (`potential_sd_fftw3.C:78-85`, `140`).
2. `backup.dat`에 `p_bc,s_bc`가 저장되므로 재시작 런마다 하부 구동 상태가 연속된다 (`WriteBackupFile`, `ReadBackupFile`).
3. `param_ambipolar>0`는 현재 코드에서 강제 abort 처리되어 실제로 사용할 수 없다 (`init.C:105-111`).

