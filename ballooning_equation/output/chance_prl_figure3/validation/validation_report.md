# Chance PRL Figure 3 재현 테스트 — 고 beta 검증 실패

논문의 Eq. (2) bean 경계와 p/q 프로파일을 native GS → direct geometry → ballooning 경로에서 실제 계산했다. 첫 안정 경계는 논문의 경계 근방이지만, 높은 beta의 안정 분류는 일치하지 않는다. 따라서 현재 workflow가 Figure 3 전체를 재현했다고 판단할 수 없다.

## 논문 입력

A=1, B=0.6, R/a=4, b/a=1.386. C로 inboard indentation d/(2a)를 조절하고 E, xbar를 다시 정한다.

```text
rho=A(1+B cos t), gamma=C sin t
R=xbar+rho cos gamma, Z=E rho sin gamma
u=psiN (axis 0, edge 1), p=p0(1-u^2)^2
q=1.03+0.84375u-1.1775u^2+3.50375u^3
beta=2 mu0 integral(p dV)/integral((Bp^2+Bphi^2) dV), dV=2 pi R dR dZ
```

이 q는 q(0)=1.03, q(1)=4.2, q_u(0)=0.84375, q_u(1)=9를 만족한다. p/q의 변수는 sqrt(psiN)이 아니다. 안정 지도에는 central beta 대신 계산한 total-field volume beta를 사용했다.

## 실제 결과

| d/(2a) | volume beta [%] | GS grid / boundary points | 계산 분류 | 논문 분류 | lambda observed | max q relative error |
|---:|---:|:---|:---|:---|---:|---:|
| 0.110 | 0.804216 | 48 x 64 / 128 | stable | stable | -4.44516894e-05 | 1.733e-04 |
| 0.110 | 4.981992 | 48 x 64 / 128 | unstable | unstable | +5.74630724e-03 | 1.165e-04 |
| 0.325 | 4.979633 | 48 x 64 / 128 | stable | stable | -2.76375022e-05 | 1.023e-03 |
| 0.250 | 0.805326 | 48 x 64 / 128 | stable | stable | -3.42894084e-05 | 9.153e-05 |
| 0.250 | 4.979705 | 48 x 64 / 128 | unstable | unstable | +2.22496481e-03 | 7.761e-05 |
| 0.250 | 14.350477 | 64 x 96 / 384 | unstable | stable | +5.77300333e-03 | 6.584e-05 |
| 0.325 | 14.315163 | 64 x 96 / 384 | unstable | stable | +1.01803076e-03 | 2.188e-04 |
| 0.250 | 14.350779 | 96 x 128 / 768 | unstable | stable | +5.85742087e-03 | 6.868e-05 |

d/(2a)=0.25의 첫 경계는 beta **[2.498114, 2.514273]%**에 bracket되었다. 논문 그림의 약 2%대 첫 경계와 가깝다. 이 경계는 고유값 부호의 bisection으로 얻었으며 reference curve에 맞춰 정하지 않았다.

stable은 요청한 표면/phase에서 안정이라는 뜻이다. 표면 0.05~0.95, 보통 8~16 phase를 조사했다. 가장 촘촘한 high-beta 점의 양의 고유값은 psiN=0.90, theta0=0에서 확인한 instability witness이며 전체 최대 고유값의 하한이다. 해당 두 표면 0.90/0.95는 모든 geometry quality gate를 통과했다.

## 높은 beta의 차이는 단순한 coarse eigenvalue 문제로 해소되지 않았다

- 고 beta에서 24 cells/turn은 음의 고유값을 주지만, 48 이상에서는 양의 값이 나온다. 768→1536 refinement에서도 양의 값이 유지된다.
- Fourier 32~48에서 high-beta 양의 고유값이 거의 유지된다. 작은 mode 수를 선택해 논문의 stable label에 맞추지 않았다.
- GS 64x96 / boundary 384에서 96x128 / boundary 768로 높여도 volume beta와 양의 고유값이 유지된다.
- finest native GS residual 2.622e-07, prescribed-q max relative error 6.868e-05. psiN=0.90에서 lambda=+5.857420873e-03.
- 동일한 native bean 행렬에서 dense와 sparse largest eigenvalue 차이 3.464e-15. largest를 nearest-to-zero mode로 잘못 고른 결과가 아니다.

GS/q residual과 표면 closure 통과만으로 논문의 ballooning 분류까지 검증된 것은 아니다. 남은 차이의 원인은 이번 테스트에서 확정하지 못했다. 논문 Eq. (3)와 독립적으로 local shear/metric/coefficient 전달을 검증하고, 인접 flux surface로 계산한 shear 및 실제 GS force balance와 대조하는 작업이 필요하다. 이 결과를 물리적인 second stability 부재로 해석할 수 없다.

## 구현 중 확인하고 수정한 문제

1. q를 .95까지만 제약한 상태에서 .95의 안정성을 평가하면 실제 dq/du가 제대로 제약되지 않을 수 있었다. paper benchmark는 .975까지 fit한다.
2. extract_flux_surface가 열린 contour를 강제로 닫힌 것으로 처리했다. 원래 endpoint로 closure를 검사하고 열린 contour를 거부하도록 수정했다.
3. contourc의 roundoff 크기 중복 점이 arclength PCHIP을 왜곡했다. 분석적인 원 테스트에서 확인하고 제거했다.
4. 압력 변화에 대한 FFprime 초기 예측과 충분한 contour grid가 high-beta native GS/q 수렴에 필요했다. 초기 예측은 최종 수렴 조건을 대체하지 않는다.
5. 모든 고유값을 구하던 dense 계산 대신 optional largest-only mode를 추가했다. 양/음의 가장 큰 고유값을 Cholesky bracket과 shift-invert로 구하며, 기존 기본 호출은 전체 dense spectrum을 반환한다. 분석해/여러 boundary condition/dense native 행렬과 비교했다.

## 실행 및 데이터

Bean geometry/q endpoint, analytic closed/open contour, stable/unstable largest eigenvalue, 여러 boundary condition과 dense 비교 테스트를 통과했다. 공통 코드 수정 후 실제 C-shape 2 surfaces x 2 phases, checkpoint resume, unknown retry, 수정 전 schema 거부 검사도 통과했다. MATLAB Code Analyzer와 git diff whitespace 검사에 문제는 없었다.

대표점, 첫 경계, theta/Fourier convergence를 수행한 테스트다. 5개 indentation x 4개 beta의 전체 map이나 second-stability 경계는 검증 완료된 결과로 제공하지 않는다.

`validation.mat`에는 native equilibria, 실제 계산점, 첫 경계 history, convergence table, dense/sparse comparison이 들어 있다. `figure3_validation.png`는 reference overlay와 theta refinement를 보여 준다.

```matlab
setup_indented_cshape_scan
addpath('verification')
test_chance_prl_geometry
test_flux_surface_closedness
test_largest_ballooning_eigenvalue
report=build_chance_validation_report;
```

geometry/profile 설정과 full-scan driver는 chance_prl_figure3_README.md를 참조한다. Scan checkpoint schema는 2이며, 수정 전 안정성 결과가 든 schema 1 checkpoint는 재사용하지 않는다.

Source: Chance, Jardin & Stix, PRL 51, 1963–1966 (1983), DOI 10.1103/PhysRevLett.51.1963, supplied local PDF. Dashed reference boundaries are approximate visual readings of Figure 3.
