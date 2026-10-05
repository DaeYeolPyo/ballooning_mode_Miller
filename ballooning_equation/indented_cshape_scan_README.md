# Indented C-shape pressure scan

현재 outboard indented C-shape에서 prescribed-p/q 평형을 다시 풀어 ideal infinite-n ballooning 안정성 지도를 만든다. 두 번째 안정 영역을 전제하지 않으며, 발견된 안정성 전이를 모두 추적한다.

2026-10-02 Chance PRL Figure 3 benchmark에서는 첫 안정 경계는 논문 근방에 나왔으나 고 beta 안정 영역은 재현되지 않았다. 따라서 현재 C-shape scan의 결과만으로 second stability의 물리적 부재를 검증했다고 해석할 수 없다. 대표점·경계·수렴 검사의 실제 결과는 [Chance benchmark report](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/ballooning_equation/output/chance_prl_figure3/validation/validation_report.md)에 기록했다. 기본 ballooning 해상도는 시작값이며, Chance test에서 24 cells/turn의 안정 판정이 refinement 후 뒤집히는 사례를 확인했다.

열린 contour 처리와 고유값 계산을 수정한 뒤 일반 C-shape checkpoint schema는 3이다. 수정 전 cached stability는 재사용하지 않으며 예제는 이전 checkpoint를 별도 파일로 보존한다.

## 실행

MATLAB에서 `ballooning_equation` 폴더를 현재 폴더로 설정한다.

```matlab
setup_indented_cshape_scan;
config = make_indented_cshape_scan_config;
scan = scan_indented_cshape_beta(config);
trace = trace_indented_cshape_marginal(scan);
plot_indented_cshape_beta(scan, trace);
```

기본 설정은 기존 fixed-q 예제의 실제 기하학적 크기를 기준으로 한다. 기본 `Indentation=[0.6 0.8 1.0]`, 중심 압력 제어값 `CentralBeta=[0.02 0.1 0.3]`이며, 총 9개 평형을 계산한다. `CentralBeta`는 그림의 체적 평균 beta가 아니다. 실제 beta는 각 평형에서 계산한다. 기본값은 시작용 해상도이며, 최종 안정 경계에는 수렴 검사가 필요하다.

checkpoint 저장·재개·그림 내보내기를 포함한 예제:

```matlab
addpath(fullfile(pwd,'examples'));
[scan,trace] = run_indented_cshape_beta_scan;
```

기본 출력은 `output/indented_cshape/beta_scan.mat`, `marginal_trace.mat`, `beta_stability.png`이다. 같은 설정의 checkpoint가 있으면 완료된 판정점을 재사용하고 `unknown` 점은 재시도한다. 이전 버전이거나 설정이 다른 checkpoint는 `beta_scan_previous_날짜시간.mat`으로 보존하고 새 계산을 시작한다. 여러 설정을 비교할 때는 별도 출력 폴더를 지정한다.

코드를 수정한 뒤 다시 실행할 때는 먼저 `clear functions`로 이미 로드된 함수를 갱신한다. 위 예제는 9개 seed 계산 이후 경계 추적까지 실행한다. 고압력의 p/q 최적화는 한 점에 수십 분이 걸릴 수 있다. 우선 평형과 seed 판정만 확인하려면 첫 번째 코드처럼 `scan_indented_cshape_beta`까지만 실행한다.

## 평형 초기화와 수렴 실패 복구

기본 q 계산 격자는 `contourNR=200`이다. 이전 값 120에서는 첫 압력점의 내부 GS가 수렴하더라도 `psiN=0.95`의 contour가 격자 밖 NaN 영역 때문에 끊겨 p-q 후보에서 탈락했다. 마지막 FFprime 후보의 GS 실패만 표시하던 오류도 수정해, 모든 초기 후보의 GS residual/update와 q contour 실패 위치를 보고한다. 격자 200의 q 계산에는 고정 mesh의 P3 값 보간을 재사용한다. canonical P3 basis와 gradient 평가는 그대로 사용한다.

ballooning 표면의 기본 해상도도 `NSurfacePoints=1024`, `FourierModes=48`로 높였다. 512점/32모드에서는 낮은 압력의 외곽면이 closure 품질 기준을 넘는 경우가 있었다. `I=0.8, psiN=0.95`에서 closure defect는 약 0.050에서 0.020으로, `I=1.0`에서는 약 0.081에서 0.036으로 감소했다. GS/q 수렴 여부와 ballooning 표면 품질 실패는 서로 구분해서 출력한다.

수렴한 이전 평형이 있으면 그 FFprime과 자속으로 목표 압력을 먼저 계산한다. 실패하면 압력 간격을 절반으로 줄이고, 중간 평형이 GS와 q 허용오차를 모두 만족한 경우에만 다음 예측값으로 사용한다. 첫 목표점이 실패하면 작은 압력의 bootstrap 평형에서 시작한다. 중간 압력은 continuation에만 쓰며 그림의 sampling 점에는 추가하지 않는다.

```matlab
config.Continuation.Enabled = true;
config.Continuation.BootstrapPressureAmplitude = 7.162e3; % Pa, BReference=3 T
config.Continuation.MinimumRelativeStep = 1e-3;
config.Continuation.MinimumAbsoluteStep = 10; % Pa
config.Continuation.StepGrowth = 1.5;
config.Continuation.MaxAttempts = 24;
config.Continuation.MaximumQIterationRestarts = 1;
config.RetryUnknown = true; % 같은 설정 checkpoint의 미판정 점만 재시도
```

허용오차는 완화하지 않았다. 실패가 계속되면 최소 압력 간격 또는 최대 시도 횟수에서 중단하며 목표점은 `unknown`으로 남는다. `scan.points{ip,ii,ia}.equilibriumAttempts`에 시도 압력, GS/q 오차, 실패 사유와 종료 조건이 저장된다. q 수렴 기준과 비교할 값은 `solverDiagnostics.qConstraintMaximumRelativeError` 및 `qConstraintRmsLogError`이다. 별도의 `qMaximumRelativeError`는 경계 extrapolation을 포함한 전체 q 진단이다.

GS는 수렴하고 q 최적화가 최대 반복 횟수에 도달한 경우에는 같은 압력의 마지막 GS 해에서 q 최적화를 한 번 더 이어간다. 이때도 q 허용오차를 충족하기 전에는 목표점이나 다른 압력의 predictor로 채택하지 않는다. 경계 추적은 저장된 bracket 평형 또는 형상별 마지막 predictor를 재사용한다.

q 반복에서 두 허용오차를 모두 만족한 trial을 발견하면 그 trial을 채택한다. 평균 오차가 더 작다는 이유로, 최대 q 오차 기준을 초과한 다른 trial로 교체하던 후보 선택도 수정했다.

고압력 계산이 오래 걸릴 때는 `config.PQConstraint.outer.verbose=true`로 각 q 반복의 오차를 출력할 수 있다. 각 q 후보마다 내부 GS를 다시 풀기 때문에, GS 수렴과 prescribed-q 수렴은 별도의 조건이다.

## 설정 예시

```matlab
% Grimm형 종횡비 비교: 실제 elongation은 동일하게 유지한다.
config = make_indented_cshape_scan_config( ...
    AspectRatio=[3 4 7 10], Indentation=[0 0.6 0.8 1.0 1.15]);

% 논문의 elongation에 맞춘 별도 비교 계열.
config = make_indented_cshape_scan_config(AspectRatio=4,Kappa=1.386);

% 실제 압력 p0 [Pa]를 직접 지정할 수도 있다.
config = make_indented_cshape_scan_config(PressureAmplitude=[1e4 1e5 5e5]);

% 경계 부근 해상도를 높이는 예시. 필드라인 길이와 turn당 해상도를
% 각각 독립적으로 증가시켜 수렴을 확인한다.
config.Ballooning.Theta0 = 2*pi*((0:15)-8)/16;
config.Ballooning.ThetaBoundary = 10*pi;
config.Ballooning.CellsPerTurn = 48;
config.Ballooning.NSurfacePoints = 1024;
config.Ballooning.FourierModes = 48;

% 가장자리 면을 추가하면 q 제약의 검사 범위도 함께 확장한다.
config.QEvaluation.interiorPsiN = [.03:.04:.95 .975 .99].';
config.PQConstraint.FFprime.fitRange = [0 .99];
config.Ballooning.RadialLabels = [.03:.04:.95 .975 .99];

config.CheckpointFile = fullfile(pwd,'output','my_scan.mat');
config.Resume = false;  % 저장 파일이 이미 있으면 덮어쓰지 않고 오류를 낸다.
config.StoreEquilibria = false; % true이면 모든 sampling 평형을 저장한다.
```

## 구현과 결과

| 함수 | 역할 |
|---|---|
| `normalize_indented_cshape_boundary` | 실제 minor radius, 종횡비, elongation을 고정하고 생성 경계의 자기접촉·방향을 검사 |
| `make_indented_cshape_scan_config` | 현재 형상에 맞춘 시작 설정과 Chance p/q 프로파일 설정 |
| `build_indented_cshape_scan_input` | dp/dpsiN 입력과 고정 Fboundary, 이전 평형의 FFprime 초기값 구성 |
| `compute_volume_beta` | P3 자속 기울기와 F로 총 자기장을 이용한 beta 체적 적분 |
| `evaluate_equilibrium_ballooning` | 실제 Direct 자속면·위상 스캔과 기하 보정 전후 진단 |
| `evaluate_indented_cshape_scan_point` | 한 pressure amplitude의 평형, beta, ballooning 평가 |
| `solve_indented_cshape_equilibrium` | 실패한 목표 압력의 bootstrap 및 압력 간격 축소, 시도 이력 저장 |
| `scan_indented_cshape_beta` | pressure continuation, prepared FEM 재사용, checkpoint 관리 |
| `trace_indented_cshape_marginal` | 인접한 완전 계산점의 모든 lambda 부호 전이를 bracket/refine |
| `plot_indented_cshape_beta` | 안정/불안정/경계/미판정 표시와 종횡비별 경계 그림 |

`solve_GS(input)` 호출은 계속 지원한다. 새 `solve_GS(input,preparedProblem)`은 동일한 경계·grid·reference·Picard·q 검사 설정의 prepared FEM을 재사용한다. 다른 형상의 prepared 객체를 전달하면 오류를 낸다. 동일 mesh의 압력 continuation에서는 이전 자속도 초기값으로 사용한다.

checkpoint에는 재개를 위해 형상당 가장 최근의 수렴한 predictor 한 개를 저장한다. `StoreEquilibria=true`이면 모든 sampling 점의 평형을 추가로 저장한다.

배열 순서는 **pressure x indentation x aspect ratio**이다. `scan.points{ip,ii,ia}`에는 실제 beta, lambda, 상태, solver diagnostics와 ballooning diagnostics가 들어간다.

- `stable`: 모든 요청 자속면·위상 계산이 완료됐고 최대 lambda < -LambdaTolerance.
- `unstable`: 요청한 면·위상 중 하나라도 lambda > LambdaTolerance. 다른 점이 실패했다면 완전한 최대값은 미정이므로 `lambdaMax=NaN`, 관측된 양의 값은 `lambdaMaxObserved`에 남는다.
- `marginal`: 완전 계산의 최대 lambda가 허용오차 안에 있다.
- `unknown`: 수렴 실패, 기하 검사 실패 등으로 안정 판정이 불가능하다.

`stable`은 요청한 sampling과 유한 field-line 영역에서의 판정이다. `complete`는 그 요청 계산의 완료 여부이며 연속적인 모든 자속면·위상에 대한 증명은 아니다. 기본 sampling은 축과 경계의 극히 가까운 면까지 포함하지 않는다.

각 자속면의 `surfaceDiagnostics`에는 자속면 오차, 법선 정렬, Bp 보정량, 보정 전후 periodicity와 FFprime이 저장된다. 보정이 설정된 품질 한계를 넘으면 해당 면은 미판정이 된다. 수렴을 확인하지 않고 품질 한계를 완화해 안정 판정을 만드는 용도로 사용하지 않는다.

`trace.roots`에는 전이 방향, p0·beta·lambda bracket, 계산점, 수렴 여부가 저장된다. 전이가 없는 경우 빈 배열을 반환한다. 미판정 sampling 구간을 넘어 연결하지 않는다. 같은 부호의 두 seed 사이에 숨어 있는 좁은 불안정 구간이나 접선 접촉은 coarse scan만으로 검출할 수 없으므로 pressure/indentation sampling을 추가한다.

그림에서 경계는 전이 방향과 압력 순위로 인접 indentation 사이를 연결한다. 여러 경계가 생기고 사라지는 접힘 근처에서는 2차원 sampling을 높여 branch 대응을 확인한다. 미판정 경계 bracket은 점선으로 표시한다. beta 자체를 계산하지 못한 실패점은 임의 위치에 표시하지 않고 개수를 제목에 명시한다.

## 검증

```matlab
addpath(fullfile(pwd,'verification'));
report = test_indented_cshape_pipeline;
report = test_indented_cshape_pipeline(RunIntegration=true);

% GS를 다시 풀지 않고 저장된 평형에서 더 촘촘한 ballooning 해상도를 검사.
% scan 계산 때 config.StoreEquilibria=true로 설정해야 한다.
config.Ballooning.NSurfacePoints = 1024;
config.Ballooning.FourierModes = 64;
report = verify_indented_cshape_equilibria(config,scan.equilibria);
```

빠른 검증은 16개 경계 조합, 해석 가능한 삼각형에서의 beta 적분, 자속 방향 반전, 고정 P3 보간의 상수·3차식·node 값 재현, 단일/다중/없는 전이, 정확히 샘플링된 0, 미판정 seed 및 refinement 구간의 비연결, 마지막 iteration에서의 bracket 수렴, 그림 생성을 검사한다. integration 검증은 실제 C-shape mode-3 평형, 2개 자속면 x 2개 위상, checkpoint 재개, 저장한 predictor로 `unknown` 재시도, 잘못된 prepared geometry 거부를 검사한다.

구현 검증에서 실제 형상 `I=0.8`, `Rgeo/ageo=4.82996`, `kappa=4.17813`에 대해 `GridSize=[24 32]`, `NBoundary=80`, `CentralBeta=[0.002 0.02]`의 두 압력점을 연결해 계산했다. 체적 beta는 각각 약 **0.068707%, 0.69465%**이고, 2개 자속면 x 2개 위상의 계산은 두 점 모두 안정이었다. 두 번째 평형에는 기본 11개 자속면 x 8개 위상도 적용했고, 모든 계산이 완료됐으며 최대 lambda는 약 **-1.06229e-5**였다. 형상 품질 기준을 만족하지 못한 계산은 `unknown`으로 판정되는 것도 확인했다.

기존 strongly indented Direct FEM 검증도 통과했다. Fourier 32/48 모드에서 lambda 상대 차이는 약 `4.90e-5`였다.

이번 수렴 실패 수정에서는 기본 `GridSize=[48 64]`, `NBoundary=128`, `contourNR=200`의 실제 압력 범위를 검사했다. 아래 8개 점은 GS/q 기준을 모두 만족했고, 1024점/48모드로 11개 자속면 x 8개 위상 계산을 완료했다. 검증한 점들의 최대 GS residual은 `4.354996e-7`, 최대 q 상대오차는 `1.612234e-3`, 최대 RMS log-q 오차는 `7.292625e-4`였다. 고정 P3 보간을 재사용한 q 계산과 원래 계산의 최대 차이는 `2.109e-15`였다.

| I | p0≈71,620 Pa: 체적 beta / 판정 | p0≈358,099 Pa | p0≈1,074,296 Pa |
|---|---|---|---|
| 0.6 | 0.71775% / stable | 3.7851% / unstable | 13.266% / unstable |
| 0.8 | 0.69446% / stable | 3.6555% / unstable | q 수렴 미확인 |
| 1.0 | 0.66544% / stable | 3.4948% / unstable | 12.090% / unstable |

`I=0.8`의 최고 압력점은 별도 고압력 초기값 진단에서도 GS residual `9.610e-8`을 얻었지만 q 최대 상대오차 `2.998e-3`으로 멈췄다. 이 점을 수렴한 평형이나 안정성 판정으로 채택하지 않았다. 이는 첫 압력점의 contour 실패와 다른 q 최적화 문제이며, 추가적인 q-response·profile 표현·mesh 수렴 검사가 필요하다. 이 실패만으로 물리적인 평형이 없다고 결론 내릴 수는 없다.

결과는 [검증 MAT](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/ballooning_equation/tmp/cshape_tests/default_validation.mat)와 [seed 판정 그림](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/ballooning_equation/tmp/cshape_tests/default_seed_validation.png)에 저장했다. 이들은 구현의 동작 검증이다. 전체 9개 seed의 확인과 marginal boundary 추적 및 최종 해상도 수렴 검증은 구분해야 한다.

전체 수식과 설계 배경은 [상세 workflow](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/ballooning_equation/indented_cshape_ballooning_workflow.md)를 참고한다.
