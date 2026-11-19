# salpha.py 사용법

원본 MATLAB 파일을 수정하지 않고 작성한 별도 Python 모듈입니다.
앞서 제시한 다음 파일과 같은 폴더에 두세요.

- eq_process.py
- local_equilibrium.py
- sturm_liouville.py

필요 패키지는 numpy, scipy, matplotlib, contourpy입니다.
직접 실행하거나 import하여 사용할 수 있습니다.

## 빠른 시작

~~~python
import numpy as np
from eq_process import Equilibrium
from salpha import (
    prepare_salpha_equilibrium,
    evaluate_salpha_point,
    scan_salpha_diagram,
    plot_salpha_diagram,
)

eq = Equilibrium.from_geqdsk("gfiles/geqdsk_PT0.6")
state = prepare_salpha_equilibrium(
    eq, 0.7, model="miller",
    n_surface_points=512, n_fit_points=256,
)

point = evaluate_salpha_point(
    state, shat=1.0, alpha=1.0,
    theta0=[-np.pi, 0.0, np.pi/2],
    n_balloon_points=151,
)
print(point.lambdaMax, point.bestTheta0)

result = scan_salpha_diagram(
    state,
    shat=np.linspace(0, 5, 11),
    alpha=np.linspace(0, 6, 11),
    theta0=2*np.pi*(np.arange(12)-6)/12,
    theta_boundary=5*np.pi,
    n_balloon_points=151,
)
figures = plot_salpha_diagram(result, show=True)
~~~

메서드 이름과 옵션은 snake_case이며 결과 필드 이름은 MATLAB과 같은
lambdaMax, lambdaTheta0, bestTheta0 등을 유지했습니다.
전처리는 MATLAB 형태 prepare_salpha_equilibrium(eq.data, eq.fields, 0.7)도 지원합니다.

## 포함 함수

- prepare_salpha_equilibrium: 평형 및 국소 형상 전처리
- prepare_fem_salpha_equilibrium: 주어진 FEM 어댑터 연결
- evaluate_salpha_point: 하나의 (shat, alpha)에서 위상별 최대 고유값
- scan_salpha_diagram: 전체 격자 스캔, 실패 위상 기록, 선택적 체크포인트
- scan_fem_salpha_diagram: FEM 연결 후 스캔; (result, state) 반환
- scan_refined_salpha_diagram: 발견 격자에서 부호 변화를 찾고 직접 재계산으로 정밀화
- trace_salpha_marginal: 지정한 stable-left/unstable-right 구간의 한 경계 추적
- plot_salpha_diagram: 안정성 및 최대 고유값 위상 플롯
- shat_to_qprime / alpha_to_pprime: 원본 s-alpha 좌표 변환
- local_volume_integrals: 국소 형상의 V, 부호 있는 Vprime, dV/dr 계산

## 모델과 부피 규약

model은 miller, circular, cshape, indented-miller, direct를 지원합니다.
indented, flux-surface 등의 별칭도 사용할 수 있습니다.

volume_source="auto"에서는 다음 규약을 사용합니다.

- Miller / Direct: 기존 prepare_salpha_equilibrium.m처럼 추출 자속면의 V, Vprime
- Circular / C-shape / Indented Miller: 해당 test_salpha*.m처럼 국소 형상의 부피 적분
- C-shape / Indented Miller의 Vprime은 물리 자속의 방향을 반영한 부호를 유지

원본 자속면 적분은 state.contourIntegrals에 별도로 남습니다.
volume_source="contour" 또는 "local"로 명시적으로 선택할 수도 있습니다.

C-shape / Indented Miller는 기본적으로 피팅에 최대 256점을 사용합니다.
더 조밀한 원본 자속면과 분리하여 n_fit_points를 지정할 수 있습니다.

~~~python
state = prepare_salpha_equilibrium(
    eq, 0.7, model="cshape",
    n_surface_points=1024, n_fit_points=256,
    radial_fit_half_width=0.03, n_radial_fit=9,
    fit_options={"max_iterations": 3000, "verbose": True},
)
print(state.diagnostics.fitConverged)
print(state.diagnostics.fitRelativeRms)
~~~

max_iterations는 반복 상한입니다. 비선형 모델은 fitConverged와 fitRelativeRms를
확인해야 합니다. 기본 라이브러리 함수는 미수렴 피팅을 자동으로 수렴했다고
표시하거나 계산에서 숨기지 않습니다.

## 결과 배열

- result.lambdaMax[shat_index, alpha_index]
- result.bestTheta0[shat_index, alpha_index]
- result.lambdaTheta0[shat_index, alpha_index, phase_index]
- result.failMessage[shat_index, alpha_index, phase_index]
- result.validPhaseCount[shat_index, alpha_index]
- result.phaseComplete[shat_index, alpha_index]
- result.FFprimeUsed / periodicityClosure[shat_index, alpha_index]

양의 lambdaMax는 불안정, 음의 값은 안정 쪽을 나타냅니다.
음의 고유값을 0으로 자르지 않습니다. 최대 고유값 하나를 계산하며,
절댓값 최대 또는 0에 가장 가까운 고유값을 대신 사용하지 않습니다.

일반 스캔은 일부 위상이 실패해도 성공한 위상의 최대값을 보관합니다.
이때 phaseComplete=False이므로 모든 위상을 포함한 안정성 판정으로 해석하면
안 됩니다. 기본 플롯은 이러한 점을 가립니다.

evaluate_salpha_point와 정밀 스캔은 위상이 하나라도 실패하면 예외를 냅니다.

## 정밀 경계 탐색

~~~python
from salpha import scan_refined_salpha_diagram

refined = scan_refined_salpha_diagram(
    state,
    shat=np.linspace(0, 5, 21),
    alpha=np.linspace(0, 6, 31),
    theta0=2*np.pi*(np.arange(12)-6)/12,
    theta_boundary=5*np.pi,
    n_balloon_points=151,
    alpha_tolerance=1e-3,
    checkpoint_file="refined_salpha.npz",
)
for s, roots in zip(refined.shat, refined.roots):
    for root in roots:
        print(s, root.alpha, root.direction, root.bracket, root.lambdaAtRoot)
~~~

direction=+1은 안정에서 불안정으로, -1은 불안정에서 안정으로 전이합니다.
발견 격자에서 정확히 0인 점은 sampled_zero로 기록합니다.
발견 격자보다 좁은 안정성 섬이나 표본 사이의 접촉근을 놓칠 수 있으므로
필요하면 alpha 발견 격자도 세분화해야 합니다.

trace_salpha_marginal은 지정한 한 구간에서 왼쪽이 안정이고 오른쪽이 불안정인
경계만 정밀화합니다. 같은 부호의 양 끝 사이에 여러 근이 있는 경우를
발견하려면 scan_refined_salpha_diagram을 사용하세요.

trace의 lambdaAtMarginal은 실제 반환한 marginalAlpha에서 다시 계산한 값입니다.
max_iterations 안에 폭 허용오차를 만족하지 못하면 iteration_limit로 표시하며
converged로 표시하지 않습니다.

## 체크포인트

일반 및 정밀 스캔 모두 checkpoint_file을 지정할 수 있습니다.
지정하지 않으면 체크포인트 파일을 쓰지 않습니다.

- 완성된 행 단위로 NPZ를 임시 파일에 저장한 뒤 교체합니다.
- resume=True가 기본이며 완료된 행을 재계산하지 않습니다.
- 평형, 스캔 격자, 위상, 허용오차, 구현 소스가 달라지면 기존 파일을 거부합니다.
- 파일은 이 Python 모듈의 전용 NPZ 형식이며 MATLAB .mat 형식이 아닙니다.
- resume=False이면 처음부터 계산합니다.
- 정밀 스캔 중 실패한 행은 완료로 표시하지 않습니다.
- 일반 스캔에서는 위상 실패도 결과로 기록하므로, 그 행을 다시 계산하려면
  resume=False 또는 새 체크포인트 파일을 사용하세요.

## FEM 연결 범위

GS solver 및 build_ballooning_adapter는 아직 Python으로 포팅하지 않았습니다.
FEM 함수는 이미 준비된 Equilibrium, (data, fields[, diagnostics]) 또는
다음 호출 규약의 어댑터를 받습니다.

~~~python
def adapter(fem_result, *, nr, nz):
    # 사용자의 FEM 결과를 앞서 정의한 data/fields 인터페이스로 변환
    return data, fields, diagnostics

state = prepare_fem_salpha_equilibrium(
    fem_result, 0.7, adapter=adapter, model="direct",
    adapter_nr=257, adapter_nz=257,
)
~~~

어댑터가 없는 raw FEM 결과를 자동으로 처리한다고 가정하지 않습니다.

## 원본에서 확인된 Direct 부호 불일치

앞서 local_equilibrium.py에서 확인한 것처럼, psi가 바깥쪽으로 감소하는 일부
Direct 사례에서는 원본 MATLAB의 scalar/vector D 및 c가 반대 부호를 가집니다.
이 모듈은 require_vector_consistency=True를 기본으로 하여 해당 계수를
검증된 안정성 결과로 사용하지 않습니다.

일반 스캔은 failMessage에 기록하고, 단일 점 및 정밀 스캔은 예외를 냅니다.
require_vector_consistency=False는 원본 수치를 대조하기 위한 옵션이며,
그 결과가 물리적으로 검증되었다는 의미는 아닙니다.

## 검증 범위

- Miller(PT/NT), Direct(Circular)의 MATLAB 전처리 및 3x3x3 스캔 비교:
  주요 수치의 최대 스케일 보정 차이 약 7e-14.
- Direct(Circular)의 첫 번째/두 번째 안정성 경계와 한계선 추적 대조.
- 실제 자료에서 다섯 모델의 단일 점 계산 및 signed Vprime 변환 확인.
- 알려진 다항식으로 복수 전이, 표본 0, 경계 포함, 반복 제한 상태 확인.
- 중단 후 재시작, 변경된 상태 거부, 완료된 행 재사용 확인.
- 일부 위상 실패, 엄격한 단일 점 실패, 비정상 입력 확인.
- FEM 어댑터/미리 변환된 인터페이스 및 CLI 실행 확인.
- 플롯을 렌더링하여 레이블과 경계 표본 표시 확인.

그리드/구간/피팅의 수렴까지 보장하는 비교는 아닙니다.
~~~text
python salpha.py gfiles/geqdsk_PT0.6 --nshat 11 --nalpha 11 --ntheta0 8 --plot
~~~
