# local_equilibrium.py 사용법

원본 MATLAB 파일을 수정하지 않은 검토용 Python 포팅입니다.
이 파일을 앞서 제시한 eq_process.py와 같은 폴더에 두세요.
필요 패키지는 numpy, scipy입니다. eq_process.py에는 contourpy도 필요합니다.

## 포함된 공개 함수

- 경계: circular_boundary, miller_boundary, cshape_boundary, indented_miller_boundary
- 피팅: fit_circular, fit_miller, fit_cshape, fit_indented_miller
- 국소 평형: circular_mercier_luc, miller_mercier_luc, cshape_mercier_luc,
  indented_miller_mercier_luc, flux_surface_mercier_luc
- 계수: evaluate_theta_psi, evaluate_gcf
- 교차점: find_r_at_z0_from_surface

## Miller 경로

~~~python
import numpy as np
from eq_process import Equilibrium
from local_equilibrium import (
    fit_miller, miller_mercier_luc, evaluate_theta_psi, evaluate_gcf,
)

eq = Equilibrium.from_geqdsk("gfiles/geqdsk_PT0.6")
psiN = 0.7
surf = eq.extract_flux_surface(psiN, npoints=1024)

param, bnd = fit_miller(
    eq, psiN, ntheta=256, dpsi=0.03, n_radial_fit=9,
)
local = miller_mercier_luc(surf, param, bnd)

Gamma, theta_psi, diag = evaluate_theta_psi(
    local,
    FFprime=float(eq.fields.FFprime(psiN)),
    pprime=float(eq.fields.pprime(psiN)),
    qprime=float(eq.fields.qprime(psiN)),
    enforce_periodicity=True,
)

bal = evaluate_gcf(
    diag, theta_bal=0.0, theta_bnd=5*np.pi,
    a_N=eq.data.rmaxis, B_N=abs(eq.data.bcentr),
)
assert bal.validation.all_finite
assert bal.validation.vector_consistent
# bal.theta, bal.g, bal.c, bal.f를 고유값 계산용 격자에 재보간합니다.
~~~

마지막 주석의 의미: bal.theta, bal.g, bal.c, bal.f를 고유값 계산용 격자에
PCHIP으로 재보간한 뒤 Sturm_Liouville에 전달하면 됩니다.

피팅 함수는 MATLAB과 비슷한 호출 형태도 지원합니다.

~~~python
param, bnd = fit_miller(eq.data, eq.fields, psiN, ntheta=256)
~~~

## 다른 형상

- Circular: fit_circular + circular_mercier_luc
- C-shape: fit_cshape + cshape_mercier_luc
- Indented Miller: fit_indented_miller + indented_miller_mercier_luc
- Direct: flux_surface_mercier_luc(surf, eq.fields)

C-shape와 Indented Miller는 다음 옵션을 사용할 수 있습니다.

~~~python
param, bnd = fit_cshape(
    eq, psiN, ntheta=256, n_radial_fit=9, dpsi=0.03,
    max_iterations=3000, n_phase_starts=32, verbose=True,
)
print(param.fit.converged, param.fit.relative_rms_distance)
~~~

Indented Miller에는 n_x_starts=33, max_abs_delta=0.98,
max_kappa=5.0, max_I0=3.0 옵션도 있습니다.
max_iterations는 반복 상한이며 수렴을 보장하지 않습니다.
비선형 피팅 후에는 param.fit.converged 및 relative_rms_distance를 확인하세요.

## 격자와 물리량

- eq_process의 surf는 중복 끝점을 제외합니다.
- Boundary와 LocalEquilibrium은 0과 2*pi를 모두 포함하는 닫힌 격자입니다.
- 피팅의 dpsi는 정규화 자속 psiN의 반폭입니다.
- qprime은 dq/dpsi이며, dq/dpsiN과 구분합니다.
- pprime, FFprime도 물리 자속 기준 원본 값을 사용합니다.
- enforce_periodicity=True는 입력 FFprime 대신 주기성 제약을 만족하는 값을
  사용합니다. diag.FFprime_input, FFprime_periodic, FFprime_used로 비교할 수 있습니다.
- theta_psi 등의 주기적인 기하는 주기 확장하지만, I의 qprime*(theta-theta0)
  항은 감지 않은 theta 좌표를 사용하여 세속적 자기 시어를 유지합니다.

## 원본에서 확인된 제한

psi가 바깥쪽으로 감소하는 일부 Direct 사례(검증한 kappa3 GEQDSK 포함)에서
MATLAB 원본의 D/c와 벡터 재구성 D_vector/c_vector가 반대 부호를 보입니다.
이 포팅은 원본 scalar 계산식을 변경하지 않았습니다. 해당 상황에서는
RuntimeWarning을 내고 bal.validation.vector_consistent=False를 반환합니다.
bal.validation.all_finite=True는 부호 일치나 물리적 타당성을 뜻하지 않습니다.
해당 Direct 사례의 c를 안정성 판정에 사용하기 전에 부호 규약을 정리해야 합니다.

## 검증 범위

- 네 가지 경계의 해석적 1차·2차 미분 및 주기성 확인.
- 실제 GEQDSK 12개에서 Direct 모델의 실행, q 재구성, theta_psi 주기성 확인.
- MATLAB과 7개 경로 비교: Circular, Miller(PT/NT), C-shape, Indented Miller,
  Direct(Circular/자속 감소 사례).
- 비선형 피팅 비교는 동일한 20회 반복 설정으로 수행했으며, 해당 비교가
  수렴한 피팅을 의미하지는 않습니다.
- 알려진 형상의 계수 회복, 두 입력 API, 비정상 입력, 주기성 설정 확인.
- GEQDSK → Miller → Mercier-Luc → theta_psi → g/c/f → Sturm-Liouville 연결 확인.
