# Indented C-shape ballooning stability workflow

작성일: 2026-10-01

현재 코드와 경계식을 기준으로 정리한 설계이며, 형상 정규화·평형 및 beta 계산·Direct ballooning 평가·스캔 및 경계 추적 driver를 구현했다. 실행 방법은 [구현 README](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/ballooning_equation/indented_cshape_scan_README.md)에 있다. 두 번째 안정 영역을 전제하지 않으며, 실제 전체 안정성 지도는 요청한 해상도와 압력 범위를 실행한 결과로 판단한다.

## 1. 가능성 판단과 재현 목표

현재 indented C-shape로 Chance Figure 3과 같은 안정/불안정 점 지도와 Grimm Figure 4와 같은 종횡비별 marginal stability 곡선을 계산할 수 있다. 이미 필요한 prescribed-p/q 평형 solver와 임의 자속면을 사용하는 Direct ballooning 경로가 존재한다.

다만 두 논문의 bean은 inboard, 작은 R 쪽의 indentation이고, 현재 코드는 outboard, 큰 R 쪽을 누르는 형상이다. 따라서 논문과 동일한 임계값이나 두 번째 안정 영역의 존재는 보장되지 않는다. 비교 대상은 우선 같은 계산 실험을 다른 경계 계열에 적용한 결과이다. 안정성은 여기서 ideal infinite-n ballooning 안정성만 뜻한다.

현재 형상은 `CShape_boundary.m`의 삼각함수 5-매개변수 경계가 아니라, `test_qsolve_indented_fixed_q.m`와 `Indented_Miller_boundary.m`에 사용된 다음 식을 뜻한다.

## 2. 현재 형상과 indentation의 정의

기존 `R0, a, kappa, delta, I0`를 각각 구성 매개변수 $R_c,a_0,\kappa_0,\delta_0,I$로 표기한다. 논문의 기하학적 치수와 혼동을 피하기 위한 표기이다.

$$
x=\arcsin\delta_0,\qquad h(\vartheta)=\left(\frac{1+\cos\vartheta}{2}\right)^3,
$$
$$
R_b(\vartheta)=R_c+a_0\cos(\vartheta+x\sin\vartheta)-Ia_0h(\vartheta),
\qquad Z_b(\vartheta)=\kappa_0a_0\sin\vartheta.
$$

여기서 $I$는 코드의 `indent`/`I0` 값이고, 실제 $d/(2a)$가 아니다.

$$
R_b(0)=R_c+a_0(1-I),\qquad R_b(\pi)=R_c-a_0.
$$

indentation 항은 $\vartheta=0$에서 최대이고 $\vartheta=\pi$에서는 0이므로, 현재 경계는 outboard indentation이다. 중앙 midplane에서의 국소 concavity 조건은 다음과 같다.

$$
R_{b,\vartheta\vartheta}(0)=a_0\left[-(1+x)^2+\frac{3I}{4}\right]>0.
$$

$\delta_0=-0.4$에서는 $I>0.46175$가 이 국소 concavity 조건이다. 이것은 안정성 임계값이 아니라 경계의 기하학적 조건이다. 다른 위치에서의 concavity나 자기축을 기준으로 한 기하를 대신 판정하는 조건도 아니다.

실제 경계에서 다음 값을 계산한다.

$$
a_g=\frac{R_{\max}-R_{\min}}2,\quad R_g=\frac{R_{\max}+R_{\min}}2,
\quad b_g=\frac{Z_{\max}-Z_{\min}}2,
$$
$$
\mathcal A=R_g/a_g,\qquad \kappa_g=b_g/a_g,
\qquad d_{\rm out}=R_{\max}-R_b(0),\qquad i_{\rm out}=d_{\rm out}/(2a_g).
$$

$i_{\rm out}$는 radial envelope를 기준으로 측정한 outboard midplane 깊이다. 국소 concavity 여부는 위의 $R_{b,\vartheta\vartheta}(0)$와 함께 기록한다. 논문의 inboard indentation과 같은 위치의 지표로 해석하지 않는다.

기존 fixed-q test의 $R_c=3, a_0=0.8, \kappa_0=3, \delta_0=-0.4, I=0.8$를 경계식으로 검산하면:

| 측정값 | 결과 |
|---|---:|
| $a_g$ | 약 0.57442 |
| $R_g$ | 약 2.77442 |
| $\mathcal A$ | 약 4.82996 |
| $\kappa_g$ | 약 4.17813 |
| $d_{\rm out}$ | 약 0.18884 |
| $i_{\rm out}$ | 약 0.16437 |

단순히 기존 $R_c,a_0,\kappa_0$를 고정하고 $I$만 바꾸면 종횡비와 elongation도 같이 바뀐다.

아래 그림의 왼쪽은 현재 경계와 선택적인 inboard 반사 대조군이며, 오른쪽은 실제 종횡비와 elongation을 고정한 C-shape 계열이다. 평형과 안정성 계산을 실행한 그림은 아니다.

![C-shape geometry check](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/ballooning_equation/indented_cshape_geometry.png)

## 3. 종횡비와 elongation을 고정하는 경계 정규화

먼저 무차원 radial shape를 계산한다.

$$
r_I(\vartheta)=\cos(\vartheta+x\sin\vartheta)-Ih(\vartheta),
$$
$$
s_I=\frac{\max r_I-\min r_I}{2},\qquad m_I=\frac{\max r_I+\min r_I}{2}.
$$

목표 minor radius $a_g=a_*$, 실제 종횡비 $\mathcal A$, 실제 elongation $\kappa_g=\kappa_*$를 정하면

$$
a_0=\frac{a_*}{s_I},\qquad \kappa_0=\kappa_*s_I,
\qquad R_c=\mathcal A a_*-a_0m_I.
$$

이 값을 원래 경계식에 넣으면 $a_g=a_*$, $R_g=\mathcal A a_*$, $b_g=\kappa_*a_*$가 유지된다. 구성 파라미터 $R_c$와 코드 reference 값 `input.dim.R0`를 구분한다. 후자는 기하학적 중심 $R_g$로 설정하고, `input.dim.aMinor=a_*`로 둔다.

$$
i_{\rm out}(I,\delta_0)
=\frac{\max r_I-(1-I)}{2s_I}.
$$

최초 스캔은 $I$를 가로축으로 수행하고 측정한 $i_{\rm out}$도 저장한다. 논문과 같은 깊이 지표로 표시하려면 $i_{\rm out}$를 가로축으로 다시 그린다. 목표 $i_{\rm out}$ 격자를 원하면 concave branch에서 $I$를 역으로 찾는다.

고정한 $\delta_0$는 구성 매개변수이다. 실제 top-point triangularity도 필요하면

$$
\delta_g=\frac{R_g-R_b(\pi/2)}{a_g}
=\frac{m_I+\delta_0+I/8}{s_I}
$$

를 별도로 저장한다. $I$를 바꾸면 $\delta_g$도 변한다. 이 실험은 고정된 $\delta_0$를 가진 C-shape 계열의 변화이다.

매 경계에서 $R_b>0$, tangent 비퇴화, 자기교차 부재, 양의 면적을 확인한다. 자기장이 정의된 $R$ 방향을 단순히 뒤집어 같은 안정성으로 간주할 수는 없다.

## 4. 논문과 비교할 압력 및 q 프로파일

$$
u=\psi_N=\frac{\psi_{\rm native}-\psi_{\rm axis}}{\psi_{\rm boundary}-\psi_{\rm axis}},\qquad 0\leq u\leq1,
$$
$$
p(u)=p_0(1-u^2)^2,\qquad \frac{dp}{du}=-4p_0u(1-u^2),\qquad p_{\rm edge}=0,
$$
$$
q(u)=1.03+0.84375u-1.1775u^2+3.50375u^3,
$$
$$
\frac{dq}{du}=0.84375-2.355u+10.51125u^2.
$$

이는 Chance의 양 끝점 조건을 만족한다. Grimm의 반올림된 값 $q_u(0)=0.84$를 정확히 적용하려면 $q(u)=1.03+0.84u-1.17u^2+3.5u^3$를 별도 설정으로 사용한다.

현재 GS 입력 `input.profile.pprime`은 $dp/d\psi_N$이다. 따라서 다음과 같이 설정한다.

```matlab
u = linspace(0,1,nProfile).';
input.profile.mode = 3;
input.profile.pprime = -4*p0*u.*(1-u.^2);
input.profile.pedge = 0;
input.profile.q = struct('psiN',u, ...
    'q',1.03 + 0.84375*u - 1.1775*u.^2 + 3.50375*u.^3);
input.profile.normalization = struct('mode','fixed_Fb', ...
    'Fb',Rgeo*Bref);
```

기존 test의 `pprime0*(1-u.^1.5).^2`는 다른 압력 프로파일을 만든다. 논문과 형상의 효과를 비교하는 계산에서는 위 식으로 교체한다.

현재 native GS 부호 규약은

$$
-\Delta^*\psi=\mu_0R^2\frac{dp}{d\psi}+F\frac{dF}{d\psi},
\qquad F=RB_\phi,
$$

이다. GS 내부에서 $dp/du$와 $F\,dF/du$를 부호가 있는 native flux span으로 나눈다. 호출하는 쪽에서 먼저 나누지 않는다. mode 3는 $F\,dF/du$의 계수를 조정하고 GS를 다시 풀어 지정한 $q$를 맞춘다. $I_p$는 결과로 결정되므로 독립적으로 고정하지 않는다.

같은 형상의 압력 continuation에서는 mesh와 FEM operator를 재사용할 수 있다. 인접한 $p_0$의 $F\,dF/du$ 등을 초기값으로 전달한다. 형상이 바뀌면 경계와 mesh를 다시 만든다.

낮은 압력부터 $p_0$를 증가시키고, `solverResult.converged`와 `constraintResult.converged`를 모두 확인한다. 현재 test는 주로 $u\leq0.95$의 q fit을 검사한다. 가장자리 ballooning을 최종 판정에 포함하려면 q evaluation과 fit range의 범위 및 해상도도 확장해 확인해야 한다.

단일 자기축과 중첩된 닫힌 자속면이 유지되는지도 확인한다. 내부 separatrix, 여러 자기축, 분리된 contour가 생기면 현재 단일 축 기반 준비 경로로 그 평형을 그대로 판정할 수 없다. 해당 경우는 별도 topology 처리 전까지 unknown으로 기록한다.

## 5. 그림의 세로축: 총 자기장을 이용한 평균 beta

$$
B_p^2=\frac{\psi_R^2+\psi_Z^2}{R^2},\qquad B_\phi^2=\frac{F(u)^2}{R^2},
\qquad dV=2\pi R\,dR\,dZ,
$$
$$
\beta_V=\frac{2\mu_0\int_\Omega p\,dV}{\int_\Omega(B_p^2+B_\phi^2)\,dV}.
$$

현재 `equilibrium.global.betaTor`는 $2\mu_0\langle p\rangle/B_0^2$이며 논문의 정의와 다르다. 총 자기장의 체적 적분을 분모로 계산하는 후처리가 필요하다.

P3 FEM과 기존 Dunavant quadrature를 사용하면 다음과 같이 계산할 수 있다.

$$
\mathcal P=\sum_{e,k}2\pi R_{ek}\,w_k\det J_e\,p(u_{ek}),
$$
$$
\mathcal M=\sum_{e,k}2\pi R_{ek}\,w_k\det J_e
\frac{\psi_{R,ek}^2+\psi_{Z,ek}^2+F(u_{ek})^2}{R_{ek}^2},
\qquad \beta_V=2\mu_0\mathcal P/\mathcal M.
$$

그림에는 $100\beta_V$를 표시한다. $p_0$가 입력이고 $\beta_V$가 출력이므로 우선 $p_0$를 스캔하고 계산된 $\beta_V$를 저장한다. 특정 beta 격자가 필요하면 같은 평형 가지에서 $\beta_V(p_0)-\beta_{\rm target}=0$을 추가로 푼다.

## 6. Direct ballooning 경로

강하게 들여쓴 내부 자속면에는 `Model="Direct"`를 사용한다. 추출한 자속면과 FEM의 자속 기울기로 기하를 구한다.

```text
solve_GS(mode=3)
  -> build_ballooning_adapter
  -> prepare_salpha_equilibrium(..., Model="Direct")
  -> state.thetaDiag
  -> evaluate_gcf(theta0)
  -> construct_matrix
  -> calculate_eigenvalue(...,'dirichlet')
```

기존 `verify_indented_direct_ballooning_point.m`의 point solve를 연결 예제로 사용할 수 있다. 해당 verification은 mode 1과 다른 압력 프로파일을 사용하므로, 이번 스캔에서는 mode 3의 출력 평형을 전달한다.

`prepare_salpha_equilibrium`은 실제 평형의 $p',q',FF'$로도 준비할 수 있다. 이번 압력 스캔에서는 `state.thetaDiag`를 직접 사용한다. `evaluate_salpha_point`는 $s,\alpha$를 자유롭게 바꾸는 국소 스캔용이므로 이번 평형 평가 driver에서는 호출할 필요가 없다.

adapter는 바깥쪽으로 증가하는 $\psi_*=D u$, $D=|\Delta\psi_{\rm native}|$를 사용한다.

$$
p'=\frac{dp}{d\psi_*}=\frac{-4p_0u(1-u^2)}D,
\qquad q'=\frac{dq}{d\psi_*}=\frac{0.84375-2.355u+10.51125u^2}D.
$$

이 변환은 adapter가 수행한다. GS 입력과 ballooning 입력에서 미분하는 자속 변수가 다르다는 점을 명시한다.

실제 평형에서 다음 값을 구한다.

$$
B_p=|\nabla\psi_*|/R,\qquad B_\phi=F/R,
\qquad q=\frac{F}{2\pi}\oint\frac{dl}{R^2B_p},
$$

법선 곡률, $B_p$, PEST angle, local shear는 실제 내부 자속면을 이용해 계산한다. 경계의 $I$를 내부 자속면의 형상 계수로 그대로 적용하지 않는다.

현재 Direct 구현에서 integrated shear에 해당하는 양은 다음과 같다.

$$
\mathscr I(\theta)=q'(\theta-\theta_0)+q\theta_\psi(\theta),
$$
$$
K=|\nabla\alpha_c|^2
=\frac1{R^2}+\frac{F^2}{R^4B_p^2}+(RB_p)^2\mathscr I^2,
\qquad \alpha_c=\zeta-q(\theta-\theta_0).
$$

이를 이용해 `evaluate_gcf`가 장력 계수 $g$, 압력 및 곡률 계수 $c$, 관성 계수 $f$를 구성한다. $\mathscr I$와 경계 indentation $I$는 서로 다른 양이다.

무차원화에서 $\nabla_N=a_N\nabla$, $B_N^*=B/B_N$, $\bar p=\mu_0p/B_N^2$, $\bar\psi=\psi_*/(a_N^2B_N)$,
$\bar H=(\mu_0p+B^2/2)/B_N^2$, $\nu=\mathbf b\cdot\nabla_N\theta$, $K_N=|\nabla_N\alpha_c|^2$로 정의하면 계수는 다음과 같다.

$$
g=\frac{\nu K_N}{B_N^*},\qquad
f=\frac{K_N}{\nu(B_N^*)^3},\qquad
c=\frac{2}{\nu(B_N^*)^4}\frac{d\bar p}{d\bar\psi}
 (\mathbf b\times\nabla_N\bar H)\cdot\nabla_N\alpha_c.
$$

스캔에서는 예를 들어 $a_N=a_*, B_N=B_{\rm ref}>0$를 일관되게 지정한다. 양의 규격화를 바꾸어도 안정성의 부호는 같지만, 고유값의 크기를 비교하려면 규격화도 기록한다.

## 7. 평형 전체의 안정성

$$
\frac{d}{d\theta}\left(g\frac{dX}{d\theta}\right)+cX=\lambda fX,
\qquad X(-\Theta)=X(\Theta)=0.
$$

현재 matrix 부호 규약은 $(C-G)x=\lambda Fx$이고 최대 고유값을 반환한다. $g,f>0$에 대해

$$
\lambda_{\max}=\max_X
\frac{\int(cX^2-gX_\theta^2)\,d\theta}{\int fX^2\,d\theta}.
$$

$\lambda_{\max}>0$이면 불안정하다. 무한 field line을 유한 구간으로 근사하므로, 특히 안정 경계 근처에서 $\Theta$와 grid 해상도의 수렴을 검사한다.

평형 전체의 대표값은 다음과 같다.

$$
\Lambda(\mathcal A,I,p_0)=\max_{u,\theta_0}\lambda_{\max}(u,\theta_0).
$$

한 자속면 또는 위상이라도 불안정하면 해당 평형을 불안정으로 판정한다. 초기 설정 제안은 다음과 같다.

- radial labels: $u=0.03:0.04:0.95$ 및 필요에 따라 $0.975,0.99$.
- phase: $\theta_0\in[-\pi,\pi)$의 8점으로 시작하고 경계 근처에서 16점 이상과 국소 최적화를 검토.
- field-line domain: $\Theta=5\pi$로 시작하고 $10\pi,20\pi$와 비교.
- field-line grid: turn당 약 24 cells로 시작하고 48 이상과 비교. 영역을 확장할 때도 turn당 cells 수를 유지.
- surface contour: 512점, FourierModes=32로 시작하고 48,64와 비교.

이는 제안값이며 논문이 정한 값이나 이번에 확인한 필요 해상도는 아니다. 가장 불안정한 자속면 주변, 축 근처, 경계 근처에서 radial sampling을 추가한다.

각 평형에서 $\Lambda$, $u_{\rm worst}$, $\theta_{0,\rm worst}$와 eigenfunction을 저장한다. 초기 스캔은 양의 고유값이 발견되면 불안정 판정을 끝내도 된다. 정확한 $\Lambda$와 가장 불안정한 자속면을 구할 때는 전체 탐색을 완료한다.

`qCorrection`에 의한 $B_p$ 보정과 periodicity에 의한 $FF'$ 보정이 오차를 가릴 수 있다. 보정 전 $q$ 적분 오차, `qCorrection-1`, `closure_input`, `FFprime_periodic-FFprime_input`를 함께 기록하고 mesh/Fourier 해상도를 높이면 개선되는지 확인한다. $q'$와 지정 프로파일의 해석적 미분이 일치하는지도 검사한다.

## 8. Chance형 및 Grimm형 스캔과 모든 경계

먼저 현재 fixed-q test의 $\delta_0=-0.4,I=0.8$ 경계를 단일 케이스로 사용한다. 이 경계를 그대로 유지하려면 $\kappa_*=4.17813$, $\mathcal A=4.82996$, $a_*=0.57442$를 기준값으로 잡는다. 이 값에서 정규화 식은 기존 구성 매개변수 $R_c=3,a_0=0.8,\kappa_0=3$을 다시 만든다. 실제 생산 코드에서는 반올림된 숫자 대신 기준 경계에서 측정한 값 자체를 사용한다.

첫 형상 스캔에서는 이 실제 elongation과 종횡비를 고정하고 $I$를 바꾼다. 이후 실제 종횡비를 4로 바꾼 비교 스캔을 추가한다. 별도로 $\kappa_g=3,\mathcal A=4$를 원하는 경우도 같은 정규화 식을 사용할 수 있지만 이는 기준 경계와 다른 elongation이다.

- $I$의 pilot grid: $0,0.25,0.45,0.6,0.8,1.0,1.15$.
- pressure control: $\Pi_0=2\mu_0p_0/B_{\rm ref}^2$. 평균 beta와는 다른 중심 압력 제어값이다.
- $\Pi_0$의 pilot grid: $0.02,0.1,0.3,0.6,1.0$ 등. 계산된 $\beta_V$를 보고 범위를 조정해 그림의 상한까지 조사한다.
- 전체 radial/phase 스캔 전에 낮은 압력의 여러 형상에서 GS/q와 Direct point solve의 수렴을 확인한다.

Chance형 그림은 $(i_{\rm out},100\beta_V)$ 또는 $(I,100\beta_V)$에 안정/불안정 점을 표시한다. 가로축 정의를 명시한다.

같은 $I$에서 압력을 증가시키며 안정 -> 불안정, 불안정 -> 안정, 고압력에서의 재불안정 전이를 모두 찾는다. 각 전이에서 $p_0$를 bracket하고 평형을 다시 풀면서 경계를 세분화한다. 유한 영역의 $\Lambda$가 0 근처에 모이는 경우 임의의 epsilon contour만으로 경계를 정하지 않고, 영역과 mesh 해상도에 따른 불확실성을 띠로 남긴다. 곡선의 접힘 근처는 $(I,p_0)$의 국소 2차원 세분화 또는 pseudo-arclength continuation을 사용한다.

Grimm형 그림은 같은 형상 계열과 프로파일로 $\mathcal A=3,4,7,10$의 모든 경계를 겹쳐 그린다. 실제 $a_g,\kappa_g$를 고정하고 $R/a=4$를 점선으로 표시한다.

논문의 elongation과 맞춘 비교 계열로 $\kappa_g=1.386$, $\delta_0=-0.4$에서도 같은 스캔을 수행한다. 현재의 높은 elongation C-shape 계산을 먼저 확인한 뒤 elongation을 바꾸면 실패 원인을 찾기 쉽다.

평형 미수렴, 기하 미수렴, ballooning solve 실패는 unknown으로 기록한다. 조사한 beta 범위 전체에서 안정하더라도 그 상한을 넘는 모든 압력에서의 안정성을 주장하지 않는다.

## 9. 추가할 driver와 저장 항목

다음 기능을 구현하고 기존 core solver에 연결했다. 시작 설정, 단일 평형 평가, 그림 생성 및 검증 함수도 추가했으며 실행 예제는 구현 README를 참고한다.

| function | role |
|---|---|
| `normalize_indented_cshape_boundary` | 경계식을 실제 $\mathcal A,\kappa_g,a_g$에 맞추고 $i_{\rm out},\delta_g$ 등을 반환 |
| `compute_volume_beta` | P3 quadrature로 총 자기장을 이용한 평균 beta 계산 |
| `evaluate_equilibrium_ballooning` | 모든 radial labels와 phase를 평가하고 최대 고유값 및 진단값 반환 |
| `scan_indented_cshape_beta` | 형상 및 압력 continuation과 checkpoint 관리 |
| `solve_indented_cshape_equilibrium` | 압력 간격 축소와 동일 압력 q 반복 재개, GS/q 수렴한 predictor만 채택 |
| `trace_indented_cshape_marginal` | 모든 전이를 세분화하고 Chance형/Grimm형 그림 생성 |

```text
for aspect ratio A
  for construction indentation I
    normalize the boundary at fixed actual a and kappa
    prepare/reuse this geometry's FEM mesh and operator
    for pressure amplitude p0 along a continuous equilibrium branch
      solve prescribed-p/q GS equilibrium (mode 3)
      reject/mark unresolved convergence failures
      compute total-field beta
      build FEM adapter once for this equilibrium
      for radial label u
        prepare actual Direct geometry and physical p', q', FF'
        for ballooning phase theta0
          build g,c,f and solve the largest eigenvalue
      store the worst eigenvalue, surface, phase and diagnostics
    locate and refine every stability transition
  save all branches for this aspect ratio
plot stable/unstable points and marginal curves
```

checkpoint의 최소 항목은 $\mathcal A,\kappa_g,\delta_0,\delta_g,I,i_{\rm out},p_0,\beta_V,I_p,\Delta\psi,R_{\rm axis},\Lambda,u_{\rm worst},\theta_{0,\rm worst}$, inner/q residual, 모든 해상도, normalization, 기하 보정량, status이다.

유용한 추가 그림은 (1) 압력에 따른 $\lambda_{\max}(u)$의 변화, (2) 가장 불안정한 자속면의 $B_p$, curvature, integrated shear, (3) 경계와 내부 자속면 및 자기축 shift이다. 두 번째 안정 영역의 유무를 평형의 변화와 연결해서 해석할 수 있다.

논문의 inboard 효과와 대조하려면 같은 기하학적 box에서 $R_b^{\rm in}=2R_g-R_b^{\rm out}$의 반사 형상을 별도 케이스로 만들어 GS부터 다시 계산할 수 있다. 반사하면 polygon orientation도 반전되므로 점 순서를 조정한다. 이 반사 형상은 선택적인 대조군이다.

## 10. 확인한 구현과 논문

- [Indented boundary](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/local_equilibrium/Indented_Miller_boundary.m)
- [Indented prescribed-q setup](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/GS_fixed/tests/integration/test_qsolve_indented_fixed_q.m)
- [Prescribed-p/q solver](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/GS_fixed/GS/solvers/solve_pq_constraint.m)
- [Pressure/standard beta postprocessing](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/GS_fixed/GS/post/finalize_GS_equilibrium.m)
- [FEM adapter](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/GS_fixed/GS/post/build_ballooning_adapter.m)
- [Direct surface preparation](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/salpha/prepare_salpha_equilibrium.m)
- [Direct geometry](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/local_equilibrium/FluxSurface_Mercier_Luc.m)
- [Ballooning coefficients](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/local_equilibrium/evaluate_gcf.m)
- [Direct point verification example](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/GS_fixed/verification/verify_indented_direct_ballooning_point.m)
- [Eigenvalue convention](C:/Users/DaeYeolPyo/Desktop/ballooning_mode/ballooning_mode/sturm_liouville_solver/calculate_eigenvalue.m)
- [Chance publisher abstract: inboard indentation](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.51.1963)
- [Chance local paper](<C:/Users/DaeYeolPyo/Desktop/CShape/[1983] M.S.Chance et al., Ballooning mode stability of bean-shaped cross sections for high-beta tokamak plasmas.pdf>) (pp. 1963-1965)
- [Grimm local paper](<C:/Users/DaeYeolPyo/Desktop/CShape/[1985] R.C.Grimm et al., MHD stability properties of bean-shaped tokmaks.pdf>) (pp. 807-809)
