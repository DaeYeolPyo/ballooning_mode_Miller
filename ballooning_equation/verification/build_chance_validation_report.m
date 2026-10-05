function report=build_chance_validation_report
%BUILD_CHANCE_VALIDATION_REPORT Publish actual results, including failed PRL match.
    root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_indented_cshape_scan;
    temp=fullfile(root,'tmp','chance_prl');
    output=fullfile(root,'output','chance_prl_figure3','validation');
    if ~isfolder(output), mkdir(output); end
    low=load(fullfile(temp,'low_results_corrected.mat'));
    first=load(fullfile(temp,'first_boundary_025.mat'));
    high25=load(fullfile(temp,'high025_mesh64_boundary384.mat'));
    high325=load(fullfile(temp,'high0325_mesh64_boundary384.mat'));
    bal25=load(fullfile(temp,'high025_mesh64_stability_refined.mat'));
    bal325=load(fullfile(temp,'high0325_mesh64_stability_refined.mat'));
    finest=load(fullfile(temp,'high025_mesh96_boundary768.mat'));
    finestBal=load(fullfile(temp,'high025_mesh96_quality_checked.mat'));
    theta=load(fullfile(temp,'theta_resolution_025.mat'));
    fourier=load(fullfile(temp,'fourier_resolution_025.mat'));
    finestValues=load(fullfile(temp,'high025_mesh96_resolved.mat'));
    eigen=load(fullfile(temp,'native_eigenvalue_check.mat'));
    points=cell(8,1); equilibria=cell(8,1);
    depths=[.11 .11 .325 .25 .25 .25 .325 .25];
    expected=["stable" "unstable" "stable" "stable" "unstable" "stable" "stable" "stable"];
    grids=[48 64;48 64;48 64;48 64;48 64;64 96;64 96;96 128];
    boundaries=[128 128 128 128 128 384 384 768];
    for k=1:3
        equilibria{k}=low.lowResults{k}.equilibrium;
        points{k}=point(equilibria{k},low.lowResults{k}.ballooning,depths(k));
    end
    for k=4:5
        equilibria{k}=first.result.scan.equilibria{k-3};
        points{k}=point(equilibria{k},first.result.scan.points{k-3}.ballooning,depths(k));
    end
    equilibria{6}=high25.result.equilibrium; points{6}=point(equilibria{6},bal25.refined,.25);
    equilibria{7}=high325.result.equilibrium; points{7}=point(equilibria{7},bal325.refined,.325);
    equilibria{8}=finest.result.equilibrium; points{8}=point(equilibria{8},finestBal.checked,.25);
    for k=1:8
        points{k}.paperStatus=expected(k); points{k}.gridSize=grids(k,:);
        points{k}.nBoundary=boundaries(k);
        points{k}.matchesPaper=points{k}.status==expected(k);
    end
    sourceConfig=make_chance_prl_figure3_config;
    report=struct('source',sourceConfig.PaperReference, ...
        'points',{points},'equilibria',{equilibria}, ...
        'firstBoundary',first.result.trace.roots, ...
        'firstBoundaryHistory',first.result.history, ...
        'thetaConvergence',theta.result,'fourierConvergence',fourier.result, ...
        'finestResolvedSurfaces',finestValues,'eigenvalueCheck',eigen, ...
        'reproductionSucceeded',false, ...
        'scope','Representative paper points, one first-boundary bracket, and independent resolution checks', ...
        'version',2);
    regression=load(fullfile(temp,'cshape_regression_after_benchmark.mat'));
    report.cshapeRegression=regression.regression;
    save(fullfile(output,'validation.mat'),'report','-v7.3');

    % Copy the existing reference overlay into a two-panel scientific figure.
    displayed=report; displayed.points=points([1:5 7:8]);
    sourceFigure=plot_chance_prl_figure3(displayed);
    sourceAxes=findobj(sourceFigure,'Type','axes');
    fig=figure('Color','w','Visible','off','Position',[100 100 1200 520]);
    layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(layout); copyobj(allchild(sourceAxes(1)),ax); hold(ax,'on');
    r=report.firstBoundary(1); center=mean(r.betaBracket);
    errorbar(ax,.25,center,center-r.betaBracket(1),r.betaBracket(2)-center,'kd', ...
        'MarkerFaceColor','k','MarkerSize',6,'DisplayName','Computed first boundary');
    xlim(ax,[.08 .35]); ylim(ax,[0 26]); grid(ax,'on');
    xlabel(ax,'Inboard indentation d/(2a)'); ylabel(ax,'Total-field volume beta [%]');
    title(ax,'PRL Fig. 3: sampled comparison'); legend(ax,'show','Location','northeast','FontSize',8);
    ax2=nexttile(layout); hold(ax2,'on');
    for k=1:numel(theta.result.psiN)
        semilogx(ax2,theta.result.cellsPerTurn,theta.result.lambda(k,:),'o-', ...
            'LineWidth',1.4,'DisplayName',sprintf('psiN = %.2f',theta.result.psiN(k)));
    end
    set(ax2,'XScale','log'); yline(ax2,0,'k--','HandleVisibility','off'); grid(ax2,'on');
    xlabel(ax2,'Ballooning cells per 2*pi'); ylabel(ax2,'Largest eigenvalue (our normalization)');
    title(ax2,'High beta 14.3505%: stability lost on refinement'); legend(ax2,'Location','southeast');
    sgtitle(fig,'First boundary close to PRL; high-beta reproduction failed','FontWeight','bold');
    exportgraphics(fig,fullfile(output,'figure3_validation.png'),'Resolution',180);
    close(sourceFigure); close(fig);

    file=fullfile(output,'validation_report.md'); fid=fopen(file,'w','n','UTF-8');
    cleanup=onCleanup(@()fclose(fid));
    line('# Chance PRL Figure 3 재현 테스트 — 고 beta 검증 실패');
    line('');
    line('논문의 Eq. (2) bean 경계와 p/q 프로파일을 native GS → direct geometry → ballooning 경로에서 실제 계산했다. 첫 안정 경계는 논문의 경계 근방이지만, 높은 beta의 안정 분류는 일치하지 않는다. 따라서 현재 workflow가 Figure 3 전체를 재현했다고 판단할 수 없다.');
    line('');
    line('## 논문 입력'); line('');
    line('A=1, B=0.6, R/a=4, b/a=1.386. C로 inboard indentation d/(2a)를 조절하고 E, xbar를 다시 정한다.');
    line(''); line('```text');
    line('rho=A(1+B cos t), gamma=C sin t');
    line('R=xbar+rho cos gamma, Z=E rho sin gamma');
    line('u=psiN (axis 0, edge 1), p=p0(1-u^2)^2');
    line('q=1.03+0.84375u-1.1775u^2+3.50375u^3');
    line('beta=2 mu0 integral(p dV)/integral((Bp^2+Bphi^2) dV), dV=2 pi R dR dZ');
    line('```'); line('');
    line('이 q는 q(0)=1.03, q(1)=4.2, q_u(0)=0.84375, q_u(1)=9를 만족한다. p/q의 변수는 sqrt(psiN)이 아니다. 안정 지도에는 central beta 대신 계산한 total-field volume beta를 사용했다.');
    line(''); line('## 실제 결과'); line('');
    line('| d/(2a) | volume beta [%] | GS grid / boundary points | 계산 분류 | 논문 분류 | lambda observed | max q relative error |');
    line('|---:|---:|:---|:---|:---|---:|---:|');
    for k=1:numel(points)
        p=points{k};
        fprintf(fid,'| %.3f | %.6f | %d x %d / %d | %s | %s | %+.8e | %.3e |\n', ...
            p.indent,p.betaPercent,p.gridSize,p.nBoundary,p.status,p.paperStatus,p.lambdaMaxObserved,p.qMaximumRelativeError);
    end
    line('');
    fprintf(fid,'d/(2a)=0.25의 첫 경계는 beta **[%.6f, %.6f]%%**에 bracket되었다. 논문 그림의 약 2%%대 첫 경계와 가깝다. 이 경계는 고유값 부호의 bisection으로 얻었으며 reference curve에 맞춰 정하지 않았다.\n',r.betaBracket);
    line('');
    line('stable은 요청한 표면/phase에서 안정이라는 뜻이다. 표면 0.05~0.95, 보통 8~16 phase를 조사했다. 가장 촘촘한 high-beta 점의 양의 고유값은 psiN=0.90, theta0=0에서 확인한 instability witness이며 전체 최대 고유값의 하한이다. 해당 두 표면 0.90/0.95는 모든 geometry quality gate를 통과했다.');
    line('');
    line('## 높은 beta의 차이는 단순한 coarse eigenvalue 문제로 해소되지 않았다'); line('');
    line('- 고 beta에서 24 cells/turn은 음의 고유값을 주지만, 48 이상에서는 양의 값이 나온다. 768→1536 refinement에서도 양의 값이 유지된다.');
    line('- Fourier 32~48에서 high-beta 양의 고유값이 거의 유지된다. 작은 mode 수를 선택해 논문의 stable label에 맞추지 않았다.');
    line('- GS 64x96 / boundary 384에서 96x128 / boundary 768로 높여도 volume beta와 양의 고유값이 유지된다.');
    fprintf(fid,'- finest native GS residual %.3e, prescribed-q max relative error %.3e. psiN=0.90에서 lambda=%+.9e.\n', ...
        points{8}.GSResidual,points{8}.qMaximumRelativeError,points{8}.lambdaMaxObserved);
    fprintf(fid,'- 동일한 native bean 행렬에서 dense와 sparse largest eigenvalue 차이 %.3e. largest를 nearest-to-zero mode로 잘못 고른 결과가 아니다.\n',abs(eigen.dense(1)-eigen.sparse));
    line('');
    line('GS/q residual과 표면 closure 통과만으로 논문의 ballooning 분류까지 검증된 것은 아니다. 남은 차이의 원인은 이번 테스트에서 확정하지 못했다. 논문 Eq. (3)와 독립적으로 local shear/metric/coefficient 전달을 검증하고, 인접 flux surface로 계산한 shear 및 실제 GS force balance와 대조하는 작업이 필요하다. 이 결과를 물리적인 second stability 부재로 해석할 수 없다.');
    line(''); line('## 구현 중 확인하고 수정한 문제'); line('');
    line('1. q를 .95까지만 제약한 상태에서 .95의 안정성을 평가하면 실제 dq/du가 제대로 제약되지 않을 수 있었다. paper benchmark는 .975까지 fit한다.');
    line('2. extract_flux_surface가 열린 contour를 강제로 닫힌 것으로 처리했다. 원래 endpoint로 closure를 검사하고 열린 contour를 거부하도록 수정했다.');
    line('3. contourc의 roundoff 크기 중복 점이 arclength PCHIP을 왜곡했다. 분석적인 원 테스트에서 확인하고 제거했다.');
    line('4. 압력 변화에 대한 FFprime 초기 예측과 충분한 contour grid가 high-beta native GS/q 수렴에 필요했다. 초기 예측은 최종 수렴 조건을 대체하지 않는다.');
    line('5. 모든 고유값을 구하던 dense 계산 대신 optional largest-only mode를 추가했다. 양/음의 가장 큰 고유값을 Cholesky bracket과 shift-invert로 구하며, 기존 기본 호출은 전체 dense spectrum을 반환한다. 분석해/여러 boundary condition/dense native 행렬과 비교했다.');
    line(''); line('## 실행 및 데이터'); line('');
    line('Bean geometry/q endpoint, analytic closed/open contour, stable/unstable largest eigenvalue, 여러 boundary condition과 dense 비교 테스트를 통과했다. 공통 코드 수정 후 실제 C-shape 2 surfaces x 2 phases, checkpoint resume, unknown retry, 수정 전 schema 거부 검사도 통과했다. MATLAB Code Analyzer와 git diff whitespace 검사에 문제는 없었다.');
    line('');
    line('대표점, 첫 경계, theta/Fourier convergence를 수행한 테스트다. 5개 indentation x 4개 beta의 전체 map이나 second-stability 경계는 검증 완료된 결과로 제공하지 않는다.');
    line(''); line('`validation.mat`에는 native equilibria, 실제 계산점, 첫 경계 history, convergence table, dense/sparse comparison이 들어 있다. `figure3_validation.png`는 reference overlay와 theta refinement를 보여 준다.');
    line(''); line('```matlab'); line('setup_indented_cshape_scan'); line("addpath('verification')");
    line('test_chance_prl_geometry'); line('test_flux_surface_closedness'); line('test_largest_ballooning_eigenvalue');
    line('report=build_chance_validation_report;'); line('```');
    line(''); line('geometry/profile 설정과 full-scan driver는 chance_prl_figure3_README.md를 참조한다. Scan checkpoint schema는 2이며, 수정 전 안정성 결과가 든 schema 1 checkpoint는 재사용하지 않는다.');
    line(''); line('Source: Chance, Jardin & Stix, PRL 51, 1963–1966 (1983), DOI 10.1103/PhysRevLett.51.1963, supplied local PDF. Dashed reference boundaries are approximate visual readings of Figure 3.');
    fprintf('Wrote %s\n',output);

    function line(value), fprintf(fid,'%s\n',value); end
end

function p=point(eq,bal,indent)
    beta=compute_volume_beta(eq);
    maximumClosure=max(cellfun(@(d)d.thetaPsi.inputClosureDefect,bal.surfaceDiagnostics));
    p=struct('indent',indent,'pressureAmplitude',eq.scanPressureAmplitude, ...
        'betaPercent',beta.betaPercent,'status',bal.status,'complete',bal.complete, ...
        'lambdaMax',bal.lambdaMax,'lambdaMaxObserved',bal.lambdaMaxObserved, ...
        'qMaximumRelativeError',eq.constraintResult.finalMaxAbsRelativeQ, ...
        'GSResidual',eq.solverResult.finalRelativeResidual,'maximumClosure',maximumClosure, ...
        'radialLabels',bal.psiN,'theta0',bal.theta0,'ballooningOptions',bal.options);
end
