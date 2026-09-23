%% test_FFprime_model.m
% Unit checks for the direct reduced-order FFprime representation.

clearvars;
close all;
clc;

thisDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fullfile(thisDirectory,'..','..');
addpath(rootDirectory);
setup_solve_GS;

psiN = linspace(0,1,101).';

%% Basis properties
basis = construct_FFprime_basis(psiN,4,true);

assert(norm(basis(:,1)-(1-psiN),inf) < 1.0e-14, ...
    'TEST_FFPRIME_MODEL:LeadingBasis', ...
    'The leading FFprime basis is not 1-psiN.');

assert(norm(basis(end,:),inf) < 1.0e-14, ...
    'TEST_FFPRIME_MODEL:BoundaryBasis', ...
    'The zero-edge FFprime basis does not vanish at psiN=1.');

%% Exact reference family
model = struct();
model.representation = 'direct-FFprime-reduced-basis';
model.psiN = psiN;
model.nBasis = 3;
model.enforceZeroEdge = true;
model.coefficients = [2;0;0];
model.Fboundary = 9;
model.signF = 1;

profile = evaluate_FFprime_model(model,psiN);
referenceFFprime = 2*(1-psiN);
referenceG = 81-2*(1-psiN).^2;

assert(norm(profile.FFprime-referenceFFprime,inf) < 1.0e-13, ...
    'TEST_FFPRIME_MODEL:FFprime', ...
    'The direct FFprime model does not reproduce the reference family.');

assert(norm(profile.G-referenceG,inf) < 2.0e-4, ...
    'TEST_FFPRIME_MODEL:GReconstruction', ...
    'G reconstructed from direct FFprime is incorrect.');

assert(abs(profile.F(end)-9) < 1.0e-14, ...
    'TEST_FFPRIME_MODEL:BoundaryF', ...
    'The reconstructed profile does not preserve boundary F.');

%% Initialization from a supplied FFprime profile
problem = struct();
problem.profile = struct('psiN',psiN);
problem.dim = struct('R0',3,'B0',3,'aMinor',0.8);

profileInput = struct();
profileInput.normalization = struct('mode','fixed_Fb','Fb',9);
profileInput.initial = struct( ...
    'type','FFprime', ...
    'psiN',psiN, ...
    'FFprime',referenceFFprime);

pprime = -1.0e6*(1-psiN.^1.5).^2;
qTarget = 1+2*psiN;
options = struct('nBasis',3,'enforceZeroEdge',true);

initialized = initialize_FFprime_model( ...
    problem,profileInput,pprime,qTarget,options);
initializedProfile = evaluate_FFprime_model(initialized,psiN);

assert(norm(initializedProfile.FFprime-referenceFFprime,inf) ...
        < 1.0e-10, ...
    'TEST_FFPRIME_MODEL:Initialization', ...
    'FFprime initialization did not recover an in-basis profile.');

fprintf('Direct FFprime model tests passed.\n');
