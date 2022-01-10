%%% Exact & Approximate Riemann Solver
%%% Developer: Mauro Rodriguez Jr. (mauro_rodriguez@brown.edu)
clear all; close all; clc
approximate = 'src/approximate_solver';
addpath(approximate);
exact_solver = 'src/exact_solver';
addpath(exact_solver);
% ratio of specific heats
g = 1.40;
% Right state
rhoR = 1.204; uR = 0; pR = 101325;
WR = [rhoR, uR, pR];
% Left state
rhoL = rhoR*1.8105824688; uL = 225.894320023; pL = pR*2.35438333333;
WL = [rhoL, uL, pL];
% initial guess for the exact Riemann solver
pguess = 0.5;
% error tolerance for the exact Riemann
tol = 1e-15;
% left distance of the domain
LX0 = 0; 
% right distance of the domain
LX1 = 1; 
% shifting length by 0.5 to be cell centered
dloc = 0.5; 
% Number of points for the exact Riemann solver
N = 1000; 
% End time of the simulation (non-dimensional)
T = 0.0003;
% CFL condition for stability of the approximate Riemann solver
cfl = 0.9; 
% number of computational points for approximate Riemann solver
NX = 200;
% flux limiter, 0: superbee, 1: minmod
flim = 0;
% root finding algorithm for p_star
ps = root_find(pguess,tol,g,WL,WR);
% plotting results
figure(1)
% exact Riemann solver evaluation and plotting
[Iex] = plot_rp(ps, g, WL, WR, LX0, LX1, dloc, N, T);
% approximate Riemann solver evaluation and plotting
Inum = plot_roe_eue(g,WL,WR,LX0,LX1,NX,T,cfl,flim);
% legend values
v1 = 'Exact';
v2 = '2nd Roe, superbee';
L = legend(v1,v2);
set(L,'Interpreter','Latex','FontSize',14);
rmpath(approximate);
rmpath(exact_solver);