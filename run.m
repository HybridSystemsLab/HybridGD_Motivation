%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clearvars
close all
clc

%% Parameters
theta = [1; 1];

gammac = 1; % 0;
gammad = 1; % inf;
lambdac = 0.1; % 0.01;
lambdad = 0.5; % 0.1;

parameters.theta = theta;
parameters.gammac = gammac;
parameters.gammad = gammad;
parameters.lambdac = lambdac;
parameters.lambdad = lambdad;
parameters.R = 0;
parameters.nu = 0;

%%
sys = hybridPE(parameters);

%% simulation horizon
tspan = [0 35];
jspan = [0 1000];
config = HybridSolverConfig('RelTol', 1e-3, 'AbsTol', 1e-4,  ...
                            'MaxStep', 0.1, 'Refine', 4);

%% simulate
x0 = [3; 6];
thetahat0 = [0; 0];
psi0 = zeros(2,2);
eta0 = -x0;
u0 = 0;
z0 = [x0;thetahat0;psi0(:);eta0;u0];

sys.R = 0;
sys.nu = 0;
solc = sys.solve(z0, tspan, jspan, config);

sys.R = 1;
sys.nu = 0;
sold = sys.solve(z0, tspan, jspan, config);

sys.R = 2;
sys.nu = 0;
solh = sys.solve(z0, tspan, jspan, config);

psi0 = zeros(2,2);
eta0 = -0.5*x0;
eps0 = x0 + eta0 - psi0*theta;
z0 = [x0;thetahat0;psi0(:);eta0;u0];
sole = sys.solve(z0, tspan, jspan, config);

sys.R = 0;
sys.nu = 1;
solnu = sys.solve(z0, tspan, jspan, config);

%%
HybridPlotBuilder.defaults.reset()

% journal
HybridPlotBuilder.defaults.set('flow line width', 3, ...
                               'jump line width', 2,...
                               'label size', 33,...
                               'tick label size', 30,...
                               't_label', '$t$ [s]')
legendFontSize = 26;

% dissertation
HybridPlotBuilder.defaults.set('flow line width', 3, ...
                               'jump line width', 2,...
                               'label size', 23,...
                               'tick label size', 23,...
                               't_label', '$t$ [s]')   
legendFontSize = 22;

orange = [1.0000    0.5984    0.2000];

figure; clf;
hpb = HybridPlotBuilder();
hpb.flowColor('blue')...
   .labels('$x$')...
   .legend({'$x_1$'},'FontSize',legendFontSize,'Orientation','horizontal')...
   .plotFlows(solh.select(1))
hold on
hpb.flowColor('green')...
   .legend({'$x_2$'},'FontSize',legendFontSize,'Orientation','horizontal')...
   .plotFlows(solh.select(2))
grid on; box on;
xlim([0, 10*pi+0.5])
xticks([0 2*pi 4*pi 6*pi 8*pi 10*pi])
xticklabels({'$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$','$10\pi$'})

pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), pos(4)])
set(gca, 'LooseInset', get(gca,'TightInset'))

figure; clf;
hpb = HybridPlotBuilder();
hpb.flowColor('blue')...
   .legend({'Continuous'},'FontSize',legendFontSize,'Orientation','horizontal')...
   .plotFlows(solc.select([3,4]),@(x) norm(theta - x))
hold on
hpb.flowColor('green')...
   .legend({'Discrete'},'FontSize',legendFontSize,'Orientation','horizontal')...
   .plotFlows(sold.select([3,4]),@(x) norm(theta - x))
hold on
hpb.flowColor(orange)...
   .legend({'Hybrid'},'FontSize',legendFontSize,'Orientation','horizontal')...
   .labels('$|\theta - \hat \theta|$')...
   .plotFlows(solh.select([3,4]),@(x) norm(theta - x))
grid on;
xlim([0, 10*pi+0.5])
ylim([0, 1.5])
xticks([0 2*pi 4*pi 6*pi 8*pi 10*pi])
xticklabels({'$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$','$10\pi$'})
yticks([0 0.5 1 1.5])
box on;

pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), pos(4)])
set(gca, 'LooseInset', get(gca,'TightInset'))


%%
eps = @(x) x(1:2) + x(9:10) - reshape(x(5:8),2,2)*theta;

figure; clf;
hpb = HybridPlotBuilder();
hpb.flowColor(orange)...
   .legend({'$\varepsilon(0,0) \not = 0$ noisy'},'FontSize',legendFontSize,'Location','NorthEast')...
   .plotFlows(solnu,@(x) norm([theta - x(3:4); eps(x)]))  
hold on
hpb.flowColor('blue')...
   .legend({'$\varepsilon(0,0) = 0$'},'FontSize',legendFontSize,'Location','NorthEast')...
   .plotFlows(solh,@(x) norm([theta - x(3:4); eps(x)])) 
hold on
hpb.flowColor('green')...
   .legend({'$\varepsilon(0,0) \not = 0$'},'FontSize',legendFontSize,'Location','NorthEast')...
   .labels('$|\xi|_{\mathcal{A}_g}$')... 
   .plotFlows(sole,@(x) norm([theta - x(3:4); eps(x)])) 
grid on; box on;
xlim([0, 10*pi+0.5])
ylim([0, 4])
xticks([0 2*pi 4*pi 6*pi 8*pi 10*pi])
xticklabels({'$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$','$10\pi$'})
yticks([0 1 2 3 4])
hpb.reorderLegendEntries([2 3 1])

pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), pos(4)])
set(gca, 'LooseInset', get(gca,'TightInset'))

%% compute \mu
idx = find(solh.j==0,1,'last');
for i = 1:idx
    temp = reshape(solh.x(i,5:8),2,2);
    psiTpsi(:,:,i) = temp'*temp;
end
temp = reshape(solh.x(idx+1,5:8),2,2);
psiplusTpsiplus = temp'*temp;

M = trapz(solh.t(1:idx),psiTpsi,3) + psiplusTpsiplus;
mu = min(eigs(M));

