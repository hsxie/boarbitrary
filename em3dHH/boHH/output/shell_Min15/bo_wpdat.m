% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)

% wpdat=[0.44,0,0.05771i; % for loss cone mirror
%     0.22,0,0.05041i;
%     0.50,0,-0.04363i;
%     0.47,0,-0.1578i;
%     0.41,0,-0.2965i;
%   ];

wpdat=[6.25,0,0.043i; % Irvine2018, theta=89.5
    7.325,0,0.098i;
    8.47,0,0.142i;
    9.52,0,0.173i;
    10.7,0,0.19i;
    11.72,0,0.191i;
    13.08,0,0.105i;
    18.2,0,0.218i;
    19.48,0,0.239i;
  ];

jselc=0; % jselc=1, alway plot the most unstable dispersion surface
% wpdat=[];

%%
run ../modules/bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
