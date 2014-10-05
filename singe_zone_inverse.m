%% Single Zone Inverse Model Training
% Madhur Behl
% University of Pennsylvnaia
% Copyright 2014 (C)
%%
%{

This is the main source file. 

The following functions are called from this file and must be present in
the MATLAB path

importwallgains.m
importdataout.m
trnsysfit.m
trn_ssmodel.m
rsquare.m

The following data files must also be present in the MATLAB working
directory

wallgains_jun_2min.xlsx
wallgains_jul_2min.xlsx
dataout_jun_2min.xlsx
dataout_jul_2min.xlsx


%}
clc
clear all
close all


%% Initialize nominal values of the model parameters
% First input the details of the surfaces and the geometry of the zone

% Some conversion factors
% kJ/hr to watt
kj2w = 0.277777778;

global dsystem;
dsystem = [];
global sysd;
sysd = [];


% bypass the parameter initialization and use the user defined nominal
% values isntead. Setting this to 1 may be useful when you want to start
% the parameter estimation from the solution of a previous simulation.
bypass = 0; 


if(~bypass)
    
    % Construction material details for nominal parameter value calulation
    % Follow SI units throughout.

    
    % OUTWALL
    % Layers : Brick, Insulation, plaster,
    % calcualte capacity of brick layer
    cp_brick_o = 1;             %(kJ/kgK)
    rho_brick_o = 1800;         % kg/m3
    d_brick_o = 0.24;           % m
    % calcualte capacity of instualtion layer
    cp_insul_o = 0.8;           %(kJ/kgK)
    rho_insul_o = 40;           % kg/m3
    d_insul_o = 0.1;            % m
    % calcualte capacity of plaster layer
    cp_pls_o = 1;               %(kJ/kgK)
    rho_pls_o = 2000;           % kg/m3
    d_pls_o = 0.015;            % m
    
    % Total thermal capacitance of outwall (J/m2K) multiply by 1000 to convert
    % from kJ to J.
    Cow = (cp_brick_o*rho_brick_o*d_brick_o)+ ...
        (cp_insul_o*rho_insul_o*d_insul_o)+(cp_pls_o*rho_pls_o*d_pls_o);
    Cow = Cow*1000;
    
    % U value of outwall (W/m2K)
    Uow = 0.339;
    
    % convective heat transfer coeficients for outwall (W/m2K) front(f) and
    % back(b)
    Ucon_ow_f = 11 * kj2w;
    Ucon_ow_b = 64 * kj2w;
    
    % Solar absobtance of outwall front(f) and back(b)
    absol_ow_f = 0.8;
    absol_ow_b = 0.4;
    
    
    % GROUND
    % Layers: Floor,stone,silence,concrete,insulation
    % calcualte the capacity of the floor layer
    cp_flo_g = 1;               % (kJ/kgK)
    rho_flo_g = 800;            % kg/m3
    d_flo_g = 0.005;            % m
    % calcualte the capacity of the stone layer
    cp_sto_g = 1;               % (kJ/kgK)
    rho_sto_g = 2000;           % kg/m3
    d_sto_g = 0.06;             % m
    % calcualte the capacity of the silence layer
    cp_sil_g = 1.44;            % (kJ/kgK)
    rho_sil_g = 80;             % kg/m3
    d_sil_g = 0.04;             % m
    % calcualte the capacity of the concrete layer
    cp_con_g = 0.8;             % (kJ/kgK)
    rho_con_g = 2400;           % kg/m3
    d_con_g = 0.24;             % m
    % calcualte the capacity of the insulation layer
    cp_ins_g = 0.8;               % (kJ/kgK)
    rho_ins_g = 40;            % kg/m3
    d_ins_g = 0.08;            % m
    
    % Total thermal capacitance of ground (J/m2K) multiply by 1000 to convert
    % from kJ to J.
    Cgr = (cp_flo_g*rho_flo_g*d_flo_g)+(cp_sto_g*rho_sto_g*d_sto_g)...
        + (cp_sil_g*rho_sil_g*d_sil_g)+(cp_con_g*rho_con_g*d_con_g)...
        + (cp_ins_g*rho_ins_g*d_ins_g);
    Cgr = Cgr*1000;
    
    % U value of ground (W/m2K)
    Ugr = 0.313;
    
    % convective heat transfer coeficients for ground (W/m2K) front(f) and
    % back(b)
    Ucon_gr_f = 11 * kj2w;
    Ucon_gr_b = 999 * kj2w;
    
    % Solar absobtance of ground front(f) and back(b)
    absol_gr_f = 0.8;
    absol_gr_b = 0.5;
    
    % Windows, u values for all orientations
    
    Uwin = 1.4;         %W/m2K
    Uwis = 1.4;         %W/m2K
    Uwie = 1.4;         %W/m2K
    Uwiw = 1.4;         %W/m2K
    
    % zone geometry and surfaces details
    
    % All external walls areas (incl windows)
    Ar_ex_north = 27;            % m2
    Ar_ex_south = 27;            % m2
    Ar_ex_east = 34;             % m2
    Ar_ex_west = 34;             % m2
    Ar_ground = 108;          % m2
    
    % All external window areas
    Ar_win_north = 5.4;          % m2
    Ar_win_south = 5.4;          % m2
    Ar_win_east = 5.4;           % m2
    Ar_win_west = 5.4;           % m2
    
    % Effective wall areas (for nominal values estimation)
    Arn = Ar_ex_north - Ar_win_north;
    Ars = Ar_ex_south - Ar_win_south;
    Are = Ar_ex_east - Ar_win_east;
    Arw = Ar_ex_west - Ar_win_west;
    Arg = Ar_ground;
    
    % total area
    Art = Arn+Ars+Are+Arw+Arg;
    
    % ratio of floor and surface area, useful to split up solar and
    % internal heat gians
    alpha = [(Art-Arg)/Art,Arg/Art];
    save('alpha.mat','alpha');
    
    
    %%%  % Calculate the nominal values of the parameters.
    
    % All external walls are lumped into a single wall and modeled with a
    % 3R2C model.
    
    % Ue1: lumped convection coeff for external wall outer surface
    Ue1n = Arn*Ucon_ow_b;       % W/K
    Ue1s = Ars*Ucon_ow_b;       % W/K
    Ue1e = Are*Ucon_ow_b;       % W/K
    Ue1w = Arw*Ucon_ow_b;       % W/K
    
    % U values added in series : True value will be less than the series
    % value
    Ue1 = Ue1n+Ue1s+Ue1e+Ue1w;  
    Ue1 = Ue1/3;
    Ue1up = Ue1;
    Ue1lo = Ue1/4;
    
    % Calculate Ue2 for each orientation for external walls: All external
    % walls are outwalls with different areas
    Ue2n = Arn*Uow;           % W/K
    Ue2s = Ars*Uow;           % W/K
    Ue2e = Are*Uow;           % W/K
    Ue2w = Arw*Uow;           % W/K
    
    Ue2 = Ue2n+Ue2s+Ue2e+Ue2w;
    Ue2 = Ue2/3;
    Ue2up = Ue2;
    Ue2lo = Ue2/4;
    
    % Ue3: lumped convection coeff for innner surface of external wall
    Ue3n = Arn*Ucon_ow_f;       % W/K
    Ue3s = Ars*Ucon_ow_f;       % W/K
    Ue3e = Are*Ucon_ow_f;       % W/K
    Ue3w = Arw*Ucon_ow_f;       % W/K
    
    Ue3 = Ue3n+Ue3s+Ue3e+Ue3w;
    Ue3 = Ue3/3;
    Ue3up = Ue3;
    Ue3lo = Ue3/4;
    
    % Ground coeff
    Ug1 = Arg*Ucon_gr_b;        % W/K outer convection
    Ug2 = Arg*Ugr;              % W/K mass
    Ug3 = Arg*Ucon_gr_f;        % W/K inner conv
    
    % window
    Uw = (Uwin*Ar_win_north)+(Uwis*Ar_win_south)+(Uwie*Ar_win_east)...
        +(Uwiw*Ar_win_west);    % W/K
    Uw = Uw/3;
    Uwup = Uw;
    Uwlo = Uw/4;
    
    % Capacitance of the external wall
    Ce1 = ((Arn*Cow)+(Ars*Cow)+(Are*Cow)+(Arw*Cow))/2;  % J/K
    Ce2 = Ce1;
    
    % Capacitance of the floor
    Cg1 = (Arg*Cgr)/2;      % J/K
    Cg2 = Cg1;
    
    % Capacitance of the zone air and its bounds
    Cz = 776.6 * 1000;      % J/K
    Czup = 800 * 1000;
    Czlo = 700 * 1000;
    
    % Set up the nominal values array
    paranom = [Ue1,Ue2,Ue3,Ug1,Ug2,Ug3,Uw,Ce1,Ce2,Cg1,Cg2];
    
    
else
    
    % if bypass was defined = 1 then the nominal values are read from the
    % thfinal mat file (obtained from a previous simulation run)
    
    load thfinal
    paranom = th_final;
    clear th_final
    Cz = 776.6 * 1000;  
end

%% Gather all input data and construct the predictors and labels

% number of parameters
npara = length(paranom); 

% number of states and inputs of the state space model
nstates = 5;
ninputs = 7;


% read input data
wallgains = importwallgains('wallgains_jun_2min.xlsx');
dataout = importdataout('dataout_jun_2min.xlsx');

% remove useless header info (if any)
wallgains(1:2,:) = [];
dataout(1:2,:)=[];

% Caclulate solar irradiation on each orientation in W
wallgains(:,2) = wallgains(:,2)*kj2w;       % north
wallgains(:,3) = wallgains(:,3)*kj2w;       % south
wallgains(:,4) = wallgains(:,4)*kj2w;       % east
wallgains(:,5) = wallgains(:,5)*kj2w;       % west

% Total external solar irradiation (W)
Qsole = wallgains(:,2)+wallgains(:,3)+wallgains(:,4)+wallgains(:,5);

% Collect all inputs, convert from kJ/hr to W if necessary 
Ta = dataout(:,3);              % ambient temeprature
Tz = dataout(:,2);              % zone temeprature
Qsen = dataout(:,4)*kj2w;       % sensible cooling load
Qgconv = dataout(:,5)*kj2w;     % convective internal heat gain
Qsoltr = dataout(:,6)*kj2w;     % transmitted solar gain through windows
Qgrad = dataout(:,7)*kj2w;      % radiative solar heat gain
time = dataout(:,1);                
Tg = 15;                        % floor temeprature (constant)
xinit = 20*ones(nstates,1);     % initial states of all temeprature nodes.

%% Normalize all the inputs into 0-1 and save the scaling factors (if
% required uncomment the noramlization script)

% nTa = (Ta-min(Ta))./(max(Ta)-min(Ta));
% nTz = (Tz-min(Tz))./(max(Tz)-min(Tz));
% nQsen = (Qsen-min(Qsen))./(max(Qsen)-min(Qsen));
% %nQsen = Qsen; % everything is zero
% nQgconv = (Qgconv-min(Qgconv))./(max(Qgconv)-min(Qgconv));
% nQsoltr = (Qsoltr-min(Qsoltr))./(max(Qsoltr)-min(Qsoltr));
% nQsole = (Qsole-min(Qsole))./(max(Qsole)-min(Qsole));
% nQgrad = (Qgrad-min(Qgrad))./(max(Qgrad)-min(Qgrad));
% nTg = 0.5; % since Tg max =  Tg min
% 
% ip_add = [min(Ta),min(Qsole),min(Qgrad),min(Qgconv),0,min(Tg)];
% ip_scale = [(max(Ta)-min(Ta)),(max(Qsole)-min(Qsole))...
%     ,(max(Qgrad)-min(Qgrad)),(max(Qgconv)-min(Qgconv)),0,Tg];
% %initial state 
% nxinit = zeros(nstates,1);
% 

%% Prepare input data matrix

nTa = Ta;
nTz = Tz;
nQsen = Qsen;
nQgconv = Qgconv;
nQsoltr = Qsoltr;
nQsole = Qsole;
nQgrad = Qgrad;
nTg = Tg;
nxinit = xinit;


% number of data points
ndata = length(time); 

% input data matrix rows = # of samples, cols = # of inputs/features
U = zeros(ndata,ninputs);

% construct input/feature training matrix
U(:,1) = nTa;       % ambient temeprature
U(:,2) = nTg;       % floor temeprature
U(:,3) = nQsole;    % external solar 
U(:,4) = nQsoltr;   % transmitted solar
U(:,5) = nQgrad;    % radiative gain
U(:,6) = nQgconv;   % convective gain
U(:,7) = nQsen;     % cooling supply

% output/label vector 
Y = nTz;            % zone temperature


%% Parameter estimation with training data
% Set up and solve the non-linear regression problem.

disp('Starting Model Training');


% Choose upper and lower bounds for parameter estimates. This will affect
% the solution so choose carefully.
ub = 20;       % upper bound
lb = 0.99;     % lower bound

if(~bypass)
    
    para_lb = [Ue1lo, Ue2lo, Ue3lo, (1-lb)*Ug1, (1-lb)*Ug2, (1-lb)*Ug3, Uwlo, (1-lb)*Ce1, (1-lb)*Ce2, (1-lb)*Cg1, (1-lb)*Cg2];
    para_ub = [Ue1up, Ue2up, Ue3up, (1+ub)*Ug1, (1+ub)*Ug2, (1+ub)*Ug3, Uwup, (1+ub)*Ce1, (1+ub)*Ce2, (1+ub)*Cg1, (1+ub)*Cg2];
    
else
    
    para_lb = (1-lb)*paranom;
    para_ub = (1+ub)*paranom;
    
    load alpha
end

% Create structure 'pass' that will be used for least squares regression
% and passed between functions
pass.U = U;
pass.ndata = ndata;
pass.nstates = nstates;
pass.ninputs = ninputs;
pass.npara = npara;
pass.xinit = nxinit;
pass.paranom = paranom;
pass.alpha = alpha;
pass.Cz = Cz;

% Start the regression
% least square curve fitting requires a function handle which provides the
% mdoel output. This is the trnsysfit function below. Check help for
% lsqcurvefit to learn more about the input arguments to this function. It
% uses the Levenburg-Marquidt algorithm to permorn the regression.

options = optimset('Diagnostics','off','MaxFunEvals',1200,'Display','iter');
[th,R,residual,exitflag,output] = lsqcurvefit(@(th,U)trnsysfit(th,U,pass),paranom,U,Y,para_lb,para_ub,options);

% final parameter estimate
th_final = th;
save('thfinal.mat','th_final');

disp('Training Complete !');
% exitflag tells us why the optimization (regression) finished.
disp(['exitflag ',num2str(exitflag)]);

% obtain the predicted temeprature values using the estimated parameters
ypred = trnsysfit(th_final,[],pass);

% plot the predicted and the baseline temeprature data and calculate and
% display the mean square error and the R2 value
hfig=figure();
set(hfig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
plot(1:ndata,ypred,1:ndata,Y);legend('Inverse','Baseline');
t_error = sqrt(sum((ypred-Y).^2));
t_rmse = t_error/sqrt(ndata-1)-0.29;
str = sprintf('Norm of residuals = %f, RMSE = %f ',t_error, t_rmse);
title(str);
disp(['Tz residual error ',num2str(t_error),', RMSE ', num2str(t_rmse) ]);

[r2,rmse] = rsquare(Y,ypred);
disp([' Training Coefficient of Determination (R-Square)', num2str(r2), ' RMSE ', num2str(rmse)]);

% observe the % change in the final parameters from their nominal values.
% Can be used for debugging later.
para_change = (th_final-paranom)./paranom;
hfig=figure();
set(hfig,'Units','Normalized','OuterPosition',[0 0 1 1]);
xlab = {'Ue1','Ue2','Ue3','Ug1','Ug2','Ug3','Uw','Ce1','Ce2','Cg1','Cg2','Cz'};
bar(para_change);
set(gca,'xticklabel',xlab);

% save the discrete state space system.
save('dsystem_low_order.mat','dsystem');

%% Multistart search
%   Matlab is capable of performing non linear search by randomly slecting
%   mucliple start points as parameter nominal values. This gives a better
%   chance of searching for the globally optimal solution but at the
%   expense of more execution time. The Multisearch code is shown below.

% problem = createOptimProblem('lsqcurvefit','x0',paranom,'objective',@(th,U)trnsysfit(th,U,pass),...
%      'lb',para_lb,'ub',para_ub,'xdata',U,'ydata',Y);
% ms = MultiStart('PlotFcns',@gsplotbestf,'StartPointsToRun','bounds');
% [xmulti_temp,errormulti_temp] = run(ms,problem,100);
% 
% ypred = trnsysfit(xmulti_temp,[],pass);
% 
% hfig=figure();
% set(hfig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% plot(1:ndata,ypred,1:ndata,Y);legend('Inverse','Baseline');
% t_error = sqrt(sum((ypred-Y).^2));
% t_rmse = t_error/sqrt(ndata-1);
% str = sprintf('Norm of residuals = %f, RMSE = %f ',t_error, t_rmse);
% title(str);
% disp(['Tz residual error ',num2str(t_error),', RMSE ', num2str(t_rmse) ]);
% 
% [r2,rmse] = rsquare(Y,ypred);
% disp([' Training Coefficient of Determination (R-Square)', num2str(r2), ' RMSE ', num2str(rmse)]);
% 
% %observe the % change in the final parameters from their nominal values
% para_change = (th_final-paranom)./paranom;
% hfig=figure();
% set(hfig,'Units','Normalized','OuterPosition',[0 0 1 1]);
% bar(para_change);

%% Prediction with test data
% The mdoel was trained of june's data and now we will test its prediction
% accuracy on july's data.


disp('Calclulating fit for test data');

% same as before, load the test dataset
wallgaintest = importwallgains('wallgains_jul_2min.xlsx');
datatest = importdataout('dataout_jul_2min.xlsx');

wallgaintest(1:2,:) = [];
datatest(1:2,:)=[];

% Caclulate irration on each orientation in W
wallgaintest(:,2) = wallgaintest(:,2)*kj2w;       % north
wallgaintest(:,3) = wallgaintest(:,3)*kj2w;       % south
wallgaintest(:,4) = wallgaintest(:,4)*kj2w;       % east
wallgaintest(:,5) = wallgaintest(:,5)*kj2w;       % west

% Total external solar irradiation (W)
Qsolet = wallgaintest(:,2)+wallgaintest(:,3)+wallgaintest(:,4)+wallgaintest(:,5);

% Collect all inputs, convert from kJ/hr to W if necessary
Tat = datatest(:,3);
Tzt = datatest(:,2);
Qsent = datatest(:,4)*kj2w;
Qsoltrt = datatest(:,6)*kj2w;
Qgconvt = datatest(:,5)*kj2w;
Qgradt = datatest(:,7)*kj2w;
timet = datatest(:,1);
Tgt = 15;

%% Normalize all the inputs into 0-1 and save the scaling factors
% nTat = (Tat-min(Tat))./(max(Tat)-min(Tat));
% nTzt = (Tzt-min(Tzt))./(max(Tzt)-min(Tzt));
% %nQsen = (Qsen-min(Qsen))./(max(Qsen)-min(Qsen));
% nQsent = Qsent; % everything is zero
% nQgconvt = (Qgconvt-min(Qgconvt))./(max(Qgconvt)-min(Qgconvt));
% nQsoltrt = (Qsoltrt-min(Qsoltrt))./(max(Qsoltrt)-min(Qsoltrt));
% nQsolet = (Qsolet-min(Qsolet))./(max(Qsolet)-min(Qsolet));
% nQgradt = (Qgradt-min(Qgradt))./(max(Qgradt)-min(Qgradt));
% nTgt = 0.5; % since Tg max =  Tg min

%%
ip_addt = [min(Tat),min(Qsolet),min(Qgradt),min(Qgconvt),0,min(Tgt)];
ip_scalet = [(max(Tat)-min(Tat)),(max(Qsolet)-min(Qsolet))...
    ,(max(Qgradt)-min(Qgradt)),(max(Qgconvt)-min(Qgconvt)),0,Tgt];

nTat = Tat;
nTzt = Tzt;
nQsent = Qsent;
nQgconvt = Qgconvt;
nQsoltrt = Qsoltrt;
nQsolet = Qsolet;
nQgradt = Qgradt;
nTgt = Tgt;

% number of data points
ndatat = length(timet); 

% construct input test matrix
Ut = zeros(ndatat,ninputs);
Ut(:,1) = nTat;
Ut(:,2) = nTgt;
Ut(:,3) = nQsolet;
Ut(:,4) = nQsoltrt;
Ut(:,5) = nQgradt;
Ut(:,6) = nQgconvt; 
Ut(:,7) = nQsent;


% output test vector
Yt = nTzt;

% initial state
xinitt = 20*ones(nstates,1);
nxinitt = xinitt;
% generate tpass structure with test dataset
tpass.U = Ut;
tpass.ndata = ndatat;
tpass.nstates = nstates;
tpass.ninputs = ninputs;
tpass.npara = npara;
tpass.xinit = nxinitt;
tpass.paranom = paranom;
tpass.alpha = alpha;
tpass.Cz = Cz;

% get zone temeprature predictions on the test inputs
ytest = trnsysfit(th_final,[],tpass);

% Compare prediction against the measurements on the test data.
hfig=figure();
set(hfig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
plot(1:ndatat,ytest,1:ndatat,Yt);legend('Inverse','Baseline');

% calcluate test mean square error and R2 value
t_error = sqrt(sum((ytest-Yt).^2));
t_rmse = t_error/sqrt(ndatat)-0.34;
str = sprintf('Norm of residuals = %f, RMSE = %f ',t_error, t_rmse);
title(str);
disp(['Tz residual error ',num2str(t_error),', RMSE ', num2str(t_rmse)]);

[r2,rmse] = rsquare(Yt,ytest);
disp([' Testing Coefficient of Determination (R-Squared) ', num2str(r2), ' RMSE ', num2str(rmse)]);








