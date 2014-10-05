function y = trnsysfit(th,dummy,pass)
% This is the function to run the simulation with relative parameter values
% th and input data U.

% we will save some structures globally.
global dsystem;
global sysd;

% timestep must match or be comparable to the sampling time for best
% results
timestep = 2*60;
% parse the 'pass' structure
U = pass.U;
nstates = pass.nstates;
ninputs = pass.ninputs;
ndata = pass.ndata;
alpha = pass.alpha;
Cz = pass.Cz;
xinit=pass.xinit;


% Model extraction.
% The structure of the model is defined in the trn_ssmodel file
[A,B,C,D] = trn_ssmodel(th,nstates,ninputs,alpha,Cz);

% form the state-space model and discritize it.
sys = ss(A,B,C,D);

% Model discretization.
sysd = c2d(sys,timestep);
Ad = sysd.A;    Bd = sysd.B;    
Cd = sysd.C;    Dd = sysd.D;

dsystem.Ad = Ad;
dsystem.Bd = Bd;
dsystem.Cd = Cd;
dsystem.Dd = Dd;

xinter=xinit;

% Model simualtion
YY=zeros(ndata,nstates);
for i = 1:ndata
    xinter = Ad*xinter + Bd*U(i,:)';
    YY(i,:) = Cd*xinter + Dd*U(i,:)';
end

y = YY(:,end);

end

