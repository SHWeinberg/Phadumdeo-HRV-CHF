function [a,r,cp,u,l,b,c,s]= qu2007_ca_apd_map_human_leak_inputT(beats_std, T, nu0,fleak0,hflag,Tscale,beta0,Nscale,Fscale,Bscale,delta)

% delta = scaling factor for SR leak

% hflag = 3; CHF restitution curve
% hflag = 2; CHF Restitution Curve, and release/Serca   
% hflag = 1; NSR Restitution Curve 
% hflag = 0; Animal Restitution Curve 


if (hflag == 2) % CHF human, Restitution and Ca
    nu = nu0*Nscale; % Gomez 2014
    fleak = fleak0*Fscale; % Gomez 2014
    beta = beta0*Bscale;
elseif (hflag == 1) % NSR  human
    nu = nu0;
    fleak = fleak0;
    beta = beta0;
elseif (hflag == 3) % CHF human, Restitution alone
    nu = nu0;
    fleak = fleak0;
    beta = beta0;
else % animal 
    nu = nu0;
    fleak = fleak0;
    beta = 5;
end 
         
n = length(T); % Number of beats

gamma = 1e-3;      % Ca->APD coupling, >0 positive, <0 negative coupling
kappa = .1;      % rate of Ca equilibration 
eta = .1;          % APD->Ca coupling
      
r = nan(1,n);       % Ca released from SR
cp = nan(1,n);      % peak cyto Ca
a = nan(1,n);       % APD 
u = nan(1,n);       % SR Ca uptake
l = nan(1,n);       % SR Ca load
b = nan(1,n);       % total Ca in cell, b_n = c_n + l_n
c = nan(1,n);       % diastolic cyto Ca
s = nan(1,n);       % leak from SR 

% initial conditions
r(1) = 100; % uM
cp(1) = 60;  % uM
a(1) = 200;  % ms
u(1) = ufun(T(1),Tscale)*hfun(cp(1),nu);   % uM
l(1) = 100;  % uM
b(1) = 200;  % uM
c(1) = b(1) - l(1);     % uM
s(1) = 0; % Ask seth

for i = 1:n-1
   r(i+1) = fleak*qfun(T(i)-a(i),Tscale)*gfun(l(i),beta);
   cp(i+1) = c(i) + r(i+1);
   a(i+1) = ffun(T(i)-a(i),hflag)/(1-pfun(cp(i+1),gamma));
   u(i+1) = ufun(T(i+1),Tscale)*hfun(cp(i+1),nu);
   s(i+1) = delta*qfun(T(i)-a(i),Tscale)*gfun(l(i),beta);
   l(i+1) = l(i) - r(i+1) + u(i+1) - s(i+1);
   b(i+1) = b(i) - kappa*(c(i)-cfun(T(i+1),Tscale)) + eta*(a(i+1)-a(i));
   c(i+1) = b(i+1) - l(i+1);
end


% Dr. Weinberg edits
a = a(beats_std:n);
r = r(beats_std:n);
cp = cp(beats_std:n);
u = u(beats_std:n);
l = l(beats_std:n);
b = b(beats_std:n);
c = c(beats_std:n);
s = s(beats_std:n);

% remove zero values
index = find(a~=0);
a = a(index);
r = r(index);
cp = cp(index);
u = u(index);
l = l(index);
b = b(index);
c = c(index);
s = s(index);


function f = ffun(d,hflag)
A1 = 0;
D1 = 0;
tau1 = 10;
% dmin = 2; 

if (hflag == 2) % CHF human, Restitution and Ca
    A0 = 372.3244; D0 = 535.1088; tau0 = 803.6337;
    dmin = 10;
    f = (A0*(1-exp(-(d+D0)./tau0))).*(d>=dmin);
elseif (hflag == 1) % NSR  human
    A0 = 275.9279 ; D0 = 437.7602; tau0 =  474.0197;
    dmin = 10;
    f = (A0*(1-exp(-(d+D0)./tau0))).*(d>=dmin);
elseif (hflag == 3) % CHF human, Restitution
    A0 = 372.3244; D0 = 535.1088; tau0 = 803.6337;
    dmin = 10;
    f = (A0*(1-exp(-(d+D0)./tau0))).*(d>=dmin);
else % animal 
    A0 = 220; D0 = 40; tau0 = 60;
    dmin = 2;
    f = (A0*(1-1./(1+exp((d-D0)/tau0))) + A1*exp(-(d-D1).^2/tau1)).*(d>=dmin);
end 

  


function g = gfun(l,beta)
alpha = 0.036;
lc = 93.5;
% beta = 5;
g = l.*(1-(1-alpha)./(1+exp((l-lc)/beta)));

function h = hfun(cp,v)
c0 = 28;
theta = 20;
h = v*cp*(1-1/(1+exp((cp-c0)/theta)));

function p = pfun(cp, gamma)
p = gamma*cp;

function c = cfun(T,Tscale)
c0 = 28;
ep = 2;
tauc = 300*Tscale;
c = c0*(1+ep*exp(-T/tauc));

function q = qfun(d,Tscale)
sigma = 0.5;
tauq = 80*Tscale;
q = 1-sigma*exp(-d/tauq);

function u = ufun(T,Tscale)
rho = 0.15;
tauu = 200*Tscale;
u = 1-rho*exp(-T/tauu);