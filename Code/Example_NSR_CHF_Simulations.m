%% Example simulations using a 5 minute RR sequence for NSR and CHF 
clear all; close all; clc

% Simulation Parameters 

% Load CHF/NSR parameters
 load('Final_NSR_CHF_params.mat','nu0','fleak0','beta0','Nscale_CHF',...
     'Fscale_CHF','Bscale_CHF','delta_CHF');
Tscale = 2;
delta_NSR = 0;
Nscale_NSR = 1;
Fscale_NSR = 1;
Bscale_NSR = 1;
hflag_NSR = 1; % NSR data 
hflag_CHF = 2; % CHF data   


% Time Window
twindowm = 5; % minutes
twindowsec = twindowm*60; % seconds


% Simulation type: NSR or CHF
sim_type = 'CHF'

switch sim_type 
    case 'NSR'
        % NSR Simulation
        % Load NSR RR file  
        NSR_file = 'nsr1.mat';
        load(NSR_file);
        
        % Select segment of RR sequence 
        RRind = ~isnan(RR);
        RRclean = RR(RRind);
        tempsum = cumsum(RRclean);
        beg_ind = 1;
        end_ind = find(tempsum<=twindowsec,1,'last');

        % Run simulation
        Tin = RRclean(beg_ind:end_ind)*1e3; % s to ms 
        beats_std_NSR = ceil(length(Tin)*.1)+1;

        [a,r,cp,u,l,b,c,s] = qu2007_ca_apd_map_human_leak_inputT(beats_std_NSR, Tin, nu0,fleak0,hflag_NSR,...
            Tscale,beta0,Nscale_NSR,Fscale_NSR,Bscale_NSR,delta_NSR);


    case 'CHF'
        % CHF simulation 
        % Load CHF RR file  
        CHF_file = 'chf201.mat';
        load(CHF_file)
        
        % Select segment of RR sequence 
        RRind = ~isnan(RR);
        RRclean = RR(RRind);
        tempsum = cumsum(RRclean);
        beg_ind = 1;
        end_ind = find(tempsum<=twindowsec,1,'last');

        % Run simulation 
        Tin_CHF = RRclean(beg_ind:end_ind)*1e3; % s to ms 
        beats_std_CHF = ceil(length(Tin_CHF)*.1)+1;

        [aCHF,rCHF,cpCHF,uCHF,lCHF,bCHF,cCHF,sCHF] = qu2007_ca_apd_map_human_leak_inputT(beats_std_CHF, Tin_CHF,...
            nu0,fleak0,hflag_CHF,Tscale,beta0,Nscale_CHF,Fscale_CHF,Bscale_CHF,delta_CHF);
end 


