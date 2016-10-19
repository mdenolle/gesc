% EXAMPLE 1:
addpath '../../sub'
clear all;
close all;
clc
mkdir('mat')


%% TIME /FREQUENCY
typ = 'f' ;          % 't': enters time series length (win= in s) and rate (rate=1/s)
                     % 'f': enters frequency series length (win=in Hz) and rate
                     % (1/Hz)
win = 5;           % (in s or Hz) window length
rate = 0.1 ;         % (in 1/s or Hz) sampling rate
FT = make_ft(win,rate,typ);
FT.Fmax=5;         % calculate until Fmax
max_mode=1;         % maximum number of surface-wave modes to ouput
freq=(0.1:0.5:1);     % plots eigenfunctions for periods 1 -> 10s



%% MEDIUM
H     = [0 5] ;        % (km) vector of interface depths 
beta  = [2];       % (km/s) S-wavespeed 
alpha = sqrt(3)*beta; % (km/s) P-wavespeed 
rho   = [1.6 ];     % (kg/dm^3) density
typ = {'cst'};          % constant within layers: 'cst', all gradients 'grad'

Dmax=10*max(beta)/FT.df;      % (km) "half space" effective thickness, matters for long period
MED = make_layers(H,Dmax,alpha,beta,rho,typ);
!mv *.ps plots


%% RESOLUTION:
% NR : resolution structure for adaptive sampling (see Section 2.3)
% typ = 'lin';       % bounded linear(ramp) adaptive sampling 
%     lambda=[2 10]; % H/lambda domains of sampling laws
%     a1=[0.5 1.5];  % ramp lower corner
%     a2=[2 8];      % ramp upper corner
% 	N1=[20 150];    % lower bound # points
% 	N2 = [40 100];  % upper bound # points
%     NR = make_npt_colloc_lin(lambda,N1,N2,a1,a2);


NR.typ='rec';       % choose recursive
NR.N1 = 15 ;        % lower bound for number of points / layer
NR.N2 = 50 ;        % maximum number of points / layers


%% GESC
bigC = gesc(MED,FT,NR,max_mode,freq);
save(['./mat/MED_' int2str(MED(1).Nl-1) 'l_' char(MED(1).typ)],'bigC','NR','MED','FT','-v7.3')
!mv *.ps plots
