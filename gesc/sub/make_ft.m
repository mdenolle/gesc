%% Create time/frequency structure FT
function FT = make_ft(win,rate,typ)
%% input
% choose between time (typ='t') or frequency ('f'):
   % enter win: window length (in s or in Hz)
   % enter rate : sampling rate (in 1/s or 1/Hz)
% choose either the frequencies to solve the eigen problem at,
% AND/OR the time domain to build seismograms
% Marine Denolle (04/2014)

FT=[];
if (strcmp(typ,'t')==1)% if TIME DOMAIN
    FT.twin = win ; % (in s) time series length
    FT.dt = rate; % (in 1/s) sampling rate
    FT.Nsps = 1/rate; % number of sample per seconds 
    FT.nwin = floor(FT.twin/FT.dt); % # points in time series
    FT.Fnyq = 1/(2*FT.dt); % Nyquist frequency
    FT.omega = 2*pi*linspace(0,FT.Fnyq,FT.nwin/2+1); % angular frequency vector
    FT.df = FT.omega(2)/(2*pi);
    FT.tt=(0:FT.dt:FT.twin-FT.dt);
elseif (strcmp(typ,'f')==1)
    FT.Fnyq = win; % (in Hz) Fmax
    FT.df = rate; % (in 1/Hz) sampling rate
    FT.dt = 1/(2*FT.Fnyq); % (in 1/s) sampling rate of hypothetical time series
    FT.nwin = 2*floor(FT.Fnyq/FT.df)-1; % # of hypothetical time series
    FT.omega=2*pi*linspace(0,FT.Fnyq,floor(FT.Fnyq/FT.df)); % angular frequency vector
    FT.Nsps = 1/FT.dt;
end
FT.Fmax=FT.Fnyq; % by default, the max frequency to look at is the Nyquist frequency
    
 return