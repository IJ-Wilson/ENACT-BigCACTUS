%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Big CACTUS case study - calculated of value-based sequential two arm trial with adaptive stopping boundaries 
%%%%%% Code modified from htadelay 
%%%%%% Laura Flight 
%%%%%% 31Aug22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all hidden; clear;
cd 'xxxxxx\Value-Based Sequential 2-Arm Design\Programs\htadelay-master'

% bring in files from current directory and from two other locations,
%   delaycore\ contains routines which are generally useful
%   delaypaper\ contains routines specific to this paper
LocalDelaySetPaths;  % Default case: Sets up PATH of Matlab to get code for running. This file may be localized.

%pause on      % set to 'on' if you want extra pauses at certain times during computation to view things, to 'off' otherwise
pause off

% CAN SET FOLLOWING VALUES DEPENDING ON ACCURACY IN PLOTS DESIRED.
doProductionRuns = true;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
doSaveMatFile = false;

%%%% the following can be customized to increase to to decrease accuracy of
%%%% both test run mode and production run mode. test run is fast but not
%%%% accurate, production mode is for generating higher resolution graphs.
PRODUCTIONREPS = 15000;
PRODUCTIONNUMBERSTD = 200;

TESTREPS = 200;             % Allow end user to configure number of simulation replications
TESTNUMPERSTD = 100;        % and fineness of grid for PDE computations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if ~exist('fignum','var'), fignum = 20; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause; close all hidden;
clear basic; clear advanced; clear mat;
tic
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetBC_basecase;

advanced.saveplot = false;    % set to true to save plots in files, false to not save files automatically
advanced.DOPLOT= false ;
advanced.keepAllOutput = true ;
if doProductionRuns
    ProductionReps = PRODUCTIONREPS;    %15000
    advanced.MinGridPerStdev = PRODUCTIONNUMBERSTD;    %200;
else
    ProductionReps = TESTREPS;
    advanced.MinGridPerStdev = TESTNUMPERSTD;    %20;
end

graphicsuffix = '';
[mat] = DoBoundaryPlotsBC(fignum, basic,advanced,graphicsuffix);
minutestocompute = toc/60

 writematrix(mat.tvec,'xxxxxx\Value-Based Sequential 2-Arm Design\Data\Derived\Tmax 95\tvec1.csv')
 writematrix(mat.bndupper,'xxxxxx\Value-Based Sequential 2-Arm Design\Data\Derived\Tmax 95\bndupper1.csv')
 writematrix(mat.bndlower,'xxxxxx\Value-Based Sequential 2-Arm Design\Data\Derived\Tmax 95\bndlower1.csv')
 writematrix(mat.Threshpoint,'xxxxxx\Value-Based Sequential 2-Arm Design\Data\Derived\Tmax 95\threspoint1.csv')
 writematrix(mat.muvec,'xxxxxx\Value-Based Sequential 2-Arm Design\Data\Derived\Tmax 95\muvec1.csv')

format long
% Expected net benefit (value of sequential trial minus value of adopting
% straight away
NB_seq=mat.B0vec(find(mat.muvec==basic.mu0))-max(0,basic.PPatients*basic.mu0)
