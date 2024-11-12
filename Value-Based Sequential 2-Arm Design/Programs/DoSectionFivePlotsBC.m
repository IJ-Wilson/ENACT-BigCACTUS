function [figout, mat] = DoSectionFivePlotsBC(fignum,basic,advanced,PathReps,ProductionReps,graphicsuffix)
% DoSectionFourPlots: This function generates plots for section which
% demonstrates sample paths and basic sensitivity analysis for the delay
% sequential sampling paper.
%
% CALLED BY: DelayExperimentsForPaper 
%
% 21 Mar 2015: Forster, Chick Pertile

% This file has been edited by Laura Flight to calculate stopping boundary
% for the CACTUS case study 

dirname = 'Figure'; % directory for putting figures
linewid = 2;

[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end; 
advanced.dt

advanced.PLOTSIMS = false ;
[~, mat] = DelayCurvesRecur(basic, advanced); % first compute the optimal stopping region for stage II (and III), from times [tau
... Tmax ] (the free boundary PDE solution
[mat] = DelayStageOne(basic, advanced, mat ); % next compute the optimal stage I experiment, from times [0 ... tau], and blend in
%the stage II solution
toneshot = (0:basic.TMax);
[mat]  = DelayOptimalBayesOneStage (basic, advanced, toneshot, mat ); % this function that Steve wrote defines the optimal fixed sample size based on a comparison of EVSI and sampling costs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% For Figure 1: Generate a few sample paths and show the general
%%%%%% features of the stopping times etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the optimal stopping rule
fignum = fignum + 1; figure(fignum); % always increment fignum and create a new figure
hold off ; 

% plot Stage 1 boundary if there exists a Stage 1
plot (mat.tvec, mat.bndupper,'--k',mat.tvec ,mat.bndlower, '--k' ,'LineWidth',linewid );
hold on ; 
if abs( mat.Threshpoint(2) - mat.Threshpoint(4) ) > 1
    scatter (mat.bestsvec(mat.Threshpoint(2):mat.Threshpoint(4)) + basic.t0, mat.muvec(mat.Threshpoint(2):mat.Threshpoint(4)),'k');
    hold on ; 
end
% plot upper Stage 1 boundary if there exists a Stage 1
if abs( mat.Threshpoint(3) - mat.Threshpoint(1)) > 1
    scatter (mat.bestsvec(mat.Threshpoint(3):mat.Threshpoint(1)) + basic.t0, mat.muvec(mat.Threshpoint(3):mat.Threshpoint(1)),'k');
    hold on ; 
end 
ylabel('Prior / Posterior Mean for E[INMB]', 'Fontsize', advanced.bigfontsize);
xlabel('n_0 + pairwise allocations', 'Fontsize', advanced.bigfontsize);
ishoriz = false; % set to true if muvec is on horizontal axis, false if on vertical axis
plot_upper = mat.muvec(mat.Threshpoint(1)) + min(mat.muvec(mat.Threshpoint(1))/4,- mat.muvec(mat.Threshpoint(2))/4) ; % set upper limit of plot region ( + /  
plot_lower =  mat.muvec(mat.Threshpoint(2)) - min(-mat.muvec(mat.Threshpoint(2))/4,mat.muvec(mat.Threshpoint(1))/4)  ; % set lower limit of plot region

plot_lower_tmp=plot_lower-500;
plot_upper_tmp=plot_upper+500;
ylim( [ plot_lower_tmp, plot_upper_tmp] ) ;
xlim( [ 0, basic.t0+basic.TMax+150] ) ;
line( [ basic.t0 + basic.tau, basic.t0 + basic.tau ], [ plot_lower_tmp, plot_upper_tmp ], 'Linestyle', ':', 'Color', 'k' ) ;
line( [ basic.t0 + basic.TMax, basic.t0 + basic.TMax ], [ plot_lower_tmp, plot_upper_tmp ], 'Linestyle', ':', 'Color', 'k' ) ;
line( [ basic.t0, basic.t0 ], [ plot_lower_tmp, plot_upper_tmp ], 'Linestyle', ':', 'Color', 'k' ) ;

tmpfontsize = round(advanced.smallfontsize * 1.1);
text(  basic.t0 + 2.5,  plot_lower_tmp + 400, 'n_0', 'Fontsize', tmpfontsize) ;
text(  basic.t0 + basic.TMax+2.5,  plot_lower_tmp + 650  , 'n_0 + T_{max}', 'Fontsize', tmpfontsize) ;
text(  basic.t0 + basic.tau + 2.5,  plot_lower_tmp + 650 , 'n_0 + \tau', 'Fontsize', tmpfontsize) ;

if length(graphicsuffix) == 0
    text(  ( basic.t0 + basic.tau ) / 5.0 , (mat.muvec(mat.Threshpoint(1)) + plot_upper_tmp)/2 , 'No trial', 'Fontsize', tmpfontsize); 
    text(  ( basic.t0 + basic.tau ) / 5.0,  (mat.muvec(mat.Threshpoint(2)) + plot_lower_tmp)/2 , 'No trial', 'Fontsize', tmpfontsize) ;
if abs(mat.Threshpoint(1) - mat.Threshpoint(3)) > 1
        text(  ( basic.t0 + basic.tau/5.0 ) , (mat.muvec(mat.Threshpoint(1)) + mat.muvec(mat.Threshpoint(3)))/2 , 'Fixed trial', 'Fontsize', tmpfontsize) ;
    end
    if abs(mat.Threshpoint(4) - mat.Threshpoint(2)) > 1
        text(  ( basic.t0 + basic.tau/5.0 ) ,  (mat.muvec(mat.Threshpoint(2)) + mat.muvec(mat.Threshpoint(4)))/2 , 'Fixed trial', 'Fontsize', tmpfontsize) ;
    end
 
    text(  ( basic.t0 + basic.tau ) / 7.0 ,0 , 'Sequential trial', 'Fontsize', tmpfontsize  ) ;
    text(  ( basic.t0 + basic.tau ) / 7.0 , -1000 , 'recruitment', 'Fontsize', tmpfontsize  ) ;
  
end
ist0 = true ; % set to true if ABCD to be plotted at t0, false if to be plotted on horizontal axis
UtilPlotABCD( basic, advanced, mat, ishoriz, ist0 );  
UtilStdizeFigure(fignum,advanced,true);
dirname = 'Figure' ;
UtilSaveFigEpsPdf(fignum,dirname,strcat('pathsBC',graphicsuffix),'-r600');

end