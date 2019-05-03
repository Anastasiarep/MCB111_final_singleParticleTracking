function [singleExpParameters, singleExpResnorm, doubleExpParameters, ...
    doubleExpResnorm] = spt_cdfDoubleExp(tLag, exposure, path, file)

% Adapted from:
% Jan 2018 Matt Holmes
% Square displacement calculation code adapted from 
% Georgia Squyres, Bisson-Filho et al 2017.
% Opens a user selected .mat track file. 
% Tracks should be stored in .mat file saved as variable "tracks"
% User should set exposure of original imaging file 
% and tLag to be used for analysis

addpath(['C:\Users\anr612\Documents\MATLAB\altmany' ...
    '-export_fig-d570645\altmany-export_fig-d570645']);
%where export_fig function lives

%% User Provided parameters
%tlag is time in seconds used to calculate 
%squared displacement distribution
% tLag = 0.02; %10ms exp + 10ms interval = 20ms.

%exposure is the exposure per frame used (in seconds)
% exposure = 0.01; %10ms exposure time

%% calculate frameLag using above parameters
if mod(tLag , exposure) ~= 0
    error('tLag is not divisible by the exposure time')
end
frameLag = tLag/exposure; %a frame every 2 exposure times

%% User selects file
% path=['C:\Users\anr612\Documents\1. Harvard18-\2. Rotations' ...
%     '\2. DenicLab\MicroscopeImages\TrackingImages\']; %where .mat files live
% condition='mock';
% n='1';
% file=['C2-23-1-19_LH_',condition,'_', n '_Tracks.mat']; %file name
% [file,path] = uigetfile('*.mat');

%% Defining exponential functions
%for particle that starts at orign at t=0
%x_0=x(1)=(r_0)^2=4Dt, where r=circle of radius r=square displacements
%sqdisp is the sqaureDisplacement values input in the equation
singleExp = @(x, sqdisp) (1 - exp(-(sqdisp/x(1))));
%sqdisp=r^2
%x(1)=(r_0)^2

%x(1) is (r_1)^2, which each equals 4(D_1)t. 
%x(2) is (r_2)^2, which each equals 4(D_2)t. 
%x(3) is alpha, the fraction of data described by x(1)
doubleExp = @(x, sqdisp) ...
    (1 - x(3)*exp(-(sqdisp/x(1))) - (1-x(3))*exp(-(sqdisp/x(2))));
%2-component system, where fast & slow mobility=
%diffusion constants D1 & D2=fractional contributions alpha & (1-alpha)

%% Load Tracks
clear tracks
load(strcat(path,'\',file))
% load([path,'/',file]);
try
    tracks;
catch
    error('Variable tracks not found in file');
end

allSqDisp = [];
    
%% Calculate all square displacements.
for i = 1:length(tracks)
    currTrack = tracks{i};
    trackLength = size(currTrack,1);
    if trackLength <= frameLag
        continue
    end

    startpos = currTrack(1:trackLength - frameLag, 2:3);
    endpos = currTrack(frameLag+1:end, 2:3);
    allSqDisp = [allSqDisp; sum((endpos-startpos).^2,2) ];
end
       
%% Make cumulative distribution function of square displacements
allSqDisp = allSqDisp';
allSqDisp = sort(allSqDisp);
%The probability is 0 at r^2 = 0;
allSqDisp = [0, allSqDisp];
%generate probability values
cdf = linspace(0, 1, length(allSqDisp)); %n values between 0 and 1
    
%% Fit to a single exponential
[singleExpParameters, singleExpResnorm, singleExpResiduals] =  ...
    lsqcurvefit(singleExp, 0.1, allSqDisp, cdf, 0, Inf);
%lscurvfit(function, x0, xdata, ydata, lower & upper bound)
%Output: solution, resnorm (squared norm of residual), residual (value of
%objective function at solution)
%Solution has size of x0

% starts at x0 and finds coefficients x to best fit the nonlinear function 
% fun(x, xdata) to the data ydata (in the least-squares sense). 
% ydata must be the same size as the vector (or matrix) F returned by fun.
    
%% Fit to a double exponential
[doubleExpParameters, doubleExpResnorm, doubleExpResiduals] = ...
    lsqcurvefit(doubleExp, [0.1 0.1 0.1], allSqDisp, cdf, ...
    [0 0 0], [Inf Inf 0.5]);

%% Plot single and double exponential curves vs. data and residuals    
% fitPlots = figure();
% p1 = plot(allSqDisp, cdf,'.');
% hold on
% singExpY = singleExp(singleExpParameters, allSqDisp);
% p2 = plot(allSqDisp, singExpY);
% doubleExpY = doubleExp(doubleExpParameters, allSqDisp);
% p3 = plot(allSqDisp, doubleExpY);
% 
% title({'Cumulative Distribution of Square Displacements',...
%     ['N = ', num2str(length(allSqDisp) - 1),'          '...
%     '\Deltat =',  num2str(tLag), 's']});
% xlabel('r^2 (µm^2)');
% ylabel('P(r^2)');
% legend([p2 p3],['Single Exponential R_0 = ',num2str(singleExpParameters(1))],...
%     ['Double Exponential R_1 = ',num2str(doubleExpParameters(1)),...
%     'R_2 = ',num2str(doubleExpParameters(2))]);
% 
% hold off 
% 
% % x_saveas=['C2-23-1-19_LH_',condition,'_', n '_fit.png'];
% % cd 'C:\Users\anr612\Documents\1. Harvard18-\2. Rotations\2. DenicLab\MicroscopeImages\TrackingImages'
% % set(gcf, 'Color', 'w');
% % saveas(gcf, x_saveas);
% 
% residualPlot = figure();
% plot(allSqDisp, singleExpResiduals, '.');   
% hold on
% plot(allSqDisp, doubleExpResiduals, '.');
% 
% title({'CDF Fit Residuals',...
%     ['N = ', num2str(length(allSqDisp) - 1),'          '...
%     '\Deltat =',  num2str(length(tLag)), 's']});
% xlabel('r^2 (µm^2)');
% ylabel('Residuals');
% legend('Single Exponential','Double Exponential');
% 
% % x_saveas=['C2-23-1-19_LH_',condition,'_', n '_res.png'];
% % cd 'C:\Users\anr612\Documents\1. Harvard18-\2. Rotations\2. DenicLab\MicroscopeImages\TrackingImages'
% % set(gcf, 'Color', 'w');
% % saveas(gcf, x_saveas);

%%
% singleExpParameters %x_0=x(1)=(r_0)^2=4Dt.
% singleExpResnorm
% doubleExpParameters
% %x(1)=(r_1)^2=4(D_1)t. 
% %x(2)=(r_2)^2=4(D_2)t. 
% %x(3)=alpha.
% doubleExpResnorm
% %resnorm (squared norm of residual)

cd C:\Users\anr612\Documents\MATLAB\SingleParticleTracking

end
