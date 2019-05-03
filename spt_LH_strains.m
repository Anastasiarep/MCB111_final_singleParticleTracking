cd C:\Users\anr612\Documents\MATLAB\SingleParticleTracking

%%
condition='mock';
nMock=string(1:10);
fileNamesMock=string(zeros(1, length(nMock)));
for i=1:length(nMock)
    fileNamesMock(i)=strcat('C2-23-1-19_LH_',...
        condition,'_', nMock(i), '_Tracks.mat'); %file name
end

condition='rap';
nRap=string(1:11);
fileNamesRap=string(zeros(1, length(nRap)));
for i=1:length(nRap)
    fileNamesRap(i)=strcat('C2-23-1-19_LH_',...
        condition,'_', nRap(i), '_Tracks.mat'); %file name
end

%%
tLag=0.02;
exposure=0.01;
path=['C:\Users\anr612\Documents\1. Harvard18-\2. Rotations' ...
    '\2. DenicLab\MicroscopeImages\23-Jan-19_TrackingImages_LH\']; 
%where .mat files live

N=length(nMock);
dSingle=zeros(1, N);
resSingle=zeros(1, N);
d1Double=zeros(1, N);
d2Double=zeros(1, N);
alphaDouble=zeros(1, N);
resDouble=zeros(1, N);

for j=1:N
    file=fileNamesMock(j);
    [singleExpParameters, singleExpResnorm, doubleExpParameters, ...
    doubleExpResnorm]=spt_cdfDoubleExp(tLag, exposure, path, file);
    dSingle(j)=singleExpParameters/(4*tLag);
    resSingle(j)=singleExpResnorm;
    d1Double(j)=doubleExpParameters(1)/(4*tLag);
    d2Double(j)=doubleExpParameters(2)/(4*tLag);
    alphaDouble(j)=doubleExpParameters(3);
    resDouble(j)=doubleExpResnorm;
end

parametersMock=horzcat(dSingle', resSingle',...
    d1Double', d2Double', alphaDouble', resDouble');

%%
tLag=0.02;
exposure=0.01;
path=['C:\Users\anr612\Documents\1. Harvard18-\2. Rotations' ...
    '\2. DenicLab\MicroscopeImages\23-Jan-19_TrackingImages_LH\']; 
%where .mat files live

N=length(nRap);
dSingle=zeros(1, N);
resSingle=zeros(1, N);
d1Double=zeros(1, N);
d2Double=zeros(1, N);
alphaDouble=zeros(1, N);
resDouble=zeros(1, N);

for j=1:N
    file=fileNamesRap(j);
    [singleExpParameters, singleExpResnorm, doubleExpParameters, ...
    doubleExpResnorm]=spt_cdfDoubleExp(tLag, exposure, path, file);
    dSingle(j)=singleExpParameters/(4*tLag);
    resSingle(j)=singleExpResnorm;
    d1Double(j)=doubleExpParameters(1)/(4*tLag);
    d2Double(j)=doubleExpParameters(2)/(4*tLag);
    alphaDouble(j)=doubleExpParameters(3);
    resDouble(j)=doubleExpResnorm;
end

parametersRap=horzcat(dSingle', resSingle',...
    d1Double', d2Double', alphaDouble', resDouble');

%%
meanparametersMock=mean(parametersMock);
meanparametersRap=mean(parametersRap);
realDiff=meanparametersMock-meanparametersRap;

STDparametersMock=std(parametersMock);
STDparametersRap=std(parametersRap);

%%
merged=vertcat(parametersMock, parametersRap);
m=1000;
bsDiff=zeros(m, 6);

%bootstrap significance
% for i=1:m
%     bsMock=datasample(merged, length(nMock));
%     bsRap=datasample(merged, length(nRap));
%     bsMockMean=mean(bsMock);
%     bsRapMean=mean(bsRap);
%     bsDiff(i, :)=bsMockMean-bsRapMean;
% end

% histogram(bsDiff(:,5))
% percentile2_5=prctile(bsDiff, 2.5); 
% percentile97_5=prctile(bsDiff, 97.5);

% [h,p]=ttest2(parametersMock,parametersRap);

dSingle=horzcat(vertcat(parametersMock(:,1), NaN), parametersRap(:,1));
resSingle=horzcat(vertcat(parametersMock(:,2), NaN), parametersRap(:,2));
d1Double=horzcat(vertcat(parametersMock(:,3), NaN), parametersRap(:,3));
d2Double=horzcat(vertcat(parametersMock(:,4), NaN), parametersRap(:,4));
alphaDouble=horzcat(vertcat(parametersMock(:,5), NaN), parametersRap(:,5));
resDouble=horzcat(vertcat(parametersMock(:,6), NaN), parametersRap(:,6));

%%
% bar(parametersMock)
% single exponential
subplot(1, 2, 1)
boxplot(dSingle, 'Labels',{'Control','Rapamycin'})
hold on
plot([1, 2], dSingle, 'o', 'Color', [204/255, 204/255, 0/255])

set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Control','Rapamycin'},...
    'FontName', 'Cambria Math', 'FontSize', 12);
ylabel('Diffusion coefficient', 'FontName', 'Cambria Math');

subplot(1, 2, 2)
boxplot(resSingle, 'Labels',{'Control','Rapamycin'})
hold on
plot([1, 2], resSingle, 'o', 'Color', [204/255, 204/255, 0/255])

set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Control','Rapamycin'},...
    'FontName', 'Cambria Math', 'FontSize', 12);
ylabel('Residuals', 'FontName', 'Cambria Math');

sgtitle('Single exponential', 'FontName', 'Cambria Math', 'FontSize', 14);
hold off

% x_saveas=['LH_singleExp.png'];
% cd 'C:\Users\anr612\Documents\1. Harvard18-\2. Rotations\2. DenicLab\MicroscopeImages\23-Jan-19_TrackingImages_LH'
% set(gcf, 'Color', 'w');
% saveas(gcf, x_saveas);

%%
% 
% boxplot([d1Double, d2Double, alphaDouble, resDouble])
% hold on
% plot([1:8], [d1Double, d2Double, alphaDouble, resDouble], 'o', 'Color', [204/255, 204/255, 0/255])

% double exponential
boxplot([d1Double, d2Double],...
    'Labels',{'Control','Rapamycin','Control','Rapamycin'})
hold on
plot([1:4], [d1Double, d2Double], 'o', 'Color', [204/255, 204/255, 0/255])

set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'Control','Rapamycin','Control','Rapamycin'},...
    'FontName', 'Cambria Math', 'FontSize', 12);
ylabel('Diffusion coefficients', 'FontName', 'Cambria Math');
text(1.4,-60,'D_{1}', 'FontName', 'Cambria Math', 'FontSize', 12)
text(3.4,-60,'D_{2}', 'FontName', 'Cambria Math', 'FontSize', 12)
sgtitle('Double exponential diffusion coefficients',...
    'FontName', 'Cambria Math', 'FontSize', 14);

hold off

% x_saveas=['LH_doubleExp_1.png'];
% cd 'C:\Users\anr612\Documents\1. Harvard18-\2. Rotations\2. DenicLab\MicroscopeImages\23-Jan-19_TrackingImages_LH'
% set(gcf, 'Color', 'w');
% saveas(gcf, x_saveas);

%%
subplot(1, 2, 1)
boxplot(alphaDouble, 'Labels',{'Control','Rapamycin'})
hold on
plot([1:2], alphaDouble, 'o', 'Color', [204/255, 204/255, 0/255])

set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Control','Rapamycin'},...
    'FontName', 'Cambria Math', 'FontSize', 12);
ylabel('\alpha', 'FontName', 'Cambria Math');

subplot(1, 2, 2)
boxplot(resDouble, 'Labels',{'Control','Rapamycin'})
hold on
plot([1:2], resDouble, 'o', 'Color', [204/255, 204/255, 0/255])

set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Control','Rapamycin'},...
    'FontName', 'Cambria Math', 'FontSize', 12);
ylabel('Residuals', 'FontName', 'Cambria Math');

sgtitle('Double exponential \alpha and residuals',...
    'FontName', 'Cambria Math', 'FontSize', 14);

hold off

% x_saveas=['LH_doubleExp_2.png'];
% cd 'C:\Users\anr612\Documents\1. Harvard18-\2. Rotations\2. DenicLab\MicroscopeImages\23-Jan-19_TrackingImages_LH'
% set(gcf, 'Color', 'w');
% saveas(gcf, x_saveas);


