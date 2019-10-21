clear
clc
dataUpper3=csvread('Upper3.csv');
dataLower3=csvread('Lower3.csv');
% Odd columns are x-coordinates
% Even columns are y-coordinates
xUpper3=dataUpper3(:,1:2:end);yUpper3=dataUpper3(:,2:2:end);
xUpper3=reshape(xUpper3,1,numel(xUpper3));
yUpper3=reshape(yUpper3,1,numel(yUpper3));
xLower3=dataLower3(:,1:2:end);yLower3=dataLower3(:,2:2:end);
xLower3=reshape(xLower3,1,numel(xLower3));
yLower3=reshape(yLower3,1,numel(yLower3));
MarkerSize=5;LineWidth=0.5;
plot(xUpper3,yUpper3,'b*',xLower3,yLower3,'b*','markersize',MarkerSize,...
    'linewidth',LineWidth);
xTip=1.0;yTip=0.0;
hold on
plot(xTip,yTip,'b*','markersize',MarkerSize,'linewidth',LineWidth);
grid on
axis equal 
axis tight





