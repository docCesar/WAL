clear;close all;clc

% Constants
h=6.6260699*10^-34;
hbar=1.0545718e-34;
me=9.10938356*10^(-31);
e=1.60217662*10^-19;

% Parameters of samples
dimension=2;
% asp=4;          % Old mask
% asp=2.5;        % Small devices of new mask
asp=24;         % Large devices of new mask
% thick=5*10^-9; % Thickness of sample layer
% thick=input('Please enter the thickness of samples in nanometers\n')
% thick=thick*10^-9; % Nanometer

fittingView=false;
fittingView=input('Do you want to check the status of fitting curve ?\nPlease enter 1 for YES or 0 for NO.\n')
fittingView=logical(fittingView)
generalView=true;

%% Import data name
fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.
number=length(fileName);
% Claim variables.

% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
dataRaw=cell(number/2,7);

hall=[];
mu=[];
taup=[];
vf=[];
dif=[];
be=[];
ne=[];
kf=[];
bi=[];
bso=[];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = importfile(filename, startRow, endRow)
%% Initialization
delimiter = '\t';
if nargin<=2
    startRow = 3;
    endRow = inf;
end

%% Format

formatSpec = '%f%f%[^\n\r]';

%% Open the data
fileID = fopen(filename,'r');

%% Read the data column.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end
%% Close the data.
fclose(fileID);

%% Name the data column.

x = dataArray{:, 1};
y = dataArray{:, 2};
x(end) = [];
y(end) = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getScatter(x,y,generalView,titleOfPlot)

if nargin < 3
    generalView=true
end

fig=figure;

if generalView==false
    set(fig,'Visible','off')
end

plot(x,y,'o','LineWidth',2.5)
set(gcf,'position',[1800 100 800 600]);
xlim([1.1*min(x)-0.1*max(x) 1.1*max(x)-0.1*min(x)]);
ylim([1.1*min(y)-0.1*max(y) 1.1*max(y)-0.1*min(y)]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)
if nargin == 4
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
end
grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getFitPlot(fitresult,xData,yData,fittingView,titleOfPlot)

if nargin < 4
    fittingView=true
end

fig=figure;
set(gcf,'position',[1800 100 800 600]);

if fittingView==false
    set(fig,'Visible','off')
end

dataPoint=plot(xData, yData ,'o');
hold on
fittingCurve=plot(fitresult,'r:');
fittingCurve.LineWidth=4;
%,'LineStyle',':', 'LineWidth',2
legend({"Data Point","Fitting Curve"});
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)

if nargin == 5
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
end

grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getPlot(x,y,generalView,titleOfPlot)

if nargin < 3
    generalView=true
end

if generalView==false
    set(fig,'Visible','off')
end

plot(x,y,'Linewidth',2);
set(gcf,'position',[1000 100 800 600]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)

if nargin == 4
    set(gcf,'Name',titleOfPlot);
    title(titleOfPlot);
end

grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
