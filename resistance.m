clear;close all;clc

%% Import data name
fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.
number=length(fileName);
% Claim variables.

rRaw=cell(number,5);       % Format: Temperature, Current, Voltage, resistance, conductance (1/¦¸).

temperatureForCheck=string;
temperatureForNum=[];

temperatureForTag=string;

%% Import data
% Format: Temperature, Current, Voltage, resistance, conductance (1/¦¸).
for i=1:number
    % Extract the thickness and magnetic field.
    patternTemperature = '(?<=Rxx0_)\d*.\d*(?=K.txt)';
    temperature=regexp(fileName(i).name,patternTemperature, 'match');

    temperatureForCheck(end+1)=temperature{1};
    temperatureForTag(end+1)=[temperature{1},' K'];
    temperatureForNum(end+1)=str2num(temperature{1});
    
    rRaw{i,1} = temperatureForNum(end);
    
    [x,y] = importfile(fileName(i).name);
    % Remove zero
    m=find(~x);
    
    if isempty(m)==0
        x(m(1)-1:end)=[];
        y(m(1)-1:end)=[];
    end    
    
    rRaw{i,2} = x./1000;
    rRaw{i,3} = y;
    
    clearvars x y m
end

rRaw = sortrows(rRaw,1); 

%% Fitting
for i = 1:number
    [rTem, cTem] = createFit(rRaw{i,2}, rRaw{i,3});
    rRaw{i,4} = rTem;
    rRaw{i,5} = 1/rTem;
    
    clearvars rTem cTem
end


%% Plot resistance
figR=figure
plot([rRaw{:,1}],[rRaw{:,4}]./1000,'o-','LineWidth',2.5)
set(gcf,'position',[800 100 800 600]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)
set(figR,'Name',"Resistance");
title("Resistance");
xlabel("T (K)")
ylabel("Resistance (k\Omega)")
grid on

%% Plot conductance
figG=figure
plot([rRaw{:,1}],[rRaw{:,5}],'o-','LineWidth',2.5)
set(gcf,'position',[800 100 800 600]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)
set(figG,'Name',"Conductance");
title("Conductance");
xlabel("T (K)")
ylabel("Conductance (1/\Omega)")
grid on

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
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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
function [r, c] = createFit(x, y)
[xData, yData] = prepareCurveData( x, y );


ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );
r = fitresult.p1;
c = fitresult.p2;

% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
% 
% xlabel x
% ylabel y
% grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%