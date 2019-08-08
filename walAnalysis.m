clear;close all;clc

% Constants
h=6.6260699*10^-34;
hbar=1.0545718e-34;
me=9.10938356*10^(-31);
e=1.60217662*10^-19;

% Parameters of samples
% asp=4;          % Old mask
% asp=2.5;        % Small devices of new mask
asp=24;         % Large devices of new mask
thick=5*10^-9; % Thickness of sample layer
% thick=input('Please enter the thickness of samples in nanometers\n')
% thick=thick*10^-9; % Nanometer

%% Import data name
fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.
number=length(fileName);
% Claim variables.
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
directionForTag=string;
temperatureForCheck=string;
temperatureForNum=[];
temperatureForTag=string;

%% Import data
% Thickness need to be added.
for i=1:number
    % Here we can add an if to separate 2 types of measurements.
    % Extract the thickness and magnetic field.
    patternDirection = '(?<=Rx)\w(?=_)';
    patternTemperature = '(?<=Rx(x|y)_)\d*(?=K_)';
    direction=regexp(fileName(i).name,patternDirection, 'match');
    % Ô¤·ÖÅäÄÚ´æ£¿
    directionForTag(end+1)=direction{1};
    temperature=regexp(fileName(i).name,patternTemperature, 'match');
    temperatureForCheck(end+1)=temperature{1};
    temperatureForNum(end+1)=str2num(temperature{1});
    temperatureForTag(end+1)=[temperature{1},' K'];
    
    [x,y] = importfile(fileName(i).name);
    % Remove zero
    m=find(~y);
    if isempty(m)==0
        x(m(1)-1:end)=[];
        y(m(1)-1:end)=[];
    end
    x=x./1000;
    if direction=="x"
        eval(['hxx_',temperature{1},'K=x;'])
        eval(['rxx_',temperature{1},'K=y;'])
        eval(['gxx_',temperature{1},'K=asp*h/(e^2)./y;'])
    elseif direction=="y"
        eval(['hxy_',temperature{1},'K=x;'])
        eval(['rxy_',temperature{1},'K=y;'])
        eval(['gxy_',temperature{1},'K=asp*h/(e^2)./y;'])
    else
        "Something wrong with file name"
    end
    clearvars x y m direction temperature;
end
temperatureForTag(1)=[];
temperatureForTag=strip(temperatureForTag,'left','0');
temperatureForCheck(1)=[];
temperatureForCheck=convertStringsToChars(temperatureForCheck);

%% rHall vs B
figure
for i=1:number/2
   eval(['plot(hxy_',temperatureForCheck{i},'K,rxy_',temperatureForCheck{i},'K,''Linewidth'',2)'])
   xlabel {H (T)}
   ylabel {R_{xy} (\Omega)}
   grid on
   hold on
    
end
title('Hall resistance')
legend(temperatureForTag(1:number/2))

%% Hall calibration
for i=1:number/2
    eval(['temR=rxy_',temperatureForCheck{i},'K;'])
    temR=flipud(temR);
    eval(['rxyCali_',temperatureForCheck{i},'K=(rxy_',temperatureForCheck{i},...
        'K-temR)./2;'])
end
clearvars temR;

%% Hall fitting
for i=1:number/2
    eval(['[xData, yData] = prepareCurveData( hxy_',temperatureForCheck{i},...
        'K, rxyCali_',temperatureForCheck{i},'K);'])
%     [xData, yData] = prepareCurveData( bH(:,i+number/2), rHN(:,i+number/2) );

    % Set up fittype and options.
    ft = fittype( 'p1*x+p2', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.9134 0.0975404049994095];
    opts.TolFun = 1e-08;
    opts.TolX = 1e-08;

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Extract fitting parameters.
    fitresult.p1
    hall(end+1)=abs(fitresult.p1)
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    plotH = plot( fitresult, xData, yData );
    legend( plotH, 'rH vs. bH', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel H
    ylabel {R_xy}
    grid on

    % Properties from Hall fitting.
    ne(i)=abs(1/e/hall(i));
    eval(['mu(i)=e^2/h*max(gxx_',temperatureForCheck{i},'K)/e/ne(i);'])
%     mu(i)=e^2/h*max(gWAL(:,i))/e/ne(i);
    taup(i)=me/e*mu(i);
    kf(i)=(2*pi*ne(i))^(1/2);  % 2D case.
    %kf(i)=(3*pi^2*ne)^(1/3);  % 3D case.
    vf(i)=hbar/me*kf(i);
    dif(i)=vf(i)^2*taup(i)/2;
    be(i)=hbar/(4*e*dif(i)*taup(i));

end
'Hall finish'
% return


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
