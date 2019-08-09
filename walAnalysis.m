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
        eval(['rxx_',temperature{1},'K=abs(y);'])
        eval(['gxx_',temperature{1},'K=asp*h/(e^2)./abs(y);'])
    elseif direction=="y"
        eval(['hxy_',temperature{1},'K=x;'])
        eval(['rxy_',temperature{1},'K=y;'])
        eval(['gxy_',temperature{1},'K=asp*h/(e^2)./y;'])
    else
        "Something wrong with file name"
        return
    end
    clearvars x y m direction temperature;
end
temperatureForTag(1)=[];
temperatureForTag=strip(temperatureForTag,'left','0');
temperatureForCheck(1)=[];
temperatureForCheck=convertStringsToChars(temperatureForCheck);

%% rHall vs H
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
    hall(end+1)=abs(fitresult.p1);
    
    % Plot fit with data.
    figure( 'Name',strcat("Hall factor fitting for ",temperatureForTag(i)) );
    plotH = plot( fitresult, xData, yData ,'o');
    legend( plotH, 'Data points', 'Fitting curve', 'Location', 'NorthEast' );
    title(strcat("Hall factor fitting for ",temperatureForTag(i)))
    % Label axes
    xlabel {H (T)}
    ylabel {R_{xy} (\Omega)}
    grid on

    % Properties from Hall fitting.
    ne(i)=abs(1/e/hall(i));
    eval(['mu(i)=e^2/h*max(gxx_',temperatureForCheck{i},'K)/e/ne(i);'])
    taup(i)=me/e*mu(i);
    kf(i)=(2*pi*ne(i))^(1/2);  % 2D case.
    %kf(i)=(3*pi^2*ne)^(1/3);  % 3D case.
    vf(i)=hbar/me*kf(i);
    dif(i)=vf(i)^2*taup(i)/2;
    be(i)=hbar/(4*e*dif(i)*taup(i));

end
"Hall fitting finished"


%% Rxx vs B
figure
for i=1:number/2
   eval(['plot(hxx_',temperatureForCheck{i},'K,rxx_',temperatureForCheck{i},...
       'K,''Linewidth'',2)'])
   xlabel {H (T)}
   ylabel {R_{xx} (\Omega)}
   grid on
   hold on
    
end
title('Longitude resistance')
legend(temperatureForTag(1:number/2),'Location','best')

%% Gxx(Normalized) vs B
figure
for i=1:number/2
   eval(['temY=gxx_',temperatureForCheck{i},'K;'])
   eval(['gxx_',temperatureForCheck{i},'K_Nor=(temY-max(temY))./max(temY);'])
   eval(['plot(hxx_',temperatureForCheck{i},'K,gxx_',temperatureForCheck{i},...
       'K_Nor,''Linewidth'',2)'])
   xlabel {H (T)}
   ylabel {G_{xx}(Normalized)}
   grid on
   hold on
   clearvars temY
end
title('Longitude resistance(Normalized)')
legend(temperatureForTag(1:number/2),'Location','best')

%% WAL fitting


for i=1:number/2
    % Pretreatment of data.
    % Get the abstract of hxx.
    eval(['xWal= abs(hxx_',temperatureForCheck{i},'K);'])
    eval(['yWal= gxx_',temperatureForCheck{i},'K-max(gxx_',temperatureForCheck{i},'K);'])
    [xData, yData] = prepareCurveData( xWal, yWal );  
    
    % Set up fittype and options.
    ft = fittype( '1/3.14159*(-psi((bso+be)/x+0.5)+log((bso+be)/x)+1.5*psi((bi+4*bso/3)/x+0.5)-1.5*log((bi+4*bso/3)/x)-0.5*psi(bi/x+0.5)+0.5*log(bi/x))+kFactor*x^2', 'independent', 'x', 'dependent', 'y' , 'problem' , 'be' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMinChange = 1e-12;
    opts.Display = 'Off'; 
    opts.Lower = [0 0];
    opts.MaxFunEvals = 3000;
    opts.MaxIter = 2000;
    opts.Robust = 'LAR';
    opts.StartPoint = [0.02 -0.05 0];
    opts.TolFun = 1e-08;
    opts.TolX = 1e-08;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts , 'problem' , be(i));
    
    bso(end+1)=fitresult.bso;
    bi(end+1)=fitresult.bi;
    
    % Plot fit with data.
    figure( 'Name', strcat("WAL fitting for ",temperatureForTag(i)) );
    plot( fitresult, xData, yData ,'o');
    legend( 'Data points', 'Fitting curve', 'Location', 'NorthEast' );
    title(strcat("WAL fitting for ",temperatureForTag(i)))
    % Label axes
    xlabel {|H| (T)}
    ylabel \DeltaG
    grid on
    
    clearvars xWal yWal xData yData
    
end

lso=sqrt(hbar./(4.*e.*bso));
lphi=sqrt(hbar./(4.*e.*bi));
ltr=(hbar/e*sqrt(2*pi)).*mu;
tauso=hbar./(4*e.*bso.*dif);
"WAL fitting finished"

%% Change the units
dif=dif.*10000;
mu=mu.*10000;

%% Ne vs temperature
figure
scatter(temperatureForNum(1:number/2),ne,'o');
xlim([0 1.1*max(temperatureForNum)])
% upLim=max(ne)+0.1*(max(ne)-min(ne));
% downLim=min(ne)-0.1*(max(ne)-min(ne));
% ylim([downLim upLim])
% clearvars upLim downLim

xlabel {T (K)}
ylabel {n_e (cm^{-3})}
grid on
box on

%% D vs temperature
figure
scatter(temperatureForNum(1:number/2),dif,'o');
xlim([0 1.1*max(temperatureForNum)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {D (cm^2/s)}
grid on
box on

%% L_SO vs temperature
figure
scatter(temperatureForNum(1:number/2),lso,'o');
xlim([0 1.1*max(temperatureForNum)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {L_{SO} (m)}
grid on
box on

%% L_Phi vs temperature
figure
scatter(temperatureForNum(1:number/2),lphi,'o');
xlim([0 1.1*max(temperatureForNum)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {L(\phi) (m)}
grid on
box on

%% L_SO vs D
figure
scatter(dif,lso,'o');
% xlim([0 1.1*max(temperatureForNum)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {D (m^2/s)}
ylabel {L(\phi) (m)}
grid on
box on

%% Tau_SO vs Tau_p
figure
scatter(taup,tauso,'o');
xlabel {\tau_p (s)}
ylabel {\tau_{SO} (s)}
grid on
box on

%% Tau_SO vs D
figure
scatter(dif,tauso,'o');
xlabel {D (cm^2/s)}
ylabel {\tau_{SO} (s)}
grid on
box on

%% mu vs temperature
figure
scatter(temperatureForNum(1:number/2),mu,'o');
xlim([0 1.1*max(temperatureForNum)])
xlabel {T (K)}
ylabel {\mu (cm^2/V/s)}
grid on
box on

%% Tau_SO vs temperature
figure
scatter(temperatureForNum(1:number/2),tauso,'o');
xlim([0 1.1*max(temperatureForNum)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {\tau_{SO} (s)}
grid on
box on

%% Tau_p vs temperature
figure
scatter(temperatureForNum(1:number/2),taup,'o');
xlim([0 1.1*max(temperatureForNum)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {\tau_p (s)}
grid on
box on

%% vf vs temperature
figure
scatter(temperatureForNum(1:number/2),vf,'o');
xlim([0 1.1*max(temperatureForNum)])
xlabel {T (K)}
ylabel {v_f }
grid on
box on

%% Bso vs temperature
figure
scatter(temperatureForNum(1:number/2),bso,'o');
xlim([0 1.1*max(temperatureForNum)])
xlabel {T (K)}
ylabel {B_{SO} (T)}
grid on
box on


"Finished"



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
