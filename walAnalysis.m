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
thick=input('Please enter the thickness of samples in nanometers\n')
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

    directionForTag(end+1)=direction{1};
    temperature=regexp(fileName(i).name,patternTemperature, 'match');
    temperatureForCheck(end+1)=temperature{1};
    temperatureForTag(end+1)=[temperature{1},' K'];
    temperatureForNum(end+1)=str2num(temperature{1});
    
    isReadTem=~isempty(find([dataRaw{:,1}]==temperatureForNum(i)));
    if isReadTem==0
        j=i;
    elseif isReadTem==1
        j=find([dataRaw{:,1}]==temperatureForNum(i));
    else
        "Wrong 1"
        return
    end
    
    [x,y] = importfile(fileName(i).name);
    % Remove zero
    m=find(~x);
    
    if isempty(m)==0
        x(m(1)-1:end)=[];
        y(m(1)-1:end)=[];
    end    
    
    x=x./1000;
    dataRaw{j,1}=temperatureForNum(i);
    if direction=="x"
        dataRaw{j,2}=x;
        dataRaw{j,3}=abs(y);
        dataRaw{j,4}=asp*h/(e^2)./abs(y);
    elseif direction=="y"
        dataRaw{j,5}=x;
        dataRaw{j,6}=y;
        temR=dataRaw{j,6};
        temR=flipud(temR);
        dataRaw{j,7}=(dataRaw{j,6}-temR)./2;
    else
        "Something wrong with file name. Wrong 2"
        return
    end
    
    clearvars x y m direction temperature j temR;
end
temperatureForTag(1)=[];
temperatureForTag=strip(temperatureForTag,'left','0');
temperatureForCheck(1)=[];
temperatureForCheck=convertStringsToChars(temperatureForCheck);

% Check data (wrong files)
for i=1:number/2
    for j=2:7
        if isempty(dataRaw{i,j})==1
          "Wrong 3"
          return
        end
    end
end
clearvars i j
%% rHall vs H
figure
for i=1:number/2
   getPlot(dataRaw{i,5},dataRaw{i,6},generalView,"R_{xy} vs H")
   xlabel {H (T)}
   ylabel {R_{xy} (\Omega)}
   hold on
    
end
title('Hall resistance')
legend(temperatureForTag(1:number/2))


%% Hall fitting
for i=1:number/2
    [xData, yData] = prepareCurveData(dataRaw{i,5},dataRaw{i,7});

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
    getFitPlot(fitresult,xData,yData,fittingView,strcat("Hall fitting result of ",temperatureForTag(i)));

    % Label axes
    xlabel {H (T)}
    ylabel {R_{xy} (\Omega)}

    % Properties from Hall fitting.
    ne(i)=abs(1/e/hall(i));
    mu(i)=e^2/h*max(dataRaw{i,4})/e/ne(i);
    taup(i)=me/e*mu(i);
    
    if dimension == 2
        kf(i)=(2*pi*ne(i))^(1/2);
    elseif dimension == 3
        kf(i)=(3*pi^2*ne)^(1/3);
    else
        "Wrong 4"
        return
    end
    
    vf(i)=hbar/me*kf(i);
    dif(i)=vf(i)^2*taup(i)/2;
    be(i)=hbar/(4*e*dif(i)*taup(i));

end
fprintf("Hall fitting finished.\nWaiting for WAL fitting.\n")


%% Rxx vs B
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
figure
for i=1:number/2
    getPlot(dataRaw{i,2},dataRaw{i,3},generalView,"R_{xx} vs B")
    xlabel {H (T)}
    ylabel {R_{xx} (\Omega)}
    hold on
    
end
title('Longitude resistance')
legend(temperatureForTag(1:number/2),'Location','northeastoutside')

%% Gxx(Normalized) vs B
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
figure('Name',"Gxx(Normalized) vs B")
gxxNor=cell(number/2,3);
% Format of gxxNor: {Temperature, hxx, gxxNor}
gxxNor(:,1:2)=dataRaw(:,1:2);
for i=1:number/2
    temY=dataRaw{i,4};
    gxxNor{i,3}=(temY-max(temY))./max(temY);
    getPlot(gxxNor{i,2},gxxNor{i,3},generalView,"Gxx(Normalized) vs B");
    xlabel {H (T)}
    ylabel {G_{xx}(Normalized)}
    hold on
    clearvars temY
end
title('Longitude resistance(Normalized)')
legend(temperatureForTag(1:number/2),'Location','northeastoutside')

%% WAL fitting
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}


for i=1:number/2
    % Pretreatment of data.
    % Get the abstract of hxx.
    xWal = abs(dataRaw{i,2});
    yWal = dataRaw{i,4}-max(dataRaw{i,4});
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
    getFitPlot(fitresult,xData,yData,fittingView,strcat("WAL fitting result of ",temperatureForTag(i)));
%     figure( 'Name', strcat("WAL fitting for ",temperatureForTag(i)) );
%     plot( fitresult, xData, yData ,'o');
%     legend( 'Data points', 'Fitting curve', 'Location', 'NorthEast' );
%     title(strcat("WAL fitting for ",temperatureForTag(i)))
    % Label axes
    xlabel {|H_{ex}| (T)}
    ylabel \DeltaG
%     grid on
    
    clearvars xWal yWal xData yData
    
end

lso=sqrt(hbar./(4.*e.*bso));
lphi=sqrt(hbar./(4.*e.*bi));
ltr=(hbar/e*sqrt(2*pi)).*mu;
tauso=hbar./(4*e.*bso.*dif);
fprintf("WAL fitting finished.\n")

%% Change the units
dif=dif.*10000;
mu=mu.*10000;

%% ne vs temperature
getScatter(temperatureForNum(1:number/2),ne,generalView,"n_e vs Temperature")
xlabel {T (K)}
ylabel {n_e (cm^{-3})}

%% D vs temperature
getScatter(temperatureForNum(1:number/2),dif,generalView,"D vs Temperature")
xlabel {T (K)}
ylabel {D (cm^2/s)}

%% L_SO vs temperature
getScatter(temperatureForNum(1:number/2),lso,generalView,"L_{SO} vs Temperature");
xlabel {T (K)}
ylabel {L_{SO} (m)}

%% L_Phi vs temperature
getScatter(temperatureForNum(1:number/2),lphi,generalView,"L_{\phi} vs Temperature");
xlabel {T (K)}
ylabel {L_\phi (m)}

%% L_SO vs D
getScatter(dif,lso,generalView,"L_{SO} vs D");
xlabel {D (m^2/s)}
ylabel {L_{SO} (m)}

%% Tau_SO vs Tau_p
getScatter(taup,tauso,generalView,"\tau_{SO} vs \tau_p");
xlabel {\tau_p (s)}
ylabel {\tau_{SO} (s)}

%% Tau_SO vs D
getScatter(dif,tauso,generalView,"\tau_{SO} vs D");
xlabel {D (cm^2/s)}
ylabel {\tau_{SO} (s)}

%% mu vs temperature
getScatter(temperatureForNum(1:number/2),mu,generalView,"\mu vs Temperature");
xlabel {T (K)}
ylabel {\mu (cm^2/V/s)}

%% Tau_SO vs temperature
getScatter(temperatureForNum(1:number/2),tauso,generalView,"\tau_{SO} vs Temperature");
xlabel {T (K)}
ylabel {\tau_{SO} (s)}

%% Tau_p vs temperature
getScatter(temperatureForNum(1:number/2),taup,generalView,"\tau_p vs Temperature");
xlabel {T (K)}
ylabel {\tau_p (s)}

%% vf vs temperature
getScatter(temperatureForNum(1:number/2),vf,generalView,"v_f vs Temperature");
xlabel {T (K)}
ylabel {v_f }

%% Bso vs temperature
getScatter(temperatureForNum(1:number/2),bso,generalView,"B_{SO} vs Temperature");
xlabel {T (K)}
ylabel {B_{SO} (T)}


fprintf("Finished\n")



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
