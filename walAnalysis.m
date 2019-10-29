clear;close all;clc

status=rmdir('temPlots','s');
status=rmdir('temFitPlots','s');
mkdir temPlots;
mkdir temFitPlots;

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
fittingView=input('Do you want to check the status of fitting curve ?\nPlease enter 1 for YES or 0 for NO.\n');
fittingView=logical(fittingView);
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
kFactor=[];
kkFactor=[];

directionForTag=string;
temperatureForCheck=string;
temperatureForNum=[];

temperatureForTag=string;

%% Import data
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
for i=1:number
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
        error('ERROR 01');
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
        error('ERROR 02: Something wrong with file name.');
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
          error('ERROR 03');
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

path=strcat("temPlots\","Hall resistance");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

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
        error('ERROR 04');
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

path=strcat("temPlots\","Longitude resistance");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% Gxx vs B
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
figure
for i=1:number/2
    dg=dataRaw{i,4}-max(dataRaw{i,4});
    getPlot(dataRaw{i,2},dg,generalView,"G_{xx} vs H")
    xlabel {H (T)}
    ylabel {\DeltaG_{xx} (e^2/h)}
    hold on
    
end
title('Longitude conductance')
legend(temperatureForTag(1:number/2),'Location','northeastoutside')
clearvars dg

path=strcat("temPlots\","Longitude conductance");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

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
title('Longitude conductivity (Normalized)')
legend(temperatureForTag(1:number/2),'Location','northeastoutside')

path=strcat("temPlots\","Longitude conductance (Normalized)");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% WAL fitting
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}


for i=1:number/2
    %%
    % Pretreatment of data.
    % Get the abstract of hxx.
    xWal = dataRaw{i,2};
    yWal = dataRaw{i,4}-max(dataRaw{i,4});
    [xData, yData] = prepareCurveData( xWal, yWal );  
    
    % Set up fittype and options.
    ft = fittype( '1/3.14159*(-psi((bso+be)/abs(x)+0.5)+log((bso+be)/abs(x))+1.5*psi((bi+4*bso/3)/abs(x)+0.5)-1.5*log((bi+4*bso/3)/abs(x))-0.5*psi(bi/abs(x)+0.5)+0.5*log(bi/abs(x)))+kFactor*x^2+kkFactor*x', 'independent', 'x', 'dependent', 'y' , 'problem' , 'be' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMinChange = 1e-12;
    opts.Display = 'Off'; 
    opts.Lower = [0 0 -Inf -Inf];
    opts.Upper = [2 80 Inf Inf];
    opts.MaxFunEvals = 50000;
    opts.MaxIter = 40000;
    opts.Robust = 'LAR';
    opts.StartPoint = [0.002 1 -0.00001 -0.00001];
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts , 'problem' , be(i));
    
    kFactor(end+1)=fitresult.kFactor;
    kkFactor(end+1)=fitresult.kkFactor;
    bso(end+1)=fitresult.bso;
    bi(end+1)=fitresult.bi;
    
    % Plot fit with data.
    getFitPlot(fitresult,xData,yData,fittingView,strcat("WAL fitting result of ",temperatureForTag(i)),256);
%     figure( 'Name', strcat("WAL fitting for ",temperatureForTag(i)) );
%     plot( fitresult, xData, yData ,'o');
%     legend( 'Data points', 'Fitting curve', 'Location', 'NorthEast' );
%     title(strcat("WAL fitting for ",temperatureForTag(i)))
    % Label axes
    xlabel {H_{ex} (T)}
    ylabel {\DeltaG (e^2/h)}
%     grid on
    
    clearvars xWal yWal xData yData
    
end

lso=sqrt(hbar./(4.*e.*bso));
lphi=sqrt(hbar./(4.*e.*bi));
ltr=sqrt(hbar./(4.*e.*be));
tauso=hbar./(4*e.*bso.*dif);
tauphi=hbar./(4*e.*bi.*dif);
gCritical=kf.*ltr;           % g¡Ô¦Ò/(e2/h)=kf*ltr=1  dimensionless conductivity,g>1 for conductor.

fprintf("WAL fitting finished.\n")

%% Change the units
dif = dif.*10000;
mu = mu.*10000;
tauso = tauso.*10^15;        % Femtosecond
tauphi = tauphi.*10^15;      % Femtosecond
taup = taup.*10^15;          % Femtosecond

lphi = lphi.*10^9;           % Nanometer
lso = lso.*10^9;           % Nanometer
ltr = ltr.*10^9;           % Nanometer

%% ne vs temperature
getScatter(temperatureForNum(1:number/2),ne,generalView,"n_e vs Temperature",{"Temperature (K)","n_e (cm^{-2})"})
% xlabel {T (K)}
% ylabel {n_e (cm^{-2})}

%% D vs temperature
getScatter(temperatureForNum(1:number/2),dif,generalView,"D vs Temperature",{"Temperature (K)","D (cm^2/s)"})
% xlabel {T (K)}
% ylabel {D (cm^2/s)}

%% L_SO vs temperature
getScatter(temperatureForNum(1:number/2),lso,generalView,"L_{SO} vs Temperature",{"Temperature (K)","L_{SO} (nm)"});
% xlabel {T (K)}
% ylabel {L_{SO} (m)}

%% L_Phi vs temperature
getScatter(temperatureForNum(1:number/2),lphi,generalView,"L_{\phi} vs Temperature",{"Temperature (K)","L_\phi (nm)"});
% xlabel {T (K)}
% ylabel {L_\phi (m)}

%% L_SO vs D
getScatter(dif,lso,generalView,"L_{SO} vs D",{"D (m^2/s)","L_{SO} (nm)"});
% xlabel {D (m^2/s)}
% ylabel {L_{SO} (m)}

%% Tau_SO vs Tau_p
getScatter(taup,tauso,generalView,"\tau_{SO} vs \tau_p",{"\tau_p (fs)","\tau_{SO} (fs)"});
% xlabel {\tau_p (fs)}
% ylabel {\tau_{SO} (fs)}

%% Tau_SO vs D
getScatter(dif,tauso,generalView,"\tau_{SO} vs D",{"D (cm^2/s)","\tau_{SO} (fs)"});
% xlabel {D (cm^2/s)}
% ylabel {\tau_{SO} (fs)}

%% mu vs temperature
getScatter(temperatureForNum(1:number/2),mu,generalView,"\mu vs Temperature",{"Temperature (K)","\mu (cm^2/V/s)"});
% xlabel {T (K)}
% ylabel {\mu (cm^2/V/s)}

%% Tau_SO vs temperature
getScatter(temperatureForNum(1:number/2),tauso,generalView,"\tau_{SO} vs Temperature",{"Temperature (K)","\tau_{SO} (fs)"});
% xlabel {T (K)}
% ylabel {\tau_{SO} (fs)}

%% Tau_phi vs thickness
getScatter(temperatureForNum(1:number/2),tauphi,generalView,"\tau_{\phi} vs Temperature",{"Temperature (K)","\tau_{\phi} (fs)"});
% xlabel {T (K)}
% ylabel {\tau_{\phi} (fs)}

%% Tau_p vs temperature
getScatter(temperatureForNum(1:number/2),taup,generalView,"\tau_p vs Temperature",{"Temperature (K)","\tau_p (fs)"});
% xlabel {T (K)}
% ylabel {\tau_p (fs)}

%% vf vs temperature
getScatter(temperatureForNum(1:number/2),vf,generalView,"v_f vs Temperature",{"Temperature (K)","v_f "});
% xlabel {T (K)}
% ylabel {v_f }

%% gCritical vs temperature
getScatter(temperatureForNum(1:number/2),gCritical,generalView,"g_{Cri} vs Temperature",{"Temperature (K)","g_{Cri} (a.u.)"});
% xlabel {T (K)}
% ylabel {v_f }

%% Bso vs temperature
getScatter(temperatureForNum(1:number/2),bso,generalView,"B_{SO} vs Temperature",{"Temperature (K)","B_{SO} (T)"});
% xlabel {T (K)}
% ylabel {B_{SO} (T)}

%% Bi vs temperature
getScatter(temperatureForNum(1:number/2),bi,generalView,"B_i vs Temperature",{"Temperature (K)","B_i (T)"});
% xlabel {T (K)}
% ylabel {B_i (T)}

fprintf("Program completed.\n")



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
function getScatter(x,y,generalView,titleOfPlot,label)

if nargin < 3
    generalView=true
end

fig=figure;

if generalView==false
    set(fig,'Visible','off')
end

plot(x,y,'o','LineWidth',2.5)
set(gcf,'position',[800 100 800 600]);
xlim([1.1*min(x)-0.1*max(x) 1.1*max(x)-0.1*min(x)]);
ylim([1.1*min(y)-0.1*max(y) 1.1*max(y)-0.1*min(y)]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)
if nargin >= 4
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
end
grid on
if nargin >= 5
    xlabel(label{1})
    ylabel(label{2})
end
titleOfPlot=strrep(titleOfPlot,'\','');
path=strcat("temPlots\",titleOfPlot);
saveas(fig,path,'meta');
saveas(fig,path,'jpeg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getFitPlot(fitresult,xData,yData,fittingView,titleOfPlot,n)

if nargin < 4
    fittingView=true
end

if nargin >= 6
    if n == 256
        fig=figure(n);
        set(fig,'Name',"WAL fitting result");
        title("WAL fitting result");
        hold on
    end
else
    fig=figure;
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
end

set(gcf,'position',[800 100 800 600]);

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


grid on
box on


titleOfPlot=strrep(titleOfPlot,'\','');
path=strcat("temFitPlots\",titleOfPlot);
saveas(fig,path,'meta');
saveas(fig,path,'jpeg');
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
