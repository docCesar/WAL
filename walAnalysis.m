clear;close all;delete(gcp('nocreate'));clc;

status=rmdir('temPlots','s');
status=rmdir('temFitPlots','s');
mkdir temPlots;
mkdir temFitPlots;
clearvars status

% Number of cores
coreNumber = 8;

% Constants
h = 6.6260699*10^-34;
hbar = 1.0545718e-34;
% me = 9.10938356*10^(-31);
me = 0.15 *  9.10938356*10^(-31);
e = 1.60217662*10^-19;

% Parameters of samples
dimension=2;
% asp=4;          % Old mask
% asp=2.5;        % Small devices of new mask
asp=24;         % Large devices of new mask
% thick=5*10^-9; % Thickness of sample layer
thick=input('Please enter the thickness of samples in nanometers\n')
thicknessForTag=strcat(thick, 'nm');
% thick=thick*10^-9; % Nanometer
fittingView=false;
fittingView=input('Do you want to check the status of fitting curve ?\nPlease enter 1 for YES or 0 for NO.\n');
fittingView=logical(fittingView);
generalView=true;

tic
%% Import data name
fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.
number=length(fileName);

% Claim variables.

% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy}
dataRaw=cell(number/2,6);
dataMod=cell(number/2,6);
sigma0=zeros(1,number/2);
hall=zeros(1,number/2);
confidenceHall=zeros(1,number/2);
mu=zeros(1,number/2);
taup=zeros(1,number/2);
vf=zeros(1,number/2);
dif=zeros(1,number/2);
be=zeros(1,number/2);
ne=zeros(1,number/2);
neBulk=zeros(1,number/2);
kf=zeros(1,number/2);
bi=zeros(1,number/2);
bso=zeros(1,number/2);
confidenceBi=zeros(1,number/2);
confidenceBso=zeros(1,number/2);
kFactor=zeros(1,number/2);
fitresultWAL=cell(3,number/2);      % Format:(xData;yData;fitresult)

directionForTag=string;
temperatureForCheck=string;
temperatureForNum=zeros(1,number/2);

temperatureForTag=string;

%% Import data
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy}
for i=1:number
    % Extract the thickness and magnetic field.
    patternDirection = '(?<=Rx)\w(?=_)';
    patternTemperature = '(?<=Rx(x|y)_)\d*(?=K_)';
    direction=regexp(fileName(i).name,patternDirection, 'match');

    directionForTag(end+1)=direction{1};
    temperature=regexp(fileName(i).name,patternTemperature, 'match');
    temperatureForCheck(i)=temperature{1};
    temperatureForTag(end+1)=[temperature{1},' K'];
    temperatureForNum(i)=str2num(temperature{1});
    
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
        if x(m(1))==x(m(1)+1)
            x(m(1)-1:end)=[];
            y(m(1)-1:end)=[];
        else
            x(m(1))=[];
            y(m(1))=[];
        end
    end    
    
    x=x./1000;
    dataRaw{j,1}=temperatureForNum(i);
    if direction=="x"
        dataRaw{j,2}=x;
        dataRaw{j,3}=abs(y);
        dataRaw{j,4}=asp*h/(e^2)./abs(y);
        sigma0(j)=max(dataRaw{j,4});
    elseif direction=="y"
        dataRaw{j,5}=x;
        dataRaw{j,6}=y;
    else
        error('ERROR 02: Something wrong with file name.');
        return
    end
    
    clearvars x y m direction temperature j temR isReadTem;
end
temperatureForTag(1)=[];
temperatureForTag=strip(temperatureForTag,'left','0');
temperatureForCheck(1)=[];
temperatureForCheck=convertStringsToChars(temperatureForCheck);

% Check data (wrong files)
for i=1:number/2
    for j=2:6
        if isempty(dataRaw{i,j})==1
          error('ERROR 03');
          return
        end
    end
end
clearvars i j

% dataMod(:,1)=dataRaw(:,1);
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
% parpool(4)
for i=1:number/2
    [xData, yData] = prepareCurveData(dataRaw{i,5},dataRaw{i,6});

    %% Remove symmetric component
    % Set up fittype and options.
    ftHallLinear = 'linearinterp';

    % Fit model to data.
    [fitresultHallLinear, gofHallLinear] = fit( xData, yData, ftHallLinear, 'Normalize', 'on' );
    
    % Substrate the linear component
    yData=(yData-fitresultHallLinear(-xData))./2;
    
    % Plot fit with data.
%     figure( 'Name', 'untitled fit 1' );
%     plot( fitresultHallLinear, xData, yData );
    
    %%
    % Set up fittype and options.
    xData = xData(50:end-50);
    yData = yData(50:end-50);
    ft = fittype( 'p1*x+p2', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [1 0.1];
    opts.TolFun = 1e-10;
    opts.TolX = 1e-10;

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Extract fitting parameters.
    hall(i)=abs(fitresult.p1);
    temConfidence=confint(fitresult,0.95);
    confidenceHall(i)=abs((temConfidence(1,1)-temConfidence(2,1))/2);
    
    % Plot fit with data.
    getFitPlot(fitresult,xData,yData,fittingView,strcat("Hall fitting result of ",temperatureForTag(i)),1,"R_{xy} (\Omega)");

    clearvars xData yData ft fitresult gof opts ftHallLinear gofHallLinear fitresultHallLinear gofHallLinear
    % Label axes
%     xlabel {H (T)}
%     ylabel {R_{xy} (\Omega)}
end
%     delete(gcp('nocreate'))
% Properties from Hall fitting.
ne=abs(1./e./hall);
neBulk=ne/(thick*1e-9);


if dimension == 2
    kf=(2.*pi.*ne).^(1/2);
    mu=(e^2/h).*sigma0./(e.*ne);
elseif dimension == 3
    kf=(3.*pi^2.*ne).^(1/3);
    mu=(e^2/h).*sigma0./(e.*neBulk);
else
    fprintf("Wrong dimension.\n")
    error('ERROR 04');
    return
end

taup=me./e*mu;
vf=hbar./me.*kf;
dif=vf.^(dimension).*taup./dimension;
be=hbar./(4.*e.*dif.*taup);

fprintf("Hall fitting finished.\nWaiting for WAL fitting.\n")


%% Rxx vs B
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy}
figure
for i=1:number/2
    getPlot(dataRaw{i,2},dataRaw{i,3},generalView,"R_{xx} vs B")
    xlabel {H (T)}
    ylabel {R_{xx} (\Omega)}
    hold on
    
end
title('Longitudinal resistance')
legend(temperatureForTag(1:number/2),'Location','northeastoutside')

path=strcat("temPlots\","Longitudinal resistance");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% Gxx vs B
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy}
figure
for i=1:number/2
    dg=dataRaw{i,4}-max(dataRaw{i,4});
    getPlot(dataRaw{i,2},dg,generalView,"G_{xx} vs H")
    xlabel {H (T)}
    ylabel {\Delta\sigma_{xx} (e^2/h)}
    hold on
    
end
title('Longitudinal conductance')
legend(temperatureForTag(1:number/2),'Location','northeastoutside')
clearvars dg

path=strcat("temPlots\","Longitudinal conductance");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% Gxx(Normalized) vs B
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy}
figure('Name',"Gxx(Normalized) vs B")
gxxNor=cell(number/2,3);
% Format of gxxNor: {Temperature, hxx, gxxNor}
gxxNor(:,1:2)=dataRaw(:,1:2);
for i=1:number/2
    temY=dataRaw{i,4};
    gxxNor{i,3}=(temY-max(temY))./max(temY);
    getPlot(gxxNor{i,2},gxxNor{i,3},generalView,"Gxx(Normalized) vs B");
    xlabel {H (T)}
    ylabel {\sigma_{xx}(Normalized)}
    hold on
    clearvars temY
end
title('Longitudinal conductivity (Normalized)')
legend(temperatureForTag(1:number/2),'Location','northeastoutside')

path=strcat("temPlots\","Longitudinal conductance (Normalized)");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% WAL fitting
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy}
dataPoint = cell(number/2,1);
fittingCurve = cell(number/2,1);
outputFitting = cell(number/2,1);
gofHLN = cell(number/2,1);


parpool('local',coreNumber);
parfor i=1:number/2
% for i=1:number/2
    %%
    % Pretreatment of data.
    % Get the abstract of hxx.
    temDataMod = cell(1, 5);
    temFitresultWAL=cell(3,1);
    xWal = dataRaw{i,2};
    yWal = dataRaw{i,3};
    [xData, yData] = prepareCurveData( xWal, yWal );  
    
    %% Remove disturbed components
    % Set up fittype and options.
    % ftInter = 'smoothingspline';
    ftInter = 'nearestinterp';
    
    % Fit model to data.
%     [fitresultInter, gof] = fit( xData, yData, ftInter);
    [fitresultInter, gof] = fit( xData, yData, ftInter, 'Normalize', 'on');
    
    
    % Remove the asymmetric component
    yData = abs( yData + fitresultInter(-xData) ) ./ 2;
%     yData = abs( yData + fitresultInter(-xData) ) ./ 2;
    
    % Select the data in high field
    temx = xData((xData > -8.9 & xData < -6) | (xData > 6 & xData < 8.9));
    temy = yData((xData > -8.9 & xData < -6) | (xData > 6 & xData < 8.9));
    weight = temx .^ 16;
    temDataMod{5}=xData;
    temDataMod{6}=yData;
    
    % Remove the classic component
    [temx, temy, weight] = prepareCurveData( temx, temy, weight );
    ft = fittype( 'p1*x^2 + p0', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.01 0];
    opts.Weights = weight;
    [fitresultRxxCali, gofRxxCali] = fit( temx, temy, ft, opts );
    
    % Plot the classic component
    getFitPlot(fitresultRxxCali,xWal,yWal,fittingView,join(["Classic component and raw data for ", temperatureForTag(i)]),temperatureForNum(i),"R_{xx} (\Omega)");
    
    %% Prepare for HLN fitting
    yData = yData - fitresultRxxCali(xData) + fitresultRxxCali.p0;
    
    % Save the modified data
    temDataMod{1} = dataRaw{i, 1};
    temDataMod{2} = xData;
    temDataMod{3} = yData;
    yDataSigma = asp.*h./(e.^2)./yData;
    yDataSigma = yDataSigma - max(yDataSigma);
    temDataMod{4} = yDataSigma;
    
    % Plot modified data.
    getScatter(xData,yDataSigma,generalView,["Modified conductance of "+temperatureForTag(i)],["H_{ex} (T)", "\Delta\sigma_{xx} (e^2/h)"]);

    %% HLN fitting
    
    % Set up fittype and options.
    ft = fittype( '1/(pi)*(-psi((bso+be)/abs(x)+0.5)+log((bso+be)/abs(x))+1.5*psi((bi+4*bso/3)/abs(x)+0.5)-1.5*log((bi+4*bso/3)/abs(x))-0.5*psi(bi/abs(x)+0.5)+0.5*log(bi/abs(x)))', 'independent', 'x', 'dependent', 'y' , 'problem' , 'be' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMinChange = 1e-18;
    opts.Display = 'Off'; 
    opts.Lower = [0 0];
%     opts.Upper = [5 80 Inf ];
    opts.MaxFunEvals = 100000;
    opts.MaxIter = 100000;
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.002 0.1];
    opts.TolFun = 1e-13;
    opts.TolX = 1e-13;
    
    % Fit model to data.
    [fitresult, gofHLN{i,1}, outputFitting{i,1}] = fit( xData, yDataSigma, ft, opts , 'problem' , be(i));
%     [fitresult, gof] = fit( xData, yData, ft, opts );
    
    bso(i)=fitresult.bso;
    bi(i)=fitresult.bi;
    temConfidence=confint(fitresult);
    confidenceBi(i)=abs((temConfidence(1,1)-temConfidence(2,1))/2);
    confidenceBso(i)=abs((temConfidence(1,2)-temConfidence(2,2))/2);
    
    temFitresultWAL = cell(3,1);
    temFitresultWAL{1} = xData;
    temFitresultWAL{2} = yDataSigma;
    temFitresultWAL{3} = fitresult;
    fitresultWAL(:,i) = temFitresultWAL;
    dataMod(i, :) = temDataMod;

    % Plot fit with data.
    [dataPoint{i,1},fittingCurve{i,1}]=getFitPlot(fitresult,xData,yDataSigma,fittingView,strcat("WAL fitting result of ",temperatureForTag(i)),256,"\Delta\sigma (e^2/h)");     


    
end
delete(gcp('nocreate'))
clearvars xWal yWal xData yData fitresult ftLinear gof temConfidence temFitresultWAL opts fitresultLinear gofLinear temx temy weight
%% Plot fitting with data
dataPoint=cell(number/2,1);
fittingCurve=cell(number/2,1);
for i=1:number/2
    [dataPoint{i,1},fittingCurve{i,1}] = getFitPlot(fitresultWAL{3,i},fitresultWAL{1,i},fitresultWAL{2,i},fittingView,strcat("WAL fitting result of ",temperatureForTag(i)),256,"\Delta\sigma (e^2/h)");
end
legend(fliplr([dataPoint{:,1},fittingCurve{1,1}]),fliplr([temperatureForTag(1:number/2),"Fitting curve"]),'Location','southeastoutside')
legend('boxoff')
saveas(gcf,"temPlots\WAL fitting results",'meta');
print(gcf,"temPlots\WAL fitting results",'-djpeg','-r600')
clearvars dataPoint fittingCurve
%%
lso=sqrt(hbar./(4.*e.*bso));
lphi=sqrt(hbar./(4.*e.*bi));
ltr=sqrt(hbar./(4.*e.*be));
tauso=hbar./(4*e.*bso.*dif);
tauphi=hbar./(4*e.*bi.*dif);
gCritical=kf.*ltr;           % g¡Ô¦Ò/(e2/h)=kf*ltr=1  dimensionless conductivity,g>>1 for conductor.

fprintf("WAL fitting finished.\n")

%% Change the units
dif = dif.*10000;
mu = mu.*10000;
ne = ne./10000;
neBulk = neBulk./1000000;
tauso = tauso.*10^15;        % Femtosecond
tauphi = tauphi.*10^15;      % Femtosecond
taup = taup.*10^15;          % Femtosecond

lphi = lphi.*10^9;           % Nanometer
lso = lso.*10^9;           % Nanometer
ltr = ltr.*10^9;           % Nanometer

%% Hall coefficients vs thickness
getScatter(temperatureForNum(1:number/2),hall,generalView,"Hall coefficient vs Temperature",{"Temperature (K)","H_e (m^2/C)"});

getErrorBar(temperatureForNum(1:number/2),hall,confidenceHall,generalView,"Hall coefficient vs Temperature with error bar",{"Temperature (K)","H_e (m^2/C)"})

%% ne (Sheet) vs temperature
getScatter(temperatureForNum(1:number/2),ne,generalView,"n_e (Sheet) vs Temperature",{"Temperature (K)","n_e (Sheet) (cm^{-2})"})
% xlabel {T (K)}
% ylabel {n_e (cm^{-2})}

%% ne (Bulk) vs temperature
getScatter(temperatureForNum(1:number/2),neBulk,generalView,"n_e (Bulk) vs Temperature",{"Temperature (K)","n_e (Bulk) (cm^{-3})"})
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
getScatter(dif,lso,generalView,"L_{SO} vs D",{"D (cm^2/s)","L_{SO} (nm)"});
% xlabel {D (cm^2/s)}
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
%%
getScatter(dif,tauso,generalView,"\tau_{SO} vs D",{"D (cm^2/s)","\tau_{SO} (fs)"});
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

getErrorBar(temperatureForNum(1:number/2),bso,confidenceBso,generalView,"B_{SO} vs Temperature with error bar",{"Temperature (K)","B_{SO} (T)"})

%% Bi vs temperature
getScatter(temperatureForNum(1:number/2),bi,generalView,"B_i vs Temperature",{"Temperature (K)","B_i (T)"});
% xlabel {T (K)}
% ylabel {B_i (T)}

getErrorBar(temperatureForNum(1:number/2),bi,confidenceBi,generalView,"B_i vs Temperature with error bar",{"Temperature (K)","B_i (T)"})

clearvars i path
fprintf("Program completed.\n")
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18);
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
function [dataPoint,fittingCurve] = getFitPlot(fitresult,xData,yData,fittingView,titleOfPlot,n,yLabel)

if nargin < 4
    fittingView=true
end

if nargin >= 6
    if n == 256
        fig=figure(n);
        set(fig,'Name',"WAL fitting result");
        title("WAL fitting result");
        hold on
    else
        fig=figure;
        set(fig,'Name',titleOfPlot);
        title(titleOfPlot);
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

dataPoint=plot(xData, yData ,'o','LineWidth',2);
hold on
fittingCurve=plot(fitresult,'r:');
fittingCurve.LineWidth=4;
%,'LineStyle',':', 'LineWidth',2
xlabel {H_{ex} (T)}
if nargin >= 7
    ylabel(yLabel)
end
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

if nargin >= 4
    set(gcf,'Name',titleOfPlot);
    title(titleOfPlot);
end

grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getErrorBar(x,y,confidence,generalView,titleOfPlot,label)

if nargin < 4
    generalView=true
end

fig=figure;

if generalView==false
    set(fig,'Visible','off')
end

errorbar(x,y,confidence,'LineWidth',2.5)
set(gcf,'position',[800 100 800 600]);
xlim([1.1*min(x)-0.1*max(x) 1.1*max(x)-0.1*min(x)]);
ylim([1.1*min(y)-0.1*max(y) 1.1*max(y)-0.1*min(y)]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)
if nargin >= 5
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
end
grid on
if nargin >= 6
    xlabel(label{1})
    ylabel(label{2})
end
titleOfPlot=strrep(titleOfPlot,'\','');
path=strcat("temPlots\",titleOfPlot);
saveas(fig,path,'meta');
saveas(fig,path,'jpeg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%