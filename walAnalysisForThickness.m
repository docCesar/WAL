clear;close all;clc

status=rmdir('temPlots','s');
status=rmdir('temFitPlots','s');
mkdir temPlots;
mkdir temFitPlots;

% Number of cores
% coreNumber=4;

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
fittingView=logical(fittingView);
generalView=true;

tic
%% Import data name
fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.
number=length(fileName);
% Claim variables.

% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
dataRaw=cell(number/2,6);
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

temperatureForCheck=string;
temperatureForNum=zeros(1,number/2);
temperatureForTag=string;

directionForTag=string;
thicknessForCheck=string;
thicknessForNum=zeros(1,number/2);
thicknessForTag=string;


%% Import data
% dataRaw format: {thickness, hxx, rxx, gxx, hxy, rxy, }

for i=1:number
    %%
    % Extract the thickness and magnetic field.
    patternDirection = '(?<=Rx)\w(?=_)';
    patternTemperature = '(?<=Rx(x|y)_)\d*(?=K_)';
    patternThickness = '(?<=K_)\d*(?=nm_)';
    direction=regexp(fileName(i).name,patternDirection, 'match');

    directionForTag(i)=direction{1};
    temperature=regexp(fileName(i).name,patternTemperature, 'match');
    temperatureForCheck(i)=temperature{1};
    temperatureForTag(i)=[temperature{1},' K'];
    temperatureForNum(i)=str2double(temperature{1});
    
    thickness=regexp(fileName(i).name,patternThickness, 'match');
    thicknessForCheck(i)=thickness{1};
    thicknessForTag(i)=[thickness{1},' nm'];
    thicknessForNum(i)=str2double(thickness{1});
    
    isReadThi=~isempty(find([dataRaw{:,1}]==thicknessForNum(i)));
    if isReadThi==0
        j=i;
    elseif isReadThi==1
        j=find([dataRaw{:,1}]==thicknessForNum(i));
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
    dataRaw{j,1}=thicknessForNum(i);
    if direction=="x"
        dataRaw{j,2}=x;
        dataRaw{j,3}=abs(y);
        dataRaw{j,4}=asp*h/(e^2)./abs(y);
        sigma0(j)=max(dataRaw{j,4});
    elseif direction=="y"
        dataRaw{j,5}=x;
        dataRaw{j,6}=y;
    else
        error('ERROR 02:Something wrong with file name.' );
        return
    end
    %%
    clearvars x y m direction temperature j temR;
end
temperatureForTag(1)=[];
temperatureForTag=strip(temperatureForTag,'left','0');
temperatureForCheck(1)=[];
temperatureForCheck=convertStringsToChars(temperatureForCheck);

% thicknessForTag(1)=[];
thicknessForTag=strip(thicknessForTag,'left','0');
% thicknessForCheck(1)=[];
thicknessForCheck=convertStringsToChars(thicknessForCheck);

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

%% rHall vs H
figure
for i=1:number/2
   getPlot(dataRaw{i,5},dataRaw{i,6},generalView,"R_{xy} vs H",temperatureForTag(i))
   xlabel {H (T)}
   ylabel {R_{xy} (\Omega)}
   hold on
    
end
title('Hall resistance')
legend(thicknessForTag(1:number/2))

path=strcat("temPlots\","Hall resistance");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% Hall fitting
for i=1:number/2
    [xData, yData] = prepareCurveData(dataRaw{i,5},dataRaw{i,6});

    %% Remove symmetric component
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
    
    
    % Set up fittype and options.
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
    getFitPlot(fitresult,xData,yData,fittingView,strcat("Hall fitting result of ",thicknessForTag(i)),temperatureForTag(i),1,"R_{xy} (\Omega)");

    % Label axes
%     xlabel {H (T)}
%     ylabel {R_{xy} (\Omega)}
    
    clearvars xData yData ft fitresult gof opts ftHallLinear gofHallLinear fitresultHallLinear gofHallLinear
end

% Properties from Hall fitting.
ne=abs(1./e./hall);
neBulk=ne./([dataRaw{:,1}].*1e-9);
mu=(e^2/h).*sigma0./(e*ne(i));
taup=me./e.*mu;

if dimension == 2
    kf=(2.*pi.*ne).^(1/2);
elseif dimension == 3
    kf=(3.*pi.^2.*ne).^(1/3);    % Here ne should be in bulk
else
    error('ERROR 04');
    return
end

vf=hbar./me.*kf;
dif=vf.^(dimension).*taup./dimension;
be=hbar./(4.*e.*dif.*taup);

fprintf("Hall fitting finished.\nWaiting for WAL fitting.\n")
% dataRaw format: {thickness, hxx, rxx, gxx, hxy, rxy, }
%% Rxx vs B

figure
for i=1:number/2
    getPlot(dataRaw{i,2},dataRaw{i,3},generalView,"R_{xx} vs H",temperatureForTag(i))
    xlabel {H (T)}
    ylabel {R_{xx} (\Omega)}
    hold on
    
end
title('Longitudinal resistance')
legend(thicknessForTag(1:number/2),'Location','northeastoutside')

path=strcat("temPlots\","Longitudinal resistance");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% Gxx vs B
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
figure
for i=1:number/2
    dg=dataRaw{i,4}-max(dataRaw{i,4});
    getPlot(dataRaw{i,2},dg,generalView,"G_{xx} vs H",temperatureForTag(i))
    xlabel {H (T)}
    ylabel {\DeltaG_{xx} (e^2/h)}
    hold on
    
end
title('Longitudinal conductance')
legend(thicknessForTag(1:number/2),'Location','northeastoutside')
clearvars dg

path=strcat("temPlots\","Longitudinal conductance");
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
    getPlot(gxxNor{i,2},gxxNor{i,3},generalView,"Gxx(Normalized) vs B",temperatureForTag(i));
    xlabel {H (T)}
    ylabel {G_{xx}(Normalized)}
    hold on
    clearvars temY
end
title('Longitudinal conductivity (Normalized)')
legend(thicknessForTag(1:number/2),'Location','northeastoutside')

path=strcat("temPlots\","Longitudinal conductance (Normalized)");
saveas(gcf,path,'meta');
saveas(gcf,path,'jpeg');

%% WAL fitting
% Format of dataRaw: {Temperature, hxx, rxx, gxx, hxy, rxy, rxyCali}
dataPoint=cell(number/2,1);
fittingCurve=cell(number/2,1);
% parpool('local',coreNumber)
for i=1:number/2
    % Pretreatment of data.
    % Get the abstract of hxx.
    % xWal = abs(dataRaw{i,2});
     temFitresultWAL=cell(3,1);
    xWal = dataRaw{i,2};
    yWal = dataRaw{i,4}-sigma0(i);
    [xData, yData] = prepareCurveData( xWal, yWal );  
    
    %% Remove asymmetric component
    % Set up fittype and options.
    ftLinear = 'linearinterp';

    % Fit model to data.
    [fitresultLinear, gofLinear] = fit( xData, yData, ftLinear, 'Normalize', 'on' );
    
    % Substrate the linear component
    yData=(yData+fitresultLinear(-xData))./2;
    
    % Plot fit with data.
    % figure( 'Name', 'untitled fit 1' );
    % plot( fitresultLinear, xData, yData );    
    
    % Set up fittype and options.
    ft = fittype( '1/3.14159*(-psi((bso+be)/abs(x)+0.5)+log((bso+be)/abs(x))+1.5*psi((bi+4*bso/3)/abs(x)+0.5)-1.5*log((bi+4*bso/3)/abs(x))-0.5*psi(bi/abs(x)+0.5)+0.5*log(bi/abs(x)))+kFactor*x^2', 'independent', 'x', 'dependent', 'y' , 'problem' , 'be' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMinChange = 1e-16;
    opts.Display = 'Off'; 
    opts.Lower = [0 0 -Inf];
    opts.MaxFunEvals = 50000;
    opts.MaxIter = 50000;
    opts.Robust = 'LAR';
    opts.StartPoint = [0.002 10 -0.00001];
    opts.TolFun = 1e-08;
    opts.TolX = 1e-08;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts , 'problem' , be(i));
    
    kFactor(i)=fitresult.kFactor;
    bso(i)=fitresult.bso;
    bi(i)=fitresult.bi;
    temConfidence=confint(fitresult);
    confidenceBi(i)=abs((temConfidence(1,1)-temConfidence(2,1))/2);
    confidenceBso(i)=abs((temConfidence(1,2)-temConfidence(2,2))/2);
    
    % Plot fit with data.
    [dataPoint{i,1},fittingCurve{i,1}]=getFitPlot(fitresult,xData,yData,fittingView,strcat("WAL fitting result of ",thicknessForTag(i)),temperatureForTag(i),256,"\DeltaG (e^2/h)");
%     figure( 'Name', strcat("WAL fitting for ",temperatureForTag(i)) );
%     plot( fitresult, xData, yData ,'o');
%     legend( 'Data points', 'Fitting curve', 'Location', 'NorthEast' );
%     title(strcat("WAL fitting for ",temperatureForTag(i)))
    % Label axes
%     xlabel {H_{ex} (T)}
%     ylabel {\DeltaG (e^2/h)}
%     grid on
%     hold on
    
    %%%%%%%%%%%%%%%%%%%%
%     if i==1
%         return
%     end
%     
    clearvars xWal yWal xData yData fitresult ftLinear gof temConfidence temFitresultWAL opts fitresultLinear gofLinear
end
% delete(gcp('nocreate'))
figure(256)
legend([dataPoint{:,1},fittingCurve{1,1}],[thicknessForTag(1:number/2),"Fitting curve"],'Location','eastoutside')
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

i=1;    % Set for tag of temperature

%% Hall coefficients vs thickness
getScatter(thicknessForNum(1:number/2),hall,generalView,"Hall coefficient vs Thickness",temperatureForTag(i),{"Thickness (nm)","H_e (m^2/C)"});
getErrorBar(thicknessForNum(1:number/2),hall,confidenceHall,generalView,"Hall coefficient vs Thickness with error bar",{"Thickness (nm)","H_e (m^2/C)"})

%% ne vs thickness
getScatter(thicknessForNum(1:number/2),ne,generalView,"n_e vs Thickness",temperatureForTag(i),{"Thickness (nm)","n_e (cm^{-2})"})
% xlabel {Thickness (nm)}
% ylabel {n_e (cm^{-2})}

%% D vs thickness
getScatter(thicknessForNum(1:number/2),dif,generalView,"D vs thickness",temperatureForTag(i),{"Thickness (nm)","D (cm^2/s)"})
% xlabel {Thickness (nm)}
% ylabel {D (cm^2/s)}

%% L_SO vs thickness
getScatter(thicknessForNum(1:number/2),lso,generalView,"L_{SO} vs thickness",temperatureForTag(i),{"Thickness (nm)","L_{SO} (nm)"});
% xlabel {Thickness (nm)}
% ylabel {L_{SO} (m)}

%% L_Phi vs thickness
getScatter(thicknessForNum(1:number/2),lphi,generalView,"L_{\phi} vs thickness",temperatureForTag(i),{"Thickness (nm)","L_\phi (nm)"});
% xlabel {Thickness (nm)}
% ylabel {L_\phi (m)}

%% L_SO vs D
getScatter(dif,lso,generalView,"L_{SO} vs D",temperatureForTag(i),{"D (cm^2/s)","L_{SO} (nm)"});
% xlabel {D (cm^2/s)}
% ylabel {L_{SO} (m)}

%% Tau_SO vs Tau_p
getScatter(taup,tauso,generalView,"\tau_{SO} vs \tau_p",temperatureForTag(i),{"\tau_p (fs)","\tau_{SO} (fs)"});
% xlabel {\tau_p (fs)}
% ylabel {\tau_{SO} (fs)}

%% Tau_SO vs D
getScatter(dif,tauso,generalView,"\tau_{SO} vs D",temperatureForTag(i),{"D (cm^2/s)","\tau_{SO} (fs)"});
% xlabel {D (cm^2/s)}
% ylabel {\tau_{SO} (fs)}

%% mu vs thickness
getScatter(thicknessForNum(1:number/2),mu,generalView,"\mu vs thickness",temperatureForTag(i),{"Thickness (nm)","\mu (cm^2/V/s)"});
% xlabel {Thickness (nm)}
% ylabel {\mu (cm^2/V/s)}

%% Tau_SO vs thickness
getScatter(thicknessForNum(1:number/2),tauso,generalView,"\tau_{SO} vs thickness",temperatureForTag(i),{"Thickness (nm)","\tau_{SO} (fs)"});
% xlabel {Thickness (nm)}
% ylabel {\tau_{SO} (fs)}

%% Tau_phi vs thickness
getScatter(thicknessForNum(1:number/2),tauphi,generalView,"\tau_{\phi} vs thickness",temperatureForTag(i),{"Thickness (nm)","\tau_{\phi} (fs)"});
% xlabel {Thickness (nm)}
% ylabel {\tau_{\phi} (fs)}

%% Tau_p vs thickness
getScatter(thicknessForNum(1:number/2),taup,generalView,"\tau_p vs thickness",temperatureForTag(i),{"Thickness (nm)","\tau_p (fs)"});
% xlabel {Thickness (nm)}
% ylabel {\tau_p (fs)}

%% vf vs thickness
getScatter(thicknessForNum(1:number/2),vf,generalView,"v_f vs thickness",temperatureForTag(i),{"Thickness (nm)","v_f "});
% xlabel {Thickness (nm)}
% ylabel {v_f }

%% gCritical vs thickness
getScatter(thicknessForNum(1:number/2),gCritical,generalView,"g_{Cri} vs thickness",temperatureForTag(i),{"Thickness (nm)","g_{Cri} (a.u.)"});
% xlabel {Thickness (nm)}
% ylabel {v_f }

%% Bso vs thickness
getScatter(thicknessForNum(1:number/2),bso,generalView,"B_{SO} vs thickness",temperatureForTag(i),{"Thickness (nm)","B_{SO} (T)"});
% xlabel {Thickness (nm)}
% ylabel {B_{SO} (T)}
getErrorBar(thicknessForNum(1:number/2),bso,confidenceBso,generalView,"B_{SO} vs Thickness with error bar",{"Thickness (nm)","B_{SO} (T)"})

%% Bi vs thickness
getScatter(thicknessForNum(1:number/2),bi,generalView,"B_i vs thickness",temperatureForTag(i),{"Thickness (nm)","B_i (T)"});

getErrorBar(thicknessForNum(1:number/2),bi,confidenceBi,generalView,"B_i vs Thickness with error bar",{"Thickness (nm)","B_i (T)"})

clearvars i path
fprintf("Program completed.\n")
toc
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
function getScatter(x,y,generalView,titleOfPlot,temperature,label)

if nargin < 3
    generalView=true;
end

fig=figure;

if generalView==false
    set(fig,'Visible','off')
end

plot(x,y,'o-','LineWidth',2.5)
set(gcf,'position',[800 100 800 600]);
xlim([1.1*min(x)-0.1*max(x) 1.1*max(x)-0.1*min(x)]);
ylim([1.1*min(y)-0.1*max(y) 1.1*max(y)-0.1*min(y)]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)

if nargin >= 4
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
end

if nargin >= 5
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
    str=strcat("T=",temperature);
    text(0.95,0.05,str,'Color','blue','FontSize',24,'Units','normalized','HorizontalAlignment','right')
end

if nargin >= 6
    xlabel(label{1})
    ylabel(label{2})
end

grid on

titleOfPlot=strrep(titleOfPlot,'\','');
path=strcat("temPlots\",titleOfPlot);
saveas(fig,path,'meta');
saveas(fig,path,'jpeg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataPoint,fittingCurve] = getFitPlot(fitresult,xData,yData,fittingView,titleOfPlot,temperature,n,yLabel)

if nargin < 4
    fittingView = true;
end
if n==256
fig=figure(n);
else
    fig=figure
end
set(gcf,'position',[800 100 800 600]);

if fittingView==false
    set(fig,'Visible','off')
end

dataPoint=plot(xData, yData ,'o','LineWidth',2);
hold on
fittingCurve=plot(fitresult,'r:');
fittingCurve.LineWidth=4;
xlabel {H_{ex} (T)}
if nargin >= 7
    ylabel(yLabel)
end
legend({"Data Point","Fitting Curve"});
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)

if nargin >= 5
    set(fig,'Name',titleOfPlot);
    title(titleOfPlot);
end

if nargin >= 6
    if n == 256
        set(fig,'Name',"WAL fitting result");
        title("WAL fitting result");
        hold on
    else
        set(fig,'Name',titleOfPlot);
        title(titleOfPlot);
        str=strcat("T=",temperature);
        text(0.95,0.05,str,'Color','blue','FontSize',24,'Units','normalized','HorizontalAlignment','right')
    end
end

grid on
% hold on
titleOfPlot=strrep(titleOfPlot,'\','');
path=strcat("temFitPlots\",titleOfPlot);
saveas(fig,path,'meta');
saveas(fig,path,'jpeg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getPlot(x,y,generalView,titleOfPlot,temperature)

if nargin < 3
    generalView = true;
end

if generalView == false
    set(fig,'Visible','off')
end

plot(x,y,'Linewidth',2);
set(gcf,'position',[800 100 800 600]);
set(gca, 'linewidth', 1.1,'fontname', 'Helvetica', 'FontSize',18)

if nargin == 4
    set(gcf,'Name',titleOfPlot);
    title(titleOfPlot);
end

if nargin == 5
    set(gcf,'Name',titleOfPlot);
    title(titleOfPlot);
    str=strcat("T=",temperature);
    text(0.95,0.05,str,'Color','blue','FontSize',24,'Units','normalized','HorizontalAlignment','right')
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