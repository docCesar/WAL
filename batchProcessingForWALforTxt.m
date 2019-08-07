clear;close all;clc

% Next update
% 0. (IMPORTANT!!!) TO solve the problem that the number of lines is
% different.
% 1. Add 3D samples cases.
% 2. Seperate Hall data and WAL data automatically.

% Constants
h=6.6260699*10^-34;
hbar=1.0545718e-34;
m=9.10938356*10^(-31);
e=1.60217662*10^-19;

% Parameters of samples
% asp=4;          % Old mask
% asp=2.5;        % Small devices of new mask
asp=24;         % Large devices of new mask
thick=5*10^-9; % Thickness of sample layer

%% Import data
fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.
number=length(fileName);
% Claim variables.
temForCheck=[];
hall=[];
bH=[];
rH=[];
gHall=[];
bWAL=[];
rWAL=[];
gWAL=[];
mu=[];
taup=[];
vf=[];
dif=[];
be=[];
ne=[];
kf=[];

bi=[];
bso=[];

for i=1:number
    name=convertCharsToStrings(fileName(i).name);
    % Get the temperature for check.
    pattern='(?<=_)\w*(?=K_Vg)';
    temReg=regexp(name,pattern,'match');
    temForCheck(end+1)=str2num(temReg{1});
    %% Initialization

    delimiter = '\t';
    startRow = 2;
    
    %% Format
    formatSpec = '%f%f%*s%[^\n\r]';
    
    %% Open the data.
    fileID = fopen(name,'r');

    %% Read the data column.

    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    [l,c]=size(dataArray{1});
    
	%% Close the data.
    fclose(fileID);
    
    %% Name the data column.

    
    if contains(name,'Rxy_')==1    % 稍后应添加一个匹配判断
        bH(:,end+1) = dataArray{:, 1}./1000;
        rH(:,end+1) = dataArray{:, 2};
        gHall(:,end+1)=asp*h/(e^2)./dataArray{:, 2};
        bWAL(1:l,end+1)=0;
        rWAL(1:l,end+1)=0;
        gWAL(1:l,end+1)=0;
    elseif contains(name,'Rxx_')==1
        bWAL(:,end+1) = dataArray{:, 1}./1000;
        rWAL(:,end+1) = dataArray{:, 2};
        gWAL(:,end+1)=asp*h/(e^2)./dataArray{:, 2};
        bH(1:l,end+1)=0;
        rH(1:l,end+1)=0;
        gHall(1:l,end+1)=0;
    else
        "Error code 01 (Wrong data name)"
        return
    end

    %% Clean the temporary variation.
    clearvars temReg pattern filename delimiter startRow formatSpec fileID dataArray ans;
end


%% rHall vs B
figure
for i=number/2+1:number
   
   plot(bH(:,i),rH(:,i),'Linewidth',2)
   xlabel {B (T)}
   ylabel {R_{xy} (\Omega)}
   grid on
   hold on
    
end
temForLegend=transpose(temForCheck(1,1:number/2));
temForLegend=num2str(temForLegend);
title('Hall resistance')
legend(temForLegend)

%% Hall calibration
tem=rH(:,number/2+1:end);
for i=1:number/2
    temR=rH(:,i+number/2);
    temR=temR';
    temR=fliplr(temR);
    temR=temR';
    tem(:,end+1)=temR;
end
rHN=(rH-tem)./2;


%% Hall fitting



for i=1:number/2
%     [xData, yData] = prepareCurveData( bH(:,i+number/2), rH(:,i+number/2) );
    [xData, yData] = prepareCurveData( bH(:,i+number/2), rHN(:,i+number/2) );

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
    hall(end+1)=abs(fitresult.p1)
%     'hall'
%     fitresult.p1
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    plotH = plot( fitresult, xData, yData );
    legend( plotH, 'rH vs. bH', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel bH
    ylabel rH
    grid on

% Properties from Hall fitting.
ne(i)=abs(1/e/hall(i));
mu(i)=e^2/h*max(gWAL(:,i))/e/ne(i);
taup(i)=m/e*mu(i);
kf(i)=(2*pi*ne(i))^(1/2);  % 2D case.
%kf(i)=(3*pi^2*ne)^(1/3);  % 3D case.
vf(i)=hbar/m*kf(i);
dif(i)=vf(i)^2*taup(i)/2;
be(i)=hbar/(4*e*dif(i)*taup(i));

end
'Hall finish'
% return

%% Rxx vs B
figure
for i=1:number/2
   
   plot(bWAL(:,i),rWAL(:,i),'Linewidth',2)
   xlabel {B (T)}
   ylabel {R_{xx} (\Omega)}
   grid on
   hold on
    
end
title('Longitude resistance')
legend(temForLegend)

%% Rxx(Normalized) vs B
figure
for i=1:number/2
   temy=(rWAL(:,i)-min(rWAL(:,i)))./min(rWAL(:,i));
   plot(bWAL(:,i),temy,'Linewidth',2)
   xlabel {B (T)}
   ylabel {R_{xx}(Normalized) (\Omega)}
   grid on
   hold on
   clearvars temy
end
title('Longitude resistance(Normalized)')
legend(temForLegend,'Location','SouthEast')


%% Gxx(Normalized) vs B
figure
for i=1:number/2
   temy=(gWAL(:,i)-min(gWAL(:,i)))./min(gWAL(:,i));
   plot(bWAL(:,i),temy,'Linewidth',2)
   xlabel {B (T)}
   ylabel {G_{xx}(Normalized)}
   grid on
   hold on
   clearvars temy
end
title('Longitude resistance(Normalized)')
legend(temForLegend,'Location','NorthEast')
%% WAL fitting
% Pretreatment of data.
bWAL=abs(bWAL);     % Get the abstract of bWAL.

% [bWALmin,bWALpos]=min(bWAL)
% bWAL=bWAL-bWALmin;
% bWAL(1800,:)=[];
% gWAL(1800,:)=[];



for i=1:number/2
%     gWAL=asp.*h./(e^2)./rWAL(:,i);
    dgWAL=gWAL(:,i)-max(gWAL(:,i));
    [xData, yData] = prepareCurveData( bWAL(:,i), dgWAL );
    

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
    i
    [fitresult, gof] = fit( xData, yData, ft, opts , 'problem' , be(i))
    
    bso(end+1)=fitresult.bso;
    bi(end+1)=fitresult.bi;
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    scatter( xData, yData,'o' );
    box on
    hold on
    plot( fitresult, xData, yData );
    legend( 'gWALhalf vs. bWALhalf', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel bWALhalf
    ylabel gWALhalf
    grid on
    
    
    
end

lso=sqrt(hbar./(4.*e.*bso));
lphi=sqrt(hbar./(4.*e.*bi));
ltr=(hbar/e*sqrt(2*pi)).*mu;
tauso=hbar./(4*e.*bso.*dif);

%% Change the units
dif=dif.*10000;
mu=mu.*10000;


%% Ne vs temperature
figure
scatter(temForCheck(1:number/2),ne,'o');
xlim([0 1.1*max(temForCheck)])
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
scatter(temForCheck(1:number/2),dif,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {D (cm^2/s)}
grid on
box on

%% L_SO vs temperature
figure
scatter(temForCheck(1:number/2),lso,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {L_{SO} (m)}
grid on
box on

%% L_Phi vs temperature
figure
scatter(temForCheck(1:number/2),lphi,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {L(\phi) (m)}
grid on
box on

%% L_SO vs D
figure
scatter(dif,lso,'o');
% xlim([0 1.1*max(temForCheck)])
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
scatter(temForCheck(1:number/2),mu,'o');
xlim([0 1.1*max(temForCheck)])
xlabel {T (K)}
ylabel {\mu (cm^2/V/s)}
grid on
box on

%% Tau_SO vs temperature
figure
scatter(temForCheck(1:number/2),tauso,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {\tau_{SO} (s)}
grid on
box on

%% Tau_p vs temperature
figure
scatter(temForCheck(1:number/2),taup,'o');
xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {T (K)}
ylabel {\tau_p (s)}
grid on
box on

%% vf vs temperature
figure
scatter(temForCheck(1:number/2),vf,'o');
xlim([0 1.1*max(temForCheck)])
xlabel {T (K)}
ylabel {v_f }
grid on
box on

%% Bso vs temperature
figure
scatter(temForCheck(1:number/2),bso,'o');
xlim([0 1.1*max(temForCheck)])
xlabel {T (K)}
ylabel {B_{SO} (T)}
grid on
box on


'Finished'