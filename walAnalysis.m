clear;close all;clc

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

%% Import data name
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
directionForTag=string;
temperatureForNum=string;
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
    
    temperatureForNum(end+1)=str2num(temperature{1});
    temperatureForTag(end+1)=temperature{1};
    
    [x,y] = importfile(fileName(i).name);
    
    
    if direction=="x"
        eval(['hxx_',temperature{1},'K=x;'])
        eval(['rxx_',temperature{1},'K=y;'])
    elseif direction=="y"
        eval(['hxy_',temperature{1},'K=x;'])
        eval(['rxy_',temperature{1},'K=y;'])
    else
        "Something wrong with file name"
    end
    
    
    
    

    
%     clearvars y isNegative magneticField thickness
end



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
