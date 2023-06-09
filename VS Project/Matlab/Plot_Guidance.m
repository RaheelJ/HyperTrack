clear
clc
close all

currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
addpath( fullfile( pathstr, '..' ) );

File1_Data = csvread('Output/Missile_States.csv', 1);
Time = File1_Data(:, 1);
Radial_Distance = File1_Data(:, 4);
Longitude = File1_Data(:, 3);
Latitude = File1_Data(:, 2);

File3_Data = csvread('Input/Target_1.csv', 1);
T_Latitude = File3_Data(ceil(1:max(Time)), 3);
T_Longitude = File3_Data(ceil(1:max(Time)), 4);
geoplot(Latitude, Longitude)

