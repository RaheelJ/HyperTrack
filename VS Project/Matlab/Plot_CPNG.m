clear
clc
close all

currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
addpath( fullfile( pathstr, '..' ) );

File0_Data = csvread('Output/CPNG_States_0.csv', 1);
Time0 = File0_Data(:, 1);
Radial_Distance0 = File0_Data(:, 4);
Longitude0 = File0_Data(:, 3);
Latitude0 = File0_Data(:, 2);

File1_Data = csvread('Output/CPNG_States_1.csv', 1);
Time1 = File1_Data(:, 1);
Radial_Distance1 = File1_Data(:, 4);
Longitude1 = File1_Data(:, 3);
Latitude1 = File1_Data(:, 2);

File2_Data = csvread('Output/CPNG_States_2.csv', 1);
Time2 = File2_Data(:, 1);
Radial_Distance2 = File2_Data(:, 4);
Longitude2 = File2_Data(:, 3);
Latitude2 = File2_Data(:, 2);

File3_Data = csvread('Output/CPNG_States_3.csv', 1);
Time3 = File3_Data(:, 1);
Radial_Distance3 = File3_Data(:, 4);
Longitude3 = File3_Data(:, 3);
Latitude3 = File3_Data(:, 2);

File4_Data = csvread('Output/CPNG_States_4.csv', 1);
Time4 = File4_Data(:, 1);
Radial_Distance4 = File4_Data(:, 4);
Longitude4 = File4_Data(:, 3);
Latitude4 = File4_Data(:, 2);

File5_Data = csvread('Output/CPNG_States_5.csv', 1);
Time5 = File5_Data(:, 1);
Radial_Distance5 = File5_Data(:, 4);
Longitude5 = File5_Data(:, 3);
Latitude5 = File5_Data(:, 2);

File6_Data = csvread('Output/CPNG_States_6.csv', 1);
Time6 = File6_Data(:, 1);
Radial_Distance6 = File6_Data(:, 4);
Longitude6 = File6_Data(:, 3);
Latitude6 = File6_Data(:, 2);

File7_Data = csvread('Output/CPNG_States_7.csv', 1);
Time7 = File7_Data(:, 1);
Radial_Distance7 = File7_Data(:, 4);
Longitude7 = File7_Data(:, 3);
Latitude7 = File7_Data(:, 2);

File8_Data = csvread('Output/CPNG_States_8.csv', 1);
Time8 = File8_Data(:, 1);
Radial_Distance8 = File8_Data(:, 4);
Longitude8 = File8_Data(:, 3);
Latitude8 = File8_Data(:, 2);
% 
% File9_Data = csvread('Output/CPNG_States_9.csv', 1);
% Time9 = File9_Data(:, 1);
% Radial_Distance9 = File9_Data(:, 4);
% Longitude9 = File9_Data(:, 3);
% Latitude9 = File9_Data(:, 2);

% File10_Data = readmatrix('Input/Target_1.csv');
% T_Latitude = File10_Data(ceil(1:max(Time0)), 3);
% T_Longitude = File10_Data(ceil(1:max(Time0)), 4);








%geoplot(Latitude1, Longitude1, 'r', T_Latitude, T_Longitude, 'b')
geoplot(Latitude1, Longitude1, 'r')
hold on
geoplot(Latitude0, Longitude0, 'g')
hold on
geoplot(Latitude2, Longitude2, 'k')
hold on
geoplot(Latitude3, Longitude3, 'r.')
hold on
geoplot(Latitude4, Longitude4, 'g.')
hold on
geoplot(Latitude5, Longitude5, 'k.')
hold on
geoplot(Latitude6, Longitude6, 'r*')
hold on
geoplot(Latitude7, Longitude7, 'g*')
hold on
geoplot(Latitude8, Longitude8, 'k*')
