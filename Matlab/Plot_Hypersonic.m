clear
clc
close all

currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
addpath( fullfile( pathstr, '..' ) );

File1_Data = csvread('Output/Hypersonic_States.csv', 1);
Time = File1_Data(:, 1);
Heading_Angle = File1_Data(:, 2);
Flight_Path_Angle = File1_Data(:, 3);
Radial_Distance = File1_Data(:, 4);
Longitude = File1_Data(:, 5);
Latitude = File1_Data(:, 6);
Lift = File1_Data(:, 7);
Bank_Angle = File1_Data(:, 8);
Line_of_Sight_Distance = File1_Data(:, 9);
Altitude_Error = File1_Data(:, 10);
Delta_V = File1_Data(:, 11);
Q_Lat = File1_Data(:, 12);
Q_Long = File1_Data(:, 13);
n_lat = File1_Data(:, 14);
n_long = File1_Data(:, 15);

% File2_Data = csvread('Output/Hypersonic_Trajectory.csv', 1);
% T_Meas = File2_Data(:, 1); 
% T_X = File2_Data(:, 2);
% T_Y = File2_Data(:, 3);
% T_Z = File2_Data(:, 4);
% I_X = File2_Data(:, 5);
% I_Y = File2_Data(:, 6);
% I_Z = File2_Data(:, 7);

File3_Data = csvread('Input/Target_2.csv', 1);
T_Latitude = File3_Data(:, 3);
T_Longitude = File3_Data(:, 4);
geoplot(Latitude, Longitude, 'r', T_Latitude, T_Longitude, 'b')

figure
plot(Time, Heading_Angle, 'b', Time, Flight_Path_Angle, 'r', Time, Q_Lat, 'b--', Time, Q_Long, 'r--')
box('on');
grid('on');
hold('all');
% Create labels
ylabel('Angle (deg)');
xlabel('Time (s)');
title('Flight Dynamics');
% Create legend
legend('Heading Angle', 'Flight Path Angle', 'Lateral LOS Angle', 'Longitudnal LOS Angle');

figure
plot(Time, n_lat, 'b', Time, n_long, 'r')
box('on');
grid('on');
hold('all');
% Create labels
ylabel('Load (g)');
xlabel('Time (s)');
title('Commands of Reaction Jet');
% Create legend
legend('Lateral Load', 'Longitudnal Load');

% figure
% plot3(I_X./1000, I_Y./1000, I_Z./1000, 'b', T_X./1000, T_Y./1000, T_Z./1000, 'r')
% box('on');
% grid('on');
% hold('all');
% % Create labels
% ylabel('North Distance (km)');
% xlabel('East Distance (km)');
% zlabel('Altitude (km)');
% title('Intercept Trajectory');
% % Create legend
% legend('Interceptor', 'Target')
% 
% figure
% plot(T_Meas, I_X, 'b', T_Meas, T_X, 'r')
% box('on');
% grid('on');
% hold('all');
% % Create labels
% ylabel('Distance (m)');
% xlabel('Time (s)');
% title('Lateral East (X) Trajectories');
% % Create legend
% legend('Interceptor', 'Target')
% 
% figure
% plot(T_Meas, I_Y, 'b', T_Meas, T_Y, 'r')
% box('on');
% grid('on');
% hold('all');
% % Create labels
% ylabel('Distance (m)');
% xlabel('Time (s)');
% title('Lateral North (Y) Trajectories');
% % Create legend
% legend('Interceptor', 'Target')
% 
% figure
% plot(T_Meas, I_Z, 'b', T_Meas, T_Z, 'r')
% box('on');
% grid('on');
% hold('all');
% % Create labels
% ylabel('Distance (m)');
% xlabel('Time (s)');
% title('Altitudes of Interceptor and Target');
% % Create legend
% legend('Interceptor', 'Target')


figure
plot(Time, Bank_Angle, 'b', Time, Lift, 'r')
box('on');
grid('on');
hold('all');
% Create labels
ylabel('Angle (Degree) / Lift (m/s^2)');
xlabel('Time (s)');
title('Commands of Interceptor');
% Create legend
legend('Bank Angle', 'Lift');

figure
plot(Time, Line_of_Sight_Distance)
box('on');
grid('on');
hold('all');
% Create labels
ylabel('Distance (m)');
xlabel('Time (s)');
title('Relative Position in Lateral Plane');
% Create legend
legend('LOS Distance');

figure
plot(Time, Altitude_Error)
box('on');
grid('on');
hold('all');
% Create labels
ylabel('Distance (m)');
xlabel('Time (s)');
title('Altitude Error');
% Create legend
legend('Altitude Error');

figure
plot(Time, Delta_V)
box('on');
grid('on');
hold('all');
% Create labels
ylabel('Speed (m/s)');
xlabel('Time (s)');
title('Relative Speed');
% Create legend
legend('Delta V');


