%% GPS_plotting.m 
% Astrodynamics Toolbox
%
% GPS constellation visualizations.
%
% Author: Josh Geiser
% Last Update: 2022-03-14



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to entire directory
addpath(genpath('..'));

% load plot parameters
pp = PLOT_PARAMETERS;



%% LOAD DATA

% GPS constellation object
load('gps_long.mat');

% chief simulation data
simdata = load_simdata('chief_simdata_long.mat');



%% POSITIONS

% ECEF positions of GPS satellites and chief spacecraft
figure(1);
background('Milky Way');
hold on;
planet3D('Earth Cloudy')
for ii = 1:32
    plot3(gps.ECEF_position(ii,:,1),gps.ECEF_position(ii,:,2),...
        gps.ECEF_position(ii,:,3),'Color',pp.matlab_light_blue,...
        'LineWidth',1);
end
plot3(simdata.r_ecef(1,:),simdata.r_ecef(2,:),simdata.r_ecef(3,:),'r',...
    'LineWidth',3);
hold off;
view(-76.5,2.5);
drawnow;

% ECI positions of GPS satellites and chief spacecraft
figure(2);
background('Milky Way');
hold on;
planet3D('Earth Cloudy')
for ii = 1:32
    plot3(gps.ECI_position(ii,:,1),gps.ECI_position(ii,:,2),...
        gps.ECI_position(ii,:,3),'Color',pp.matlab_light_blue,...
        'LineWidth',1);
end
plot3(simdata.r_eci(1,:),simdata.r_eci(2,:),simdata.r_eci(3,:),'r',...
    'LineWidth',3);
hold off;
view(-76.5,2.5);
drawnow;



%% ANIMATION

% initializes figure
f = figure('Position',pp.Position);

% creates plot for every point in time
for k = 1:20:length(simdata.t)
    
    % clears figure
    clf(f);
    hold on;

    % animated line for satellite orbit
    h = animatedline('Color','r','LineWidth',2);
    
    % animated line for GPS satellites
    hs1 = arrayfun(@(x)animatedline('MaximumNumPoints',1,'LineWidth',1,...
        'Marker','o','MarkerFaceColor',[1,1,1]),1:32);
    
    % animated line for communication with closest 8 GPS satellites
    hs2 = arrayfun(@(x)animatedline('MaximumNumPoints',2,'Color','g'),1:8);

    % background and grid formatting
    background('black');
    grid on;
    ax = gca;
    ax.GridColor = [1,1,1];
    ax.GridAlpha = 0.25;

    % rotates Earth
    opts.RotAngle = rad2deg(gmst_linear(simdata.MJD_UT1(k)));
    planet3D('Earth Cloudy',opts);
    
    % counter used to create the 8 lines for the 8 closest GPS satellites
    j = 1;
    
    % SVIDs of closest four satellites
    SVID = gps.get_closest_SVIDs(simdata.r_ecef(:,k),simdata.MJD_GPS(k),8);

    % loops over all 32 GPS satellites
    for i = 1:32

        % red GPS satellite if one of closest four, otherwise white
        if ismember(i,SVID)
            hs1(i).Color = 'r';
            hs1(i).MarkerFaceColor = 'r';
        else
            hs1(i).Color = 'w';
            hs1(i).MarkerFaceColor = 'w';
        end
        
        % adds line between our satellite ith GPS satellite if it is one of
        % the closest 4 satellites
        if ismember(i,SVID)
            addpoints(hs2(j),[gps.ECI_position(i,k,1),...
                simdata.r_eci(1,k)],[gps.ECI_position(i,k,2),...
                simdata.r_eci(2,k)],[gps.ECI_position(i,k,3),...
                simdata.r_eci(3,k)]);
            j = j+1;
        end

        % plot ith GPS satellite
        addpoints(hs1(i),gps.ECI_position(i,k,1),...
            gps.ECI_position(i,k,2),gps.ECI_position(i,k,3));
    end

    % traces our satellite's ECI orbit
    plot3(simdata.r_eci(1,1:k),simdata.r_eci(2,1:k),simdata.r_eci(3,...
        1:k),'LineWidth',2,'Color','r');

    % plots our satellite
    plot3(simdata.r_eci(1,k),simdata.r_eci(2,k),simdata.r_eci(3,k),...
        'o','LineWidth',1,'Color','r','MarkerFaceColor','y');

    % axis labels and title
    xlabel('$I\;[\mathrm{m}]$','Interpreter','latex','FontSize',...
        pp.FontSize_axis);
    ylabel('$J\;[\mathrm{m}]$','Interpreter','latex','FontSize',...
        pp.FontSize_axis);
    zlabel('$K\;[\mathrm{m}]$','Interpreter','latex','FontSize',...
        pp.FontSize_axis);
    title('\textbf{ECI Positions of Spacecraft and GPS Satellites}',...
        'Interpreter','latex','FontSize',pp.FontSize_title);

    % updates plot (limits to 20 fps)
    drawnow limitrate
    
    % view point
    view(-45,10);
    
    % axis limits
    xlim([-25000000,25000000]);
    ylim([-25000000,25000000]);
    zlim([-25000000,25000000]);
    
    % captures the plot as an image 
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % writes to GIF
    if k == 1
        imwrite(imind,cm,'gps_measurement.gif','gif','Loopcount',inf); 
    else
        imwrite(imind,cm,'gps_measurement.gif','gif','WriteMode','append'); 
    end
    
end