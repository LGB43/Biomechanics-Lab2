# Biomechanics-Lab2

%BIOENG 1631 
%Lab 2 
%Finley Siegel and Nick Koontz

clc; clear; close all;

%%  User / file inputs 
markerfile_t1 = "Trial01_Markers.txt";
markerFiles = ["Trial02_Markers.txt" "Trial03_Markers.txt"];
forceFiles  = ["Trial02_Forces.txt" "Trial03_Forces.txt"];

%%  Constants / Parameters 
sampleFrequency = 120; % Hz
dt = 1 / sampleFrequency;
mass = 67.8;      % kg
height = 174.5;   % cm 
footMassFraction = 0.0145;
g = 9.81;
m_foot = mass * footMassFraction;

%%  Load baseline data (for avg foot length) 
markerdata_t1 = readmatrix(markerfile_t1);
avgFootLength = mean(vecnorm(markerdata_t1(:,17:19) - markerdata_t1(:,11:13), 2, 2));

%%  Loop over trials 
nTrials = numel(markerFiles);
results = struct('markerfile', [], 'forcefile', [], 'time_rel', [], 'ankleAngle', [], 'ankleAcceleration', [], 'ankleNetMoment', []);
results = repmat(results, nTrials, 1);

for i = 1:nTrials
    %%  Load trial files 
    markerfile = markerFiles(i);
    forcefile  = forceFiles(i);
    markerdata = readmatrix(markerfile);
    forcedata  = readmatrix(forcefile);
    
    %%  Detect stance intervals 
    time_all = forcedata(:,1);
    Fz_all = forcedata(:,10);
    contactLogical = Fz_all > 0;
    
    % find start/end indices of each contact episode
    d = diff([0; contactLogical; 0]);
    stanceStart_all = find(d==1);
    stanceEnd_all   = find(d==-1) - 1;
    
    % If multiple episodes exist, use the first one
    stanceStart = stanceStart_all(1);
    stanceEnd   = stanceEnd_all(1);
    
    % Extract stance arrays aligned to the same indices for marker/force data
    stance_forces = forcedata(stanceStart:stanceEnd, :);
    stance_time = stance_forces(:,1);
    
    % Copy columns to named variables
    COPY = stance_forces(:,3); % mm
    ForceY = stance_forces(:,9); % N
    ForceZ = stance_forces(:,10); % N
    MomentX = stance_forces(:,11); % Nmm
    
    % Get ankle marker positions from markerdata for the stance frames
    ankle_xyz = markerdata(stanceStart:stanceEnd, 11:13); % mm
    ankleY = ankle_xyz(:,2); % mm
    ankleZ = ankle_xyz(:,3); % mm
    
    %%  Compute ankle direction & angle 
    foot_vector = markerdata(stanceStart:stanceEnd,17:19) - markerdata(stanceStart:stanceEnd,11:13);
    foot_dir = foot_vector ./ norm(foot_vector);
    
    % compute ankle angle in sagittal plane
    nFrames = size(foot_dir,1);
    ankleAngle = nan(nFrames,1);
    for k = 1:nFrames
        v = foot_dir(k, [2 3]); % [0, y, z]
        ankleAngle(k) = atan2d(v(1), v(2)) - 90; % allows for negative angles below horizontal
    end
    
    %%  Angular acceleration 
    ankleAcceleration = nan(nFrames,1);
    for k = 2:(nFrames-1)
        ankleAcceleration(k) = (ankleAngle(k+1) - 2*ankleAngle(k) + ankleAngle(k-1)) / dt^2;
    end
    
    %%  External moment from GRF about ankle 
    r_y = COPY - ankleY; % mm
    r_z = 0 - ankleZ; % mm
    extMoment_GRF = r_y .* ForceZ - r_z .* ForceY; % Nmm
    
    %%  External moment from foot weight about ankle 
    delta_y = 0.5 * avgFootLength .* cosd(ankleAngle); % mm
    Fz_foot = - m_foot * g; % N
    extMoment_foot = delta_y .* Fz_foot; % Nmm
    
    %%  Rotational inertia 
    footLength_m = avgFootLength / 1000; % m
    radius_of_gyration = 0.69 * footLength_m; % m
    I_foot = m_foot * radius_of_gyration^2; % kg·m^2
    alpha_rad_s2 = deg2rad(ankleAcceleration); % rad/s^2
    inertia_Nm = I_foot .* alpha_rad_s2; % Nm
    
    %%  Net muscle moment 
    extMoment_total = extMoment_GRF + extMoment_foot; % Nmm
    ankleNetMoment = extMoment_total/1000 - inertia_Nm; % Nm
    
    %%  Store results 
    t = stance_time - stance_time(1);
    results(i).markerfile = markerfile;
    results(i).forcefile  = forcefile;
    results(i).time_rel   = t;
    results(i).ankleAngle = ankleAngle;
    results(i).ankleAcceleration = ankleAcceleration;
    results(i).ankleNetMoment = ankleNetMoment;
end

%%  Final Plot 
figure('Name', 'Trial Data', 'Position',[150 150 900 450]);

ax(1) = subplot(2, 2, 1);
plot(results(1).time_rel/sampleFrequency, results(1).ankleAngle, 'LineWidth', 1.2);
ylabel('Ankle angle (deg)'); xlabel('Time (s)');
grid on;
title('Ankle Angle - Trial 2 - Walking');

ax(2) = subplot(2, 2, 2);
plot(results(2).time_rel/sampleFrequency, results(2).ankleAngle, 'LineWidth', 1.2);
ylabel('Ankle angle (deg)'); xlabel('Time (s)');
grid on;
title('Ankle Angle - Trial 3 - Running');

ax(3) = subplot(2, 2, 3);
plot(results(1).time_rel/sampleFrequency, results(1).ankleNetMoment, 'LineWidth', 1.2);
ylabel('Net Muscle Moment (Nm)'); xlabel('Time (s)');
grid on;
title('Net Muscle Moment - Trial 2 - Walking');

ax(4) = subplot(2, 2, 4);
plot(results(2).time_rel/sampleFrequency, results(2).ankleNetMoment, 'LineWidth', 1.2);
ylabel('Net Muscle Moment (Nm)'); xlabel('Time (s)');
grid on;
title('Net Muscle Moment - Trial 3 - Running');
linkaxes(ax([1,2]),'y');
linkaxes(ax([3,4]),'y');
linkaxes(ax([1,3]),'x');
linkaxes(ax([4,2]),'x');
saveas(gcf,'Trial Data.png')

%%  Results Table 
meanAnkleAngle = [mean(results(1).ankleAngle,'omitnan') mean(results(2).ankleAngle,'omitnan')]
maxAnkleAngle = [max(results(1).ankleAngle) max(results(2).ankleAngle)]
minAnkleAngle = [min(results(1).ankleAngle) min(results(2).ankleAngle)]

meanAnkleMoment = [mean(results(1).ankleNetMoment,'omitnan') mean(results(2).ankleNetMoment,'omitnan')]
maxAnkleMoment = [max(results(1).ankleNetMoment) max(results(2).ankleNetMoment)]
minAnkleMoment = [min(results(1).ankleNetMoment) min(results(2).ankleNetMoment)]
