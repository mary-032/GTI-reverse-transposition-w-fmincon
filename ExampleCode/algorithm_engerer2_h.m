%% Example script of reverse transposition algorithm (RUN the entire script)
% ATTENTION: the following code expects a good quality of input data, please perform suitable 
% data cleaning and quality control of your input data.
% ------INPUTS------ %
%     GTI = global tilted irradiance [W/m2] (array)
%     Eexth = extraterrestrial irradiance to an horizontal plane [W/m2] (array)
%     Eextn = extraterrestrial irradiance to a normal plane [W/m2] (array)
%     Incidence_angle = solar incidence angle [deg] (array)
%     Zenith_angle = solar zenith angle [deg] (array)
%     Rb = cos(Incidence_angle)/cos(Zenith_angle) (array)
%     Rr = (1-cos(Tilt))/2; (value)
%     Albedo = surface albedo [0-1] (array)
%     ast = apparent solar time [0-24h] (array)
%     Gcs = clear-sky global horizontal irradiance [W/m2] (array)
%     Solar_azimuth = solar azimuth angle [deg] (array) - Convention: South = 0
%       deg, North = -180 deg
%     mode = "1-min" or "1-h", temporal resolution (string)
%     Tilt = surface tilt angle [deg] (value)
%     Surface_Azimuth = surface azimuth angle [deg] (value) - Convention: South = 0
%       deg, North = -180 deg

% ------OUTPUTS------ %
%     GHI_E2_P90 = global horizontal irradiance [W/m2] (array) calculated
%       from Engerer2 [1] and Perez1990 [2]
%     DHI_E2_P90 = diffuse horizontal irradiance [W/m2] (array) calculated
%       from Engerer2 [1] and Perez1990 [2]
%     GHI_E2_P90 = global horizontal irradiance [W/m2] (array) calculated
%       from Engerer2 [1] and Hay & Davies [3]
%     DHI_E2_P90 = diffuse horizontal irradiance [W/m2] (array) calculated
%       from Engerer2 [1] and Hay & Davies [3]

% References
% [1] Engerer, N.A., 2015. Minute resolution estimates of the diffuse fraction of
%   global irradiance for southeastern australia. Solar Energy 116, 215–237.
%   doi:10.1016/j.solener.2015.04.012
% [2] Perez, R., Ineichen, P., Seals, R., Michalsky, J., Stewart, R., 1990. Modeling
%   daylight availability and irradiance components from direct and global
%   irradiance. Solar Energy 44, 271–289. doi:10.1016/0038-092X(90)90055-H.
% [3] Hay, J., Davies, J., 1980. Calculations of the solar radiation incident on an
%   inclinedsurface, in: Hay, J., Won, T.(Eds.), ProceedingsofFirstCanadian
%   Solar Radiation Data Workshop, Ministry of Supply and Services, Canada.
%   p. 59
% [4] Bright, J.M., Engerer, N.A., 2019. Engerer2: Global re-parameterisation,
%   update, and validation of an irradiance separation model at different
%   temporal resolutions. Journal of Renewable and Sustainable Energy 11.
%   doi:10.1063/1.5097014.
% [5] Sandia Labs. MATLAB_PV_LIB. https://github.com/sandialabs/MATLAB_PV_LIB

% Load example data 1-h
clear all
load("example_input_1h.mat")

%% Reverse Transposition Algorithm with Engerer2 + P90
% Calculate GHI from GTI, minimizing RMSE of GTI predicted and GTI measured.
GHI_E2_P90 = arrayfun(@(gti, eexth, rb, alb, zen, ast, gcs, saz, extn) ...
    getGHI_engerer2_P90(gti, eexth, rb, Rr, alb, zen, ast, gcs, Tilt, Surface_Azimuth, saz, extn, mode), ...
    GTI, Eexth, Rb, Albedo, Zenith_angle, ast, Gcs, Solar_azimuth, Eextn);

% Calculate DHI from GHI with Engerer2
DHI_E2_P90 = from_engerer2(GHI_E2_P90, Eexth, Zenith_angle, ast, Gcs, mode);

%% Reverse Transposition Algorithm with Engerer2 + HD
% Calculate GHI from GTI, minimizing RMSE of GTI predicted and GTI measured.
GHI_E2_HD = arrayfun(@(gti, eexth, rb, alb, zen, ast, gcs) ...
    getGHI_engerer2_HD(gti, eexth, rb, Rr, alb, zen, ast, gcs, Tilt, mode), ...
    GTI, Eexth, Rb, Albedo, Zenith_angle, ast, Gcs);

% Calculate DHI from GHI with Engerer2
DHI_E2_HD = from_engerer2(GHI_E2_HD, Eexth, Zenith_angle, ast, Gcs, mode);

%% Plot Data
TT_E2_HD = timetable (time, GTI, DHI_E2_HD, GHI_E2_HD);
TT_E2_P90 = timetable (time, GTI, DHI_E2_P90, GHI_E2_P90);

figure;

% Plot TT_E2_HD
subplot(2,1,1);
hold on;
plot(TT_E2_HD.time, TT_E2_HD.GTI, 'r', 'DisplayName', 'GTI');
plot(TT_E2_HD.time, TT_E2_HD.DHI_E2_HD, 'b', 'DisplayName', 'DHI');
plot(TT_E2_HD.time, TT_E2_HD.GHI_E2_HD, 'k', 'DisplayName', 'GHI');
hold off;
legend('Location', 'best');
xlabel('Time');
ylabel('Irradiance (W/m²)');
title('E2 with HD');
grid on;

% Plot TT_E2_P90
subplot(2,1,2);
hold on;
plot(TT_E2_P90.time, TT_E2_P90.GTI, 'r', 'DisplayName', 'GTI');
plot(TT_E2_P90.time, TT_E2_P90.DHI_E2_P90, 'b', 'DisplayName', 'DHI');
plot(TT_E2_P90.time, TT_E2_P90.GHI_E2_P90, 'k', 'DisplayName', 'GHI');
hold off;
legend('Location', 'best');
xlabel('Time');
ylabel('Irradiance (W/m²)');
title('E2 with P90');
grid on;


%% Function Definitions
function GHI = getGHI_engerer2_P90(GTI, Eexth, Rb, Rr, Albedo, Zenith_angle, ast, Gcs, Tilt, Surface_Azimuth, Solar_azimuth, Eextn, mode)
    FitnessFunction = @(x) transp_engerer2_P90(x, GTI, Eexth, Rb, Rr, Albedo, Zenith_angle, ast, Gcs, Tilt, Surface_Azimuth, Solar_azimuth, Eextn, mode);
    opts = optimoptions('fmincon', 'MaxIterations', 1000, 'Display', 'off');
    GHI = fmincon(FitnessFunction, 0, [], [], [], [], 0, 1200, [], opts);
end

function GHI = getGHI_engerer2_HD(GTI, Eexth, Rb, Rr, Albedo, Zenith_angle, ast, Gcs, Tilt, mode)
    FitnessFunction = @(x) transp_engerer2_HD(x, GTI, Eexth, Rb, Rr, Albedo, Zenith_angle, ast, Gcs, Tilt, mode);
    opts = optimoptions('fmincon', 'MaxIterations', 1000, 'Display', 'off');
    GHI = fmincon(FitnessFunction, 0, [], [], [], [], 0, 1200, [], opts);
end

function error = transp_engerer2_P90(GHI, GTI, Eexth, Rb, Rr, Albedo, Zenith_angle, ast, Gcs, Tilt, Surface_Azimuth, Solar_azimuth, Eextn, mode)
    kt = GHI ./ Eexth;
    deltaKtc = max(0, Gcs ./ Eexth - kt);
    Kde = max(0, 1 - max(0, Gcs ./ GHI));

    % Engerer2 decomposition model
    kd_engerer2 = calculate_kd_engerer2(kt, ast, Zenith_angle, deltaKtc, Kde, mode);

    % Diffuse and beam components
    DHI = kd_engerer2 .* GHI;
    DHI = max(0, DHI); % Ensure DHI is non-negative for pvl_perez
    BHI = GHI - DHI;

    % Perez model taken from PVLib [2],[5]
    AM = pvl_relativeairmass(Zenith_angle, 'kastenyoung1989');
    DNI = max(0, BHI ./ cosd(Zenith_angle));

    [GTI_d, ~, ~, ~] = pvl_perez(Tilt, Surface_Azimuth + 180, DHI, DNI, Eextn, Zenith_angle, Solar_azimuth + 180, AM, '1990');
    GTI_engerer2_P90 = BHI .* Rb + GTI_d + Albedo .* GHI .* Rr;

    error = rmse(GTI_engerer2_P90, GTI, 'omitnan');
end

function error = transp_engerer2_HD(GHI, GTI, Eexth, Rb, Rr, Albedo, Zenith_angle, ast, Gcs, Tilt, mode)
    kt = GHI ./ Eexth;
    deltaKtc = max(0, Gcs ./ Eexth - kt);
    Kde = max(0, 1 - max(0, Gcs ./ GHI));

    % Engerer2 decomposition model
    kd_engerer2 = calculate_kd_engerer2(kt, ast, Zenith_angle, deltaKtc, Kde, mode);

    % Diffuse and beam components
    DHI = kd_engerer2 .* GHI;
    BHI = GHI - DHI;

    %Transposition model (Hay&Davies) [3]
    Ai=BHI./Eexth;
    Rd=(1-Ai).*cos(deg2rad(Tilt)/2)^2+Ai*Rb;

    GTI_b = BHI .* Rb;
    GTI_d = DHI .* Rd;
    GTI_r = Albedo .* GHI .* Rr;

    GTI_engerer2_HD = GTI_b + GTI_d + GTI_r;
    error = rmse(GTI_engerer2_HD, GTI, 'omitnan');
end

function kd = calculate_kd_engerer2(kt, ast, Zenith_angle, deltaKtc, Kde, mode)
    if strcmp(mode, '1-min')
        % Global parameterization of Engerer2 1-min (new) [4]
        C = 0.10562;
        b = [-4.1332, 8.2578, 0.010087, 0.00088801, -4.9302, 0.44378];
    elseif strcmp(mode, '1-h')
        % Global parameterization of Engerer2 1-h [4]
        C = -0.0097539;
        b = [-5.3169, 8.5084, 0.013241, 0.0074356, -3.0329, 0.56403];
    else
        error('Invalid mode. Choose either "1-min" or "1-h"');
    end
    % Engerer2 [1]
    kd = (C + (1 - C) ./ (1 + exp(b(1) + b(2) * kt + b(3) * ast + b(4) * Zenith_angle + b(5) * deltaKtc)) + b(6) * Kde);
end

function DHI = from_engerer2(GHI, Eexth, Zenith_angle, ast, Gcs, mode)
    kt = GHI ./ Eexth;
    deltaKtc = max(0, Gcs ./ Eexth - kt);
    Kde = max(0, 1 - max(0, Gcs ./ GHI));
    kd_engerer2 = calculate_kd_engerer2(kt, ast, Zenith_angle, deltaKtc, Kde, mode);
    DHI = kd_engerer2 .* GHI;
end
