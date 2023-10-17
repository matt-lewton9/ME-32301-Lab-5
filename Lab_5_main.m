clear 
clc

strains = readmatrix("CSV Data.csv",'Range',2); % [e-06] input strains
P = strains(:,8); % [lbf]
Cy = 1; % [in] dist from midspan for X

%% Flexural Stress
% Izz calc
h = 4; % [in] depth
b = 4; % [in] Flange Width
t_w = 0.3675; % [in] Web thickness
t_f = 0.3675; % [in] Flange thickness

Izz_vert = (t_w * (h - 2* t_f)^3) / 12; % [in^4] vertical rectangle
Izz_horiz = 2 * ((b * (t_f^3)/12) + (b * t_f * (((h - t_f)/2)^2))); % [in^4] two horizontal rectangles
Izz = Izz_vert + Izz_horiz; %[in^4]

% First area moments 
s = h/2; % dist from neutral axis to top of beam
y = Cy ; % dist from neutral axis SG
Qf = s * t_f *((h/2) - (t_f/2));
Qw = ( (b/2) * (((h/2)^2) - (((h-(2*t_f))/2)^2))) + (((t_w/2) * ( ((h-(2*t_f)/2)^2) - (y^2) )));

% Flexural Stress
LR = 40; %[in]
LP = 24; %[in]

Mx = P .* (LR-LP) ./4; % [in-lbs]
sigma_xx = Mx .* Cy ./ Izz; % [psi] flexural stress at every load for X mounted beams.
sigma_xz = Mx .* (h/2) ./ Izz; % [psi] flexural stress at every load for Z mounted beams.

txy = 2*(P./2) .* (Qf ./ (Izz .* t_f)); %pure bending
txz = 2*(P./2) .* (Qw ./ (Izz .* t_w)); %Pure pure bending

%% Strain Directional Transformations
for i = 1:8
    SG911(i,:) = strain_transform(-112.5, -67.5, -22.5, [strains(i,9) strains(i,10) strains(i,11)]); % [e-06] Strain transform of SG 9-10-11 to ex, ey, gxy
    SG1214(i,:) = strain_transform(-45, 0, 45, [strains(i,12) strains(i,13) strains(i,14)]); % [e-06] Strain transform of SG 12-13-14 to ez, ex, gxz
    SG1517(i,:) = strain_transform(-90, -45, 0, [strains(i,15) strains(i,16) strains(i,17)]); % [e-06] Strain transform of SG 15-16-17 to ex, ey, gxy
end

E = [sigma_xx ./ SG911(:, 1), sigma_xz ./ SG1214(:, 2), sigma_xx ./ SG1517(:, 1)];
SEM = std(E,1)./sqrt(8); 
 ts = 2.365;
E_Bar= mean(E,1)
CI_E = [E_Bar - (ts*SEM); E_Bar + (ts*SEM)]

v = [SG911(:, 2) ./ SG911(:, 1), SG1214(:, 1) ./ SG1214(:, 2), SG1517(:,2) ./ SG1517(:, 1)];
v_Bar= mean(v,1)
SEMv = std(v,1)/sqrt(8); 
CI_v = [v_Bar - (ts*SEMv); v_Bar + (ts*SEMv)]


   %beam_out911 = beam_model(sigma_xx, txy, SG911(:,1), SG911(:,2), SG911(:,3));
   %beam_out1214 = beam_model(sigma_xz, txz, SG1214(:,2), SG1214(:,1), SG1214(:,3));
   %beam_out1517 = beam_model(sigma_xx, txy, SG1517(:,1), SG1517(:,2), SG1517(:,3));
% 
%     boot_out911 = beam_model_bootstrap(sigma_xx, txy, SG911(:,1), SG911(:,2), SG911(:,3));
%     boot_out1214 = beam_model_bootstrap(sigma_xz, txz, SG1214(:,1), SG1214(:,2), SG1214(:,3));
%     boot_out1517 = beam_model_bootstrap(sigma_xx, txy, SG1517(:,1), SG1517(:,2), SG1517(:,3));
