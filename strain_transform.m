function [strain_components] = strain_transform(thetaA, thetaB, thetaC, gauge_strain)
    % Transformation matrix
    transform_matrix = [cosd(thetaA)^2, sind(thetaA)^2, 2*cosd(thetaA)*sind(thetaA);
                        cosd(thetaB)^2, sind(thetaB)^2, 2*cosd(thetaB)*sind(thetaB);
                        cosd(thetaC)^2, sind(thetaC)^2, 2*cosd(thetaC)*sind(thetaC)];
    % Strain components epsilon_x, epsilon_y and gamma_yx
    %strain_components = inv(transform_matrix) * ((gauge_strain.*(10e-06))');
    strain_components = transform_matrix \ ((gauge_strain.*(10^-6))');


    end