function [E, nu] = beam_model(sigma_x, tau_xy, epsilon_x, epsilon_y, gamma_xy)

    epsilon_x = epsilon_x';
    epsilon_y = epsilon_y';
    gamma_xy = gamma_xy';
    sigma_x = sigma_x';
    tau_xy = tau_xy';

    % Initial guess for E and nu
    E0 = 10e6;
    nu0 = 0.3;
    % Define objective function

% % ORIGINAL
%     objfun = @(x)...
% sum((sigma_x-x(1)*(epsilon_x + x(2)*epsilon_y)./(1-x(2)^2)).^2 ...
% + (-x(1)*(x(2)*epsilon_x + epsilon_y)./(1-x(2)^2)).^2 ...
% + (tau_xy-x(1)*gamma_xy/2./(1+x(2))).^2);

 % FIXED
    objfun = @(x)...
    sum( (  (sigma_x - (x(1)*((epsilon_x + (x(2)*epsilon_y)))./(1-(x(2)^2)) )  ).^2 ) ...
+ ( (-x(1)*((x(2)*epsilon_x) + epsilon_y)./(1-(x(2).^2))).^2) ...
+ ((tau_xy - (x(1)*gamma_xy/(2.*(1+x(2)))  ) ).^2) );


    % Define constraints for nu (-1<nu<0.5 for most materials)
    options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
    lb = [0,-1]; % lower bounds
    ub = [Inf,0.5]; % upper bounds

    % Call fmincon to minimize objective function
    x = fmincon(objfun,[E0;nu0],[],[],[],[],lb,ub,[],options);
    % Display the optimal values of E and nu
    E = x(1);
    nu = x(2);

    disp(['Optimal value of E: ', num2str(E)])
    disp(['Optimal value of nu: ', num2str(nu)])
    fprintf("\n")
    end