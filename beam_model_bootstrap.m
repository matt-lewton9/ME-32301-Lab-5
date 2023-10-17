function [E_conf_int, nu_conf_int] = beam_model_bootstrap(sigma_x, tau_xy, epsilon_x, epsilon_y, gamma_xy)


n = length(sigma_x); % number of data points
numBootstraps = 1000; % number of bootstrap samples

% Initial guess for E and nu
E0 = 1;
nu0 = 0.3;

% Define objective function
objfun = @(x,sigma_x,tau_xy,epsilon_x,epsilon_y,gamma_xy)...
sum((sigma_x-x(1)*(epsilon_x + x(2)*epsilon_y)./(1-x(2)^2)).^2 ...
+ (-x(1)*(x(2)*epsilon_x + epsilon_y)./(1-x(2)^2)).^2 ...
+ (tau_xy-x(1)*gamma_xy/2./(1+x(2))).^2);

% Define constraints for nu (-1<nu<0.5 for most materials)
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
lb = [0,-1]; % lower bounds
ub = [Inf,0.5]; % upper bounds

% Initialize array to hold bootstrap estimates
bootstrap_estimates = zeros(numBootstraps,2);

% Perform bootstrap estimation
for i = 1:numBootstraps
    % Generate bootstrap sample
    ind = randi(n,n,1);
    sigma_x_b = sigma_x(ind);
    tau_xy_b = tau_xy(ind);
    epsilon_x_b = epsilon_x(ind);
    epsilon_y_b = epsilon_y(ind);
    gamma_xy_b = gamma_xy(ind);
    % Call fmincon to minimize objective function for bootstrap sample
    x_b = fmincon(@(x)...
    objfun(x,sigma_x_b,tau_xy_b,epsilon_x_b,epsilon_y_b,gamma_xy_b),...
    [E0;nu0],[],[],[],[],lb,ub,[],options);
    bootstrap_estimates(i,:) = x_b;
end

% Compute confidence intervals
alpha = 0.05; % 95% confidence interval
E_conf_int = quantile(bootstrap_estimates(:,1),[alpha/2, 1-alpha/2]);
nu_conf_int = quantile(bootstrap_estimates(:,2),[alpha/2, 1-alpha/2]);

% Display the confidence intervals
disp(['95% confidence interval for E:', num2str(E_conf_int(1)), ', ', num2str(E_conf_int(2)), ']'])
disp(['95% confidence interval for nu: ', num2str(nu_conf_int(1)), ', ', num2str(nu_conf_int(2)), ']'])
fprintf("\n")

end