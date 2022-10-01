% Name: Anuradha Agarwal 
% Advanced numerical methods : Optimization 
% Last edited September 25, 2021

% This script solves the function using newtons method and steepest descent method with initial guesses [1.2,1.2] and [-1.2,1]
% We are trying the find the minimum value of this function —> where the tolerance of the gradient function value is 1.0e-08

alpha_ini = 1; % initial alpha
rho = 0.5; % rho value 
c = 0.0001; % value of c  
x_ini = [1.2;1.2]; % intial guess 1
x_ini2 = [-1.2;1]; % initial guess 2

%% Calling all the functions
[iter, x_p, f_val, alpha, xval, pnew] = backtracking_newton(x_ini, alpha_ini, rho, c); % newton method with first guess
[iters, x_ps, f_vals, alphas, x_vals, psd] = backtracking_sd(x_ini, alpha_ini, rho, c); % sd method with first guess
[iter2, x_p2, f_val2, alpha2, xval2,pnew2] = backtracking_newton(x_ini2, alpha_ini, rho, c); % newton method with second guess
[iters2, x_ps2, f_vals2, alphas2, x_vals2, psd2] = backtracking_sd(x_ini2, alpha_ini, rho, c); % sd method with second guess

%% Plotting


new1 = 1:iter-1;
new2 = 1:iter2-1;
sd1 = 1:iters-1;
sd2 = 1:iters2-1;

figure()
semilogy(new1,f_val(2:end),'LineWidth',5)
title('Semilog plot of Newtons method with initial guess [1.2,1.2]')
xlabel('Iterations')
ylabel('f(x)')

figure()
semilogy(new2,f_val2(2:end),'LineWidth',5)
title('Semilog plot of Newtons method with initial guess [-1.2,1]')
xlabel('Iterations')
ylabel('f(x)')

figure()
semilogy(sd1,f_vals(2:end),'LineWidth',5)
title('Semilog plot of SD method with initial guess [1.2,1.2]')
xlabel('Iterations')
ylabel('f(x)')

figure()
semilogy(sd2,f_vals2(2:end),'LineWidth',5)
title('Semilog plot of SD method with initial guess [-1.2,1]')
xlabel('Iterations')
ylabel('f(x)')

%% functions for backtracking newton 
function [iteration, x, f_val, alpha, x_val, pnew] = backtracking_newton(x_ini, alpha_ini, rho, c)
% The inputs of this function are:
	% x_ini : the initial guess
	% alpha_ini: the initial value of alpha
	% rho : value of rho
	% cv: value of c
% The outputs of this function are:
	% iteration: the number of iterations
	% x: gives the final x value where the minimum of the function is located 
	% f_val: 1d array of function values at each iteration
	% alpha: 1d array of alpha values calculated using opt_alpha function 
	% x_val: 2d array of all the x_values that this function goes through
	% pew: 2d array of p values
iteration = 1; % initializing the iteration 
tolerance = 1e-8; % function tolerance
max_iterations = 1e5; % maximum iterations not to be exceeded
% initializing the vectors to save the outputs needed and allocating space 
f_val = zeros(max_iterations, 1);
x_val = zeros(max_iterations, 2);
pnew = zeros(max_iterations, 2);
alpha = alpha_ini*ones(max_iterations,1);
x = x_ini;
f_val(iteration) = rb_func(x);
x_val(iteration, :) = x_ini;

% checking the condition
while norm(rb_grad(x)) > tolerance && iteration < max_iterations % not to exceed the max iterations and reaching the tolerance level
    iteration = iteration + 1; % increasing the counter
    opt_new = newtonopt(x); % newtons method
    [alpha(iteration), f_val(iteration)] = opt_alpha(x, opt_new, alpha(iteration), rho, c); % updating alpha and function value
    x_val(iteration,:) = transpose(x + alpha(iteration)*opt_new); % updating value of x_val
    x = x + alpha(iteration)*opt_new;
    pnew(iteration,:) = opt_new; % p_newton
end 
f_val = f_val(1:iteration); % final f_val
alpha = alpha(1:iteration); % final alpha values
x_val = x_val(1:iteration, :); % final x_val
pnew = pnew(1:iteration, :); % final p newton
end 

% Backtracking sd
function [iteration, x, f_val, alpha, x_val, psd] = backtracking_sd(x_ini, alpha_ini, rho, c)
% The inputs of this function are:
	% x_ini : the initial guess
	% alpha_ini: the initial value of alpha
	% rho : value of rho
	% cv: value of c
% The outputs of this function are:
	% iteration: the number of iterations
	% x: gives the final x value where the minimum of the function is located 
	% f_val: 1d array of function values at each iteration
	% alpha: 1d array of alpha values calculated using opt_alpha function 
	% x_val: 2d array of all the x_values that this function goes through
	% pew: 2d array of p values
iteration = 1; % initializing interation 
tolerance = 1e-8; % tolerance
max_iterations = 1e5; % maximum iterations
% initializing the vectors to save the outputs needed and allocating space 
f_val = zeros(max_iterations, 1);
x_val = zeros(max_iterations, 2);
psd = zeros(max_iterations, 2);
alpha = alpha_ini*ones(max_iterations,1);
alpha(1) = 1;
x = x_ini;
f_val(iteration) = rb_func(x);
x_val(iteration, :) = x_ini;
% checking the condition 
while norm(rb_grad(x)) > tolerance && iteration < max_iterations % not to exceed the max iterations and reaching the tolerance level
    iteration = iteration + 1; % increasing the counter
    opt_sd = sdopt(x); % sd method
    [alpha(iteration), f_val(iteration)] = opt_alpha(x, opt_sd, alpha(iteration), rho, c); % updating alpha and function value
    x_val(iteration,:) = transpose(x + alpha(iteration)*opt_sd); % updating value of x_val
    x = x + alpha(iteration)*opt_sd;
    psd(iteration,:) = opt_sd; % p_sd
end 
f_val = f_val(1:iteration); % final f_val
alpha = alpha(1:iteration); % final alpha values
x_val = x_val(1:iteration, :); % final x_val
psd = psd(1:iteration, :); % final p sd
end 

function [alpha, f_x_k] = opt_alpha(x_k, p_k, alpha_ini, rho, c)
% This function return the function value at x and updated optimized alpha
% The inputs are x_k which is the x value, p_k which is the p value, alpha_ini which is 1, rho and c which are 0.5 and 0.0001 respectively
% The function returns alpha and f_x_k
% Return the step length based on first Wolfe condition
alpha = alpha_ini;
f_x_k = rb_func(x_k); % finding the function value 
while (rb_func(x_k + (alpha * p_k)) > f_x_k + c * alpha * transpose(p_k) * rb_grad(x_k) ) % condition
  alpha = rho * alpha; % updating alpha
end
end

%% Newtons method
function opt_new = newtonopt(x_k)
opt_new = (-rb_hess(x_k))^-1 * rb_grad(x_k); % formula from slides
end 

%% Stepest descent 
function opt_sd = sdopt(x_k)
opt_sd = -rb_grad(x_k);
opt_sd = opt_sd/norm(opt_sd); % formula from slides 

end 

%% Below are the Rosenbrock function, the gradient of the function and the hessian 
function func = rb_func(x)
func = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
end

% The gradient and the hessian are solved by hand
function gradfunc = rb_grad(x)
gradfunc = [2 * x(1) - 400 * x(1) * (- x(1)^2 + x(2)) - 2;
      200 * (x(2) - x(1)^2)];
end

function hessfunc = rb_hess(x)
hessfunc = [2 + 1200 * x(1)^2 - 400 * x(2), -400*x(1);
      -400 * x(1), 200];
end