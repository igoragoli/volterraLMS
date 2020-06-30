%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volterraLMSAlgorithm as seen on      %
% Paulo Diniz - Adaptive Filtering     %
% Igor Oliveira                        % 
% 30/06/2020                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements the LMS Volterra Filter  for a second-order series
% and a N-th order filter.
clc;
clear all;
close all;

%% Defining initial parameters

numberOfSamples = 500 + 1;
N = 10; % Filter order

mu1 = 0.001; % 1st order kernels convergence coefficient 
mu2 = 0.001; % 2nd order kernels convergence coefficient

m = (N+1) + (N+1)^2; % LMS Volterra Filter vector size (the size of the
                     % input vector x, of the weight vector d, etc.)
                     
% Initializing the input vector
x = zeros(m, 1);

% Initializing the weight vector
% The weight vector in the LMS Volterra filter is composed of 
% Volterra kernels, "w_{oi}".
w = zeros(m, 1);

% Initializing the input information, which we will call u (gaussian distribution with sigma^2 = 1)
u = randn(numberOfSamples, 1);

% Initializing the desired plant vector (Example 11.1, page 459 on Diniz's book)
t = linspace(0, 10, numberOfSamples);
d = zeros(numberOfSamples, 1);
d(1) = 0;
for i = 2:numberOfSamples
    d(i) = u(i) + 0.5*u(i - 1);
end
for i = 1:numberOfSamples
    d(i) = d(i) + 0.2*d(i)^2 + 0.1*d(i)^3;
end

% Let's add some gaussian noise to d
sig_v2 = 1e-6;
v = randn(numberOfSamples, 1)*sig_v2;
d = d + v;

L = 10; % Number of runs
numberOfIterations = numberOfSamples - N; % Number of iterations per run

% ==============================
% Running the algorithm
% ==============================

% Create the convergence coefficient matrix
mu = muMatrix(N, mu1, mu2);
% Initialize mean square error;
mse = zeros(numberOfIterations,1);
for i = 1:L
    for k = 0:numberOfIterations-1
        x = inputVector(u, N, k);
        y = x'*w;
        e = d(k + 1) - x'*w;
        se = e^2;
        mse(k+1) = mse(k+1) + se/L;
        w = w + 2*mu*e*x;
    end
end

% ==============================
% Plotting graphs
% ==============================

fig = figure(1);
t =  0 : 1 : numberOfIterations-1;
plot(t, 10*log10(mse));

% ==============================
% Function definitions
% ==============================

% newIndex(N, k):
% "samples" vector ranges from 1 to numberOfSamples. But we would like it to range from
% -N to numberOfSamples - (N + 1). newIndex(N, k) takes an 
% index "k" in this desired range and transforms it to "sample"'s real range.
% Arguments:
%   - k : index in [-N, numberOfSamples - (N + 1)]
function i = newIndex(N, k)
    i = k + (N + 1);
end

% inputVector(samples, N, k):
% Constructs the LMS Volterra filter input vector "x", which contains all
% the x(k - l_1)x(k - l_2) ... x(k - l_i) in the volterra series expansion.
% Arguments:
%   - samples : column vector of all samples.
%   - k : index in [-N, numberOfSamples - (N + 1)].
function x = inputVector(u, N, k)
    % Initialize the inputVector "x"
    x = zeros((N+1) + (N+1)^2, 1);
    % Correct the index
    k = newIndex(N, k); 
    % Construct the first part, with only one term
    for i = 0:N
        x(i+1) = u(k-i);
    end
    % Construct the second part, with two terms
    for i = 0:N
        for j = 0:N
            x(N+2 + i*N + i+j) = u(k-i)*u(k-j);
        end
    end
end

% muMatrix(N, mu1, mu2)
% Creates the convergence coefficient matrix
% Arguments:
%   - N : filter order
%   - mu1, mu2 : convergence coefficients
function M = muMatrix(N, mu1, mu2)
    m = (N+1) + (N+1)^2;
    M = zeros(m);
    for i = 1:N+1
        M(i,i) = mu1;
    end
    for i = (N+2):m
        M(i,i) = mu2;
    end
end 