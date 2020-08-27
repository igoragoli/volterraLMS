%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DVM-LMS Algorithm Implementation     %
% Igor Oliveira                        % 
% 26/08/2020                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements a parametrized Decomposable Volterra Model LMS Filter for a
% K-th order series and a M-th order filter.

clc;
clear;
close;

% Filter parameters
numberOfSamples = 1000;
numberOfRuns = 50;
K = 2; % Series order 
M = 8; % Filter order
eps = 1e-3;

% Convergence parameters
mu0 = 0.01; % Convergence coefficient

% Input signal correlation parameters
sig_u2 = 1; % Input signal variance
a = 0; % Regulates the correlation
b = sqrt(1 - a^2) * sig_u2; % Corrects the input signal variance to sig_u2

% Desired signal noise parameters
sig_v2 = 1e-6; % Corruption noise variance

for l = 1:numberOfRuns
    % Initializing the input information, u (gaussian distribution with variance
    % sig_u2)
    
    u = rand(numberOfSamples, 1) * sig_u2;
   
    % Correlating the input information
    
    whiteNoise = rand(numberOfSamples, 1) * 1; % White noise with variance = 1
    for i = 2:numberOfSamples
        u(i) = a * u(i - 1) + b * whiteNoise(i);        
    end
   
    % Initializing the delayed input vector
    ui = zeros(1, M);
    
    % Defining the desired plant vector (Example 11.1, pg. 459 on Diniz's book)
    d = zeros(numberOfSamples, 1);
    p = zeros(numberOfSamples, 1);
    for i = 2:numberOfSamples
        p(i) = u(i) + 0.5 * u(i - 1);
    end
    for i = 1:numberOfSamples
        d(i) = p(i) + 0.2 * p(i)^2 + 0.1 * p(i)^3;
    end

    % Let's add some gaussian noise to d
    v = randn(numberOfSamples, 1) * sig_v2;
    dNoise = d + v;
    
    % Creating error vector
    e = zeros(1, numberOfSamples);
    
    % Initialization
    w = zeros(M, K);
    for i = 1:K - 1
        w(1, i) = 2^(1 - i); 
    end 
    
    % Iteration
    for i = 1:numberOfSamples
        % Creating row input regressor
        ui = [u(i) ui(1:M - 1)];
        
        % Computing y_{os}
        yo = zeros(1, K);
        for s = 1:K
            yo(s) = ui * w(:, s);
        end
        
        % Computing y_s
        y = zeros(1, K);
        for s = 1:K
            y(s) = prod(yo(1:s - 1))*prod(yo(s + 1:K));
        end
        
        % Computing the output and the error
        yout = y(K) * yo(K);
        e(i) = d(i) - yout;
        
        % Updating weight matrix
        mu = zeros(1, K); 
        for s = 1:K
            mu(s) = mu0 / (abs(y(s))^2 + eps);
            w(:, s) = w(:, s) + mu(s) * e(i) * conj(y(s)) * ui';
        end
    end
    squaredErrorMatrix(l, :) = e.^2;
end

mse = mean(squaredErrorMatrix);

% Plotting Graphs
figure;
plot(10 * log(mse));
xlabel('i');
ylabel('Mean Squared Error [dB])');
ylim([-100, 20])
grid;