%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DVM-LMS-1 Algorithm Implementation   %
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
numberOfIterations = 10;
K = 2; % Series order 
M = 10; % Filter order
eps = 1e-3; 

% Convergence parameters
mu0 = 0.04; % Convergence coefficient

% Input signal correlation parameters
sig_u2 = 1; % Input signal variance
a = 0; % Regulates the correlation
b = sqrt(1 - a^2) * sig_u2; % Corrects the input signal variance to sig_u2

% Desired signal noise parameters
sig_v2 = 1e-6; % Corruption noise variance

% The desired plant will be a Ko-th series order and Mo-th filter order 
% decomposable Volterra plant. To impose the decomposability condition, the
% vector containing the kernels of the Volterra plant (wo) has to be 
% equal to the kronecker product of wAux(1), wAux(2), ..., wAux(Ko).
Ko = 2;
Mo = 10;
wAux = rand(Mo, Ko) * 0.01;
wo = 1;
for k = 1:Ko
    wo = kron(wo, wAux(:, k));
end

for l = 1:numberOfIterations
    % Initializing the input information, u (gaussian distribution with variance
    % sig_u2)
    
    input = rand(numberOfSamples, 1) * sig_u2;
   
    % Correlating the input information
    
    whiteNoise = rand(numberOfSamples, 1) * 1; % White noise with variance = 1
    for i = 2:numberOfSamples
        input(i) = a * input(i - 1) + b * whiteNoise(i);        
    end
   
    % Initializing row input regressors
    ui = zeros(1, M);
    uio = zeros(1, Mo); % Desired Volterra plant
    % Let's add some gaussian noise to d
    
    % Initializing error, desired and output vectors
    e = zeros(1, numberOfSamples);
    d = e;
    output = e;
    
    % Initializing the DVM-LMS-1 Algorithm weights
    w = zeros(M, K);
    for i = 1:K - 1
        w(1, i) = 2^(1 - i); 
    end 
    
    % Iteration
    for i = 1:numberOfSamples
        % Creating row input regressor
        ui = [input(i) ui(1:M - 1)];
        uio = [input(i) ui(1:Mo - 1)];
        
        % Defining desired output
        uioK = 1;
        for k = 1:Ko
            uioK = kron(uioK, uio); 
        end
        d(i) = uioK * wo;
        v = randn() * sig_v2;
        d(i) = d(i) + v;
        
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
        output(i) = y(K) * yo(K);
        e(i) = d(i) - output(i);
        
        % Updating weight matrix
        mu = zeros(1, K); 
        for s = 1:K
            mu(s) = mu0 / (abs(y(s))^2 + eps);
            w(:, s) = w(:, s) + mu(s) * e(i) * conj(y(s)) * ui';
        end
    end
    squaredErrorMatrix(l, :) = e.^2;
    outputMatrix(l, :) = output;
end

mse = mean(squaredErrorMatrix);
outputMean = mean(outputMatrix);

% Plotting Graphs
figure;
subplot(1,2,1)
plot(10 * log(mse));
xlabel('i');
ylabel('Mean Squared Error [dB])');
ylim([-100, 20])
grid;
subplot(1,2,2)
hold on;
plot(outputMean, 'r');
%plot(d, 'r');
xlabel('i');
ylabel('Outputs');
grid;
hold off;
