%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volterraLMSAlgorithmFull             %
% Igor Oliveira                        % 
% 30/06/2020                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements a parametrized LMS Volterra Filter for a K-th order series
% and a N-th order filter.

clc;
clear;
close;

%% Parameter definition
% Filter Parameters
numberOfRuns = 50;
numberOfSamples = 1000;
K = 2; % Series order 
N = 10; % Filter order

% Convergence parameters
mu = 0.01; % Convergence coefficient

% Input signal correlation parameters
sig_u2 = 1; % Input signal variance
a = 0.05; % Regulates the correlation
b = sqrt(1 - a^2) * sig_u2; % Corrects the input signal variance to sig_u2

% Desired signal noise parameters
sig_v2 = 1e-6; % Corruption noise variance

% Calculating the Volterra LMS Filter vector size (the size of the
% input vector x, of the weight vector w etc.)
m = 0;
for k = 1:K
  m = m + N^k;
end

%% Running the algorithm

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
    z = zeros(N, 1);
    
    % Initializing weight vector w. The weight vector in the LMS Volterra
    % filter is composed of Volterra kernels "w_{oi}"
    w = zeros(m, 1);
    
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
    
    for i = 1:numberOfSamples
        % Delaying input vector samples
        z = [u(i); z(1:N - 1)];
        
        % Creating input vector
        x = z;
        xAux = x;
        for j = 1:K - 1
            xAux = kron(z, xAux);
            x = [x; xAux];    
        end
        
        % Calculating the output
        y = x' * w;
        
        % Calculating the error
        a = dNoise(i) - y;
        e(i) = a;
       
        % Updating the weight vector
        w = w + 2 * mu * e(i) * x;
        
    end
    squaredErrorMatrix(l, :) = e.^2;
end

% Calculating MSE
mse = mean(squaredErrorMatrix);

% Plotting Graphs
figure;
plot(10 * log(mse));
xlabel("i")
ylabel("Mean Squared Error [dB])");
ylim([-100, 20])
grid;