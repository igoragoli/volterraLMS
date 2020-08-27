%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volterraLMSAlgorithmKron as seen on  %
% Paulo Diniz - Adaptive Filtering,    %
% now using the Kronecker product      %
% Igor Oliveira                        % 
% 30/06/2020                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements the LMS Volterra Filter  for a second-order series
% and a N-th order filter.

clc;
clear;
close;

%% Defining initial parameters

numberOfRuns = 50;
numberOfSamples = 1000;
N = 10; % Filter order

% Calculating the Volterra LMS Filter vector size (the size of the
% input vector x, of the weight vector w etc.)
m = 2 * N + factorial(N) / (factorial(N - 2) * 2);

% Defining convergence coefficient
mu = 0.01;

%% Running the algorithm

for l = 1:numberOfRuns
    % Initializing the input information, u (gaussian distribution with 
    % sigma^2 = 1)
    u = rand(numberOfSamples, 1);
   
    % Initializing the delayed input vector
    z = zeros(N, 1);
    
    % Initializing weight vector w. The weight vector in the LMS Volterra
    % filter is composed of Volterra kernels "w_{oi}"
    w = zeros(m, 1);
    
    % Defining the desired plant vector (Example 11.1, pg. 459 on Diniz's book)
    d = zeros(numberOfSamples, 1);
    p = zeros(numberOfSamples, 1);
    d(1) = 0;
    
    % Initializing error vector
    for i = 2:numberOfSamples
        p(i) = u(i) + 0.5 * u(i - 1);
    end
    for i = 1:numberOfSamples
        d(i) = p(i) + 0.2 * p(i)^2 + 0.1 * p(i)^3;
    end

    % Let's add some gaussian noise to d
    sig_v2 = 1e-6;
    v = randn(numberOfSamples, 1) * sig_v2;
    dNoise = d + v;
    
    for i = 1:numberOfSamples
        % Delaying input vector samples
        z = [u(i); z(1:N - 1)];
        
        % Creating input vector
        zAux = [1; z];
        x = z;
        for j = 1:N
            x = [x; kron(z(j), z(j:N))];            
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