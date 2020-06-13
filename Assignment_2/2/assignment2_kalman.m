%{
    Author: Sai Himal Allu
    Course: ECN-614
    Instructor: Dr Debashis Ghosh
    Email: sallu@ec.iitr.ac.in

    The following module is intended to be a part of the solution set to an ECN-614 
    assignement that was undertaken by the author as a part of his regular coursework 
    for Spring 2020 semester. Duplication of any part of this assignment without the 
    prior permission of the author is strictly prohibited.
%}

%{
    Question:
        Given the measured positions s(k) at time instants k = 1, 2, ..., a Kalman filter is used
        to recursively estimate the position and velocity of the truck at all times. Write a
        MATLAB program for the same.
        Generate an arbitrary test sequence of 500 samples of measured truck positions as
        s(k) = s(k − 1) + r(k), 1 ≤ k ≤ 500, where r(k) is a random number distributed
        uniformly in the range 0 to 0.5, and s(0) = 0.
        Assuming that the initial state distribution is Gaussian with zero mean and identity
        covariance matrix, i.e. N (0, I), run the above program to estimate the position of the
        truck at all times till k = 501. Hence, plot the measured path of the truck, estimated
        path and the error in estimating the position of the truck.

%}

%{
    Solution Idea:
        Generate the observed positions of the truck by using the relation given in the question.
        Use the Kalman filter principles to estimate the actual position of the truck. These positions 
        are then plotted. The Kalman filter principles which are used in this solution set uses the 
        relations defined in Simon Haykin's Adaptive Filter Theory
   
%}

clear all;
clc;
clf;

%%define the meta-variables
duration = 500;                                                  %number of samples
dt = 1;                                                         %sampling interval

F = [1 dt; 0 1] ;                                               %state transition matrix
C = [1 0];                                                      %measurement matrix

%%initial distribution of states is over a N(0, I) 
mu = [0 0];
sigma = [1 0; 0 1];
rng('default')                                                  %set seed for reproducibility
Q_estimate = mvnrnd(mu, sigma, 1);                              %simulate a initial state from N(0, I)
Q_estimate = Q_estimate';

G = zeros(2, 1);                                                %initialize kalman gain
K = eye(2);                                                     %initialize the K matrix as defined in the question
Ez = 0.25;                                                      %covariance matrix for gps error
Ex = 1 * [dt^4/4 dt^3/2; dt^3/2 dt^2];                          %covariance error matrix for acceleration  


%Generate the measured values for the positions
rng(0, 'twister')                                               %set seed for reproducability
r = 0.5 .* rand(duration, 1);                                   %uniform distribution between 0 and 0.5
S_measured = zeros(duration, 1); 

for k=1:duration
    if k == 1
        S_measured(k) = r(k);
    else
        S_measured(k) = S_measured(k-1) + r(k);
    end
end

disp("Run the Kalman filtering on the measured positions");

S_estimate = [];                                                %array for storing kalman estimated values of truck position
vel_estimate = [];                                              %array for stroing kalman estimated values of truck velocity

for t = 1: dt: duration
    
    R = C * K * C' + Ez;                                        %
    
    G = F * K * C' * inv(R);                                    %compute the kalman gain

    alpha_t = S_measured(t) - C *  Q_estimate;                  %calculate the innovation
    
    Q_estimate = F * Q_estimate + G * alpha_t;                  %update the state estimate
    
    S_estimate = [S_estimate; Q_estimate(1)];                   %store position estimate for plotting
    
    vel_estimate = [vel_estimate, Q_estimate(2)];               %store velocity estimate for plotting
    
    k_hat = K - F * G * C * K;                                  

    K = F * k_hat * F' + Ex;                                    %update the covariance estimate

end


error_sig = S_estimate - S_measured;                            %calculate the error signal between the estimated and the measured values 

% Plot the results
figure(1);
plot(S_measured, '-r.');

figure(2);
plot(S_estimate, '-g.');

figure(3);
plot(error_sig,  '-k.');  