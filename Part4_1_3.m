% Code for Part IV, Step 1 with Plot
beta0 = 0.3;      % Base transmission rate
A = 5;            % Amplitude of variation
omega= 2 * pi / 365; % Angular frequency for daily periodicity
gamma = 0.1;      % Recovery rate
h = 0.1;          % Time step (days)
T = 30 ;           % Total simulation time (days)
S0 = 990;         % Initial susceptible population
I0 = 10;          % Initial infected population
R0 = 0;           % Initial recovered population
N= S0 +I0 +R0; % Total population

% Time vector
time =0 :h:T;

% Initialize arrays for S, I, and R
S = zeros( size( time ));
I = zeros(size(time));
R = zeros( size(time ));

% Set initial conditions
S(1) = S0;
I(1) = I0;
R(1) = R0;

for t = 1 :length( time )  - 1
    % Compute beta(t) with periodic variation
    beta_t = beta0 * (1 + A * sin(omega * time(t)));
    
    k1_S = -beta_t * S(t) * I(t) / N;
    k1_I =beta_t * S(t) * I(t) / N - gamma * I(t);
    k1_R= gamma * I ( t);
    
    k2_S = -beta_t * (S(t) + h/2 * k1_S) * (I(t) + h/2 * k1_I) / N;
    k2_I = beta_t * (S(t) + h/2 * k1_S) * (I(t) + h/2 * k1_I) / N - gamma * (I( t) + h/2 * k1_I);
    k2_R = gamma * (I(t) + h/2 * k1_I);
    
    k3_S = -beta_t * (S(t) + h/2 * k2_S) * (I(t) + h/2 * k2_I )/ N;
    k3_I = beta_t * (S(t) + h/2 * k2_S) * (I(t) + h/2 * k2_I) / N - gamma * (I(t) + h/2 * k2_I);
    k3_R = gamma * (I(t) + h/2 * k2_I);
    
    k4_S = -beta_t * (S(t) + h * k3_S) * (I(t) + h * k3_I)/ N;
    k4_I = beta_t * (S(t) + h * k3_S) * (I(t) + h * k3_I) / N - gamma * (I(t) + h * k3_I);
    k4_R = gamma * (I(t) + h * k3_I);
    
    S(t+1) = S(t) + h/ 6 * (k1_S + 2*k2_S + 2*k3_S + k4_S);
    I(t+1) = I(t) + h/6 * (k1_I + 2 *k2_I + 2*k3_I + k4_I);
    R(t+1) = R(t) + h/6 * (k1_R + 2* k2_R + 2*k3_R + k4_R);
end

% Plot
figure;
plot(time, S, 'b', 'DisplayName', 'S(t)'); hold on;
plot(time, I, 'r', 'DisplayName', 'I(t)');
plot(time, R, 'g', 'DisplayName', 'R(t)');
xlabel('Time (days)');
ylabel('Population');
legend;
title('SIR Model with Periodic Transmission Rate without B equation');
grid on;
