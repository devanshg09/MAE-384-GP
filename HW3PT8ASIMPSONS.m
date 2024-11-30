%% Homework 3 

% Problem 8
% Part A Part II
% Simpson's 1/3

a = 0.7; % lower limit of integration
b = 6.3; % upper limit integration
c = 7; % chord length in meters
n = 8; % number of subintervals
rho = 0.302; % density in kg/(m^3)
S = 100; % surface area of the wing in m^2
U = 240; % velocity in m/s

h = (b-a)/n; % width of subintervals

x = a:h:b; % coordinates of subintervals

y = [0.7249, 0.6952, 0.7409, 0.6648, 0.5757, 0.5896, 0.5959, 0.2807, 0.4429];   % Cp,l - Cp,u

% calculations

SimpInt = (h/3)*((y(1) + ((4)*sum(y(:,2:2:end))) +((2)*sum(y(:,3:2:7))) + y(9))) % Simpsons 1/3 approximation of the integral

Cl = SimpInt * (1/c) % lift coefficient dimensionless

L = (1/2)*rho*(U)^2*S*Cl % lift force in Newtons
