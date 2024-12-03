% Givens
S(1) = 990;
I(1) = 10;
R(1) = 0;
h = 1;
Time = 100;
%% Seasonal Influenza
% Beta = 0.3 Gamma = 0.4
N = @(S, I, R) S + I + R;
Beta = 0.3;
gamma = 0.1;
dSdt = @(N, S, I) (-Beta/N)*S*I;
dIdt = @(N, S, I) (Beta/N)*S*I - gamma*I;
dRdt = @(I) gamma*I;
Population(1) = S(1) + I(1);
Population2(1) = 0;
Population3(1) =  0;
for i=1:Time
    Population(i+1) = N(S(i),I(i),R(i));
    K1Susceptible = dSdt(Population(i), S(i), I(i));
    K1Infected = dIdt(Population(i), S(i), I(i));
    K1Recovered = dRdt(I(i));
    K2StepsizeS = S(i)+1/2*K1Susceptible*h; %K1 Seasonal Influenza
    K2StepsizeI = I(i)+1/2*K1Infected*h; %K1 Infected
    K2StepsizeR = R(i) + 1/2*K1Recovered*h; %K1 Recovery
    Population2(i+1) = N(K2StepsizeS,K2StepsizeI,K2StepsizeR);
    K2SusceptibleS = dSdt(Population2(i+1), K2StepsizeS, K2StepsizeI); %K2 Seasonal Influenza
    K2Infected = dIdt(Population2(i+1), K2StepsizeS, K2StepsizeI); %K2 Infected
    K2Recovered = dRdt(K2StepsizeI); %K2 Recovery
    K3StepsizeS = S(i)+1/2*K2SusceptibleS*h; %K3 Seasonal Influenza
    K3StepsizeI = I(i) + 1/2*K2Infected*h; %K3 Infected
    K3StepsizeR = R(i) + 1/2*K2Recovered*h; %K3 Recovery
    Population3(i+1) = N(K3StepsizeS,K3StepsizeI,K3StepsizeR);
    K3SusceptibleS = dSdt(Population3(i+1), K3StepsizeS, K3StepsizeI);
    K3Infected = dIdt(Population3(i+1), K3StepsizeS, K3StepsizeI);
    K3Recovered = dRdt(K3StepsizeI);
    K4StepsizeS = S(i) + K3SusceptibleS; %K4 Seasonal Influenza
    K4StepsizeI = I(i) + K3Infected; %K4 Infected
    K4StepsizeR = R(i) + K3Recovered; %K4 Recovery
    Population4(i+1) = N(K4StepsizeS,K4StepsizeI,K4StepsizeR);
    K4SusceptibleS = dSdt(Population4(i+1),K4StepsizeS, K4StepsizeI);
    K4Infected = dIdt(Population4(i+1), K4StepsizeS,K4StepsizeI);
    K4Recovered = dRdt(K4StepsizeI);
    % Adding all the K values up
    S(i+1) = S(i) + 1/6*(K1Susceptible + 2*K2SusceptibleS + 2*K3SusceptibleS + K4SusceptibleS);  
    I(i+1) = I(i) + 1/6*(K1Infected + 2*K2Infected + 2*K3Infected + K4Infected);
    R(i+1) = R(i) + 1/6*(K1Recovered + 2*K2Recovered + 2*K3Recovered + K4Recovered);
end

%% Covid
Sc(1) = 990;
Ic(1) = 10;
Rc(1) = 0;
N = @(S, I, R) S + I + R;
Beta = 1;
gamma = 0.1;
dSdt = @(N, S, I) (-Beta/N)*S*I;
dIdt = @(N, S, I) (Beta/N)*S*I - gamma*I;
dRdt = @(I) gamma*I;
Population(1) = S(1) + I(1);
Population2(1) = 0;
Population3(1) = 0;
Population4(1) = 0;
for i=1:Time
    Population(i+1) = N(Sc(i),Ic(i),Rc(i));
    K1Susceptible = dSdt(Population(i), Sc(i), Ic(i)); %K1 Seasonal Influenza
    K1Infected = dIdt(Population(i), Sc(i), Ic(i)); %K1 Infected
    K1Recovered = dRdt(Ic(i)); %K1 Recovery
    K2StepsizeS = Sc(i)+1/2*K1Susceptible*h;
    K2StepsizeI = Ic(i)+1/2*K1Infected*h;
    K2StepsizeR = Rc(i) + 1/2*K1Recovered*h;
    Population2(i+1) = N(K2StepsizeS,K2StepsizeI,K2StepsizeR);
    K2SusceptibleS = dSdt(Population2(i+1), K2StepsizeS, K2StepsizeI);  %K2 Seasonal Influenza
    K2Infected = dIdt(Population2(i+1), K2StepsizeS, K2StepsizeI); %K2 Infected
    K2Recovered = dRdt(K2StepsizeI); %K2 Recovery
    K3StepsizeS = Sc(i) + 1/2*K2SusceptibleS*h;
    K3StepsizeI = Ic(i) + 1/2*K2Infected*h;
    K3StepsizeR = Rc(i) + 1/2*K2Recovered*h;
    Population3(i+1) = N(K3StepsizeS,K3StepsizeI,K3StepsizeR);
    K3SusceptibleS = dSdt(Population3(i+1), K3StepsizeS, K3StepsizeI); %K3 Seasonal Influenza
    K3Infected = dIdt(Population3(i+1), K3StepsizeS, K3StepsizeI); %K3 Infected
    K3Recovered = dRdt(K3StepsizeI); %K3 Recovery
    K4StepsizeS = Sc(i) + K3SusceptibleS;
    K4StepsizeI = Ic(i) + K3Infected;
    K4StepsizeR = Rc(i) + K3Recovered;
    Population4(i+1) = N(K4StepsizeS,K4StepsizeI,K4StepsizeR);
    K4SusceptibleS = dSdt(Population4(i+1),K4StepsizeS, K4StepsizeI); %K4 Seasonal Influenza
    K4Infected = dIdt(Population4(i+1), K4StepsizeS,K4StepsizeI); %K4 Infected
    K4Recovered = dRdt(K4StepsizeI); %K4 Recovery
    % Adding all the K values up
    Sc(i+1) = Sc(i) + 1/6*(K1Susceptible + 2*K2SusceptibleS + 2*K3SusceptibleS + K4SusceptibleS);
    Ic(i+1) = Ic(i) + 1/6*(K1Infected + 2*K2Infected + 2*K3Infected + K4Infected);
    Rc(i+1) = Rc(i) + 1/6*(K1Recovered + 2*K2Recovered + 2*K3Recovered + K4Recovered);
end
%% Measles
Sm(1) = 990;
Im(1) = 10;
Rm(1) = 0;
N = @(S, I, R) S + I + R;
Beta = 2;
gamma = 0.2;
dSdt = @(N, S, I) (-Beta/N)*S*I;
dIdt = @(N, S, I) (Beta/N)*S*I - gamma*I;
dRdt = @(I) gamma*I;
Population(1) = S(1) + I(1);
Population2(1) = 0;
Population3(1) =  0;
Population4(1) = 0;
for i=1:Time
    Population(i+1) = N(Sm(i),Im(i),Rm(i));
    K1Susceptible = dSdt(Population(i), Sm(i), Im(i)); %K1 Seasonal Influenza
    K1Infected = dIdt(Population(i), Sm(i), Im(i)); %K1 Infected
    K1Recovered = dRdt(Im(i)); %K1 Recovery
    K2StepsizeS = Sm(i)+1/2*K1Susceptible*h;
    K2StepsizeI = Im(i)+1/2*K1Infected*h;
    K2StepsizeR = Rm(i) + 1/2*K1Recovered*h;
    Population2(i+1) = N(K2StepsizeS,K2StepsizeI,K2StepsizeR);
    K2SusceptibleS = dSdt(Population2(i+1), K2StepsizeS, K2StepsizeI); %K2 Seasonal Influenza 
    K2Infected = dIdt(Population2(i+1), K2StepsizeS, K2StepsizeI); %K2 Infected
    K2Recovered = dRdt(K2StepsizeI); %K2 Recovery
    K3StepsizeS = Sm(i)+1/2*K2SusceptibleS*h;
    K3StepsizeI = Im(i) + 1/2*K2Infected*h;
    K3StepsizeR = Rm(i) + 1/2*K2Recovered*h;
    Population3(i+1) = N(K3StepsizeS,K3StepsizeI,K3StepsizeR);
    K3SusceptibleS = dSdt(Population3(i+1), K3StepsizeS, K3StepsizeI); %K3 Seasonal Influenza
    K3Infected = dIdt(Population3(i+1), K3StepsizeS, K3StepsizeI); %K3 Infected
    K3Recovered = dRdt(K3StepsizeI); %K3 Recovery
    K4StepsizeS = Sm(i) + K3SusceptibleS;
    K4StepsizeI = Im(i) + K3Infected;
    K4StepsizeR = Rm(i) + K3Recovered;
    Population4(i+1) = N(K4StepsizeS,K4StepsizeI,K4StepsizeR);
    K4SusceptibleS = dSdt(Population4(i+1),K4StepsizeS, K4StepsizeI); %K4 Seasonal Influenza
    K4Infected = dIdt(Population4(i+1), K4StepsizeS,K4StepsizeI); %K4 Infected
    K4Recovered = dRdt(K4StepsizeI); %K4 Recovery
    %Summing them all up
    Sm(i+1) = Sm(i) + 1/6*(K1Susceptible + 2*K2SusceptibleS + 2*K3SusceptibleS + K4SusceptibleS);
    Im(i+1) = Im(i) + 1/6*(K1Infected + 2*K2Infected + 2*K3Infected + K4Infected);
    Rm(i+1) = Rm(i) + 1/6*(K1Recovered + 2*K2Recovered + 2*K3Recovered + K4Recovered);
end
%% Plots
T = 0:h:100;
figure(1)
plot(T,S, 'r', T,I, 'k', T,R, 'b')
grid on
title('Seasonal Infuenza')
legend('S(t) Susceptible Individuals', 'I(t) Infected Individuals', 'R(t) Recovered Individuals')
figure(2)
plot(T,Sc, 'r', T,Ic, 'k', T,Rc, 'b')
title('Covid-19')
legend('S(t)Susceptible Individuals', 'I(t) Infected Individuals', 'R(t) Recovered Individuals')
figure(3)
plot(T,Sm, 'r', T,Im, 'k', T,Rm, 'b')
title('Measles')
legend('S(t) Susceptible Individuals', 'I(t) Infected Individuals', 'R(t) Recovered Individuals')
