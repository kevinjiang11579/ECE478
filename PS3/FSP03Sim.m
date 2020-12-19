%% Simulating Stochastic Differential Equations, Kevin Jiang
N = 250;
del = 0.01;
drift = 0.1;
volat = 0.2;
r = 0.05;
% (a)
%Under risk-netural measure, dS(t) = r*S(t)*dt + sigma*S(t)*dW~(t).
    %r = 0.05, sigma = 0.2
    %dW~(t) = theta*dt + dW(t), theta = (alpha-r)/sigma = 0.25
% (b)
S = zeros(1000,N);
S(:,1) = ones(1000,1);
for i = 2:N % Generate 1000 paths
    dW = randn(1000,1)*sqrt(del);
    dS = drift*S(:,i-1)*del + volat*S(:,i-1).*dW;
    S(:,i) = dS + S(:,i-1);
end

% (c)
ESN2 = mean(S(:,floor(N/2)));
ESN = mean(S(:,N));

% (d);
K = ESN;
tau = del*N - del*N/2;
VN2 = zeros(1000, 1);
SN2values = sort(S(:,floor(N/2)));
for i = 1:1000
    VN2(i,1) = BSM(tau, SN2values(i,1), K, r, volat);
    if VN2(i,1) < 0
        VN2(i,1) = 0;
    end
end
figure;
plot(SN2values, VN2);
title("V[N/2] with BSM");
ylabel("Option value V(t)");
xlabel("Stock price S(t)");

% (e)
Si = S(1:10,:);
SiN2 = Si(:,floor(N/2));
Ke = mean(Si(:,N));
Vi = zeros(10,1);
for i = 1:10
    Vi(i,1) = BSM(tau, SiN2(i,1), K, r, volat);
    if Vi(i,1) < 0
        Vi(i,1) = 0;
    end
end
ViMC = zeros(10,1);
for i = 1:10
    ViMC(i,1) = montecarlo(SiN2(i,1), floor(N/2), N, r, del, drift, volat, K);
    ViMC(i,1) = ViMC(i,1)/exp(-r*del*N/2);
end
figure;
scatter(SiN2, Vi, 'r');
hold on;
scatter(SiN2, ViMC, 'b');
title("V[N/2] two ways: BSM (red) vs Monte Carlo (Blue)");
ylabel("Option value V(t)");
xlabel("Stock price S(t)");

function [c] = BSM(tau, S, K, r, volat)
    dplus = (log10((S/K) + (r+(volat^2)/2)*tau))/(volat*sqrt(tau));
    dminus = (log10((S/K) + (r-(volat^2)/2)*tau))/(volat*sqrt(tau));
    c = S*normcdf(dplus) - K*exp(-r*tau)*normcdf(dminus);
end

function [DV] = montecarlo(S_init, start, N, r, del, drift, volat, K)
    Spaths = zeros(1000, N-start);
    Spaths(:,1) = S_init;
    for i = (start+2):N
        pI = i - start;
        dW = randn(1000,1)*sqrt(del);
        dWt = ((drift-r)/volat)*del + dW;
        dS = r*Spaths(:,pI-1)*del + volat*Spaths(:,pI-1).*dW;
        Spaths(:,pI) = dS + Spaths(:,pI-1);
    end
    V = Spaths(:,N-start)-K;
    V(find(V < 0)) = 0;
    DV = exp(-r*del*N)*V;
    DV = mean(DV);
end