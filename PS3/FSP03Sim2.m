%% CIR Interest Rate Model, Kevin Jiang
del = 0.01;
t = 0:del:10;
N = length(t);
beta = 1;
alpha = 0.1*beta;
sig = 0.5;
r = 0.05;
%(a)
R = zeros(1000, N);
R(:,1) = r*ones(1000,1);
validPaths = 1:1:1000; % vector keeping track of valid paths with row indices
for i = 2:N
    dW = randn(1000,1)*sqrt(del);
    dR = (alpha - beta*R(validPaths,i-1))*del + sig*(R(validPaths,i-1).^0.5).*dW(validPaths);
    R(validPaths,i) = dR + R(validPaths,i-1);
    lastSize = size(validPaths,1);
    validPaths = find(R(:,i) > 0); % check for paths that are now <= 0
    if size(validPaths,1) < lastSize % More paths became <= 0
        fprintf("Stopped %d paths, %d paths remain\n", lastSize - size(validPaths,1), size(validPaths,1));
    end
end
if size(validPaths,1) < 1000
    disp("R <= 0 occured in at least one path");
end
% (b)
figure;
plot(t, R(validPaths(1),:));
for i = 2:10
    hold on;
    plot(t, R(validPaths(i),:));
end
title("First 10 valid paths for R(t)");
ylabel("R(t)");
xlabel("t");
%(c)
ER1 = mean(R(validPaths, 101));
VarR1 = mean(R(validPaths, 101).^2) - ER1^2;
ER1calc = getE(alpha,beta,r,t(101));
VarR1calc = getVar(alpha,beta,sig,r,t(101));
ER10 = mean(R(validPaths, N));
VarR10 = mean(R(validPaths, N).^2) - ER10^2;
ER10calc = getE(alpha,beta,r,t(N));
VarR10calc = getVar(alpha,beta,sig,r,t(N));
fprintf("Estimated E[R(1)]: %1.4f, Estimated Var(R(1)): %1.4f\n", ER1, VarR1);
fprintf("Caluclated E[R(1)]: %1.4f, Calculated Var(R(1)): %1.4f\n", ER1calc, VarR1calc);
fprintf("Estimated E[R(10)]: %1.4f, Estimated Var(R(10)): %1.4f\n", ER10, VarR10);
fprintf("Caluclated E[R(10)]: %1.4f, Calculated Var(R(10)): %1.4f\n", ER10calc, VarR10calc);
function [Var] = getVar(alpha, beta, sig, r, t)
    Var = (sig^2)*r*(exp(-beta*t)-exp(-2*beta*t))/beta+(alpha*sig^2)*(1-2*exp(-beta*t)-exp(-2*beta*t))/(2*beta^2);
end

function [E] = getE(alpha, beta, r, t)
    E = exp(-beta*t)*r + (alpha/beta)*(1-exp(-beta*t));
end