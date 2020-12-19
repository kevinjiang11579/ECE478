clear;
% Q4 Arima with s_t
% Same as before, much closer to AR model.
N = 250;
M = 10;
b = [1, 0.7, 0.1];
a = [1, -2.49, 2.045, -0.5544];
borig = [1, 0.7, 0.1];
aorig = [1, -1.5, 0.56];
v = trnd(5,N,1);
v = v/sqrt(var(v));
r = filter(b,a,v); %get r using the operator notation for ARMA model
s = r - [0;r(1:end-1)];
rorig = filter(borig,aorig,v); %get r using the operator notation for ARMA model
% (a)
[sauto,lags] = xcorr(s, M, 'normalized');
sauto = sauto(M+1:end); %only take positive time lag
figure;
stem(lags(M+1:end),sauto);
% (b)
C = toeplitz(sauto);
% (c) using Cholesky because ldl() messes up rows
Lchol = transpose(chol(C));
D = (eye(M+1).*Lchol).^2;
L = Lchol*D^(-0.5);
Linv = L^(-1);
% (d)
F = eye(M+1); 
refMat = zeros(M); % Holds all reflection coeffs
P = zeros(M,1); % Holds all prediction error powers
for i = 1:M
    [fpc,pep,ref] = levinson(sauto, i);
    F(i+1,1:i) = fliplr(fpc(2:end));
    refMat(1:i,i) = ref;
    P(i) = pep;
end
% (e)
FCF = F*C*transpose(F);
compPD = P - diag(D(2:end,2:end)); %compare P_M and D
compLinvF = F - Linv; %compare L^-1 and F
fprintf("Elements greater than 10^12 in P_M and D difference: %d\n", length(find(compPD > 10^(-12))));
fprintf("Elements greater than 10^12 in L^-1 and F difference: %d\n", length(find(compLinvF > 10^(-12))));

% (f)
A = ones(N-M, M);
y = s(M+1:end);
for i = 1:M
    A(:,i) = s(i:N-M-1+i); % Each column is r with i time lag
end
A = fliplr(A);
Apseudo = ((transpose(A)*A)^(-1))*transpose(A); %pseudo inverse
w = Apseudo*y; %solve for w
compWfpc = w + transpose((fpc(2:end))) % the error vector v, varies each time

figure;
plot(s);
hold on;
plot(rorig, 'r');
legend("s_t", "Original");
% You can tell there is unit root nonstationarity by the autocorrelation
% coefficients, there are no peaks, it stays high. The autocorrelation
% coefficents here vary greatly from the original.