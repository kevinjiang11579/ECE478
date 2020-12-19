clear;
close all;
%Q5 with normal v
% (a)
% Multiplying both sides of the first equation by A
% Get delta' = A*delta, G' = A*G*A^-1, since [a;b] = (A^-1)*[r1,r2]
d = [0.3;0.5];
G = [0.99,0;0.3,0.3];
A = [3,4;2,3];
d_tild = A*d;
G_tild = A*G*(A^(-1));
% (b)
N = 250;
hb = [1, 0.5, 0.06];
ha = [1, -0.7, 0.1];
a = zeros(N,1);
wn = randn(N,1);
u = filter(hb,ha,wn);
a(1) = d(1) + u(1);
for i = 2:N
    a(i) = d(1) + G(1)*a(i-1) + u(i);
end
var_0 = mean(a.^2)/16;
v = randn(N,1);
v = v*sqrt(var_0/var(v));
b = zeros(N,1);
b(1) = d(2) + v(1);
for i = 2:N
    temp = G*[a(i-1);b(i-1)];
    b(i) = d(2) + temp(2) + v(i);
end
ab = transpose([a,b]);
r1r2 = A*ab;
r1 = r1r2(1,:);
r2 = r1r2(2,:);
% (c)
figure;
subplot(2,1,1);
plot(u);
title("u_t");
subplot(2,1,2);
plot(v);
title("v_t");
figure;
subplot(2,1,1);
plot(1:1:N, a, 1:1:N, b);
title("a_t and b_t");
legend("a_t", "b_t");
subplot(2,1,2);
plot(1:1:N, r1, 1:1:N, r2);
title("r1_t and r2_t");
legend("r1_t", "r2_t");
% (d)
% bt is the only one that is stationary
% We know [a;b] = A^(-1)*[r1;r2]
% So coefficients to get stationary bt is bottom row of A^(-1);
Ainv = A^(-1);
stat_coeff = Ainv(2,:);
% (e)
M = 20;
r1corr = xcorr(r1, M, 'normalized');
r2corr = xcorr(r2, M, 'normalized');
[r1r2xc,r1r2lags] = xcorr(r1,r2,M,'normalized');
figure;
subplot(2,1,1);
stem(r1corr(21:end));
title("Autocorrelation coefficents of r1");
subplot(2,1,2);
stem(r2corr(21:end));
title("Autocorrelation coefficents of r2");
figure;
stem(r1r2lags, r1r2xc);
title("Cross-correlation coefficents of r1 and r2");
% (f)
% The correlation decays very fast but also goes negative...
s1 = r1 - [0, r1(1:end-1)];
s2 = r2 - [0, r2(1:end-1)];
s1corr = xcorr(s1, M, 'normalized');
s2corr = xcorr(s2, M, 'normalized');
[s1s2xc,s1s2lags] = xcorr(s1,s2,M,'normalized');
figure;
subplot(2,1,1);
stem(s1corr(21:end));
title("Autocorrelation coefficents of s1");
subplot(2,1,2);
stem(s2corr(21:end));
title("Autocorrelation coefficents of s2");
figure;
stem(s1s2lags, s1s2xc);
title("Cross-correlation coefficents of s1 and s2");
% (g)
% magnitude of cross correlation is signficantly less, however overall
% slope varies
wn2 = randn(1,N);
c = zeros(1,N);
c(1) = wn2(1);
for i = 2:N
    c(i) = 0.99*c(i-1) + wn2(i);
end
ccorr = xcorr(c,M,'normalized');
cr1corr = xcorr(c,r1,M,'normalized');
cr2corr = xcorr(c,r2,M,'normalized');
figure;
stem(ccorr(21:end));
title("Autocorrelation coefficients of c");
figure;
subplot(2,1,1);
stem(s1s2lags, cr1corr);
title("Cross-correlation coefficents of c and r1");
subplot(2,1,2);
stem(s1s2lags, cr2corr);
title("Cross-correlation coefficents of c and r2");