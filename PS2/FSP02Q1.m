% Q1   Kevin Jiang
clear;
FarmaFrench = readmatrix('FarmaFrench.csv');
libor = readmatrix('LIBOR.csv');
libor = libor(:,2);
FF2011 = FarmaFrench(253:504,2:end); %use entries in year 2011, discard the date
K = transpose(FF2011);
m = (mean(K,2));
C = cov(transpose(K));
R2011 = libor(2)/(100*360);
[varMVP, muMVP, muM, varM, muAll48, varAll48, varMVPns, muMVPns,...
    varMns, muMns, weightMVP, weightM, weightMVPns, weightMns] = efrontier(m, C, R2011);
figure;
plot(varAll48.^0.5, muAll48);
hold on;
plot([0, varM].^0.5, [R2011, muM]);
hold on;
plot(varMVPns, muMVPns, '*b', varMns, muMns, '*r');
K2 = transpose([FF2011(:, 10), FF2011(:, 30)]);
m2 = (mean(K2,2));
C2 = cov(transpose(K2));
[varMVP2, muMVP2, muM2, varM2, muAll2, varAll2, ~, ~, ~, ~, weightMVP2, weightM2, ~, ~] = efrontier(m2, C2, R2011);
hold on;
plot(varAll2.^0.5, muAll2);
legend("48 portfolios","Capital Market line","MVP approx.","M approx.","2 portfolios");
title("Efficient frontier and CML");
K_MVP = transpose(weightMVP)*K;
K_MVPns = transpose(weightMVPns)*K; %MVP approximation
K_Mns = transpose(weightMns)*K; %M approximation
K_MVP2 = transpose(weightMVP2)*K2;
K_M2 = transpose(weightM2)*K2;
K_M = transpose(weightM)*K;
covMandMVP = cov(K_M, K_MVP);
covMandMVPns = cov(K_M, K_MVPns);
covMandMns = cov(K_M, K_Mns);
covMandMVP2 = cov(K_M, K_MVP2);
covMandM2 = cov(K_M, K_M2);
corrMandMVP = getcorr(covMandMVP);
corrMandMVPns = getcorr(covMandMVPns);
corrMandMns = getcorr(covMandMns);
corrMandMVP2 = getcorr(covMandMVP2);
corrMandM2 = getcorr(covMandM2);
betaMVP = getbeta(covMandMVP);
betaMVPns = getbeta(covMandMVPns);
betaMns = getbeta(covMandMns);
betaMVP2 = getbeta(covMandMVP2);
betaM2 = getbeta(covMandM2);
[divMVP, sysMVP] = risks(K_MVP, K_M, betaMVP, muMVP, muM, covMandMVP);
[divMVPns, sysMVPns] = risks(K_MVPns, K_M, betaMVPns, muMVPns, muM, covMandMVPns);
[divMns, sysMns] = risks(K_Mns, K_M, betaMns, muMns, muM, covMandMns);
[divMVP2, sysMVP2] = risks(K_MVP2, K_M, betaMVP2, muMVP2, muM, covMandMVP2);
[divM2, sysM2] = risks(K_M2, K_M, betaM2, muM2, muM, covMandM2);
% The point that acheives maximum beta should be the market portfolio.
% Covariance between K_m and K_m should be the max, and beta is 1 in this
% case.
% a)
eigVals = eig(C); % All positive
% b)
cond = max(eigVals)/min(eigVals);
% c)
figure;
semilogy(1:48, flip(eigVals));
title("Log scale plot of eigenvalues of C");
%There is a sudden drop from the largest eigenvalue to the second largest,
%which matches with the large condition number


function [div, sys] = risks(Kv, Km, beta, muV, muM, covMV)
    eV = (Kv-muV)-beta*(Km-muM);
    div = (beta^2)*covMV(1,1);
    sys = var(eV);
end
function [beta] = getbeta(covM)
    beta = covM(1,2)/covM(1,1);
end

function [corr] = getcorr(covM)
    corr = covM(1,2)/sqrt(covM(1,1)*covM(2,2));
end

function [varMVP, muMVP, muM, varM, muAll, varAll, varMVPns, muMVPns,...
    varMns, muMns, weightMVP, weightM, weightMVPns, weightMns] = efrontier(m, C, R)
    % Find weightMVP and weightM
    onevector = ones(length(m),1);
    varMVP = 1/(transpose(onevector)*((C^-1)*onevector));
    weightMVP = varMVP*((C^-1)*onevector);
    muMVP = transpose(m)*weightMVP;
    mex = m - R.*onevector;
    weightM = (1/(transpose(onevector)*((C^-1)*mex)))*((C^-1)*mex);
    muM = transpose(m)*weightM;
    varM = transpose(weightM)*C*weightM;
    % Get efficient frontier
    m_tilda = [m, onevector];
    B = transpose(m_tilda)*((C^-1)*m_tilda);
    u = 0.01:0.01:0.99;
    muV = zeros(1, length(u));
    varV = zeros(1, length(u));
    for i = 1:length(u)
        weight = u(i)*weightMVP + (1-u(i))*weightM;
        mu_tilda = transpose(m_tilda)*weight;
        muV(i) = mu_tilda(1);
        varV(i) = transpose(mu_tilda)*(B^-1)*transpose(m_tilda)*(C^-1)*m_tilda*(B^-1)*mu_tilda;
    end
    muAll = [muM, muV, muMVP];
    varAll = [varM, varV, varMVP];
    weightMVPns = weightMVP;
    weightMVPns(find(weightMVPns<0)) = 0;
    weightMVPns = weightMVPns./sum(weightMVPns);
    weightMns = weightM;
    weightMns(find(weightMns<0)) = 0;
    weightMns = weightMns./sum(weightMns);
    muMVPns = transpose(m)*weightMVPns;
    varMVPns = transpose(weightMVPns)*C*weightMVPns;
    muMns = transpose(m)*weightMns;
    varMns = transpose(weightMns)*C*weightMns;
end