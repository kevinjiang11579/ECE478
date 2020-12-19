% Q2 monthly  Kevin Jiang
clear;
Stonks = readmatrix('Stonks.csv');
SP500 = readmatrix('S&P500Data.csv');
libor = readmatrix('LIBOR.csv');
dates = readmatrix('Dates.csv');
monthEnds = zeros(1,11);
monthCount = 1;
for i = 2:length(dates)
    if (mod(dates(i),10000)-mod(dates(i-1),10000)) > 50
        monthEnds(monthCount) = i-1;
        monthCount = monthCount + 1;
    end
end
libor = libor(:,2);
Kstonksdaily = zeros(size(Stonks,1), size(Stonks,2));
Ksp500daily = zeros(size(SP500,1), size(SP500,2));
for i = 1:length(Stonks)
    Kstonksdaily(i,:) = (Stonks(i,:)-Stonks(1,:))./(Stonks(1,:));
    Ksp500daily(i,:) = (SP500(i,:)-SP500(1,:))./(SP500(1,:));
end
Kstonksdaily = transpose(Kstonksdaily(253:504,:));
Ksp500daily = transpose(Ksp500daily(253:504,:));
R2011 = libor(2)/(100*12);
Kstonks = zeros(5,12);
Kstonks(:, 1) = mean(Kstonksdaily(:,1:monthEnds(1)), 2);
Ksp500 = zeros(1,12);
Ksp500(:, 1) = mean(Ksp500daily(:,1:monthEnds(1)));
for i = 2:11
    Kstonks(:, i) = mean(Kstonksdaily(:,(monthEnds(i-1)+1):monthEnds(i)), 2);
    Ksp500(:, i) = mean(Ksp500daily(:,(monthEnds(i-1)+1):monthEnds(i)));
end
Kstonks(:, 12) = mean(Kstonksdaily(:,(monthEnds(11)+1):end), 2);
Ksp500(:, 12) = mean(Ksp500daily(:,(monthEnds(11)+1):end));
muSP = mean(Ksp500); % mean and variance for S&P500
varSP = var(Ksp500);
muStonks = zeros(1,5);
varStonks = zeros(1,5);
for i = 1:5
    muStonks(i) = mean(Kstonks(i,:)); %mean and variance for individual stocks
    varStonks(i) = var(Kstonks(i,:));
end
mStonks = (mean(Kstonks,2));
CStonks = cov(transpose(Kstonks));
[varMVPSt, muMVPSt, muMSt, varMSt, muAllSt, varAllSt, weightMVPSt, weightMSt] = efrontier(mStonks, CStonks, R2011);
figure;
plot(varAllSt.^0.5, muAllSt);
hold on;
plot(varSP^0.5, muSP, '*r', varStonks(1)^0.5, muStonks(1), '*b',...
    varStonks(2)^0.5, muStonks(2), '*g', varStonks(3)^0.5, muStonks(3), '*bl',...
    varStonks(4)^0.5, muStonks(4), '*c',varStonks(4)^0.5, muStonks(4), '*m');
legend("Efficient Frontier", "S&P 500", "AAPL", "AMZN", "ATVI", "LRCX", "TXN");
betas = zeros(1,5); %betas
divs = zeros(1,5); %diversifiable risk
syses = zeros(1,5); %systematic risk
for i = 1:5
    [divs(i), syses(i), betas(i)] = riskandbeta(Kstonks(i,:),Ksp500, muSP, muStonks(i));
end


function [div, sys, beta] = riskandbeta(Kv, Km, muM, muV)
    covM = cov(Km,Kv);
    beta = covM(1,2)/covM(1,1);
    eV = (Kv-muV)-beta*(Km-muM);
    div = (beta^2)*covM(1,1);
    sys = var(eV);
end

function [varMVP, muMVP, muM, varM, muAll, varAll, weightMVP, weightM] = efrontier(m, C, R)
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

end