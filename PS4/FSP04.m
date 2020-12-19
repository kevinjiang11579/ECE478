%1
N = 1e6;
nsamples = randn(1,N);
csamples = trnd(1,1,N);
tfsamples = trnd(5,1,N);
ttsamples = trnd(10,1,N);
tfsamples = tfsamples/sqrt(var(tfsamples));
ttsamples = ttsamples/sqrt(var(ttsamples));
ngt4 = length(find(abs(nsamples)>4))/N;
cgt4 = length(find(abs(csamples)>4))/N;
tfgt4 = length(find(abs(tfsamples)>4))/N;
ttgt4 = length(find(abs(ttsamples)>4))/N;
figure;
plot(1:1:N, csamples);
title("Cauchy Random Samples");
