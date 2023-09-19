function [SqI, PE, TS] = Sequentiality_Index(FR)
% FR es la matriz de respuesta, cada rengñon es una neurona y cada columna
% es un bin de tiempo.
[N, M] = size(FR);
FR(FR<1e-10) = 1e-10;

pj = mean(FR); % promedio de activación entre neuronas por bin
pj = pj / sum(pj); % normalización para asegurar pdf
PE = sum(-pj .* log(pj)) / log(M); % Peak entropy

rit = FR ./ repmat(sum(FR), N, 1); % normalización de respuesta por bin de tiempo
TS = 1-mean(sum(-rit .* log(rit)) / log(N)); % temporal sparcity


SqI = sqrt(PE * TS);