function [info,stderr] = kernelvineinfo_a(vines,copula,pcond,alpha,erreps,cases)
% KERNELVINEINFO Mixed copula vine mutual information estimate.
%   [INFO,STDERR] = KERNELVINEINFO(VINES,PCOND,ALPHA,ERREPS,CASES) returns
%   an estimate INFO of the mutual information between a discrete random
%   variable and a set of conditional mixed vines specified in the cell
%   VINES. Each element of this cell specifies as mixed vine as returned by
%   MIXEDVINEFIT. The probability of occurence of each vine is specified in
%   the vector PCOND which must have the same size as VINES.
%   ALPHA is the significance level of the estimate (default ALPHA = 0.05).
%   ERREPS is the maximum standard error of the mutual information estimate
%   as a stopping criterion (default ERREPS = 1e-3).
%   CASES is the number of samples that are drawn in each iteration of the
%   Monte Carlo estimation (default CASES = 1000).
%
%   INFO and STDERR are scalars in unit bit (base 2 logarithm).

% Argument checks
if nargin < 2
    error('kernelvineinfo: Usage [info,stderr] = kernelvineinfo(vines,pcond,alpha,erreps,cases)');
end
if ~iscell(vines)
    error('kernelvineinfo: Argument "vines" must be a cell array');
end
if ~isvector(pcond) || any(pcond<0) || any(pcond>1)
    error('kernelvineinfo: pcond must be a vector of probabilities');
end
if nargin < 3 || isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha)
    error('kernelvineinfo: Argument "alpha" must be a scalar');
end
if nargin < 4 || isempty(erreps)
    erreps = 1e-3;
elseif ~isscalar(erreps)
    error('kernelvineinfo: Argument "erreps" must be a scalar');
end
if nargin < 5
    cases = 1000;
elseif ~isscalar(cases)
    error('kernelvineinfo: Argument "cases" must be a scalar');
end

% Number of conditions
ncond = length(pcond);
% Probabilities of conditions
pcond = repmat(pcond(:)',cases,1);

% Gaussian confidence interval for erreps and level alpha
conf = norminv(1 - alpha,0,1);

stderr = inf;
h = 0;
varsum = 0;
k = 0;
x = zeros(cases,length(vines{1}.margins));
pcum = cumsum(pcond(:,end:-1:1),2);
pcum = pcum(:,end:-1:1);
while stderr >= erreps
    % Generate samples
    c = repmat(rand(cases,1),1,ncond);
    c = sum(c <= pcum,2);
    for i = 1:ncond
        sel = c==i;
        cases_sel = sum(sel);
        x(sel,:) = distrnd(vines{i},copula{i},cases_sel);
    end
    p = zeros(cases,1);
    for i = 1:ncond
        [kp,~,~] = Fit_vCopula(vines{i},x,size(copula{1},2),[],0,copula{i},'rand');
        p = p + pcond(1,i) * kp;
    end
    logp = log(p);
    log2p = logp(~isinf(logp)) / log(2);
    k = k + 1;
    % Monte-Carlo estimate of entropy
    h = h + (-mean(log2p) - h) / k;
    % Estimate standard error
    varsum = varsum + sum((-log2p - h) .^ 2);
    stderr = conf * sqrt(varsum / (k * cases * (k * cases - 1)));
    
    
    
    
    figure(100)
    hold on
    subplot(1,2,1);plot(k,h,'O')
    title('Info')
    hold on
    subplot(1,2,2);plot(k,stderr,'*')
    title('stderr')
    hold on
    drawnow
    
end

% Subtract conditional entropies
info = h;
for i = 1:ncond
    errepscond = erreps;%pcond(1,i) * erreps / 2;
    [hcond,stderrcond] = kernelvineentropy(vines{i},copula{i},alpha,errepscond,cases);
    info = info - pcond(1,i) * hcond;
    stderr = stderr + pcond(1,i) * stderrcond; 
end

end

function x = distrnd(vine,copula,cases)
y = zeros(numel(vine.margins{1}.ker),numel(vine.margins));
for hh = 1:numel(vine.margins)
    y(:,hh) = vine.margins{hh}.ker;
end
x = kerncoprnd(copula,y,cases);
end

function [h,stderr] = kernelvineentropy(vine,copula,alpha,erreps,cases)
% Gaussian confidence interval for erreps and level alpha
conf = norminv(1 - alpha,0,1);
stderr = inf;
h = 0;
varsum = 0;
k = 0;
while stderr >= erreps
    % Generate samples
    x = distrnd(vine,copula,cases);
    [kp,~,~] = Fit_vCopula(vine,x,size(copula,2),[],0,copula,'rand');
    logp = log(kp);
    log2p = logp(~isinf(logp)) / log(2);
    k = k + 1;
    % Monte-Carlo estimate of entropy
    h = h + (-mean(log2p) - h) / k;
    % Estimate standard error
    varsum = varsum + sum((-log2p - h) .^ 2);
    stderr = conf * sqrt(varsum / (k * cases * (k * cases - 1)));
    
    figure(100)
    hold on
    subplot(1,2,1);plot(k,h,'*')
    title('Info')
    hold on
    subplot(1,2,2);plot(k,stderr,'o')
    title('stderr')
    hold on
    drawnow    
    
end
end
