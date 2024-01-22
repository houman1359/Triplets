function bw=bw_tll(data,deg)

    n = size(data,1);
%     bw = 3 * n^(-1 / (4 * deg + 2)) * (chol(cov(qnorm(data))))';
    bw = 5 * n^(-1 / (4 * deg + 2)) * (chol(cov(data)))';   %%%% for the pc transformed data
    
    
