

function [intg]=densitynorm(vines,copula,erreps,cases,vineest)


intg=1;
return

V=prod(abs([min(vines.range')-max(vines.range')]))

alpha=0.05;

if isempty(erreps)
    erreps = 1e-3;
end

conf = norminv(1 - alpha,0,1);

stderr=inf;
intg=0;
varsum = 0;
k=0;


while stderr >= erreps && k<1000
    
    
    k = k + 1;
    clear X pp cc
    %     for hh=1:numel(vines.margins)
    %         X(:,hh)=vines.margins{hh}.ker;
    %         x=linspace(min(vines.range(hh,:)),max(vines.range(hh,:)),5000);
    %         cc(:,hh)=datasample(x,cases);
    %     end
    
    %     cc=kerncoprnd(copula,X,cases);
    
    cc = mixedvinernd(vineest,cases);
    
    
    [pp,~,~] = Fit_vCopula(vines,cc,size(copula,2),[],3,copula,'rand',[]);
%     p=p*V;
 
    q = mixedvinepdf(vineest,cc);

    q=(q+0.05/V)/1.05;
    
    p=pp./q;
    
    intg = intg + (mean(p)- intg)/k;
    
    varsum = varsum + sum(((p - intg)) .^ 2);
    stderr = conf * sqrt(varsum / (k * cases * (k * cases - 1)));
    
    
    figure(100)
    hold on
    subplot(1,2,1);plot(k,intg,'O')
    hold on
    subplot(1,2,2);plot(k,stderr,'*')
    drawnow
    
end


end