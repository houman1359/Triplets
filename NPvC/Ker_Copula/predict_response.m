function [y_mle,y_em,LL]=predict_response(PP,ygv,y)




% q1=sum(PP(:,:)*mean((diff(ygv))));
% q11=repmat(q1,numel(ygv),1);
% q2=PP./q11;

for ii=1:size(PP,2)
    
    [~,j]=max(smooth(PP(1:end,ii),4));
%     [~,j]=max(PP(1:end,ii));
    
    y_mle(ii)=ygv(j);
    
%     q=q2(:,ii);%this is f(n|B)
    dX=diff(ygv)';dX(numel(ygv))=0;
    q1=sum(dX.*PP(:,ii));
%     q1=trapz(ygv,PP(:,ii));

    if sum(q1)~=0        
    q=PP(:,ii)/q1;
%   y_em(ii)=nansum(q.*ygv'*mean(diff(ygv))); %%%% expectation value

    y_em(ii)=sum(dX.*q.*ygv'); %%%% expectation value
%     y_em(ii)=trapz(ygv,q.*ygv'); %%%% expectation value
    else  %%% I changed this on May2019 because i don;t think it is correct to assume uniformc distirbution when the distribution is zero everywhere
    q=PP(:,ii);%1/(max(ygv)-min(ygv))*ones(size(PP(:,ii)));
    y_em(ii)=y_mle(ii);
    end


%     LL(ii)=nansum(q'.*(ygv-y_em(ii)).^2*mean(diff(ygv))); %%%% VARIANCE of the distribution

LL(ii)=sum(dX'.*q'.*(ygv-y_em(ii)).^2); %%%% VARIANCE of the distribution
% LL(ii)=trapz(ygv,q'.*(ygv-y_em(ii)).^2); %%%% VARIANCE of the distribution
% 
%     if size(PP,2)==numel(y)
%         [~,h]=min(abs(y(ii)-ygv));
%         LL(ii)=log2(q(h));
%     else
%         LL=NaN;
%     end
    
end

% y_em(isnan(y_em))=0;
