function [C,Mar]=kernelcdf(x,y,par)


if par.fit~=0
u=x*0;
for j=1:numel(x)
    u(j)=sum(x<=x(j))/(numel(x)+1);
    %u(j)=sum(x<=x(j))/(numel(x));
end

[uu nn]=unique(u); 
Mar.s=x(nn);
Mar.p=u(nn);
s=x(nn);
pp=u(nn);
else
    
s=par.s;
pp=par.p;
Mar=NaN; 


end

% pp=(pp-min(pp))*(par.max-par.min)/(max(pp)-min(pp))+par.min;
% 
% y(abs(y-max(s))<1e-10)=par.max-eps;%max(s);
% y(abs(y-min(s))<1e-10)=par.min+eps;%min(s);

% C=interp1(s,pp,y,'nearest');
% C=interp1(s,pp,y,'linear','extrap');


%%%%%%%%%%%%%%%%%%%%%%% removed this on 30 may 2018
% if par.fit==-1
% y=y+rand(size(y))*1e-10;
% end

% C=interp1(s,pp,y,'nearest');
C=interp1(s,pp,y,'linear');

rng(1234)
C(isnan(C) & y>max(s))=par.max*(1-1e-10*abs(rand(size(C(isnan(C) & y>max(s))))));%1-1/(numel(s)+1);
C(isnan(C) & y<min(s))=par.min*(1+1e-10*abs(rand(size(C(isnan(C) & y<min(s))))));

% C(isnan(C) & y<min(s))=par.min+eps;%1/(numel(s)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS IS USED FOR ENTROPY PAPER
% C(C>=1)=par.max-1e-10;%1-1/(numel(s)+1);
% C(C<=0)=par.min+1e-10;%1/(numel(s)+1);
%%%%%%%%%%%%%%%%%% I changed it to this on 7 july 2018
C(C>=par.max)=par.max*(1-1e-10*abs(rand(size(C(C>=par.max)))));%1-1/(numel(s)+1);
C(C<=par.min)=par.min*(1+1e-10*abs(rand(size(C(C<=par.min)))));%1-1/(numel(s)+1);

rng('shuffle')




