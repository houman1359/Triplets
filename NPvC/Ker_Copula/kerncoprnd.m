
function [ura u] = kerncoprnd(copula,data,cases,vine,tim)

%%%%%%% this is from 5.15 and Algorithm 5.3 of the following book
% % % Simulating Copulas: Stochastic Models, Sampling Algorithms and Applications


rng('shuffle', 'v5uniform');
if size(copula,2)>1
[pts,GRID_u]= mk_grid(size(copula{1,2}.C_grid{1},1),vine.METH{1,2});
mag=max(GRID_u(:));
mig=min(GRID_u(:));

d = length(copula);
w = rand(cases,d);
% mign=mig+1e-10;
% magn=mag-1e-10;
% w=(magn-mign)*(w-min(w))./(max(w)-min(w))+mign;
w=(mag-mig)*((w-min(w))./(max(w)-min(w)))+mig;

copula1=copula;
copula2=copula;


%%%%%%%%%%% for r-vine gen look at 2019_Book_AnalyzingDependentDataWithVine.pdf

v = zeros(cases,d,d);
v(:,1,1) = reshape(w(:,1),[cases 1 1]);
v2 = zeros(cases,d,d);
v2(:,1,1) = reshape(w(:,1),[cases 1 1]);
for i = 2:d
    v(:,i,i) = reshape(w(:,i),[cases 1 1]);
%     v2(:,i,i) = reshape(w(:,i),[cases 1 1]);
    for k = (i-1):-1:1
        
        tr=k;
        str=i-k+1;
        method=vine.METH{tr,str};
        
        if ~strcmp(vine.families{tr,str},'ind')
            
            if strcmp(vine.type,'c-vine')
                copula1{tr,str}.C_grid=copula{tr,str}.C_grid{1};
                v(:,k,i)=kerncopccdfinv(copula1{tr,str},[v(:,k,k) v(:,k+1,i)],method);%copula{tr,str}.fit.F.ccdf([v(:,tr,1) v(:,tr,str)]);
            elseif strcmp(vine.type,'d-vine')
                copula1{tr,str}.C_grid=copula{tr,str}.C_grid{1};
                v(:,k,i)=kerncopccdfinv(copula1{tr,str},[v2(:,k,i-1) v(:,k+1,i)],method);%copula{tr,str}.fit.F.ccdf([v(:,tr,1) v(:,tr,str)]);
                if i<d
                copula2{tr,str}.C_grid=copula{tr,str}.C_grid{2}';
                v2(:,k+1,i)=kerncopccdf(copula2{tr,str},[v(:,k,i) v2(:,k,i-1)],method);%copula{tr,str}.fit.F.ccdf([v(:,tr,1) v(:,tr,str)]);
                end
            end
            
        else
%             wi = rand(cases,1);
%             wi=(magn-mign)*(wi-min(wi))./(max(wi)-min(wi))+mign;
            wi = rand(cases,1);
            wi=(mag-mig)*((wi-min(wi))./(max(wi)-min(wi)))+mig;

            v(:,k,i)=wi;
            
        end
        %         v(:,i,1) = copulaccdfinv(vine.families{k,i},[v(:,k,k) v(:,i,1)],vine.theta{k,i},1);
    end
    
    v2(:,1,i)=v(:,1,i);
    
    %     if i < d
    %         for j = 1:(i-1)
    %             v(:,i,j+1) = copulaccdf(vine.families{j,i},[v(:,j,j) v(:,i,j)],vine.theta{j,i},1);
    %         end
    %     end
end

u=squeeze(v(:,1,:));



%%%%%%%% this part is to make the sampling continous, specially for time
%%%%%%% domain
uu=u*0;
U=sort(unique(GRID_u(:)));
dU=diff(U);
for i=1:size(u,2)
    % [~,m]=histc(u(:,i)+1e-100,U);
    % m(m==numel(U))=numel(U)-1;
    % m(m==0 & u(:,i)<=mig)=1;
    % m(m==0 & u(:,i)>=mag)=numel(U)-1;
    % uu(:,i)=U(m)+dU(m).*rand(size(dU(m)));
    
    
    % [~,m]=histc(u(:,i)+1e-100,U);
    % m(m==numel(U))=numel(U)-1;
    % m(m==0 & u(:,i)<=mig)=1;
    % m(m==0 & u(:,i)>=mag)=numel(U)-1;
    % uu(:,i)=U(m)+dU(m).*rand(size(dU(m)));
    
    
    
    % U=sort(unique(u(:,i)));
    % dU=diff(U);
    % for y=1:size(u,1)
    % % uu(y,i)=u(y,i)+dU(m(y))*rand;
    % uu(y,i)=U(m(y))+dU(m(y))*rand;
    % end
    
    U=sort(unique(u(:,i)));
    for y=2:numel(U)
        G=find(u(:,i)==U(y));
        %     D=U(y+1)-U(y);
        D=U(y)-U(y-1);
        %     uu(G,i)=u(G,i)-D*rand(numel(G),1);
        uu(G,i)=U(y-1)+D*rand(numel(G),1);    %%%% changed 19 feb 2019
    end
    
end

if nargin>4
    if tim~=-1
        u(:,tim)=uu(:,tim);   %%%%%%%%%%%%% only to be used if the second dimension is time. I did not use this for the entropy paper ******** INPORTANT
    else
        u=uu;
    end
    
    % u=uu;  %%%%% may 2018, this is something to be considered
end
u=uu;
%%%%%% from CDF to samples

else
    
d = length(copula);
u = rand(cases,d);

end

for i=1:size(u,2)
    
    %%%%% below is the old, ~ used in entropy paper
    u(u(:,i)>=max(copula{i,1}.MarginG.p),i)=max(copula{i,1}.MarginG.p)-1e-10.*abs(rand(sum(u(:,i)>=max(copula{i,1}.MarginG.p)),1));
    u(u(:,i)<=min(copula{i,1}.MarginG.p),i)=min(copula{i,1}.MarginG.p)+1e-10.*abs(rand(sum(u(:,i)<=min(copula{i,1}.MarginG.p)),1));
    
    if strcmp(vine.families{i,1},'discrete')
        ura(:,i)=interp1(copula{i,1}.MarginG.p,copula{i,1}.MarginG.s,u(:,i),'linear');
    else
        ura(:,i)=interp1(copula{i,1}.MarginG.p,copula{i,1}.MarginG.s,u(:,i),'linear');
    end
    
end


end



