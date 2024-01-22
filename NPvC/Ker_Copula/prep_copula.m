



function [vine,X_new,ord]=prep_copula(X,margins,families,iscons,type_c,sort_n,range,tim)

%%%% tim is the number of the time variable, if no time exist it can be nan

% X  t x d matrix, with the first column being the neuron
% margins    1 x d     in this version can be only 'kernel' but easily we
% add other parmeteric

% families 'ker' or others;
% iscont    for 'ker'-->%  1 = binary , 2 = positive or zero
% it is different from continous or descrete one


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d= numel(iscons);
vine.margins = cell(d,1);
vine.type = type_c;%'c-vine';'nc-vine';'d-vine'; % Canonical vine type

vine.families = cell(d);
vine.theta = cell(d);

if strcmp(sort_n,'sort')
for i=1:d
    for j=1:d
        C(i,j)=abs(corr(X(:,i),X(:,j),'type','kendall'));
    end
end

clear ord
% [~,f2]=max(mean(C));
ord(1)=1;
for i=2:d-1
    ss=setdiff(1:d,ord);
    [~,n]=max(C(ord(i-1),ss));
    ord(i)=ss(n);
end
ord(d)=setdiff(1:d,ord);
else
ord=1:size(X,2);
end


X_new=X(:,ord);


for i=1:d-1
    for j=2:d-i+1
        vine.families{i,j} = families{i,j};
        vine.families{j,i} = families{i,j};
%         if i==1 && j==tim
%         vine.METH{i,j}=[1 0];
%         else
%         vine.METH{i,j}=[1 1];
%         end
    end
end



for dd=1:d

    RA=rand(size(X_new(:,1)));

    vine.margins{dd}.dist = margins{dd};
    vine.margins{dd}.theta = [];
    vine.margins{dd}.iscont = iscons(dd); %  1 = binary , 2 = positive or zero
        
    
    if strcmp(families{2,dd},'kercop') | isempty(families{2,dd})
    
%     mm=max(X_new(:,dd));    
%     vine.margins{dd}.ker = X_new(:,dd)+abs(1e-10*rand(size(X_new(:,dd))));
%     vine.margins{dd}.ker(vine.margins{dd}.ker>mm)=X_new(vine.margins{dd}.ker>mm,dd)+abs(1e-10*rand(size(X_new(vine.margins{dd}.ker>mm,dd))));

    [nn xx]=unique(X_new(:,dd));
    a=sort(unique(nn));
    ad=diff(a);
    ad(numel(ad)+1)=ad(end);
    Di=X_new(:,dd);
    
    for i=a'
        Di(X_new(:,dd)==i)=ad(find(a==i));
    end    
    vine.margins{dd}.ker = X_new(:,dd)+Di.*rand(size(X_new(:,dd)))*1e-10;%    CHANGES ON MAY 2019
    
    %%%%%%%%%% I copied the following part on july 2022 from the 2d package
    


  %  vine.margins{dd}.ker = X_new(:,dd)+Di.*rand(size(X_new(:,dd))); 
    
    
    
    elseif strcmp(families{2,dd},'discrete')
        
    [nn xx]=unique(X_new(:,dd));    
    a=sort(unique(nn));ad=diff(a);ad(numel(ad)+1)=1e-10;Di=X_new(:,dd);
    for i=a'
    Di(X_new(:,dd)==i)=ad(find(a==i));
    end
    vine.margins{dd}.ker = X_new(:,dd)+Di.* RA;
    end
    
%     vine.margins{dd}.ker = X_new(:,dd)+abs(1e-10*rand(size(X_new(:,dd))));
    vine.theta{dd,dd} = vine.margins{dd}.ker;

%     vine.range(1,dd)=max(vine.margins{dd}.ker);range(2,dd)=max(vine.margins{dd}.ker);

end
vine.range=range;

