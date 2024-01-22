function [Dt dd_nonzero]=vinecopula_sample(vine,copula,tvec,NSAM,cases,nn)

% NSAM=NSAM*2;

% NSAM=50;
% Generate samples

Nz=ones(size(tvec));
ddd=-inf;h22=[];cc=[];dd0=0;


NW=0;
NN=0;
while ddd<NSAM
    clear X
    
    for hh=1:numel(vine.margins)
        X(:,hh)=vine.margins{hh}.ker;
    end
    
    cc=[];
    while isempty(cc)
    try
%     cc=kerncoprnd(copula,X,cases,vine,nn);
    cc=kerncoprnd(copula,X,cases,vine,-1);
    end
    end
    
    if dd0==0
        D=cc;
    else
        D=cat(1,D,cc);
    end
    
    [~,h22]=histc(D(:,nn),[tvec]);
    if min(hh)==0
        h22=h22+1;
    end
    
    %                                 h22{i}=h22{i}+1;
    for tt=1:numel(tvec)
        dd(tt)=sum((h22==tt));
    end
    
    NN=NN+1;
    
    if NW==0 
    dd_nonzero=1:numel(dd);%find(dd~=0);%
%     dd_nonzero=find(dd<20);%
    dd=dd(dd_nonzero);
    NW=1;
    tvec=tvec(dd_nonzero);
    end
    

%     if NN>10%NW==0 
% %     dd_nonzero=1:numel(dd);%find(dd~=0);%
%     dd_nonzero=find(dd>20);%
%     dd=dd(dd_nonzero);
%     NW=1;
%     tvec=tvec(dd_nonzero);
%     end    
    
    ddd=min(dd(1:end-1))
    dd0=1;
end



Dt=D;

% return  %%%% remove the following part on 7 august 2018. I think this is not sampling properly the time



[~,h22]=histc(D(:,nn),tvec);

for tt=1:numel(tvec)-1
    ffi=find(h22==tt);
    if tt==1
        Dt=D(randsample(ffi,NSAM),:);
    else
        Dt=cat(1,Dt,D(randsample(ffi,NSAM),:));
    end
end

