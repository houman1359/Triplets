function [p,Mar]=kernelpdf(x,TYPE,y,par)


if par.fit~=0
    
    if TYPE==2 %| sum(x>-1e-10)==numel(x)
        % [bw,pp,s,c1]=kde(x,100,min([x;y]),max([x;y])+eps);
        X0=find(x<1e-10);
        X1=find(x>=1e-10);
        
        [~,in]=min(x(X0));
        X00=setdiff(X0,X0(in));
        if isfield(par,'MAX')
            if ~isempty(par.MAX)
                x(X00)=(par.MAX-2e-10)*rand(size(X00));
            end
        end
        
        %         x(x<1e-6)=1e-6;
        %         y(y<1e-6)=1e-6;
        
        if sum(x<0)==0 & numel(X0)>1 %& numel(X0)>1
            
            %             x(X0)=x(X0)*1e5;
            %             y(y<1e-9)=y(y<1e-9)*1e5;
            
            %             x=[x;par.MAX];
            
            %             [bw,pp,s,c]=kde((x),256,min((x)),max((x))+eps);
            %
            % %             [bw,pp,s,c]=kde(log2(x),128,min(log2(x)),max(log2(x))+eps);
            % %             s=2.^(s);
            %
            %             s0=[];s1=[];
            %             if numel(X0)>1
            %             s0=linspace(min(x(X0)),max(x(X0)),10);
            %             end
            %             if numel(X1)>1
            %             s1=s(s>max(x(X0)));
            % %             s1=s1(1:end-1);
            %             end
            %
            %             pp0=ones(size(s0));pp1=[];
            %             if numel(X0)>1
            %             pp0=pp0/trapz(s0,pp0)*numel(X0)/numel(x);
            %             end
            %             if numel(X1)>1
            %             pp1=pp(s>max(x(X0)));
            % %             pp1=pp1(1:end-1);
            %             pp1=pp1/trapz(s1,pp1)*numel(X1)/numel(x);
            %             end
            %
            %             SM=[];PPM=[];
            %             if numel(X1)>1 & numel(X0)>1
            %                 SM=[max(s0)+1e-6 par.MAX-1e-10];
            %                 PPM=zeros(size(SM));
            %             end
            %
            %
            %             s=[s0 SM s1];
            %             pp=[pp0';PPM';pp1]';
            %             dX=diff(s);dX(numel(s))=0;
            %             pp=pp/sum(pp.*dX);
            % %             pp=pp/trapz(s,pp);
            %
            
            x1=x(X1);
            x0=x(X0);
            s0=[];pp0=[];s1=[];pp1=[];
            if numel(X0)>1
                [bw,pp0,s0,c1]=kde(x0,32,min([x0]),max([x0])+eps);
                s0=s0(1:end);pp0=mean(pp0)*ones(size(pp0));
                [mi ni]=sort(s0);
                pp0=pp0-linspace(0,1e-9,numel(s0))';
                %             pp0=pp0/trapz(s0,pp0);                
            end
            SM=[];
            if numel(X0)>1 %& numel(X1)>1 
                %             SM=linspace(max(s0)+1e-7,min(s1)-1e-7,100);
                SM=[par.MAX-2e-10 par.MAX-1e-10];%[max(s0)+1e-6];% par.MAX-1e-10];
                if ~isempty(SM)
                PPM=[pp0(end) 0]';
                pp0=[pp0;PPM];
                end
            end
            
            s0=[s0 SM];
            if numel(X0)>1
                dX=diff(s0);dX(numel(s0))=0;
                pp0=pp0/sum(pp0(:).*dX(:));
            end
      
                  
            if numel(X1)>1
                [bw,pp1,s1,c1]=kde(x1,128,min([x1]),max([x1])+eps);
                %             pp1=pp1/trapz(s1,pp1);
                dX=diff(s1);dX(numel(s1))=0;
                pp1=pp1/sum(pp1.*dX');
            end
            
            
            s=([s0 s1]);
            pp=[pp0*numel(X0)/numel(x);pp1*numel(X1)/numel(x)]';
            dX=diff(s);dX(numel(s))=0;
            pp=pp/sum(pp.*dX);
            %             pp=pp/trapz(s,pp);
            
            
            %             s=([s0 s1]);
            %             pp=[pp0*numel(X0)/numel(x);pp1*numel(X1)/numel(x)]';
            %             pp=pp/trapz(s,pp);
            
            %             sv=linspace(min(s),max(s),1e8);
            %             pv=interp1(s,pp,sv,'linear');
            %
            %             pv=pv/trapz(sv,pv);
            %             s=sv;
            %             pp=pv;
            
        else
            [bw,pp,s,c1]=kde(x,128,min([x]),max([x])+eps); %%128
            dX=diff(s);dX(numel(s))=0;
            pp=pp/sum(pp.*dX');
            %             pp=pp/trapz(s,pp);
        end
        
    else
%         [bw,pp,s,c1]=kde(x,128,min([x])-eps,max([x])+eps); %64
        [bw,pp,s,c1]=kde(x,200,min([x]),max([x])+eps); %%changed august 2021 to have the same number of bins as the f(n|BC) from copula for the deviance analysis

        dX=diff(s);dX(numel(s))=0;
        pp=pp/sum(pp.*dX');
        %         pp=pp/trapz(s,pp);
    end
    
    % X=linspace(min(x)-eps,max(x)+eps,1000);
    % dX=mean(diff(X));
    % px=interp1(s,pp,X,'linear','extrap');
    
    
    %     dX=mean(abs(diff(s)));
    % px=interp1(s,pp,X,'linear','extrap');
    
    if numel(unique(s))~=1 & numel(unique(s))~=2
        % pp=pp/sum(pp*dX);
        % pp=pp/traps(s*pp);
    else
        % pp=pp/trapz(s,pp);    %%% this should be ok for later but i keep it like this to be consistent with Alice analysis. I will remove the sum later.
    end
    
    
    Mar.s=s;
    Mar.p=pp;
    
else
    
    s=par.s;
    pp=par.p;
    Mar=NaN;
    
end

y(abs(y-max(s))<1e-7)=max(s)-eps;  %%%% it used to be 1e-9 i changed it to -15 on 26 August 2019
y(abs(y-min(s))<1e-7)=min(s)+eps;


if numel(unique(s))~=1 & numel(unique(s))~=2
%     p=interp1(s,pp,y,'linear');
        p=interp1(s,pp,y,'nearest'); %%changed august 2021
elseif s(end)==y
    p=pp;
elseif s(end)~=y
    p=mean(pp)*ones(size(y));
end

 

p(y>max(s) | y<min(s))=0;

if TYPE==2
    if isfield(par,'MAX')
        if ~isempty(par.MAX)
            p(y<=par.MAX-2e-10)=pp(1)-linspace(0,1e-9,sum(y<=par.MAX-2e-10));
            p(y>par.MAX-2e-10 & y<=par.MAX-1e-10)=0;
        end
    end
end
 

return %%%% june 2017