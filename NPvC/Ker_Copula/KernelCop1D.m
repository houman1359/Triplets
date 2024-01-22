

function [bw]=KernelCop1D(knots,dat,poi,method,Lfit,fitt)


%%%%%    _u --> uv space
%%%%%    _S --> \Phi space
%%%%%    _X --> PC spacec

% method='fix','LL1','LL2',  then 'nn';



n=size(dat,1);
d=size(dat,2);
%%%%%%%%%
data.u=dat;
data.S=norminv(data.u,0,1);    %%%%norminv
% qrs=linsolve(bw,zdata');
% [COEFF,data.X,~] = pca(data.S);
data.X=data.S;
%%%%%%%%

% knots
bw=bw_tll(data.X,2);
bw(bw<1e-3)=1e-3;
bw=bw/10;
%%%%%%% fit the model

points.u=poi;
points.S=norminv(points.u);    %%%%norminv, qnorm
% points.X=points.S*(COEFF);
points.X=points.S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Grid=GRID_Bands(dat,knots,1);
GRID_X=Grid.X;
GRID_S=Grid.S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.u=double(data.u);
data.S=double(data.S);
data.X=double(data.X);

points.u=double(points.u);
points.S=double(points.S);
points.X=double(points.X);

xx=GRID_S(1:round(sqrt(size(GRID_X,1))),1);
xxd=[diff(xx)];xxd=[xxd;xxd(end)];
xd1=mean(xxd);

xd=xd1;

rng(1234);

samp=1:size(data.X,1);

afk = cvpartition(numel(samp),'KFold',5); % Stratified cross-validation

bw0=double(bw);
    
    if strcmp(method,'LL1') | strcmp(method,'nn') | strcmp(method,'LL2')
        for i=1:numel(afk.TestSize)
            dd{i}.u=data.u(samp(~afk.test(i)==1),:);
            dd{i}.S=data.S(samp(~afk.test(i)==1),:);
            dd{i}.X=data.X(samp(~afk.test(i)==1),:);
            ddT{i}.u=data.u(samp(afk.test(i)==1),:);
            ddT{i}.S=data.S(samp(afk.test(i)==1),:);
            ddT{i}.X=data.X(samp(afk.test(i)==1),:);
        end
        
        options=optimset('MaxIter',100,'Display','off','TolX',1e-3,'MaxFunEvals',100);
          
        if strcmp(method,'LL1')
            bw1=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0),[0],[5],options);
            bw=abs(bw0*bw1);
            disp('Opt=LL1')
        end
        
        if strcmp(method,'nn')            
            alpha=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,1),[0],[1],options);
            bw=dist_nn(alpha,data,Grid,bw0);            
            disp('Opt=nn')
        end
        
    else
        disp('Opt=fix')
    end
    
    
   
    


return

%%%%%% we can use the matlab package locfit by Brian Lau   https://github.com/brian-lau/locfit

%%%%%%
