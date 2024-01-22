

function [ccdf_points,pd_points,ccdf_grid,pd_grid,ccdf_data,pd_data,COPULA]=KernelCop3D(knots,dat,poi,method,Lfit,fitt,bw2)


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
[COEFF,data.X,~] = pca(data.S(:,1:2));
data.X(:,3)=data.S(:,3);
%%%%%%%%

if fitt==1
bw=bw_tll(data.X(:,1:2),2);
bw(bw<1e-3)=1e-3;
bw=[bw(1,1) bw(2,2)];
bw=bw/10;
else
% bw=Lfit.bw;
end
%%%%%%% fit the model
% bw(1:2)=bw2;
% bw(3)=0.5;
bw=bw2;

% bw(3)=bw(3)*(size(data.u,1)^(1/90));  %%%%% MODIFY it later
%%%%%% i have to transform the bandwidth according to the dimensional
%%%%%% change 2-->3 and 1-->3.








points.u=poi;
points.S=norminv(points.u);    %%%%norminv, qnorm
points.X(:,1:2)=points.S(:,1:2)*(COEFF);
points.X(:,3)=points.S(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(method)
Grid=GRID_Bands(dat,knots,[1 1]);
else
Grid=GRID_Bands(dat,knots,method);
method='LL1';
end
GRID_X=Grid.X;
GRID_S=Grid.S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_lim=[min(GRID_X(:,1)) max(GRID_X(:,1))];
Y_lim=[min(GRID_X(:,2)) max(GRID_X(:,2))];

[G1,G2]=meshgrid(linspace(Y_lim(1),Y_lim(2),knots),linspace(X_lim(1),X_lim(2),knots));
GRID_X_mesh=[G2(:) G1(:)];

data.u=double(data.u);
data.S=double(data.S);
data.X=double(data.X);

points.u=double(points.u);
points.S=double(points.S);
points.X=double(points.X);


xx=GRID_S(1:round(sqrt(size(GRID_X,1))),1);
xxd=[diff(xx)];xxd=[xxd;xxd(end)];
xd1=mean(xxd);
xx=GRID_S(1:round(sqrt(size(GRID_X,1))):end,2);
xxd=[diff(xx)];xxd=[xxd;xxd(end)];
xd2=mean(xxd);


xd=xd1*xd2;


rng(1234);

% NN=min([10000,size(data.X,1)]);
% samp=randsample(1:size(data.X,1),NN);

samp=1:size(data.X,1);

afk = cvpartition(numel(samp),'KFold',5); % Stratified cross-validation


bw0=double(bw);

if fitt==1
    
    if strcmp(method,'LL1') | strcmp(method,'nn') | strcmp(method,'LL2')
        for i=1:numel(afk.TestSize)
            dd{i}.u=data.u(samp(~afk.test(i)==1),:);
            dd{i}.S=data.S(samp(~afk.test(i)==1),:);
            dd{i}.X=data.X(samp(~afk.test(i)==1),:);
            ddT{i}.u=data.u(samp(afk.test(i)==1),:);
            ddT{i}.S=data.S(samp(afk.test(i)==1),:);
            ddT{i}.X=data.X(samp(afk.test(i)==1),:);
        end
        
        options=optimset('MaxIter',50,'Display','off','TolX',1e-3,'MaxFunEvals',100);
          
        if strcmp(method,'LL2')
            bw1=fminsearchbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0),bw0,[0 0],[10 10],options);
            % bw1=fminsearch(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0),bw0,options);
            bw=abs(bw1);
            disp('Opt=LL2')
        end
        
        if strcmp(method,'LL1')
            bw1=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0),[0],[5],options);
            % bw1=fminsearch(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0),1,options);
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
    
    
    lfit = loclik_fit(bw,data,Grid);
    [pd_data,ccdf_data,pd_grid,ccdf_grid,pd_points,ccdf_points,F,norma]=func_tll(lfit,Grid,points.u,data,1,0);
    
    fit.F=F;
    fit.bw=bw;
    fit.norm=norma;
    COPULA.data=data;
    COPULA.fit=fit;                  
    
else
    
    lfit = loclik_fit(bw,data,Grid);    
    [pd_data,ccdf_data,pd_grid,ccdf_grid,pd_points,ccdf_points,F,norma]=func_tll(lfit,Grid,points.u,data,fitt,0);
    
    fit.F=F;
    fit.bw=bw;
    fit.norm=norma;
    COPULA.data=data;
    COPULA.fit=fit;                  
    
end

return

%%%%%% we can use the matlab package locfit by Brian Lau   https://github.com/brian-lau/locfit

%%%%%%
