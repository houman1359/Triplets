
function [CCDF_points,pd_points,CCDF_grid,pd_grid,CCDF_data,pd_data,COPULA]=KernelCop(knots,dat,poi,method,method_fit,Lfit,fitt)


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
[COEFF,data.X,~,~,~,mu] = pca(data.S);
data.X=data.S * COEFF - repmat(mu,size(data.S,1),1) * COEFF;
if (abs(det(COEFF))-1)<1e-6
    data.X=data.X+1e-6 * rand(size(data.X));
end
%%%%%%%%

% knots
if fitt==1
    bw=bw_tll(data.X,2);
    bw(bw<1e-3)=1e-3;
    bw=[bw(1,1) bw(2,2)]/10;    
else
    if isfield(Lfit,'fitted')==1
        bw0=bw_tll(data.X,2);
        bw0(bw0<1e-3)=1e-3;
        bw0=[bw0(1,1) bw0(2,2)];
        Lfit.bw=bw0*Lfit.bw;
    end
    bw=Lfit.bw;
end
%%%%%%% fit the model

points.u=poi;
points.S=norminv(points.u,0,1);    %%%%norminv, qnorm
points.X=points.S * COEFF - repmat(mu,size(points.S,1),1) * COEFF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(method)
    Grid=GRID_Bands(dat,knots,method);
else
    Grid=GRID_Bands(dat,knots,method);
%     method=method_fit;
end
GRID_X=Grid.X;
GRID_S=Grid.S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.u=double(data.u);
data.S=double(data.S);
data.X=double(data.X);

points.u=double(points.u);
points.S=double(points.S);
points.X=double(points.X);

[x1,~,s1]=unique(sort(Grid.S(:,1)));
xd11=[diff(x1)];xd11=[xd11;xd11(1)];
[x2,~,s2]=unique(sort(Grid.S(:,2)));
xd22=[diff(x2)];xd22=[xd22;xd22(1)];
xd1=nan(size(Grid.u,1),1);
xd2=nan(size(Grid.u,1),1);
xd=xd1;
for k=1:numel(x1)
    xd1(Grid.S(:,1)==x1(k))=xd11(k);
    xd2(Grid.S(:,2)==x2(k))=xd22(k);
end

for k1=1:numel(x1)
    for k2=1:numel(x2)
        xd(Grid.S(:,1)==x1(k1) & Grid.S(:,2)==x2(k2))=xd11(k1)*xd22(k2);
    end
end

if numel(uniquetol(xd(:),1e-3))==1
    xd=uniquetol(xd(:),1e-3);
end

P1=normpdf(x1);
P2=normpdf(x2);
NORM=(P1*P2')';

% NN=min([20000,size(data.X,1)]);
% samp=randsample(1:size(data.X,1),NN);
% if size(data.X,1)<20000
    samp=1:size(data.X,1);
% end

afk = cvpartition(numel(samp),'KFold',5); % 5 was used, changed 11 may 2018  % Stratified cross-validation

bw0=double(bw);

if fitt==1
    
    if isfield(Lfit,'fitted')==0
        if strcmp(method_fit,'LL1') | strcmp(method_fit,'nn') | strcmp(method_fit,'LL2')
            for i=1:numel(afk.TestSize)
                dd{i}.u=data.u(samp(~afk.test(i)==1),:);
                dd{i}.S=data.S(samp(~afk.test(i)==1),:);
                dd{i}.X=data.X(samp(~afk.test(i)==1),:);
                ddT{i}.u=data.u(samp(afk.test(i)==1),:);
                ddT{i}.S=data.S(samp(afk.test(i)==1),:);
                ddT{i}.X=data.X(samp(afk.test(i)==1),:);
            end
            
            %         options=optimset('MaxIter',150,'Display','off','TolX',1e-4,'MaxFunEvals',150);
            options=optimset('MaxIter',100,'Display','off','TolX',1e-1,'TolFun',1e-1,'MaxFunEvals',100);
            
            if strcmp(method_fit,'LL2')
                %             bw1=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM),[1e-10],[5],options);
                %             bw0=abs(bw0*bw1);
                options=optimset('MaxIter',50,'Display','off','TolX',1e-1,'TolFun',1e-1,'MaxFunEvals',50);
                bw1=fminsearchbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,1),bw0,[1e-4 1e-4],[2 2],options);
                options=optimset('MaxIter',50,'Display','off','TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50);
                bw1=fminsearchbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,0),bw0,bw1*0.8,bw1*1.2,options);
                bw=abs(bw1);
                disp('Opt=LL2')
            end
            
            if strcmp(method_fit,'LL1')
                %             options=optimset('MaxIter',50,'Display','off','TolX',1e-1,'TolFun',1e-1,'MaxFunEvals',50);
                %             bw1=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,1),[1e-4],[2],options);
                %             options=optimset('MaxIter',50,'Display','off','TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50);
                %             bw1=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,0),[bw1*0.8],[bw1*1.2],options);
                %              bw=abs(bw0*bw1);
                %             disp('Opt=LL1')
                options=optimset('MaxIter',100,'Display','off','TolX',1e-1,'TolFun',1e-1,'MaxFunEvals',100);
                bw1=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,1),[1e-6],[5],options);
                options=optimset('MaxIter',100,'Display','off','TolX',1e-3,'TolFun',1e-2,'MaxFunEvals',100);
                bw1=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,0),[bw1*0.8],[bw1*1.2],options);
                bw=abs(bw0*bw1);
                %             disp('Opt=LL1')
            end
            
            if strcmp(method_fit,'nn')
                
                error('the analytic LL works only for constant bw. for this we have to use non-analytic version if we want')
                alpha=fminbnd(@(x) MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,1,NORM),[0],[1],options);
                bw=dist_nn(alpha,data,Grid,bw0);
                
                disp('Opt=nn')
            end
            
        else
            disp('Opt=fix')
        end
    else
        
                bw=bw0*Lfit.bw*10;
        
    end
%     
%     tic
%     nn=linspace(1e-4,1,100);
%     for ii=1:100
%         ii
%     err(ii,:)=MISE(nn(ii),data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM);
%     end
%     toc
%     figure;plot(nn*bw0(1),err,'-O')
%     [tt t]=min(err);bw=abs(bw0*nn(t))
%      
     
    lfit = loclik_fit(bw,data,Grid);
    [pd_data,ccdf_data,ccdf2_data,pd_grid,ccdf_grid,ccdf2_grid,pd_points,ccdf_points,ccdf2_points,F,norma]=func_tll(lfit,Grid,points.S,data,1,0,NORM);

    fit.F=F;
    fit.bw=bw;
    fit.norm=norma;
    COPULA.data=data;
    COPULA.fit=fit;

else
    
    if fitt==-1        
        lfit = loclik_fit(bw,data,Grid);
        fit.lfit=[];
    end
    [pd_data,ccdf_data,ccdf2_data,pd_grid,ccdf_grid,ccdf2_grid,pd_points,ccdf_points,ccdf2_points,Fn,~]=func_tll(Lfit,Grid,points.S,data,fitt,0,NORM);
    
    if fitt==0
        COPULA.fit=0;
        COPULA.Grid=0;
        COPULA.data=0;
        COPULA.points=0;
    elseif fitt==-1
        COPULA=Lfit;
        fit.F=[];
        fit.F=Fn;
        fit.bw=Lfit.bw;
        COPULA.data=data;
        COPULA.fit=fit;
    end
    
end

CCDF_data{1}=ccdf_data;
CCDF_data{2}=ccdf2_data;
CCDF_points{1}=ccdf_points;
CCDF_points{2}=ccdf2_points;
CCDF_grid{1}=ccdf_grid;
CCDF_grid{2}=ccdf2_grid;



return
