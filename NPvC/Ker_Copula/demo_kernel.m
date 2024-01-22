
clear all

addpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Selmaan_Switch/glmnet'))
addpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno'))


system('module load stats/R/3.3.1')

disp('Constructing mixed copula vine...');
d = 3; % Dimension

global vine;
vine.type = 'c-vine'; % Canonical vine type

%% Test probability density function

% Set margins
% vine.margins = cell(d,1);
% Standard normal margin
% vine.margins{1}.dist = 'norm';
% vine.margins{1}.theta = [0;1];
% vine.margins{1}.iscont = true; % Continuous margin
% Gamma margin
% vine.margins{2}.dist = 'norm';'gam';
% vine.margins{2}.theta = [0;1];[2;4];
% vine.margins{2}.iscont = true; % Continuous margin% Set margins
% vine.margins = cell(d,1);

% Standard normal margin
% vine.margins{3}.dist = 'norm';
% vine.margins{3}.theta = [0;1];
% vine.margins{3}.iscont = true; % Continuous margin
% Gamma margin
% vine.margins{4}.dist = 'norm';'gam';
% vine.margins{4}.theta = [0;1];[2;4];
% vine.margins{4}.iscont = true; % Continuous margin

% Poisson margin
% vine.margins{3}.dist = 'norm';'poiss';
% vine.margins{3}.theta = [0;1];10;
% vine.margins{3}.iscont = true; % Discrete margin

clear x
x(:,1)=rand(2000,1);%[0:0.01:1 0:0.01:1 0:0.01:1 0:0.01:1];
x(1:1000,2)=x(1:1000,1)+rand(1000,1)/20;%+rand(100,1)/4;
x(1001:2000,2)=1-x(1001:2000,1)+rand(1000,1)/20;%+rand(100,1)/4;

DDD=500*abs(randn(2000,1));
X(:,1) = DDD+randn(2000,1)*10;
X(1:1000,2) = 1*X(1:1000,1)+randn(1000,1)*10;
% X(1001:2000,2) = 1500-1*(X(1001:2000,1)+randn(1000,1)*10);
X(1001:2000,2) = 1*(X(1001:2000,1)+randn(1000,1)*10);
X=X+1;
x=X/500;


% x(:,1)=x(:,2);
f=x;

% x(:,2)=real(x(:,2).^(2))/1.5;
x(:,1)=real(x(:,1)*1);


X1 = x2fx(f(1:2000,[1 2])/1,'interaction');
X1(:,1)=1;

X2 = x2fx(x(1:2000,[1 2])/1,'interaction');
X2(:,1)=1;

%%% mu = exp(X2(:,[1 2 3 4])*[1;0.1;0.1;0.1] + 1);     %%%%%%%%% the x case
mu(1,find(x(:,1)>median(x(:,1)))) = (X2(find(x(:,1)>median(x(:,1))),[1 2 3 4])*[1;1;0;0.0] + 1)+(X2(find(x(:,1)>median(x(:,1))),[1 2 3 4])*[1;0;2;0.0] + 1);    
mu(1,find(x(:,1)<=median(x(:,1)))) = (X2(find(x(:,1)<=median(x(:,1))),[1 2 3 4])*[1;-2;0;0.0] + 1)+ (X2(find(x(:,1)<=median(x(:,1))),[1 2 3 4])*[1;0;-2;0.0] + 1);    




% M(:,1)=exp(X1(:,[2])*[0.1] + 1);   
% M(:,2)=exp(X1(:,[3])*[0.1] + 1);   

y = poissrnd(mu);


% y(find(x(:,1)>median(x(:,1))))=round(exp(2-x(find(x(:,1)>median(x(:,1))),1)))+round(exp(2+x(find(x(:,1)>median(x(:,1))),2)));
% y(find(x(:,1)<=median(x(:,1))))=round(exp(2+x(find(x(:,1)<=median(x(:,1))),1)))+round(exp(2-x(find(x(:,1)>median(x(:,1))),2)));

% y(find(x(:,1)>median(x(:,1))))=round(exp(2-x(find(x(:,1)>median(x(:,1))),1)));
% y(find(x(:,1)<=median(x(:,1))))=round(exp(2+x(find(x(:,1)<=median(x(:,1))),1)));

% y(150)=700; 

% 
% clear A1 A2
% aa=0;
% for alp=eps:0.1:1-eps
% aa=aa+1;
% [B FitInfo] = lassoglm(X1(:,[1 2 3 4]),y,'poisson','CV',10,'Alpha',alp);
% A1(:,aa)=B(:,FitInfo.Index1SE);
% A2(:,aa)=B(:,FitInfo.IndexMinDeviance);
% end
% 
% figure;plot(A1(1,:),'b');hold on;plot(A1(2,:),'r');hold on;plot(A1(3,:),'g');hold on;plot(A1(4,:),'k')
% figure;plot(A2(1,:),'b');hold on;plot(A2(2,:),'r');hold on;plot(A2(3,:),'g');hold on;plot(A2(4,:),'k')


opt = glmnetSet;
opt.alpha = 1/2;
opt.thresh = 1e-6;
GG = cvglmnet(X1(:,[1 2 3 4]),y,'poisson',opt,'deviance',[],[]);
            

ind1se = find(GG.lambda == GG.lambda_1se);
beta = GG.glmnet_fit.beta(:,ind1se);
prediction1se = cvglmnetPredict(GG,X1(:,[1 2 3 4]),GG.lambda(ind1se),'response');
% prediction = cvglmnetPredict(GG,X1(:,[1 2 3 4]),[],'response');% is the same as above
predictionnull = cvglmnetPredict(GG,X1(:,[1 2 3 4]),GG.lambda(1),'response');

nullDev = GG.cvm(1);
cvDev = GG.cvm(ind1se);
% rrr=corr(prediction,y)
rrr1=1-sum(abs(prediction1se'-y).^2)/sum(abs(predictionnull'-y).^2)
rrr2 = 1 - [cvDev]./[nullDev]


% x(:,3)=x(:,1);%rand(1000,1)/20;%[0:0.01:1 0:0.01:1 0:0.01:1 0:0.01:1];%rand(1000,1);
% x(1:500,4)=x(1:500,2);%rand(1000,1)/20;%[0:0.01:1 0:0.01:1 0:0.01:1 0:0.01:1];%rand(1000,1);
% x(501:1000,4)=rand(500,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x(:,1)=((x(:,1)-min(x(:,1)))/(max(x(:,1)-min(x(:,1)))));
% x(:,2)=((x(:,2)-min(x(:,2)))/(max(x(:,2)-min(x(:,2)))));


% x(:,3)=((x(:,3)-min(x(:,3)))/(max(x(:,3)-min(x(:,3)))));
% x(:,4)=((x(:,4)-min(x(:,4)))/(max(x(:,4)-min(x(:,4)))));
% x(:,3)=((x(:,3)-min(x(:,3)))/(max(x(:,3)-min(x(:,3)))));

% x(x>1)=1;
% x(x<0)=0;
% y=x;
% x=M;


% y1=x(:,2)+rand(size(x(:,2)))/5;

% [q1,b1]=sort(y1);
% [q2,b2]=sort(x(:,2));

% for i=1:numel(y1)
% y(i)=sum(y1>=y1(i))/(numel(y1)+1);
% z(i)=sum(x(:,2)>=x(i,2))/(numel(y1)+1);
% end

x=f;

% x(:,1)=x(:,2);
% y=x(:,1);
z=x(:,2)';%.^(-0.3);

% ddd=find(x(:,1)>1.5);
% y=y(ddd);
% z=z(ddd);
% x=x(ddd,1);

% y=x(:,1)+x(:,2);

% y=sqrt(2)*1/2*erfinv(2*x-1);
% x=y;
% x(x>1)=1;
% x(x<-1)=-1;
% 
% x=(x/2+0.5);
% x=u;


%%%%%%%%%%%%%%%%%%%%%%%%
%%%% transform into equipop data

% Nb=100;
% Mb=round(numel(y)/Nb);
% 
% for t=1:numel(y)
% b1(t)=sum(x(:,1)<=x(t,1));
% c1(t)=sum(z<=z(t));
% d1(t)=sum(y<=y(t));
% end
% b2=floor(b1/Mb)+1;
% c2=floor(c1/Mb)+1;
% d2=floor(d1/Mb)+1;
% 
% % clear x1
% % for t=1:Nb
% % x1(t)=mean(x(b2==t,1));
% % end    

b2=x(:,1);
c2=z';
d2=y';


%%%%%%%%%%%%%%%%%%%%%%%%%

NK=128;
% Set margins
vine.margins = cell(d,1);
% Standard normal margin
if d==3
vine.margins{3}.dist = 'kernel';
vine.margins{3}.theta = NK;
vine.margins{3}.iscont = 2; %  1 = binary , 2 = positive or zero
vine.margins{3}.ker=b2;
end

% Gamma margin
vine.margins{2}.dist = 'kernel';
vine.margins{2}.theta = NK;
vine.margins{2}.iscont = 2; 
vine.margins{2}.ker=c2;
% Gamma margin
vine.margins{1}.dist = 'kernel';
vine.margins{1}.theta = NK;
vine.margins{1}.iscont = 2; 
vine.margins{1}.ker=d2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set copula families
vine.families = cell(d);
vine.theta = cell(d);
% Gaussian copula family
vine.families{1,2} = 'kercop'; %%kercop
vine.families{2,1} = 'kercop'; %%kercop

if d==3
vine.families{1,3} = 'kercop'; %%kercop
vine.families{3,1} = 'kercop'; %%kercop
vine.families{2,3} = 'kercop'; %%kercop
vine.families{3,2} = 'kercop'; %%kercop
vine.theta{3,3} = b2;
end

vine.theta{2,2} = c2;
vine.theta{1,1} = d2;





fprintf('\n');
system('module load stats/R/3.3.1')

%% Test probability density function
disp('Calculating probability density function on a grid...');
% Calculate probability density function on lattice
vineW=vine;

x1gv = linspace(min(x(:,1))-0.05*abs(min(x(:,1))),max(x(:,1))-0.05*abs(max(x(:,1))),100);

x2gv = linspace(min(z)-0.05*abs(min(z)),max(z)+0.05*abs(max(z)),100);%linspace(0.5,25,50);

ygv = linspace(min(y)-0.05*abs(min(y)),max(y)+0.05*abs(max(y)),100);%min(y):max(y);%min(y):max(y);%min(y)-1:max(y)+1;% 


if d==3
    [x1,x2,x3] = ndgrid(ygv,x2gv,x1gv);
%     [p logp] = mixedvinepdf_HS(vineW,[x1(:),x2(:),x3(:)]);
    [p logp] = mixedvinepdf_HS_New(vineW,[x1(:),x2(:),x3(:)]);
    pp=reshape(p,numel(ygv),numel(x2gv),numel(x1gv));
    
    pp=pp(3:end-3,3:end-3,3:end-3);
    ygv=ygv(3:end-3);x2gv=x2gv(3:end-3);x1gv=x1gv(3:end-3);
    
    ppp=squeeze(sum(pp));
    figure;imagesc(x1gv,x2gv,ppp')

    ppp(ppp>1e-3)=1;
    figure;imagesc(x2gv,ygv,reshape(ppp,numel(x2gv),numel(x1gv)));

for ii=1:numel(z)
    [~,j2]=min(abs(x2gv(:)-z(ii)));
    [~,j3]=min(abs(x1gv(:)-x(ii,1)));
    q=squeeze(pp(:,j2,j3));

    [~,j]=max(smooth(q(1:end),1));
    y2(ii)=ygv(j);
%         y2(ii)=sum(q.*ygv'); %%%% expectation value

end
    
elseif d==2
    [x1,x2] = ndgrid(ygv,x2gv);
%     p = mixedvinepdf_HS(vineW,[x1(:),x2(:)]);
    p = mixedvinepdf_HS_New(vineW,[x1(:),x2(:)]);
    pp=reshape(p,numel(ygv),numel(x2gv));
    figure;imagesc(ygv(5:end),x2gv(5:end),pp(5:end,5:end)')
    

for ii=1:numel(z)
    [~,j2]=min(abs(x2gv(:)-z(ii)));
    q=squeeze(pp(:,j2));

    [~,j]=max(smooth(q(1:end),1));
    y2(ii)=ygv(j);    %%%% maximum likelihood

    
%     y2(ii)=sum(q.*ygv'); %%%% expectation value


%     y2(ii)=ygv(discretesample(q/sum(q), 1));%%%%%%% Sampling from q distribution
end

end

% y=smooth(y)';y2=smooth(y2)';prediction1se=smooth(prediction1se)';

figure;plot(y2);hold on;plot(y)
hold on;
plot(prediction1se,'g')
[corr(y2',y','type','Spearman') corr(prediction1se,y','type','Spearman')]

rrr3=1-sum(abs(y2-y).^2)/sum(abs(mean(y)-y).^2)
rrr4=1-sum(abs(prediction1se'-y).^2)/sum(abs(mean(y)-y).^2)

{'' 'cop' 'glm';'dev-full matrix' rrr3 rrr4;'Spear-corr' corr(y2',y,'type','Spearman') corr(prediction1se,y,'type','Spearman');'corr' corr(y2',y) corr(prediction1se,y)}


q=p;q(q>1e-4)=1;qq=reshape(q,100,200);figure;imagesc(x2gv,ygv,qq);



for u=1:20
cf(u)=corr(y((u-1)*100+1:u*100)',y2((u-1)*100+1:u*100)');
gf(u)=corr(y((u-1)*100+1:u*100)',prediction1se((u-1)*100+1:u*100));
end

figure;plot(cf);hold on;plot(gf,'g')

return

%%%%% using the data itself
tic
ygv = linspace(min(y),max(y),100);unique(y);%unique(y);linspace(min(y),max(y),100);min(y)-1:max(y)+1;%

clear XXXX PP y3
V=[z' x(:,1)];
% V=[z];
n=0;
XXXX=zeros(size(V,1)*numel(ygv),size(V,2)+1);

for j=1:size(V,1)
for i=1:numel(ygv)
    n=n+1;
XXXX(n,:)=[ygv(i) V(j,1) V(j,2)];
% XXXX(n,:)=[ygv(i) V(j,1)];
% XXXX(n,2)=V(j,1);
% XXXX(n,3)=V(j,2);
end
end

[pd logp] = mixedvinepdf_HS_New(vineW,XXXX);

n=0;
for j=1:size(V,1)
for i=1:numel(ygv)
    n=n+1;
PP(i,j)=pd(n);
end
end


for ii=1:size(V,1)
    [~,j]=max(smooth(PP(:,ii),1));
    y3(ii)=ygv(j);
end

figure;plot(y3);hold on;plot(y)
hold on;
plot(prediction1se,'g')
[corr(y3',y') corr(prediction1se,y')]

rrr5=1-sum(abs(y3-y).^2)/sum(abs(mean(y)-y).^2)
rrr6=1-sum(abs(prediction1se'-y).^2)/sum(abs(mean(y)-y).^2)

{'' 'cop' 'glm';'dev-full matrix' rrr3 rrr4;'dev-data' rrr5 rrr6;'corr' corr(y3',y) corr(prediction1se,y)}

toc




%%%%% using the data itself
tic
ygv = linspace(min(y),max(y),100);unique(y);%unique(y);linspace(min(y),max(y),100);min(y)-1:max(y)+1;%

clear XXXX PP y3
V=[z' x(:,1)];
% V=[z];
n=0;
XXXX=zeros(size(V,1)*numel(ygv),size(V,2)+1);

for j=1:size(V,1)
for i=1:numel(ygv)
    n=n+1;
XXXX(n,:)=[ygv(i) V(j,1) V(j,2)];
% XXXX(n,:)=[ygv(i) V(j,1)];
% XXXX(n,2)=V(j,1);
% XXXX(n,3)=V(j,2);
end
end

[pd logp] = mixedvinepdf_HS_New(vineW,XXXX);

n=0;
for j=1:size(V,1)
for i=1:numel(ygv)
    n=n+1;
PP(i,j)=pd(n);
end
end


for ii=1:size(V,1)
    [~,j]=max(smooth(PP(:,ii),1));
    y3(ii)=ygv(j);
end

figure;plot(y3);hold on;plot(y)
hold on;
plot(prediction1se,'g')
[corr(y3',y') corr(prediction1se,y')]

rrr5=1-sum(abs(y3-y).^2)/sum(abs(mean(y)-y).^2)
rrr6=1-sum(abs(prediction1se'-y).^2)/sum(abs(mean(y)-y).^2)

{'' 'cop' 'glm';'dev-full matrix' rrr3 rrr4;'dev-data' rrr5 rrr6;'corr' corr(y3',y) corr(prediction1se,y)}

toc

iscont = false(d,1);
for i = 1:d
    iscont(i) = true;%vine.margins{i}.iscont;
end

vineest = mixedvinefit(XXXX,vine.type,iscont)






return


% see=[2 3];
% n=0;
% for s=see
%     
%     n=n+1;
%     vine_b.type=vine.type;
%     vine_b.margins{n}=vine.margins{s};
%     m=0;
%     for s1=see
%         m=m+1;
%     vine_b.families{n,m}=vine.families{n,m};
%     vine_b.theta{n,m}=vine.theta{n,m};
%     end
%     
% end
% 
% p12 = mixedvinepdf(vine_b,[x2(:),x3(:)]);


pp=reshape(p,numel(x1gv),numel(x1gv),numel(x1gv));
pp12=squeeze(sum(pp));%reshape(p12,numel(x1gv),numel(x1gv));


for ii=1:size(x,1)
    [~,j2]=min(abs(x2gv(:)-x(ii,2)));
    [~,j3]=min(abs(x1gv(:)-x(ii,1)));
    q=squeeze(pp(:,j2,j3)/pp12(j2,j3));

    [~,j]=max(smooth(q(:),1));
    y2(ii)=ygv(j);
end

figure;plot(y2);hold on;plot(y)

% figure;plot(sum(p12./p1)/sum(sum(p12./p1)));hold on;plot(sum(p2)/sum(sum(p2)))

% figure;plot(sum(p12'./p2')/sum(sum(p12'./p2')));hold on;plot(sum(p1')/sum(sum(p1')))

rrr3=sum(abs(y2-y').^2)/sum(abs(mean(y)-y).^2)






return

x1gv = linspace(min(x(:,1)),max(x(:,1)),200);
vineW=vine;

ygv = linspace(min(y),max(y),100);
[x1,x2] = ndgrid(ygv,x1gv);

pd = mixedvinepdf(vineW,[x1(:),x2(:)]);

[o1 v1]=hist(y,100);
[o2 v2]=hist(x(:,2),200);

figure
plot(o2/sum(o2))
hold on
plot(sum(pdd)/(sum(sum(pdd))))

figure
plot(o1/sum(o1))
hold on
plot(sum(pdd')/(sum(sum(pdd'))))

for ii=1:size(x,1)
    [~,j2]=min(abs(x1gv(:)-x(ii,2)));
    q=squeeze(pdd(:,j2));

    [~,j]=max(smooth(q(:),1));
    y2(ii)=ygv(j);
end

figure;plot(y2);hold on;plot(y)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:2
switch i
    case 1
        FF{i}=marginpdf(vine.margins{i},x1(:));
        FFc{i}=margincdf(vine.margins{i},x1(:));
    case 2
        FF{i}=marginpdf(vine.margins{i},x2(:));
        FFc{i}=margincdf(vine.margins{i},x2(:));
end
end

%     FF{3} = copulaccdf(vine.families{1,2},[FF{1} FF{2}],vine.theta{2,2},2);
    FF{3} = copulaccdf(vine.families{1,2},[FFc{1} FFc{2}],[FFc{1} FFc{2}]);
%     FF{3} = copulapdf(vine.families{1,2},[x1(:),x2(:)],[x1(:),x2(:)]);
    
    FFF=FF{1}.*FF{2}.*FF{3};

    f=reshape(FFF,100,200);
    figure
    imagesc(f)
    
    
    
    
    

%%%%%%%%


% pt = mixedvinepdf(vineW,[x1(5,5),x2(10,10)]);

% Plot 2D margins
figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);
subplot(2,3,1);
margin12 = reshape(sum(sum(p,4),3),[length(x1gv),length(x2gv)]);
imagesc(x1gv,x2gv,margin12');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 2');
subplot(2,3,2);
margin13 = reshape(sum(sum(p,4),2),[length(x1gv),length(x3gv)]);
imagesc(x1gv,x3gv,margin13');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 3');
subplot(2,3,3);
margin14 = reshape(sum(sum(p,3),2),[length(x1gv),length(x4gv)]);
imagesc(x1gv,x4gv,margin14');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 4');
subplot(2,3,4);
margin23 = reshape(sum(sum(p,4),1),[length(x2gv),length(x3gv)]);
imagesc(x2gv,x3gv,margin23');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 3');
subplot(2,3,5);
margin24 = reshape(sum(sum(p,3),1),[length(x2gv),length(x4gv)]);
imagesc(x2gv,x4gv,margin24');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 4');
subplot(2,3,6);
margin34 = reshape(sum(sum(p,2),1),[length(x3gv),length(x4gv)]);
imagesc(x3gv,x4gv,margin34');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 3');
ylabel('Margin 4');
fprintf('\n');

%% Test sampling
disp('Sampling from mixed copula vine...');
% Draw samples
cases = 1000;
x = mixedvinernd(vine,cases);
% Plot samples in 2D
figure('Name','Mixed vine copula samples in 2D','Position',[0,0,1600,1000]);
subplot(2,3,1);
scatter(x(:,1),x(:,2),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
subplot(2,3,2);
scatter(x(:,1),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 3');
subplot(2,3,3);
scatter(x(:,1),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 4');
subplot(2,3,4);
scatter(x(:,2),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 3');
subplot(2,3,5);
scatter(x(:,2),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 4');
subplot(2,3,6);
scatter(x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 3');
ylabel('Margin 4');
% Plot samples in 3D
figure('Name','Mixed vine copula samples in 3D','Position',[0,0,1000,1000]);
subplot(2,2,1);
scatter3(x(:,1),x(:,2),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
zlabel('Margin 3');
subplot(2,2,2);
scatter3(x(:,1),x(:,2),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
zlabel('Margin 4');
subplot(2,2,3);
scatter3(x(:,1),x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 3');
zlabel('Margin 4');
subplot(2,2,4);
scatter3(x(:,2),x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 3');
zlabel('Margin 4');
fprintf('\n');

%% Test vine fit
disp('Fitting parameters to samples...');
% Construct argument to specify which margins are continuous
iscont = false(d,1);
for i = 1:d
    iscont(i) = true;%vine.margins{i}.iscont;
end
vineest = mixedvinefit(x,vine.type,iscont);

% Compare ground-truth and estimated mixed vine
fprintf('\nGround-truth mixed vine:\n');
for i = 1:d
    fprintf(' Margin %d:     %s with parameters\t',i,vine.margins{i}.dist);
    for j = 1:length(vine.margins{i}.theta)
        fprintf('\t%.2f ',vine.margins{i}.theta(j));
    end
    fprintf('\n');
end
for i = 1:d
    for j = (i+1):d
        fprintf(' Copula (%d,%d): %s with parameters\t',i,j,vine.families{i,j});
        for k = 1:length(vine.theta{i,j})
            fprintf('\t%.2f ',vine.theta{i,j}(k));
        end
        fprintf('\n');
    end
end
fprintf('\n');
disp('Estimated mixed vine:');
for i = 1:d
    fprintf(' Margin %d:     %s with parameters\t',i,vineest.margins{i}.dist);
    for j = 1:length(vineest.margins{i}.theta)
        fprintf('\t%.2f ',vineest.margins{i}.theta(j));
    end
    fprintf('\n');
end
for i = 1:d
    for j = (i+1):d
        fprintf(' Copula (%d,%d): %s with parameters\t',i,j,vineest.families{i,j});
        for k = 1:length(vineest.theta{i,j})
            fprintf('\t%.2f ',vineest.theta{i,j}(k));
        end
        fprintf('\n');
    end
end
fprintf('\n');

%% Estimate entropy
disp('Estimating entropy of mixed copula vine...');
alpha = 0.05; % Significance level of estimate (95% confidence)
erreps = 1e-1; % Maximum standard error
[h,stderr] = mixedvineentropy(vine,alpha,erreps);
% Display results
disp([' Estimated entropy:          ' num2str(h) ' bit']);
disp([' Standard error of estimate: ' num2str(stderr) ' bit']);
fprintf('\n');


% save('/home/hs258/Data_Folder/TEMP/tmp/test.mat')


%% Estimate mutual information
disp('Estimating mutual information for mixed vine and estimated mixed vine...');
% We calculate the mutual information between a binary variable and the
% vine and estimated vine as conditional multivariate distributions
vines = cell(2,1);
vines{1} = vine;
vines{2} = vineest;
pcond = [0.5;0.5]; % Equal probability of observing either conditional vine
alpha = 0.05; % Significance level of estimate (95% confidence)
erreps = 1e-1; % Maximum standard error
[info,stderr] = mixedvineinfo(vines,pcond,alpha,erreps);
% Display results
disp([' Estimated mutual information: ' num2str(info) ' bit']);
disp([' Standard error of estimate:   ' num2str(stderr) ' bit']);
fprintf('\n');
