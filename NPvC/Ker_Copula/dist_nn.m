

function bw_nn=dist_nn(alpha,data,grid,bw)

bw=bw/norm(bw);

dat=data.X;

alp=0.05:0.05:5;
nn=zeros(numel(alp),size(dat,1));
kn=sqrt(size(grid.X,1));

for m=alp
    mm=find(alp==m);
    for i1=1:size(dat,1)
%         D=((dat(:,1)-dat(i1,1))/(m*bw(1))).^2+((dat(:,2)-dat(i1,2))/(m*bw(2))).^2;
        D=((dat(:,1)-dat(i1,1))/(m*bw(1))).^2+((dat(:,2)-dat(i1,2))/(m*bw(2))).^2;
        
        nn(mm,i1)=sum(D<1)/numel(D);
        
    end
end


% per=[0.05:0.5:2.05 2.5:5:50];per=per/100;
bw_nn=zeros(2,size(dat,1));

p=alpha;
N=abs(nn-p*ones(size(nn)));
[~,U]=min(N);
bw_nn(1,:)=alp(U)*bw(1);
bw_nn(2,:)=alp(U)*bw(2);







