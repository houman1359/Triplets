

function err=MISE(a,data,Grid,points,xd,afk,dd,ddT,bw0,nn,NORM,uspace)

if numel(a)==1 && nn==0
    bw=abs(a*bw0);
    
else
    if nn==0
        bw=abs(a);
        
    else
        bw=dist_nn(abs(a),data,Grid,bw0);
    end
end



% uspace=0;



if uspace==0
    x1=unique(Grid.S(:,1));
    x2=unique(Grid.S(:,2));
    
    lfit = loclik_fit(bw,data,Grid);
%         pd_grid=lfit.Kergrid/(sum(lfit.Kergrid(:).*xd));
    [~,~,pd_grid,~,~,~,~]=func_tll(lfit,Grid,points,data,0,2,NORM);
    
    for i=1:numel(afk.TestSize)   %%% used to be parfo
        if nn==0
            lfitCV = loclik_fit(bw(:),dd{i},Grid);
        else
            lfitCV = loclik_fit(bw(:,~afk.test(i)),dd{i},Grid);
        end
        [kkk{i},~,pp,~,~,~,~]=func_tll(lfitCV,Grid,points,ddT{i},0,2,NORM);
%         kkk{i}=kkk{i}/(sum(pp(:).*xd));
%         ppp(i,:,:)=pp;
    end
    for i=1:numel(afk.TestSize)
        pd_dataCV(afk.test(i)==1)=kkk{i}/sum(pd_grid(:)*xd);
    end
    
%     err=trapz(x2,trapz(x1,pd_grid.^2,1))-2*mean(pd_dataCV);
    e=pd_grid.^2;
    err=sum(e(:).*xd)-2*mean(pd_dataCV);    

elseif uspace==1
    
    %%%%%% THIS IS THE OLD ORIGINAL 
    
for i=1:numel(afk.TestSize)
    if nn==0
    lfitCV = loclik_fit(bw(:),dd{i},Grid);
    else
    lfitCV = loclik_fit(bw(:,~afk.test(i)),dd{i},Grid);
    end
    [kkk{i},~,~,~,~,~,~]=func_tll(lfitCV,Grid,points,ddT{i},0,1,NORM);
    kkk{i}=kkk{i}/(sum(lfitCV.Kergrid(:).*xd));
end
for i=1:numel(afk.TestSize)
    pd_dataCV(afk.test(i)==1)=kkk{i};
end

lfit = loclik_fit(bw,data,Grid);
pd_grid=lfit.Kergrid/(sum(lfit.Kergrid(:).*xd));

% [~,~,pd_grid,~,~,~,~]=func_tll(lfit,Grid,points,data,0,1);
%  pd_grid=pd_grid/sum(lfit.Kergrid.*xd);

e=pd_grid.^2;

% err=sum(e(:).*xd)-2/numel(pd_dataCV)*sum(pd_dataCV);    
err=sum(e(:).*xd)-2*mean(pd_dataCV);    
    
elseif uspace==2
    
    x1=unique(Grid.u(:,1));
    x2=unique(Grid.u(:,2));
    xd=x1*x2';
    
    lfit = loclik_fit(bw,data,Grid);
    [~,~,pd_grid,~,~,~,~]=func_tll(lfit,Grid,points,data,-1,0,NORM);
    
    for i=1:numel(afk.TestSize)
        if nn==0
            lfitCV = loclik_fit(bw(:),dd{i},Grid);
        else
            lfitCV = loclik_fit(bw(:,~afk.test(i)),dd{i},Grid);
        end
        [kkk{i},~,~,~,~,~,~]=func_tll(lfitCV,Grid,points,ddT{i},-1,0,NORM); %%%-1,0
    end
    for i=1:numel(afk.TestSize)
        pd_dataCV(afk.test(i)==1)=kkk{i};
    end
    
    err=trapz(x2,trapz(x1,pd_grid.^2,1))-2/numel(pd_dataCV)*sum(pd_dataCV);
    
end