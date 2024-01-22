function [I_dec perf]=kernelvinedecoder(vines,copula,data,nFolds)



ncond=2;
D=cat(1,data{1},data{2});
REAL=cat(1,ones(size(data{1},1),1),0*ones(size(data{2},1),1));
PREDICT=REAL*NaN;

if nFolds~=1
for nc=1:2
trialIDs = 1:size(data{nc},1);
trialPartition = cvpartition(length(trialIDs), 'KFold', nFolds);
foldIDs{nc} = nan(size(data{nc},1),1);
for nFold = 1:nFolds
    foldTrials = trialIDs(find(test(trialPartition,nFold)));
    foldFrames = find(ismember(1:size(data{nc},1),foldTrials));
    foldIDs{nc}(foldFrames) = nFold;
end
end

else
for nc=1:2
foldIDs{nc} = ones(size(data{nc},1),1);
end
end

lik1=nan(size(D,1),2);

for FOL=1:nFolds
    

    if nFolds~=1
    D0=cat(1,data{1}(foldIDs{1}==FOL,:),data{2}(foldIDs{2}==FOL,:));
    else
    D0=cat(1,data{1},data{2});
    end
    
    vines_F=vines;
    if nFolds~=1
    for nc=1:ncond
        for j=1:size(vines{nc}.margins,1)
    vines_F{nc}.margins{j}.ker=vines{nc}.margins{j}.ker(foldIDs{nc}~=FOL);
    vines_F{nc}.theta{j,j} = vines{nc}.margins{j}.ker(foldIDs{nc}~=FOL);
        end
    end
    end
    
    dims=[1,2];
    pdf=NaN*ones(size(D0,1),ncond);
    pdfn=NaN*ones(size(D0,1),ncond);
    qq=NaN*ones(size(D0,1),ncond);
    qqc=NaN*ones(size(D0,1),ncond);
    pp=NaN*ones(size(D0,1),ncond);
    qqcc=NaN*ones(size(D0,1),ncond);
    pm0=NaN*ones(ncond,1000,size(D0,1));
    ppT0=NaN*ones(ncond,1000,size(D0,1));
    BC=NaN*ones(ncond,1000,size(D0,1));
    Bf=NaN*ones(ncond,size(D0,1));
    
    y_vectorT=linspace(vines{1}.range(1,1),vines{1}.range(1,2),1000);
    
    for j = 1:ncond
        
        [p,~,~,~,~] = Fit_vCopula(vines_F{j},D0,size(copula{1},2),[],-1,copula{j},'rand',[]);
        
        X=p;
        qq(:,j)=X;   %%% f(B,n|  C)
                
        % Marginalizing the n
        mm=1;
        for sor=1:5
            
            disp(['fold = ', num2str(FOL),', condition = ', num2str(j),' , sor = ',num2str(sor)])
            y_vector=y_vectorT((sor-1)*200+1:sor*200);
            ne=0;
            points=zeros(numel(y_vector)*size(D0,1),size(D0,2));
            
            for i2=1:size(D0,1)
                for i1=1:numel(y_vector)
                    ne=ne+1;
                    Z(mm)=y_vector(i1);
                    Z(setdiff(1:size(D0,2),mm))=D0(i2,setdiff(1:size(D0,2),mm));
                    points(ne,:)=Z;
                end
            end
            
            [pT,~,~,~,~] = Fit_vCopula(vines_F{j},points,size(copula{1},2),[],-1,copula{j},'rand',[]);
            
            ppp=reshape(pT,numel(y_vector),[]);
            
            if sor==1
                pm=ppp;
            else
                pm=cat(1,pm,ppp);
            end
        end
        
        Bf(j,:)=sum(pm,1)';
        
        X=Bf(j,:);
        qqcc(:,j)=X*mean(diff(y_vectorT));      %f(B|C)
    end
    
    %%% likelihood f(n|CB)
        
    lik1([find(foldIDs{1}==FOL);find(foldIDs{2}==FOL)+size(data{1},1)],1:2)=qq./qqcc;
            
    PREDICT(find(foldIDs{1}==FOL))=lik1(find(foldIDs{1}==FOL),1)>lik1(find(foldIDs{1}==FOL),2);
    PREDICT(find(foldIDs{2}==FOL)+size(data{1},1))=lik1(find(foldIDs{2}==FOL)+size(data{1},1),1)>lik1(find(foldIDs{2}==FOL)+size(data{1},1),2);
end

I_dec=CONF_INFO(REAL,PREDICT);
perf=(sum(lik1(1:size(data{1},1),1)>lik1(1:size(data{1},1),2))+sum(lik1(size(data{1},1)+1:end,1)<lik1(size(data{1},1)+1:end,2)))/(size(lik1,1));


