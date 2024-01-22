function [info,stderr_tot,in,En] = kernelvineinfo(vines,copula,cases,KMIN,SET_CO,pcond,neuro_N,time_N,tvec,par)

if nargin>6
    mm=neuro_N;
    nn=time_N;
end
INTE=1;
if nargin==10
    INTE=par.INT;
end

alpha=0.05;
erreps = 5e-3;
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% code to compute mutual information given two vine copulas or between two vine copulas
% copyright Houman Safaai March 2018
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Argument checks
if nargin < 2
    error('kernelvineinfo: Usage [info,stderr] = mixedvineinfo(vines,pcond,alpha,erreps,cases)');
end
if nargin < 3 || isempty(alpha)
    alpha = 0.01;
elseif ~isscalar(alpha)
    error('kernelvineinfo: Argument "alpha" must be a scalar');
end
if nargin < 4 || isempty(erreps)
    erreps = 1e-3;
elseif ~isscalar(erreps)
    error('kernelvineinfo: Argument "erreps" must be a scalar');
end

if nargin>5
    pcond0=pcond;
    % Number of conditions
    ncond = length(cell2mat(pcond));
    % Probabilities of conditions
else
    ncond=1;
end

% Gaussian confidence interval for erreps and level alpha
conf = norminv(1 - alpha,0,1);

h=0;
hcond=0;


kmin=KMIN;%    this is used for copula info paper 40;
k = 0;

if nargin==9
    
    for tt=1:numel(tvec)-1
        info(tt)=0;
        N=zeros(ncond,numel(tvec));
        
        for w=1:18
            En(w,tt) = 0;
            varsum(w,tt) = 0;
            stder(w,tt) = 0;
        end
        
        
        for w=1:ncond
            En1(w+2,tt) = 0;
            En2(w+12,tt) = 0;
            En3(w+6,tt) = 0;
            En4(w+9,tt) = 0;
            En5(w+16,tt) = 0;
        end
    end
    
else
    varsum1 = 0;
    varsum2 = 0;
    varsum_tot = 0;
    info = 0;
    infoc1 = 0;
    infoc2 = 0;
    infoc = 0;
    infot = 0;
end

for w=1:18
    En(w,1) = 0;
    varsum(w,1) = 0;
    stder(w,1) = 0;
end

if nargin==10 & strcmp(par,'par')==1
    param=1;
else
    param=0;
end

stderr1 = inf;
stderr2 = inf;
stderr_tot = inf;


while (stderr1 >= erreps | stderr2 >= erreps | stderr_tot >= erreps) & k<kmin
    
    if param==0
        if length(vines)==1
            
            
            infok = zeros(cases,1);
            infokc1 = zeros(cases,1);
            infokc = zeros(cases,1);
            
            k = k + 1;
            % Generate samples
            
            clear X pp qq pp_par qq_par
            for hh=1:numel(vines.margins)
                X(:,hh)=vines.margins{hh}.ker;
            end
            
            if k==1
                %                                         [~,~,copS,~,~] = Fit_vCopula(vines,X(1,:),size(copula,2),vines.METH,-1,copula,'rand',[],500);  %%%% we generate a 500 x 500 template copula to use in sampling only. This make the sampling more uniform in the u plane.
            end
            copS=copula;
            
            [D du]=kerncoprnd(copS,X,cases,vines);
            for dd=1:size(D,2)
                %                         D(:,dd)= D(:,dd)+abs(1e-10*rand(size(D(:,dd))));
                 %p.fit=1;p.max=max(copula{1,1}.MarginG.p);p.min=min(copula{1,1}.MarginG.p);
                 %[pc(dd,:),~]=kernelcdf(vines.margins{dd}.ker,D(:,dd),p);
                 %[pcD(dd,:),~]=kernelcdf(D(:,dd),D(:,dd),p);
            end
            
%             DD=[];
%             xx=floor(vines.margins{2}.ker);
%             xu=unique(xx);
%             for xc=1:numel(xu)
%                 n_x(xc)=round(cases*sum(xx==xu(xc))/numel(xx));
%                 DD=[DD;D(randsample(find(floor(D(:,2))==xu(xc)),n_x(xc)),:)];
%             end            
%             D=DD;
            
            %[~,~,~,~,p_cop] = Fit_vCopula(vines,D,size(copula,2),[],-1,copula,'rand',[]);
            [~,~,~,~,p_cop] = Fit_vCopula(vines,du,size(copula,2),[],-4,copula,'rand',[]); %-4 is when the input is the copula samples
                

            
            log2pp=log2(p_cop);%.*prod(pcD)'./prod(pc)';
            log2pp(p_cop==0)=0;
            I= log2pp  ;
            infokc1 =  I;
            
            infoc1 = infoc1 + ( nanmean(infokc1(~isinf(infokc1))) - infoc1) / k;
            
            varsum1 = varsum1 + nansum((infokc1(~isinf(infokc1)) - infoc1) .^ 2);
            stderr1 = conf * sqrt(varsum1 / (k * cases * (k * cases - 1)));
            info = (infoc1);
            
            in.info(k)=info;
            in.stderr(k)=stderr1;
            
            figure(100)
            hold on
            subplot(1,3,1);plot(k,info,'O')
            title('Info')
            hold on
            subplot(1,3,2);plot(k,stderr1,'*');hold on;plot(k,stderr_tot,'O')
            title('stderr')
            hold on
            subplot(1,3,3);plot(k,mean(infokc1),'*')%;hold on;plot(k,mean(infokt),'O')
            title('E(infok)')
            drawnow
            
            En=NaN;
            
        else
            clear infok infokc1 infokc2 infokc infokt pp qq pp_par qq_par D cc
            
            SET_CON=cell2mat(SET_CO);
            pcon=cell2mat(pcond)/sum(cell2mat(pcond));
            
            SETS=[];
            for i=1:numel(SET_CO)
                SETS=[SETS i*ones(1,numel(SET_CO{i}))];
            end
            conditions=unique(SETS);
            
            NSAM=100;
            k = k + 1;
            % Generate samples
            
            clear Dt
            if k==1
                for i=1:ncond
                    D{i}=[];
                    clear DA
                    for j=1:numel(vines{SET_CON(i)}.margins)
                        DA(:,j)=vines{SET_CON(i)}.margins{j}.ker;
                    end
                    [~,~,copS{SET_CON(i)},~,~] = Fit_vCopula(vines{SET_CON(i)},DA(1,:),size(copula{SET_CON(i)},2),vines{SET_CON(i)}.METH,-1,copula{SET_CON(i)},'rand',[],50); %200
                    copula{SET_CON(i)}=copS{SET_CON(i)};
                    %                     vines{SET_CON(i)}.copulaB=copSB{SET_CON(i)};
                end
            end
            
            for i=1:ncond
                %                 [Samp{i} TV{i}]=vinecopula_sample(vines{SET_CON(i)},copS{SET_CON(i)},tvec,round(NSAM*pcon(i)/max(pcon)),cases,nn);    %%%% nn is time
                [Samp{i} TV{i}]=vinecopula_sample(vines{SET_CON(i)},copS{SET_CON(i)},tvec,NSAM,cases,nn);    %%%% nn is time
                if i==1
                    tind=TV{1};
                    SampSHUF=Samp{i};
                    %                     clear N
                    %                     for j=1:size(copula{SET_CON(1)},2)
                    %                         N(:,j)=vines{SET_CON(i)}.margins{j}.ker;
                    %                     end
                    %                     SampSHUF=N;
                else
                    tind= intersect(tind,TV{i});
                    SampSHUF=cat(1,SampSHUF,Samp{i});
                    %                     clear N
                    %                     for j=1:size(copula{SET_CON(1)},2)
                    %                         N(:,j)=vines{SET_CON(i)}.margins{j}.ker;
                    %                     end
                    %                     SampSHUF=cat(1,SampSHUF,N);
                end
            end
            tvec=tvec(tind);            

            Dtsh=SampSHUF(randsample(1:size(SampSHUF,1),size(SampSHUF,1)),:);
            for i=1:ncond
                if i==1
                    Dt=Samp{1};
                    %                     clear N
                    %                     for j=1:size(copula{SET_CON(1)},2)
                    %                         N(:,j)=vines{SET_CON(i)}.margins{j}.ker;
                    %                     end
                    N=Samp{i};Ds=Dt;
                    Ds=N;
                    SampSH{i}=Dtsh(1:size(Ds,1),:);
                else
                    %                     clear N
                    %                     for j=1:size(copula{SET_CON(1)},2)
                    %                         N(:,j)=vines{SET_CON(i)}.margins{j}.ker;
                    %                     end
                    N=Samp{i};Ds=Dt;
                    SampSH{i}=Dtsh(size(Ds,1)+1:size(Ds,1)+size(N,1),:);
                    Ds=cat(1,Ds,N);
                    Dt=cat(1,Dt,Samp{i});
                end
            end
            
            vinesSH=vines;
            for j=1:ncond
                for i=1:size(copula{SET_CON(1)},2)
%                 vines{SET_CON(j)}.margins{i}.ker=Samp{j}(:,i);
%                 vines{SET_CON(j)}.theta{i,i}=Samp{j}(:,i);
                vinesSH{SET_CON(j)}.margins{i}.ker=SampSH{j}(:,i);
                vinesSH{SET_CON(j)}.theta{i,i}=SampSH{j}(:,i);
                end
            end
            
            clear D Q
            D=Dt;
            
            [~,h2]=histc(D(:,nn),tvec);
            
            if min(h2)==0
                h2=h2+1;
            end
            
            DS=D(:,2:end);
            
            QB=[];
            clear Q Qsh
            for j=1:ncond
                
                [Q(:,j),~,CC{SET_CON(j)},~,~] = Fit_vCopula(vines{SET_CON(j)},D,size(copula{SET_CON(1)},2),vines{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
                [Qsh(:,j),~,CCsh{SET_CON(j)},~,~] = Fit_vCopula(vinesSH{SET_CON(j)},D,size(copula{SET_CON(1)},2),vinesSH{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
                
                if isfield(vines{SET_CON(j)},'B')
                    [QB(:,j),~,~,~,~] = Fit_vCopula(vines{SET_CON(j)}.B,DS,size(vines{SET_CON(j)}.copulaB,2),vines{SET_CON(j)}.B.METH,-1,vines{SET_CON(j)}.copulaB,'rand',[]);
                end
                
%                 copula{SET_CON(j)}=CC{SET_CON(j)};
%                 [Q(:,j),~,~,~,~] = Fit_vCopula(vines{SET_CON(j)},D,size(copula{SET_CON(1)},2),vines{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
%                 [Qsh(:,j),~,~,~,~] = Fit_vCopula(vinesSH{SET_CON(j)},D,size(copula{SET_CON(1)},2),vinesSH{SET_CON(j)}.METH,-1,CCsh{SET_CON(j)},'rand',[]);
                
            end
            
            
            for i=1:ncond
                if i==1
                    LAB=ones(size(Samp{1},1),1);
                else
                    LAB=cat(1,LAB,SETS(i)*ones(size(Samp{i},1),1));
                end
            end
            
            
            NNN=200
            dims=[1,2];
            Qnb=Q;
            
            ncond=numel(SET_CON);
            clear pMar pdf pdfTV pdfn pp qq qqc qqcc pq pmmm fm pmmmT fmT pppp
            for j = 1:ncond
                
                y_vectorT=linspace(min(vines{SET_CON(j)}.margins{1}.ker),max(vines{SET_CON(j)}.margins{1}.ker),200);%200;

%                 y_vectorT=uniquetol(vines{SET_CON(j)}.margins{1}.ker,1e-4)';
                
                %             y_vectorT=sort(unique(yy));                
                %             y_vectorT=sort(unique([y_vectorT;vines{1}.y_vector']));
                NBI=floor(numel(y_vectorT)/NNN);
                if NBI*NNN<numel(y_vectorT)
                    NBI=NBI+1;
                end
                
                Q0=Qnb;
                if size(vines{SET_CON(1)}.range,1)>2
                    [vineMar{j},copulaMar{j}]=kernelmarginalize(dims,vines{SET_CON(j)},copula{SET_CON(j)});
                    [vineMarSH{j},copulaMarSH{j}]=kernelmarginalize(dims,vinesSH{SET_CON(j)},CCsh{SET_CON(j)});
                    [pMar(:,j),~,~,~,~] = Fit_vCopula(vineMar{(j)},D(:,dims),size(copulaMar{1},2),vineMar{(j)}.METH,-1,copulaMar{(j)},'rand',[]); %% f(n,t|C)
                    [pdfSH(:,j),~,~,~,~] = Fit_vCopula(vineMarSH{(j)},D(:,dims),size(copulaMarSH{1},2),vineMarSH{(j)}.METH,-1,copulaMarSH{(j)},'rand',[]); %% f(n,t|C)
                else
                    vineMar{j}=vines{SET_CON(j)};
                    copulaMar{j}=copula{SET_CON(j)};
                    pMar(:,j)=Q(:,j);
                end
                
                %%%%%% LATER I SHOULD ADD I(n;b_i| B_{-i} ) for
                %%%%%% each component of B
                %                             for x=1:size(D,2)-2
                %                             [pSingle{x}(:,j),~,~,~,~] = Fit_vCopula(vineMar{j},D(:,[1 x+2]),size(copulaMar{1},2),[],-1,copulaMar{j},'rand',[]); %% f(n,t|C)cc
                %                             end
                
                p=Qnb(:,j);
                
                par.fit=0;par.s=copula{SET_CON(j)}{nn,1}.MarginS.s;par.p=copula{SET_CON(j)}{nn,1}.MarginS.p;
                [pdf(:,j),~]=kernelpdf(D(:,nn),2,D(:,nn),par);
                par.fit=0;par.s=CCsh{SET_CON(j)}{nn,1}.MarginS.s;par.p=CCsh{SET_CON(j)}{nn,1}.MarginS.p;
                [pdfSH(:,j),~]=kernelpdf(D(:,nn),2,D(:,nn),par);
                
                par.fit=0;par.s=copula{SET_CON(j)}{mm,1}.MarginS.s;par.p=copula{SET_CON(j)}{mm,1}.MarginS.p;
                [pdfn0V(:,j),~]=kernelpdf(D(:,mm),2,D(:,mm),par);
%                 [pdfSH(:,j),~]=kernelpdf(Dtsh(:,mm),2,Dtsh(:,mm),par);
                
                X=p./pdf(:,j);
                qq0V(:,j)=X;   %%% f(B,n|t,C)
                
                X=Qsh(:,j)./pdfSH(:,j);
                qq0Vsh(:,j)=X;   %%% f(B,n|t,C)
                
                
                
                X=p./pMar(:,j); %%% f(B|t,n,C)
                qqcV(:,j)=X;
                
                X=pMar(:,j)./pdf(:,j);
                pp0V(:,j)=X;   %%% f(n|t,C)
                
                X=pdf(:,j);
                pppp(:,j)=X;   %%% f(t|C)
                
                % Marginalizing the n
                if nn~=1 & isempty(QB) & INTE==1
                    
                    for sor=1:NBI
                        
                        disp(['condition = ', num2str(j),' , sor = ',num2str(sor),'/',num2str(NBI)])
                        
                        if sor*NNN<numel(y_vectorT)
                            y_vector=y_vectorT((sor-1)*NNN+1:sor*NNN);
                        else
                            y_vector=y_vectorT((sor-1)*NNN+1:end);
                        end
                        ne=0;
                        points=zeros(numel(y_vector)*size(D,1),size(D,2));
                        mmm=setdiff(1:size(D,2),mm);
                        
                        for i2=1:size(D,1)
                            for i1=1:numel(y_vector)
                                ne=ne+1;
                                Z(mm)=y_vector(i1);
                                Z(mmm)=D(i2,mmm);
                                points(ne,:)=Z;
                            end
                        end
                        
                        [pT,~,~,~,~] = Fit_vCopula(vines{SET_CON(j)},points,size(copula{SET_CON(j)},2),vines{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
                        
                        ppp=reshape(pT,numel(y_vector),[]);
                        
                        if sor==1
                            pm=ppp;
                        else
                            pm=cat(1,pm,ppp);
                        end
                        
                    end
                    
                    clear pmm
                    for jj=1:size(pm,2)
                        pmm(jj,1)=trapz(y_vectorT,pm(:,jj));
                    end
                    %                     pmm=sum(pm(1:end-1,:).*diff(y_vectorT),1)';
                    X=pmm./pdf(:,j);
                    qqccV(:,j)=X;
                    
                else
                    
                    if ~isempty(QB)
                        qqcc0(:,j)=QB(:,j)./pdf(:,j);
                    else
                        qqcc0=qq0;
                    end
                    
                end
                
            end
            
            if k==1
                for tt=1:numel(tvec)-1
                    info(tt)=0;
                    N=zeros(ncond,numel(tvec));
                    for w=1:18
                        En(w,tt) = 0;
                        varsum(w,tt) = 0;
                        stder(w,tt) = 0;
                    end
                    for w=1:ncond
                        En1(w+2,tt) = 0;
                        En2(w+12,tt) = 0;
                        En3(w+6,tt) = 0;
                        En4(w+9,tt) = 0;
                        En5(w+16,tt) = 0;
                    end
                end
            end
            
                        Vnbt=1/1e6;%prod(diff(vines{1}.range(:,:)'));
                        Vbt=1/1e6;%prod(diff(vines{5}.range(setdiff(1:size(Dt,2),1),:)'));
                        Vn=1/1;%prod(diff(vines{5}.range(1,:)'));
            
            %                         qqcc(qqcc<1e-1/Vbt)=0;
            %                         qq(qq<1e-1/Vnbt)=0;
            %                         pp(pp<1e-1/Vn)=0;
            
            %%%% convert to condtions, summing over sub_conditions
            %             clear qq pp qqcc Q pdfn
            
            
            
            %             pcon(:)=0.25;
            PO=repmat(pcon,size(qq0V,1),1);
            clear Q qqsh qq qqcc pp pdfn
            for i=1:numel(conditions)
                Q(:,i)=sum(PO(:,find(SETS==i)).*Q0(:,find(SETS==i)),2);
                qq(:,i)=sum(PO(:,find(SETS==i)).*qq0V(:,find(SETS==i)),2);
                qqsh(:,i)=sum(PO(:,find(SETS==i)).*qq0Vsh(:,find(SETS==i)),2);
                qqcc(:,i)=sum(PO(:,find(SETS==i)).*qqccV(:,find(SETS==i)),2);
                pp(:,i)=sum(PO(:,find(SETS==i)).*pdf(:,find(SETS==i)),2);
                pdfn(:,i)=sum(PO(:,find(SETS==i)).*pp0V(:,find(SETS==i)),2);
%                 pdfn(:,i)=sum(PO(:,find(SETS==i)).*pdfSH(:,find(SETS==i)),2);
            end
%             qq=pp;
            qq0=qqsh;qqcc0=qq;pdfn0=pdfn;
            for i=1:ncond
%             qqsh(max(qq0')<1e2,i)=0;%mean(qq0(mean(qq0,2)<1/Vnbt,:),2);
%             qqcc(mean(qqcc0,2)<1e8,i)=0;%mean(qqcc0(mean(qqcc0,2)<1/Vbt,:),2);
%             pdfn(mean(pdfn0,2)<1e1,i)=0;%mean(pdfn0(mean(pdfn0,2)<1/Vn,:),2);
            end
            
            clear In
            for tt=1:numel(tvec)-1
                ff=find(h2==tt);
                clear PR
                for i=conditions
                    PR(1,i)=sum(LAB(ff)==i)/numel(ff);
                end
                PR=repmat(PR,numel(ff),1);
                
                clear posnb posb posn posnb_SH
                for i=1:numel(conditions)
                    posnb(:,i)=qq(ff,i)*PR(1,i)./nansum(PR.*qq(ff,:),2);
                    posnb_SH(:,i)=posnb(:,i);%qqsh(ff,i)*PR(1,i)./nansum(PR.*qqsh(ff,:),2);
                    posb(:,i)=qqcc(ff,i)*PR(1,i)./nansum(PR.*qqcc(ff,:),2);
                    posn(:,i)=pdfn(ff,i)*PR(1,i)./nansum(PR.*pdfn(ff,:),2);
                    posnb(isnan(sum(posnb(:,i),2)),i)=PR(1,i);
                    posnb_SH(isnan(sum(posnb_SH(:,i),2)),i)=PR(1,i);
                    posb(isnan(sum(posb(:,i),2)),i)=PR(1,i);
                    posn(isnan(sum(posn(:,i),2)),i)=PR(1,i);
                end
                                
                HC=-sum(PR(1,:).*log2(PR(1,:)));
                Inb(tt)=HC+nanmean(nansum(posnb.*log2(posnb),2),1);
                Inb_sh(tt)=HC+nanmean(nansum(posnb_SH.*log2(posnb_SH),2),1);
                Ib(tt)=HC+nanmean(nansum(posb.*log2(posb),2),1);
                In(tt)=HC+nanmean(nansum(posn.*log2(posn),2),1);                
                
                fn=nansum(PR.*pdfn(ff,:),2); %%% f(n)
                fb=nansum(PR.*qqcc(ff,:),2); %%% f(b)
                fnb=nansum(PR.*qq(ff,:),2); %%% f(b)
                fnCb=fnb./fb;
                fbCn=fnb./fn;
                
                Hn=-mean(log2(fn(fnCb~=0 & fb~=0)));
                HnCb=-mean(log2(fnCb(fnCb~=0 & fb~=0)));
                InAb(tt)=Hn-HnCb;
                
                Hb=-mean(log2(fb(fbCn~=0 & fn~=0)));
                HbCn=-mean(log2(fbCn(fbCn~=0 & fn~=0)));
                IbAn(tt)=Hb-HbCn;
            end
            
            figure(5)
            subplot(1,3,1)
            hold on;plot(Inb-Inb_sh*0);
            %hold on;plot(Ib);
            axis square;title('I(nB;C),I(B;C)')
            subplot(1,3,2)
            hold on;plot(Inb-Ib);axis square;title('I(n;C|B)')
            subplot(1,3,3)
            hold on;plot(In);axis square;title('I(n;C)')
            
            figure;plot(smooth(Inb));hold on;plot(smooth(Inb_sh))
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            for tt=1:numel(tvec)-1
                
                ff=find(h2==tt);
                qw=Q(ff,:);
                yy = quantile(qw(:),0.05);
                %                 ff=find(h2==tt & (Q(:,1)>yy | Q(:,2)>yy)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% this threshhold SHOULD BE CHECKED
                
                clear PR
                for i=conditions
                    PR(1,i)=sum(LAB(ff)==i)/numel(ff);
                end
                
                %                                 PR=[sum(LAB(ff)==1) sum(LAB(ff)==2)]/numel(ff);%[0.5 0.5];%
                
                PR=repmat(PR,numel(ff),1);
                %                                             PR=(PR.*pdf(ff,:))./repmat(sum(PR.*pdf(ff,:),2),1,2);
                
                pdfnc=nansum(PR .* qqcc(ff,:),2);
                pdfncMeasure0=nansum(PR .* qq(ff,:),2);
                pdfncMeasure00=nansum(PR .* pp(ff,:),2);
                pdfncMeasure00n=nansum(PR .* pdfn(ff,:),2);
                
                
                for co=conditions
                    G=-(qqcc(ff,co).*log2(qqcc(ff,co))./pdfnc);%G(qq(ff,co)==0 | pdfncMeasure0==0)=0;
                    In= PR(:,co).*G;
                    
                    %                                                                 G=-log2(qqcc(ff(LAB(ff)==co),co));%G(qqcc(ff(LAB(ff)==co),co)==0)=0;
                    %                                                                 In= PR(LAB(ff)==co,co).*G;
                    In(isinf(abs(In)) | isnan(In))=0;
                    H{co+2,tt}=In;
                    H1{co+2,tt}=In;
                end
                
                G=-log2(pdfnc);%G(pdfncMeasure0==0)=0;
                In= G;
                In(isinf(abs(In)) | isnan(In))=0;
                H{15,tt}=In;                   %%H(B|C)
                
                %%%%%%%%%%%
                
                G=-log2(pdfncMeasure0);
                I= G;
                I(isinf(abs(I)) | isnan(I))=0;
                H{5,tt}=I;
                
                G=-log2(pdfnc.*pdfncMeasure00);
                I= G;
                I(isinf(abs(I)) | isnan(I))=0;
                H{6,tt}=I;
                
                for co=conditions
                    G=log2(qq(ff(LAB(ff)==co),co)./(pp(ff(LAB(ff)==co),co).*qqcc(ff(LAB(ff)==co),co)));
                    In= PR(LAB(ff)==co,co).*G;
                    In(isinf(abs(In)) | isnan(In))=0;
                    H{12+co,tt}=In;
                    H2{12+co,tt}=In;
                end
                
                %%%%%%%%%%%
                
                
                for co=conditions
                    G=-pp(ff,co).*log2(pp(ff,co))./pdfncMeasure00; %pp(ff(sum(isnan(pp(ff,:)),2)~=0),:)=0;
                    I= PR(:,co).*G;
                    I(isinf(abs(I)) | isnan(I))=0;
                    H{6+co,tt}=I;
                    H3{6+co,tt}=I;
                end
                
                G=-log2(pdfncMeasure00);
                In= G;
                In(isinf(abs(In)) | isnan(In))=0;
                H{9,tt}=In;
                
                
                for co=conditions
                    G=-pdfn(ff,co).*log2(pdfn(ff,co))./pdfncMeasure00n; %pp(ff(sum(isnan(pp(ff,:)),2)~=0),:)=0;
                    I= PR(:,co).*G;
                    I(isinf(abs(I)) | isnan(I))=0;
                    H{16+co,tt}=I;
                    H5{16+co,tt}=I;
                end
                
                G=-log2(pdfncMeasure00n);
                In= G;
                In(isinf(abs(In)) | isnan(In))=0;
                H{16,tt}=In;
                
                %%%%%%%%%%
                
                for co=conditions
                    
                    G=-(qq(ff,co).*log2(qq(ff,co))./pdfncMeasure0);  %G(qqcc(ff,co)==0 | pdfnc==0)=0;
                    I= PR(:,co).*G;
                    
                    %                                                                 G=-log2(qq(ff(LAB(ff)==co),co));%G(qq(ff(LAB(ff)==co),co)<1e-6)=NaN;
                    %                                                                 I= PR(LAB(ff)==co,co).*G;
                    
                    I(isinf(abs(I)) | isnan(I))=0;
                    H{9+co,tt}=I;
                    H4{9+co,tt}=I;
                end
                
                G=-log2(pdfncMeasure0);%G(pdfnc==0)=0;
                G(isinf(abs(G)) | isnan(G))=0;
                H{12,tt}=G;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ff~=0
                    for w=1:18
                        En(w,tt) = En(w,tt)+  (nanmean(H{w,tt}(~isinf(H{w,tt}))) - En(w,tt)) /k;
                        varsum(w,tt) = varsum(w,tt) + nansum((H{w,tt}(~isinf(H{w,tt})) - En(w,tt)).^ 2);
                        NN=numel(H{w,tt});
                        stder(w,tt) = conf * sqrt(varsum(w,tt) / (NN*k * (NN*k - 1)));
                    end
                    
                    
                    for w=conditions
                        En1(w+2,tt) = En1(w+2,tt)+  (nanmean(H1{w+2,tt}(~isinf(H1{w+2,tt}))) - En1(w+2,tt)) /k;
                        En5(w+16,tt) = En5(w+16,tt)+  (nanmean(H5{w+16,tt}(~isinf(H5{w+16,tt}))) - En5(w+16,tt)) /k;
                        En2(w+12,tt) = En2(w+12,tt)+  (nanmean(H2{w+12,tt}(~isinf(H2{w+12,tt}))) - En2(w+12,tt)) /k;
                        En3(w+6,tt) = En3(w+6,tt)+  (nanmean(H3{w+6,tt}(~isinf(H3{w+6,tt}))) - En3(w+6,tt)) /k;
                        En4(w+9,tt) = En4(w+9,tt)+  (nanmean(H4{w+9,tt}(~isinf(H4{w+9,tt}))) - En4(w+9,tt)) /k;
                    end
                end
                
            end
            
            
            clear in stde
            ncond=numel(conditions);
            
            for tt=1:numel(tvec)-1
                
                %%%   I(n;C)
                
                %                 in.info(1,tt)=En(9,tt)-En(8,tt)-En(7,tt);
                in.info(1,tt)=En(9,tt)-sum(En3(7:6+ncond,tt));
                
                %%%   I(n;B)
                
                in.info(2,tt)=En(6,tt)-En(5,tt);
                
                %%%   I(B;C)
                
                %                 in.info(3,tt)=En(15,tt)-En(3,tt)-En(4,tt);
                in.info(3,tt)=En(15,tt)-sum(En1(3:2+ncond,tt));
                
                %%% I(n,B;C)
                
                %                 in.info(4,tt)=En(12,tt)-En(10,tt)-En(11,tt);
                in.info(4,tt)=En(12,tt)-sum(En4(10:9+ncond,tt));
                
                %%% I(n;C | B)
                
                %                             in.info(in.info<0)=0;
                
                in.info(5,tt)=in.info(4,tt) - in.info(3,tt);
                
                %%% I(n;B | C)
                
                %                 in.info(6,tt)=En(13,tt)+En(14,tt);
                in.info(6,tt)=sum(En2(13:12+ncond,tt));
                
                %%% I(n;B , C)
                
                in.info(7,tt)=in.info(5,tt) + in.info(2,tt);
                
                
                %%% I(n;B , C)
                
                in.info(8,tt)=En(16,tt)-sum(En(17:16+ncond,tt));
                
                
                stde(1,tt)=stder(3,tt);
                stde(2,tt)=stder(4,tt);
                stde(3,tt)=stder(15,tt);
                stde(4,tt)=stder(10,tt);
                stde(5,tt)=stder(11,tt);
                stde(6,tt)=stder(12,tt);
                
            end
            
            
            %             in.info(in.info<0)=0;
            
            
            stderr1=0;%norm(stderr1T(~isnan(stderr1T) & ~isinf(stderr1T)));
            stderr2=0;%norm(stderr2T(~isnan(stderr2T) & ~isinf(stderr2T)));
            stderr_tot=max(stder(:));
            
            in.stderr(k)=stderr_tot;
            
            
            
            figure(100)
            subplot(2,3,1)
            plot(tvec(1:end-1),in.info(1,:))
            hold on
            plot(tvec(1:end-1),in.info(8,:))
            hold on
            %                         xlim([5 290])
            title('I(n;C)')
            
            subplot(2,3,2)
            plot(tvec(1:end-1),in.info(4,:))
            hold on
            %                                                 plot(tvec(1:end-1),squeeze(mean(HH(:,12,:)-HH(:,10,:)-HH(:,11,:),1)),'.-')
            title('I(n,B;C)')
            
            subplot(2,3,3)
            plot(tvec(1:end-1),in.info(5,:))
            hold on
            %                                                 xlim([5 290])
            title('I(n;C|B)')
            
            subplot(2,3,4)
            plot(tvec(1:end-1),in.info(3,:))
            hold on
            %                         plot(tvec(1:end-1),inSH.info(3,:),'r')
            %                         xlim([5 290])
            title('I(B;C)')
            
            subplot(2,3,5)
            plot(tvec(1:end-1),in.info(2,:))
            hold on
            %                         xlim([5 290])
            title('I(n;B)')
            
            subplot(2,3,6)
            plot(tvec(1:end-1),in.info(6,:))
            hold on
            %                         plot(tvec(1:end-1),in.info(4,:)-in.info(1,:))
            %                         xlim([5 290])
            title('I(n;B|C)')
            drawnow
            
            
            
            
            
            %          close all
            
        end
    end
    
    
    if param==1
        V=1;
        
        infok = zeros(cases,1);
        infokc1 = zeros(cases,1);
        
        k = k + 1;
        % Generate samples
        
        D = mixedvinernd(vines,cases);
        
        [p,~,pdf] = mixedvinepdf(vines,D);%Fit_vCopula(vines{j},cc{j},size(copula{1},2),[],0,copula{j},'rand');
        pcop=p./pdf;
        
        log2pp=log2(pcop);log2pp(pcop==0)=0;
        I=log2pp  ;
        I(isinf(I) | isnan(I))=0;
        infokc1 =  I;
        
        infoc1 = infoc1 + ( nanmean(infokc1(~isinf(infokc1))) - infoc1) / k;
        varsum1 = varsum1 + nansum((infokc1(~isinf(infokc1)) - infoc1) .^ 2);
        stderr1 = conf * sqrt(varsum1 / (k * cases * (k * cases - 1)));
        
        info = (infoc1);
        
        in.info(k)=info;
        in.stderr(k)=stderr1;
        
        figure(100)
        hold on
        subplot(1,3,1);plot(k,info,'O')
        title('Info')
        hold on
        subplot(1,3,2);plot(k,stderr1,'*');hold on;plot(k,stderr2,'*');hold on;plot(k,stderr_tot,'O')
        title('stderr')
        hold on
        subplot(1,3,3);plot(k,mean(infokc1),'*')
        title('E(infok)')
        drawnow
    end
    
    
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% OLD OPTIONS BEFORE AUGUST 2018
switch 'condition'
    case 'uni'
        
        switch C
            case 0
                
                infok = zeros(cases,1);
                infokc = zeros(cases,1);
                infokt = zeros(cases,1);
                
                k = k + 1;
                % Generate samples
                
                for i = 1%:ncond
                    
                    clear X pp cc
                    for hh=1:numel(vines.margins)
                        X(:,hh)=vines.margins{hh}.ker;
                        Ma(i,hh)=max(X(:,hh));
                        Mi(i,hh)=min(X(:,hh));
                    end
                end
                
                for hh=1:numel(vines.margins)
                    x=linspace(min(Mi(:,hh)),max(Ma(:,hh)),cases);
                    cc(:,hh)=datasample(x,cases);
                end
                
                V=prod(abs([min(vines.range')-max(vines.range')]));
                
                for j = 1%:ncond
                    [p,~,~,~,pcop] = Fit_vCopula(vines,cc,size(copula,2),[],-1,copula,'rand',[]);
                    pp(:,j)=(pcop);
                end
                
                % Monte-Carlo estimate of information
                
                log2pp=pcop.*log2(pcop)*V;log2pp(p==0)=0;
                I=log2pp  ;
                I(isinf(I) | isnan(I))=NaN;
                infokc1 =  I;
                
                infoc1 = infoc1 + ( nanmean(infokc1(~isinf(infokc1))) - infoc1) / k;
                varsum1 = varsum1 + nansum((infokc1(~isinf(infokc1)) - infoc1) .^ 2);
                stderr1 = conf * sqrt(varsum1 / (k * cases * (k * cases - 1)));
                
                info = (infoc1);
                
                in.info(k)=info;
                in.stderr(k)=stderr1;
                
                figure(100)
                hold on
                subplot(1,3,1);plot(k,info,'O')
                title('Info')
                hold on
                subplot(1,3,2);plot(k,stderr1,'*');hold on;plot(k,stderr1,'*');hold on;plot(k,stderr_tot,'O')
                title('stderr')
                hold on
                subplot(1,3,3);plot(k,mean(infok),'*')
                title('E(infok)')
                drawnow
                
            case 2
                
                
                if length(vines)~=1
                    
                    V=1;
                    clear infok infokc1 infokc2 infokc infokt
                    
                    k = k + 1;
                    % Generate samples
                    
                    
                    
                    
                    
                    for i = 1:ncond
                        ddd=-inf;h22{i}=[];cc{i}=[];dd0=0;
                        while ddd<100
                            
                            clear X pp cc
                            for hh=1:numel(vines{i}.margins)
                                X(:,hh)=vines{i}.margins{hh}.ker;
                                Ma(i,hh)=vines{1}.range(hh,1);%max(X(:,hh));
                                Mi(i,hh)=vines{1}.range(hh,2);%min(X(:,hh));
                            end
                            
                            for hh=1:numel(vines{i}.margins)
                                
                                rng('shuffle', 'v5uniform')
                                
                                x=linspace(min(Mi(:,hh)),max(Ma(:,hh)),1e7);
                                cc{i}(:,hh)=datasample(x,cases);
                            end
                            
                            [h1 h22{i}]=histc(cc{i}(:,nn),tvec);
                            h22{i}=h22{i}+1;
                            
                            if dd0==0
                                D=cc{i};
                            else
                                D=cat(1,D,cc{i});
                            end
                            
                            [h1 h22{i}]=histc(D(:,nn),tvec);
                            h22{i}=h22{i}+1;
                            for tt=2:numel(tvec)-1
                                dd(tt-1)=sum((h22{i}==tt));
                            end
                            ddd=min(dd)
                            dd0=1;
                        end
                        
                        
                        [h1 h22{i}]=histc(D(:,nn),tvec);
                        h22{i}=h22{i}+1;
                        
                        for tt=2:numel(tvec)-1
                            ffi=find(h22{i}==tt);
                            if i==1 & tt==2
                                Dt=D(randsample(ffi,100),:);
                            else
                                Dt=cat(1,Dt,D(randsample(ffi,100),:));
                            end
                        end
                        
                    end
                    D=Dt;
                    [h1 h2]=histc(D(:,nn),tvec);
                    h2=h2+1;
                    
                    
                    
                    VT=prod(abs(vines{1}.range(2,1)'-vines{1}.range(2,2)'));
                    VF=prod(abs(vines{1}.range([1 3:size(vines{1}.range,1)],1)'-vines{1}.range([1 3:size(vines{1}.range,1)],2)'))'*mean(diff(tvec));
                    VB=prod(abs(vines{1}.range([3:size(vines{1}.range,1)],1)'-vines{1}.range([3:size(vines{1}.range,1)],2)'))*mean(diff(tvec));
                    Vn=prod(abs(vines{1}.range(1,1)'-vines{1}.range(1,2)'))*mean(diff(tvec));
                    
                    dims=[1,2];
                    
                    for j = 1:ncond
                        
                        [vineMar{j},copulaMar{j}]=kernelmarginalize(dims,vines{j},copula{j});
                        [pMar,~,~,~,~] = Fit_vCopula(vineMar{j},D(:,dims),size(copulaMar{1},2),[],0,copulaMar{j},'rand',[]);
                        
                        
                        [p,~,~,~,~] = Fit_vCopula(vines{j},D,size(copula{1},2),[],0,copula{j},'rand',[]);
                        
                        par.fit=0;par.s=copula{j}{nn,1}.MarginS.s;par.p=copula{j}{nn,1}.MarginS.p;
                        [pdf(:,j),~]=kernelpdf(D(:,nn),2,D(:,nn),par);
                        
                        par.fit=0;par.s=copula{j}{1,1}.MarginS.s;par.p=copula{j}{1,1}.MarginS.p;
                        [pdfn(:,j),~]=kernelpdf(D(:,1),2,D(:,1),par);
                        
                        
                        %                             pdf(:,j)=1;
                        
                        
                        X=p./pdf(:,j);%./pdf;
                        qq(:,j)=X;   %%% f(B,n|t  C)
                        
                        X=p./pMar;    %%% f(B|t n C)
                        qqc(:,j)=X;
                        
                        X=pMar./pdf(:,j);%./pdf;
                        pp(:,j)=X;   %%% f(n|t C)
                        
                        % Marginalizing the n
                        if mm~=0
                            
                            y_vectorT=linspace(vines{1}.range(mm,1),vines{1}.range(mm,2),200);
                            for sor=1:4
                                
                                disp(['condition = ', num2str(j),' , sor = ',num2str(sor)])
                                y_vector=y_vectorT((sor-1)*50+1:sor*50);
                                ne=0;
                                points=zeros(numel(y_vector)*size(D,1),size(D,2));
                                
                                for i2=1:size(D,1)
                                    for i1=1:numel(y_vector)
                                        ne=ne+1;
                                        Z(mm)=y_vector(i1);
                                        Z(setdiff(1:size(D,2),mm))=D(i2,setdiff(1:size(D,2),mm));
                                        points(ne,:)=Z;
                                    end
                                end
                                
                                [pT,~,~,~,~] = Fit_vCopula(vines{j},points,size(copula{1},2),[],0,copula{j},'rand',[]);
                                
                                ppp=reshape(pT,numel(y_vector),[]);
                                
                                if sor==1
                                    pm=ppp;
                                else
                                    pm=cat(1,pm,ppp);
                                end
                            end
                            
                            pm=sum(pm,1)'*mean(diff(y_vectorT));
                            
                            
                            X=pm./pdf(:,j);    %%% f(B|t C)
                            qqcc(:,j)=X;
                        end
                        mmm=mm;
                    end
                    
                    
                    
                    for tt=2:numel(tvec)-1
                        ff=find(h2==tt);
                        PR=[0.5 0.5];%[sum(h22{1}==tt) sum(h22{2}==tt)]/(sum(h22{1}==tt)+sum(h22{2}==tt));
                        
                        
                        
                        
                        Q=sum(PR .* qq(ff,:),2);
                        
                        pC=sum(PR .* qqc(ff,:),2);
                        pdfnc=sum(PR .* qqcc(ff,:),2);
                        
                        
                        for co=1:ncond
                            if co==1
                                fff=1:200;
                            else
                                fff=1:200;
                            end
                            G=qqc(ff(fff),co);%-qqc(ff(fff),co).*log2(qqc(ff(fff),co));G(qqc(ff(fff),co)==0)=0;
                            I= 0.5*G;
                            I(isinf(I) | isnan(I))=0;
                            H{co,tt}=I * VB;
                            
                            G=qqcc(ff(fff),co);%-qqcc(ff(fff),co).*log2(qqcc(ff(fff),co));G(qqcc(ff(fff),co)==0)=0;
                            In= 0.5*G;
                            In(isinf(In) | isnan(In))=0;
                            H{co+2,tt}=In * VB;
                            
                        end
                        
                        %%%%%%%%%%%
                        
                        G=pC;%-pC.*log2(pC);G(pC==0)=0;
                        I= G;
                        I(isinf(I) | isnan(I))=0;
                        H{5,tt}=I * VB;
                        
                        
                        G=pdfnc;%-pdfnc.*log2(pdfnc);G(pdfnc==0)=0;
                        In= G;
                        In(isinf(In) | isnan(In))=0;
                        H{6,tt}=In * VB;
                        
                        
                        
                        %%%%%%%%%%%
                        
                        
                        pC=sum(PR .* pp(ff,:),2);
                        
                        for co=1:ncond
                            if co==1
                                fff=1:200;
                            else
                                fff=1:200;
                            end
                            G=pp(ff(fff),co);%-pp(ff(fff),co).*log2(pp(ff(fff),co));G(pp(ff(fff),co)==0)=0;
                            I= 0.5*G;
                            I(isinf(I) | isnan(I))=0;
                            H{6+co,tt}=I * Vn;
                            
                        end
                        
                        G=pC;%-pC.*log2(pC);G(pC==0)=0;
                        In= G;
                        In(isinf(In) | isnan(In))=0;
                        H{9,tt}=In * Vn;
                        
                        
                        %%%%%%%%%%
                        
                        pC=sum(PR .* qq(ff,:),2);
                        
                        for co=1:ncond
                            if co==1
                                fff=1:200;
                            else
                                fff=1:200;
                            end
                            G=-qq(ff(fff),co).*log2(qq(ff(fff),co));G(qq(ff(fff),co)==0)=0;
                            I= 0.5*G;
                            I(isinf(I) | isnan(I))=0;
                            H{9+co,tt}=I * VF;
                            
                        end
                        
                        G=-pC.*log2(pC);G(pC==0)=0;
                        G(isinf(G) | isnan(G))=0;
                        H{12,tt}=G * VF;
                        
                        
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        if ff~=0
                            for w=1:12
                                En(w,tt) = En(w,tt)+  (nanmean(H{w,tt}(~isinf(H{w,tt}))) - En(w,tt)) /k;
                                varsum(w,tt) = varsum(w,tt) + nansum((H{w,tt}(~isinf(H{w,tt})) - En(w,tt)).^ 2);
                                stder(w,tt) = conf * sqrt(varsum(w,tt) / (200*k * (200*k - 1)));
                            end
                            
                        end
                        
                        
                    end
                    
                    
                    
                    for tt=2:numel(tvec)-1
                        
                        %%%   I(n;C)
                        
                        in.info(1,tt)=En(9,tt)-1*En(8,tt)-1*En(7,tt);
                        
                        %%%   I(n;B)
                        
                        in.info(2,tt)=En(6,tt)-En(5,tt);
                        
                        %%%   I(B;C)
                        
                        in.info(3,tt)=En(6,tt)-1*En(3,tt)-1*En(4,tt);
                        
                        %%% I(n,B;C)
                        
                        in.info(4,tt)=En(12,tt)-1*En(10,tt)-1*En(11,tt);
                        
                        %%% I(n;C | B)
                        
                        in.info(5,tt)=in.info(4,tt) - in.info(3,tt);
                        
                        %%% I(n;B | C)
                        
                        in.info(6,tt)=in.info(5,tt) + in.info(2,tt) - in.info(1,tt);
                        
                        %%% I(n;B , C)
                        
                        in.info(7,tt)=in.info(5,tt) + in.info(2,tt) ;
                        
                        
                    end
                    
                    
                    stderr1=0;%norm(stderr1T(~isnan(stderr1T) & ~isinf(stderr1T)));
                    stderr2=0;%norm(stderr2T(~isnan(stderr2T) & ~isinf(stderr2T)));
                    stderr_tot=max(stder(:));
                    
                    
                    in.stderr(k)=stderr_tot;
                    
                    figure(100)
                    subplot(1,3,[1])
                    hold on;
                    plot(En(1,:))
                    
                    %                             plot(tvec(2:end),in.info(1,:));hold on;plot(tvec(2:end),in.info(4,:))
                    
                    axis tight
                    title('Info')
                    
                    subplot(1,3,2);hold on;plot(k,stderr_tot,'O')
                    title('stderr')
                    
                    subplot(1,3,3)
                    hold on
                    plot(tvec(2:end),in.info(2,:));hold on;plot(tvec(2:end),in.info(6,:))
                    %                         plot(tvec(2:end),cell2mat(infot)-cell2mat(infotc))
                    hold on
                    %                         plot(tvec(2:end),cell2mat(infotc))
                    axis tight
                    
                    title(num2str(sum(N(1,:),2)))
                    drawnow
                    
                end
                
                %
                %
                %                     t=tvec(2:end);
                %                     t=t(t<250);
                %                     i1=in.info(1,t<250);
                %                     i2=in.info(5,t<250);
                % %
                % %                     ii=linspace(0,0.05,6);
                % %                     i1(1:6)=(i1(1:6)+ii)/2;
                % %                     i2(1:6)=(i2(1:6)+ii)/2;
                %
                %                     figure
                %                     plot(t,smooth(i1,3))
                %                     hold on
                %                     plot(t,smooth(i2+i1/5,3))
                %                     axis tight
                %
                %
                
                
                
        end
        
        
        
        
        
        
        
        
    case 'coppar'
        
        V=1;
        
        infok = zeros(cases,1);
        infokc1 = zeros(cases,1);
        infokc2 = zeros(cases,1);
        infokc = zeros(cases,1);
        infokt = zeros(2*cases,1);
        
        k = k + 1;
        % Generate samples
        
        D = mixedvinernd(vineest,cases);
        
        
        
        [pc,~,~,~,pcor] = Fit_vCopula(vines,D,size(copula,2),[],-1,copula,'rand',[]);
        [p,~,pdf] = mixedvinepdf(vineest,D);
        pcop=p./pdf;
        % Monte-Carlo estimate of information
        
        V=1;
        
        %             for i=1:ncond
        %                 log2pp=log2(pp(:,i));log2pp(pp(:,i)==0)=0;
        %                 I=pcond(1,i)  .* log2pp  ;
        %                 I(isinf(I) | isnan(I))=0;
        %                 infokc = infokc + I;
        %             end
        log2pp=pcor.*log2(pcor)./pcop;log2pp(pcor==0)=0;
        I= log2pp  ;
        I(isinf(I) | isnan(I))=0;
        infokc1 =  I;
        
        
        infoc1 = infoc1 + ( nanmean(infokc1(~isinf(infokc1))) - infoc1) / k;
        
        
        varsum1 = varsum1 + nansum((infokc1(~isinf(infokc1)) - infoc1) .^ 2);
        stderr1 = conf * sqrt(varsum1 / (k * cases * (k * cases - 1)));
        
        
        info = (infoc1);
        
        in.info(k)=info;
        in.stderr(k)=stderr1;
        
        
        figure(100)
        hold on
        subplot(1,3,1);plot(k,info,'O')
        %     hold on;subplot(1,3,1);plot(k,infoc,'*')
        %     hold on;subplot(1,3,1);plot(k,infot,'*')
        title('Info')
        hold on
        subplot(1,3,2);plot(k,stderr1,'*');hold on;plot(k,stderr2,'*');hold on;plot(k,stderr_tot,'O')
        title('stderr')
        hold on
        subplot(1,3,3);plot(k,mean(infokc),'*');%hold on;plot(k,mean(infokt),'O')
        title('E(infok)')
        drawnow
        
    case 'par'
        
        V=1;
        
        infok = zeros(cases,1);
        infokc1 = zeros(cases,1);
        
        k = k + 1;
        % Generate samples
        
        D = mixedvinernd(vines,cases);
        
        [p,~,pdf] = mixedvinepdf(vines,D);%Fit_vCopula(vines{j},cc{j},size(copula{1},2),[],0,copula{j},'rand');
        pcop=p./pdf;
        
        log2pp=log2(pcop);log2pp(pcop==0)=0;
        I=log2pp  ;
        I(isinf(I) | isnan(I))=0;
        infokc1 =  I;
        
        infoc1 = infoc1 + ( nanmean(infokc1(~isinf(infokc1))) - infoc1) / k;
        varsum1 = varsum1 + nansum((infokc1(~isinf(infokc1)) - infoc1) .^ 2);
        stderr1 = conf * sqrt(varsum1 / (k * cases * (k * cases - 1)));
        
        info = (infoc1);
        
        in.info(k)=info;
        in.stderr(k)=stderr1;
        
        figure(100)
        hold on
        subplot(1,3,1);plot(k,info,'O')
        title('Info')
        hold on
        subplot(1,3,2);plot(k,stderr1,'*');hold on;plot(k,stderr2,'*');hold on;plot(k,stderr_tot,'O')
        title('stderr')
        hold on
        subplot(1,3,3);plot(k,mean(infokc1),'*')
        title('E(infok)')
        drawnow
        
        
end