function Simulation_Pairwise_git(N_or,NL,alpha_N)

%%% October 2023
%%%% the code first simulats activities of three neurons for 2 stimuli and
%%%% compute the conditional information and interaciton informations using
%%%% a decoder, NPvC, and GLM


%%% inputs
%%% N_or:  simmmulation index
%%% NL : 1=linear model, 2=nonlinear
%%% alpha_N: nonlinearity weight which is between 0 and 2

%%% outputs
%     %%% conditional informations computed using the decoder, GLM, and the
%     %%% copula
%     dec_inf{alpha_N}
%     glm_inf{alpha_N}
%     cop_inf{alpha_N}
% 
%     %%% interactions informations computed using the decoder, GLM, and the
%     %%% copula
%     dec_interaction{alpha_N}
%     glm_interaction{alpha_N}
%     cop_interaction{alpha_N}





addpath(genpath('../NPvC/glmnet'))
addpath(genpath('../NPvC/Ker_Copula'))
rng(N_or)


    N_sample = 1000;
    mu=[5 6]+5;
    for S=1:2
        TH=0;
        while TH==0
            muu = [0 0 0];
            ro=-0.1;
            Sigma = [1 ro ro;ro 1 ro;ro ro 1];
            x_s{S}=mvnrnd(muu,Sigma,N_sample);
            c(1)=corr(x_s{S}(:,1),x_s{S}(:,2));
            c(2)=corr(x_s{S}(:,1),x_s{S}(:,3));
            c(3)=corr(x_s{S}(:,3),x_s{S}(:,2));
            if sum(abs(c-ro)<0.001)==3
                TH=1;
            end
        end

        muuu = [1 1];
        ro=0.8;
        Sigma = [1 ro;ro 1]/2;
        B{S}=mvnrnd(muuu,Sigma,N_sample);
        xx=B{S}(:,1)';
        yy=B{S}(:,2)';        
            
        alpha_LIST=linspace(0,2,20);
        alpha=alpha_LIST(alpha_N);

        B{S}=(cat(2,xx',yy'))+1;

        if NL==1 
            INP{S}=cat(2,xx',yy',xx'+yy');
        elseif NL==2
            INP{S}=cat(2,xx'.^2,yy'.^2,xx'.*yy');
        end

        INP{S}=alpha*(INP{S}-mean(INP{S}));
        x_samp{S}=x_s{S}+INP{S};
        x_samp{S}=x_samp{S}-mean(x_samp{S})+mu(S);
        x_s{S}=x_s{S}-mean(x_s{S})+mu(S);

        if N_or==1 & S==2
            DAT{alpha_N}.x_sampe=x_samp;
            DAT{alpha_N}.x_s=x_s;
            DAT{alpha_N}.INP=INP;
            DAT{alpha_N}.B=B;
        end

    end

    pll=0;
    if pll==1
        i=0;s1=[2 3 6];s2=[4 7 8];
        figure(50)
        for n1=1:2
            for n2=n1+1:3
                i=i+1;
                subplot(5,3,i+6)
                plot(x_samp{1}(:,n1),x_samp{1}(:,n2),'.')
                hold on
                plot(x_samp{2}(:,n1),x_samp{2}(:,n2),'.r')
                axis square
                title('data')
                subplot(5,3,i)
                plot(x_s{1}(:,n1),x_s{1}(:,n2),'.')
                hold on
                plot(x_s{2}(:,n1),x_s{2}(:,n2),'.r')
                title('conditioned')
                axis square
                subplot(5,3,i+3)
                plot(INP{1}(:,n1),INP{1}(:,n2),'.')
                hold on
                plot(INP{2}(:,n1),INP{2}(:,n2),'.r')
                title('B-input')
                axis square
            end
        end
    end
figure;plot(xx,yy,'.');axis equal;axis square



    for n1=1:2
        for n2=n1+1:3
            ne1=cat(1,x_s{1}(:,n1),x_s{2}(:,n1));
            ne2=cat(1,x_s{1}(:,n2),x_s{2}(:,n2));
            Y=cat(2,ne1,ne2);
            stimulus=cat(1,ones(size(x_s{1}(:,n1))),zeros(size(x_s{2}(:,n1))));
            Mdl = fitcdiscr(cat(2,ne1,ne2),stimulus);
            decod = predict(Mdl,Y);
            info_int_dec(n1,n2)=CONF_INFO(stimulus,decod);

            for sh=1:1000
                ne2=cat(1,x_s{1}(randsample(1:size(x_s{1},1),size(x_s{1},1)),n2),x_s{2}(randsample(1:size(x_s{2},1),size(x_s{1},1)),n2));
                Mdl = fitcdiscr(cat(2,ne1,ne2),stimulus);
                decod = predict(Mdl,cat(2,ne1,ne2));
                infoDecsh(sh)=CONF_INFO(stimulus,decod);
            end
            info_ind_dec(n1,n2)=mean(infoDecsh);
        end
    end

    for n1=1:3
        ne1=cat(1,x_s{1}(:,n1),x_s{2}(:,n1));
        stimulus=cat(1,ones(size(x_s{1}(:,n1))),zeros(size(x_s{2}(:,n1))));
        Mdl = fitcdiscr(ne1,stimulus);
        decod = predict(Mdl,ne1);
        I=CONF_INFO(stimulus,decod);

        for sh=1:1000
            ne1=cat(1,x_s{1}(:,n1),x_s{2}(:,n1));
            Mdl = fitcdiscr(ne1,stimulus(randsample(1:numel(stimulus),numel(stimulus))));
            decod = predict(Mdl,ne1);
            infoDecsh(sh)=CONF_INFO(stimulus(randsample(1:numel(stimulus),numel(stimulus))),decod);
        end

        InfoDec(n1)=I-1*mean(infoDecsh);
    end

    clear y
    for neu=1:3
        for stim=1:2

            X=cat(2,x_samp{stim}(:,neu),B{stim});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% design the vine

            clear range iscont iscons families Z
            for d=1:size(X,2)
                margins{d}='kernel';
                if d==1
                    iscons(d)=0;
                else
                    iscons(d)=0;
                end
            end
            range(:,1)=min(X);
            range(:,2)=max(X);
            for i=1:size(X,2)
                for j=1:size(X,2)
                    families{i,j}='kercop';
                end
            end
            vine{neu,stim}=prep_copula(X,margins,families,iscons,'c-vine','rand',range);
            mm=1;
            y{neu,stim}=X(:,mm);
            z{neu,stim}=X;
            Xpred=X(:,setdiff(1:size(X,2),mm));

            %%%%%%%% points: the set of points which we want to compute the
            %%%%%%%% density and use to integrate to compute conditional CDF
            %%%%%%%% and also density over data points

            knots=100;
            grid_size = 200;
            y_vector{neu} =linspace(min(cat(1,x_samp{1}(:,neu),x_samp{2}(:,neu)))-eps,max(cat(1,x_samp{1}(:,neu),x_samp{2}(:,neu)))+eps,grid_size);
            dy=mean(diff(y_vector{neu}));
            points{neu,stim}=zeros(numel(y_vector{neu})*size(Xpred,2),size(Xpred,2)+1);
            set_train=1:size(Xpred,1);
            ne=0;
            for j=1:numel(set_train)
                for i=1:numel(y_vector{neu})
                    ne=ne+1;
                    Z(mm)=y_vector{neu}(i);
                    Z(setdiff(1:size(Xpred,2)+1,mm))=X(set_train(j),setdiff(1:size(Xpred,2)+1,mm));
                    points{neu,stim}(ne,:)=Z;
                end
            end
            for i = 1:d
                iscont(i) = true;
            end
            vine{neu,stim}.condition=0;
            for i=1:d
                for j=1:d
                    vine{neu,stim}.METH{i,j}='LL1';
                end
            end

            %%%%% fiting copula over data and buiding the vine copula trained
            %%%%% on the data to be used later to estimate density function on
            %%%%% data vectors
            [f_pointR{neu,stim},f_data1{neu,stim},copula{neu,stim},~] = Fit_vCopula(vine{neu,stim},points{neu,stim}(1,:),20,'LL1',1,0,'rand',[],knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package
            [f_pointR{neu,stim},f_data1{neu,stim},copula{neu,stim},~] = Fit_vCopula(vine{neu,stim},points{neu,stim}(1,:),20,'LL1',-1,copula{neu,stim},'rand',[],knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package

            disp(['Vine for neuron ',num2str(neu),' ,Stimulus ',num2str(stim),' fitted' ])
        end
    end

    for neu=1:3
        for stim=1:2

            %%% estimating densities using the pre-trained vine copula
            %%% bivariates
            [f_points_train{neu,1}((stim-1)*N_sample*numel(y_vector{neu})+1:stim*N_sample*numel(y_vector{neu})),~,~,~] = Fit_vCopula(vine{neu,1},points{neu,stim},20,[],-1,copula{neu,1},'rand',[],knots);
            [f_points_train{neu,2}((stim-1)*N_sample*numel(y_vector{neu})+1:stim*N_sample*numel(y_vector{neu})),~,~,~] = Fit_vCopula(vine{neu,2},points{neu,stim},20,[],-1,copula{neu,2},'rand',[],knots);

            f_points{neu,stim}=f_points_train{neu,stim}((stim-1)*N_sample*numel(y_vector{neu})+1:stim*N_sample*numel(y_vector{neu}));
            f_vecXB=reshape(f_points{neu,stim},numel(y_vector{neu}),[])+1e-16;
            f_B=sum(f_vecXB*dy)+1e-16;
            f_vecXconB=f_vecXB./repmat(f_B,size(f_vecXB,1),1);
            CDF_vecXconB=cumsum(f_vecXconB)*dy;
            for i=1:numel(y{neu,stim})
                f_XconB{neu}(i+N_sample*(stim-1)) = interp1(y_vector{neu} , f_vecXconB(:,i) , y{neu,stim}(i),'linear');
                CDF_XconB{neu}(i+N_sample*(stim-1)) = interp1(y_vector{neu} , CDF_vecXconB(:,i) , y{neu,stim}(i),'nearest');%,'extrap'
            end

            [y_mle{neu,stim},y_em{neu,stim},LL]=predict_response(f_vecXB,y_vector{neu},y{neu,stim});
            yt{neu,stim}=y{neu,stim};
            CORR(stim,neu)=corr(yt{neu,stim}(~isnan(y_em{neu,stim}(:))),y_em{neu,stim}');
        end


        f_vecXB_1=reshape(f_points_train{neu,1},numel(y_vector{neu}),[])+1e-16*rand(size(reshape(f_points_train{neu,1},numel(y_vector{neu}),[])));
        f_vecXB_2=reshape(f_points_train{neu,2},numel(y_vector{neu}),[])+1e-16*rand(size(reshape(f_points_train{neu,1},numel(y_vector{neu}),[])));
        f_B_1=sum(f_vecXB_1*dy)+1e-16*rand(size(sum(f_vecXB_1)));
        f_B_2=sum(f_vecXB_2*dy)+1e-16*rand(size(sum(f_vecXB_1)));
        f_vecXconB_1=f_vecXB_1./repmat(f_B_1,size(f_vecXB_1,1),1);
        f_vecXconB_2=f_vecXB_2./repmat(f_B_2,size(f_vecXB_2,1),1);
        CDF_vecXconB_1=cumsum(f_vecXconB_1)*dy;
        CDF_vecXconB_2=cumsum(f_vecXconB_2)*dy;

        for stim=1:2
            [glikliCop_1{neu}((stim-1)*N_sample+1:stim*N_sample),~,~,~] = Fit_vCopula(vine{neu,1},z{neu,stim},20,[],-1,copula{neu,1},'rand',[],knots);
            [glikliCop_2{neu}((stim-1)*N_sample+1:stim*N_sample),~,~,~] = Fit_vCopula(vine{neu,2},z{neu,stim},20,[],-1,copula{neu,2},'rand',[],knots);

            for i=1:numel(y{neu,stim})
                likliCop_1{neu}(i+N_sample*(stim-1)) = interp1(y_vector{neu} , f_vecXconB_1(:,i+(stim-1)*numel(y{neu,stim})) , y{neu,stim}(i),'linear');
                CDF_XconB_1{neu}(i+N_sample*(stim-1)) = interp1(y_vector{neu} , CDF_vecXconB_1(:,i+(stim-1)*numel(y{neu,stim})) , y{neu,stim}(i),'nearest');%,'extrap'
                likliCop_2{neu}(i+N_sample*(stim-1)) = interp1(y_vector{neu} , f_vecXconB_2(:,i+(stim-1)*numel(y{neu,stim})) , y{neu,stim}(i),'linear');
                CDF_XconB_2{neu}(i+N_sample*(stim-1)) = interp1(y_vector{neu} , CDF_vecXconB_2(:,i+(stim-1)*numel(y{neu,stim})) , y{neu,stim}(i),'nearest');%,,'extrap'
            end
        end

        %%%% computing informations by decoding the trial type using the
        %%%% likelihoods

        glikliCop_1{neu}=glikliCop_1{neu}+1e-16*rand(size(glikliCop_1{neu}));
        glikliCop_2{neu}=glikliCop_2{neu}+1e-16*rand(size(glikliCop_2{neu}));
        decoded = double(glikliCop_1{neu}>glikliCop_2{neu});
        decoded(decoded==0)= -1;
        stimulus = cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1));
        InfoCopnB(neu)=CONF_INFO(stimulus,decoded);
        for sh=1:1000
        InfoCopnBsh(neu)=CONF_INFO(randsample(stimulus,numel(stimulus)),decoded);          
        end

        decoded = double(f_B_1>f_B_2);
        decoded(decoded==0)= -1;
        stimulus = cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1));
        InfoCopB(neu)=CONF_INFO(stimulus,decoded);
        for sh=1:1000
        InfoCopBsh(neu)=CONF_INFO(randsample(stimulus,numel(stimulus)),decoded);          
        end
        InfoCop(neu)=(InfoCopnB(neu)-InfoCopB(neu))-(InfoCopnBsh(neu)-InfoCopBsh(neu))



    end

    if pll==1
        figure(50)
        i=0;
        for n1=1:2
            for n2=n1+1:3
                i=i+1;
                subplot(5,3,i+9)
                plot(y_em{n1,1},y_em{n2,1},'.')
                hold on
                plot(y_em{n1,2},y_em{n2,2},'.r')
                title('copula fits')
                axis square
            end
        end
    end

%%%%% computing copula interaction for three pairs

    for n1=1:2
        for n2=n1+1:3

            for sh=1:1000
                trs1 = cat(2,randsample(1:N_sample,N_sample),randsample(N_sample+1:2*N_sample,N_sample));
                trs2 = cat(2,randsample(1:N_sample,N_sample),randsample(N_sample+1:2*N_sample,N_sample));
                li_1 = likliCop_1{n1}(trs1) .* likliCop_1{n2}(trs2);
                li_2 = likliCop_2{n1}(trs1) .* likliCop_2{n2}(trs2);
                li_1 = li_1 ;%+ 1e-16 * rand(size(li_1));
                li_2 = li_2 ;%+ 1e-16 * rand(size(li_2));

                decoded = double(li_1>li_2);
                decoded(decoded==0)= -1;
                stimulus = cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1));
                info_s(sh)=CONF_INFO(stimulus,decoded);
            end

            info_sh_cop(n1,n2)=mean(info_s);

            for stim=1:2
                X=cat(1,CDF_XconB{n1}(1+N_sample*(stim-1):N_sample*stim),CDF_XconB{n2}(1+N_sample*(stim-1):N_sample*stim))';
                clear vine iscons iscont families range
                for d=1:size(X,2)
                    margins{d}='kernel';
                    if d==1
                        iscons(d)=0;
                    else
                        iscons(d)=0;
                    end
                end
                range(:,1)=min(X);
                range(:,2)=max(X);
                for i=1:size(X,2)
                    for j=1:size(X,2)
                        families{i,j}='kercop';
                    end
                end
                vine=prep_copula(X,margins,families,iscons,'c-vine','rand',range);
                for i = 1:d
                    iscont(i) = true;
                end
                vine.condition=0;
                for i=1:d
                    for j=1:d
                        vine.METH{i,j}='LL1';
                    end
                end
                [f_p,f_pair,copula_pair,~] = Fit_vCopula(vine,X(1,:),20,'LL1',1,0,'rand',[],knots);
                if stim==1
                    X=cat(1,CDF_XconB_1{n1},CDF_XconB_1{n2})';
                    [f_p_1,f_pair,~,~,copula_pair_1] = Fit_vCopula(vine,X,20,'LL1',-1,copula_pair,'rand',[],knots);
                else
                    X=cat(1,CDF_XconB_2{n1},CDF_XconB_2{n2})';
                    [f_p_2,f_pair,~,~,copula_pair_2] = Fit_vCopula(vine,X,20,'LL1',-1,copula_pair,'rand',[],knots);
                end
            end

            li_1 = likliCop_1{n1} .* likliCop_1{n2} .* copula_pair_1';
            li_2 = likliCop_2{n1} .* likliCop_2{n2} .* copula_pair_2';
            li_1 = li_1 ;%+ 1e-16 * rand(size(li_1));
            li_2 = li_2 ;%+ 1e-16 * rand(size(li_2));

            decoded = double(li_1>li_2);
            decoded(decoded==0)= -1;
            stimulus = cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1));
            info_int_cop(n1,n2)=CONF_INFO(stimulus,decoded);

            disp(['Copula interaction for neuron pairs ',num2str(n1),' & ',num2str(n2),' computed' ])

        end
    end
    clear vine y


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   GLM GAUSSIAN

    for neu=1:3
        X=cat(2,cat(1,x_samp{1}(:,neu),x_samp{2}(:,neu)),cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1)),cat(1,B{1},B{2}));
        y=X(:,1);
        Xpred=X(:,2:end);

        nFolds = 10;
        trialIDs = 1:size(Xpred,1);
        trialPartition = cvpartition(length(trialIDs), 'KFold', nFolds);
        foldIDs = nan(size(Xpred,1),1);
        for nFold = 1:nFolds
            foldTrials = trialIDs(find(test(trialPartition,nFold)));
            foldFrames = find(ismember(1:size(Xpred,1),foldTrials));
            foldIDs(foldFrames) = nFold;
        end

        opt = glmnetSet;
        opt.alpha = 1/2;
        opt.thresh = 1e-6;

        GL{neu} = cvglmnet(Xpred,y,'gaussian',opt,'deviance',[],foldIDs);
        GLM_predictionF{neu} = cvglmnetPredict(GL{neu},Xpred,[],'response');

        CORR_F{neu}=[corr(y(1:N_sample),GLM_predictionF{neu}(1:N_sample)) corr(y(N_sample+1:2*N_sample),GLM_predictionF{neu}(N_sample+1:2*N_sample))]

        in_ls=find(GL{neu}.lambda==GL{neu}.lambda_1se);

        be=GL{neu}.glmnet_fit.beta(:,in_ls);
        a0=GL{neu}.glmnet_fit.a0(in_ls);

        %%%%%% Gaussian likelihood
        %%% compute the sigma2 for S=1
        X=cat(2,x_samp{1}(:,neu),ones(size(B{1},1),1),B{1});
        y=X(:,1);
        Xpred=X(:,2:end);
        mu=Xpred*be+a0;
        sigma2_1 = 1/(numel(mu)-size(Xpred,2))*sum((mu-y).^2);
        
        %%% f(n_i | b_i S=1)
        X=cat(2,cat(1,x_samp{1}(:,neu),x_samp{2}(:,neu)),cat(1,ones(size(B{1},1),1),1*ones(size(B{2},1),1)),cat(1,B{1},B{2}));        
        y=X(:,1);
        Xpred=X(:,2:end);
        mu=Xpred*be+a0;
        likli_1{neu} = 1/(sqrt(sigma2_1 * 2 * pi)) * exp(-(mu-y).^2/(2 * sigma2_1));

        
        %%% compute the sigma2 for S=2
        X=cat(2,x_samp{2}(:,neu),-1 * ones(size(B{2},1),1),B{2});
        y=X(:,1);
        Xpred=X(:,2:end);
        mu=Xpred*be+a0;
        sigma2_2 = 1/(numel(mu)-size(Xpred,2))*sum((mu-y).^2);

        %%% f(n_i | b_i S=1)
        X=cat(2,cat(1,x_samp{1}(:,neu),x_samp{2}(:,neu)),cat(1,-1*ones(size(B{1},1),1),-1*ones(size(B{2},1),1)),cat(1,B{1},B{2}));
        y=X(:,1);
        Xpred=X(:,2:end);
        mu=Xpred*be+a0;
        likli_2{neu} = 1/(sqrt(sigma2_2 * 2 * pi)) * exp(-(mu-y).^2/(2 * sigma2_2));



        decoded = double(likli_1{neu}>likli_2{neu});
        decoded(decoded==0)= -1;
        stimulus = cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1));
        I=CONF_INFO(stimulus,decoded);
        for sh=1:1000
            decoded = double(likli_1{neu}(randsample(1:numel(likli_1{neu}),numel(likli_1{neu})))>likli_2{neu}(randsample(1:numel(likli_2{neu}),numel(likli_2{neu}))));
            decoded(decoded==0)= -1;
            InfoGlmsh(sh)=CONF_INFO(stimulus,decoded);
        end

        InfoGlm(neu)=I-1*mean(InfoGlmsh)

    end


    
    if pll==1
        figure(50)
        i=0;
        for n1=1:2
            for n2=n1+1:3
                i=i+1;
                subplot(5,3,i+12)
                plot(GLM_predictionF{n1}(1:N_sample),GLM_predictionF{n2}(1:N_sample),'.')
                hold on
                plot(GLM_predictionF{n1}(1+N_sample:2*N_sample),GLM_predictionF{n2}(1+N_sample:2*N_sample),'.r')
                title('copula fits')
                axis square
            end
        end
    end


    for n1=1:2
        for n2=n1+1:3

            li_1 = likli_1{n1} .* likli_1{n2};
            li_2 = likli_2{n1} .* likli_2{n2};
            li_1 = li_1 ;%+ 1e-16 * rand(size(li_1));
            li_2 = li_2 ;%+ 1e-16 * rand(size(li_2));

            decoded = double(li_1>li_2);
            decoded(decoded==0)= -1;
            stimulus = cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1));
            info_int_glm(n1,n2)=CONF_INFO(stimulus,decoded);

            for sh=1:1000
                trs1 = cat(2,randsample(1:N_sample,N_sample),randsample(N_sample+1:2*N_sample,N_sample));
                trs2 = cat(2,randsample(1:N_sample,N_sample),randsample(N_sample+1:2*N_sample,N_sample));
                li_1 = likli_1{n1}(trs1) .* likli_1{n2}(trs2);
                li_2 = likli_2{n1}(trs1) .* likli_2{n2}(trs2);
                li_1 = li_1 ;%+ 1e-16 * rand(size(li_1));
                li_2 = li_2 ;%+ 1e-16 * rand(size(li_2));

                decoded = double(li_1>li_2);
                
                decoded(decoded==0)= -1;
                stimulus = cat(1,ones(size(B{1},1),1),-1*ones(size(B{2},1),1));
                info_s(sh)=CONF_INFO(stimulus,decoded);
            end

            info_ind_glm(n1,n2)=mean(info_s);

        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%% conditional informations computed using the decoder, GLM, and the
    %%% copula
    dec_inf{alpha_N} = InfoDec;
    glm_inf{alpha_N} = InfoGlm;
    cop_inf{alpha_N} = InfoCop;

    %%% interactions informations computed using the decoder, GLM, and the
    %%% copula
    dec_interaction{alpha_N} = info_int_dec-info_ind_dec;
    glm_interaction{alpha_N} = info_int_glm-info_ind_glm;
    cop_interaction{alpha_N} = info_int_cop-info_sh_cop;


if N_or==1 & NL==1 & alpha_N==1
    save(['results/simN_',num2str(N_or),'_alpha_',num2str(alpha_N),'_NL_',num2str(NL),'.mat'],'dec_inf','glm_inf','cop_inf','dec_interaction','cop_interaction','glm_interaction','GLM_predictionF','y_em','x_samp','x_s','INP','xx','yy')
else

    save(['results/simN_',num2str(N_or),'_alpha_',num2str(alpha_N),'_NL_',num2str(NL),'.mat'],'dec_inf','glm_inf','cop_inf','dec_interaction','cop_interaction','glm_interaction')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
