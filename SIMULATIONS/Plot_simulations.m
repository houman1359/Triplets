

    clear all
    close all
    
    di = dir(['./results/s*']);

    for al=1:4
        for ab=1:20
            sdecINF{al,ab}=[];sglmINF{al,ab}=[];scopINF{al,ab}=[];
            n(al,ab)=0;
        end
    end
    for i=3:numel(di)
        Da=load(['./results/',di(i).name]);
        str = di(i).name;
        numbers = regexp(str, '\d+', 'match');
        numbers = cellfun(@str2double, numbers);
        alpha_NN=numbers(3);
        alpha_N=numbers(2);
        alpha=alpha_N;

        n(alpha_NN,alpha)=n(alpha_NN,alpha)+1;

        sdecINF{alpha_NN,alpha}=[sdecINF{alpha_NN,alpha} Da.dec_inf{alpha_N}];
        sglmINF{alpha_NN,alpha}=[sglmINF{alpha_NN,alpha} Da.glm_inf{alpha_N}];
        scopINF{alpha_NN,alpha}=[scopINF{alpha_NN,alpha} Da.cop_inf{alpha_N}];

        sdec{alpha_NN,alpha}(n(alpha_NN,alpha),1:3)=[Da.dec_interaction{alpha_N}(1,2);Da.dec_interaction{alpha_N}(1,3);Da.dec_interaction{alpha_N}(2,3)];
        scop{alpha_NN,alpha}(n(alpha_NN,alpha),1:3)=[Da.cop_interaction{alpha_N}(1,2);Da.cop_interaction{alpha_N}(1,3);Da.cop_interaction{alpha_N}(2,3)];
        sglm{alpha_NN,alpha}(n(alpha_NN,alpha),1:3)=[Da.glm_interaction{alpha_N}(1,2);Da.glm_interaction{alpha_N}(1,3);Da.glm_interaction{alpha_N}(2,3)];

        strip_dec{alpha_NN,alpha}(n(alpha_NN,alpha))=sum([Da.dec_interaction{alpha_N}(1,2);Da.dec_interaction{alpha_N}(1,3);Da.dec_interaction{alpha_N}(2,3)]>0);
        strip_cop{alpha_NN,alpha}(n(alpha_NN,alpha))=sum([Da.cop_interaction{alpha_N}(1,2);Da.cop_interaction{alpha_N}(1,3);Da.cop_interaction{alpha_N}(2,3)]>0);
        strip_glm{alpha_NN,alpha}(n(alpha_NN,alpha))=sum([Da.glm_interaction{alpha_N}(1,2);Da.glm_interaction{alpha_N}(1,3);Da.glm_interaction{alpha_N}(2,3)]>0);

    end


    alpha_LIST=linspace(0,2,20);
    for nonlinear=1:2
        for i=1:20
            decInt(i)=mean(sdec{nonlinear,i}(:));
            copInt(i)=mean(scop{nonlinear,i}(:));
            glmInt(i)=mean(sglm{nonlinear,i}(:));
            sdecInt(i)=std(sdec{nonlinear,i}(:))/sqrt(numel(sdec{nonlinear,i}))*2;
            scopInt(i)=std(scop{nonlinear,i}(:))/sqrt(numel(sdec{nonlinear,i}))*2;
            sglmInt(i)=std(sglm{nonlinear,i}(:))/sqrt(numel(sdec{nonlinear,i}))*2;
        end

        hk=figure(66)
        subplot(2,4,nonlinear+2)
        errorbar(alpha_LIST,smooth(decInt),smooth(sdecInt),'k','linewidth',3)
        hold on
        errorbar(alpha_LIST,smooth(glmInt),smooth(sglmInt),'r','linewidth',3)
        hold on
        errorbar(alpha_LIST,smooth(copInt),smooth(scopInt),'b','linewidth',3)
        yline(0,'linewidth',3)
        if nonlinear ==1
            title('f_i(b)=\{x,y,x+y\}')
            ylabel('Pairwise interaction information (bits)')
            legend({'I^{int}(r^0_1,r^0_2;S)','I^{int}_{GLM}(r_1,r_2;S|B)','I^{int}_{NPvC}(r_1,r_2;S|B)'},'Location','SouthWest')
        else
            title('f_i(b)=\{x^2,y^2,xy\}')
        end
        ylim([-0.024 0.022])
        xlabel('\alpha weight')
    end

    for nonlinear=1:2
        for i=1:20
            decInt(i)=mean(sdecINF{nonlinear,i});
            copInt(i)=mean(scopINF{nonlinear,i});
            glmInt(i)=mean(sglmINF{nonlinear,i});
            sdecInt(i)=std(sdecINF{nonlinear,i})/sqrt(numel(sdecINF{nonlinear,i}))*2;
            scopInt(i)=std(scopINF{nonlinear,i})/sqrt(numel(sdecINF{nonlinear,i}))*2;
            sglmInt(i)=std(sglmINF{nonlinear,i})/sqrt(numel(sdecINF{nonlinear,i}))*2;
        end

        figure(66)
        subplot(2,4,nonlinear)
        errorbar(alpha_LIST,smooth(decInt),smooth(sdecInt),'k','linewidth',3)
        hold on
        errorbar(alpha_LIST,smooth(glmInt),smooth(sglmInt),'r','linewidth',3)
        hold on
        errorbar(alpha_LIST,smooth(copInt),smooth(scopInt),'b','linewidth',3)
        if nonlinear ==1
            title('f_i(b)=\{x,y,x+y\}')
            ylabel('Single neuron information (bits)')
            legend({'I(r^0;S)','I_{GLM}(r;S|B)','I_{NPvC}(r;S|B)'},'Location','SouthWest')

        else
            title('f_i(b)=\{x^2,y^2,xy\}')
        end
        ylim([0 0.12])
        xlabel('\alpha weight')
    end



    for nl=1:2
        for al=1:20

            %clear dec_1 dec_2 dec_3 dec_4 cop_1 cop_2 cop_3 cop_4 glm_1 glm_2 glm_3 glm_4
            for i=1:1000
                F=randsample(1:numel(strip_dec{nl,al}),round(0.5*numel(strip_dec{nl,al})));
                dec_1(i,nl,al)=mean(strip_dec{nl,al}(F)==1);
                dec_2(i,nl,al)=mean(strip_dec{nl,al}(F)==3);
                dec_3(i,nl,al)=mean(strip_dec{nl,al}(F)==0);
                dec_4(i,nl,al)=mean(strip_dec{nl,al}(F)==2);
                glm_1(i,nl,al)=mean(strip_glm{nl,al}(F)==1);
                glm_2(i,nl,al)=mean(strip_glm{nl,al}(F)==3);
                glm_3(i,nl,al)=mean(strip_glm{nl,al}(F)==0);
                glm_4(i,nl,al)=mean(strip_glm{nl,al}(F)==2);
                cop_1(i,nl,al)=mean(strip_cop{nl,al}(F)==1);
                cop_2(i,nl,al)=mean(strip_cop{nl,al}(F)==3);
                cop_3(i,nl,al)=mean(strip_cop{nl,al}(F)==0);
                cop_4(i,nl,al)=mean(strip_cop{nl,al}(F)==2);


                dec_trip(i,nl,al)=mean(strip_dec{nl,al}(F)==3);
                glm_trip(i,nl,al)=mean(strip_glm{nl,al}(F)==3);
                cop_trip(i,nl,al)=mean(strip_cop{nl,al}(F)==3);
            end

        end
    end
    %set(gcf,'PaperPositionMode','Auto','Units','normalized','Position',[0 0 1 1]/2);

    hk=figure(66)
    nj=0;
    for i=2
        for j=[1 19]
            nj=nj+1;
            subplot(2,4,nj+6)
            ba=bar([mean(dec_1(:,i,j)) mean(dec_2(:,i,j)) mean(dec_3(:,i,j)) mean(dec_4(:,i,j));mean(glm_1(:,i,j)) mean(glm_2(:,i,j)) mean(glm_3(:,i,j)) mean(glm_4(:,i,j));mean(cop_1(:,i,j)) mean(cop_2(:,i,j)) mean(cop_3(:,i,j)) mean(cop_4(:,i,j))]', 'FaceColor','flat');
            ba(1).CData = [0 0 0];
            ba(2).CData = [1 0 0];
            ba(3).CData = [0 0 1];
            if j==1
                ylabel('Triplet probability')
                legend({'Groun truth','GLM','NPvC'},'Location','SouthWest')
                title('f_i(b)=\{x^2,y^2,xy\}, \alpha=0')
            else
                title('f_i(b)=\{x^2,y^2,xy\}, \alpha=2')
            end
            xticklabels({'(-,-,+)','(+,+,+)','(-,-,-)','(-,+,+)'})
            h=xline(1,'lineStyle','--','color',[1 1 1]/3,'linewidth',3,'alpha',0.3)
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h=xline(2,'lineStyle','--','color',[1 1 1]/3,'linewidth',3,'alpha',0.3)
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h=xline(3,'lineStyle','--','color',[1 1 1]/3,'linewidth',3,'alpha',0.3)
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h=xline(4,'lineStyle','--','color',[1 1 1]/3,'linewidth',3,'alpha',0.3)
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end


    for i=1:2
        subplot(2,4,4+i)
        errorbar(alpha_LIST,smooth(squeeze(mean(dec_trip(:,i,:)))),smooth(squeeze(std(dec_trip(:,i,:)))),'k','linewidth',3)
        hold on
        errorbar(alpha_LIST,smooth(squeeze(mean(glm_trip(:,i,:)))),smooth(squeeze(std(glm_trip(:,i,:)))),'r','linewidth',3)
        hold on
        errorbar(alpha_LIST,smooth(squeeze(mean(cop_trip(:,i,:)))),smooth(squeeze(std(cop_trip(:,i,:)))),'b','linewidth',3)
        ylim([0 1])
        if i==1

        end
        if i ==1
            title('f_i(b)=\{x,y,x+y\}')
            legend({'Groun truth','GLM','NPvC'},'Location','SouthWest')
            ylabel('Probability of (+,+,+) triplet')
        else
            title('f_i(b)=\{x^2,y^2,xy\}')
        end
        xlabel('\alpha weight')
    end


    set(gcf,'PaperPositionMode','Auto','Units','normalized','Position',[0 0 1 1]/2);

