%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
load('Pool_Data.mat')

%%%%% The POOL_Data.mat file contains three variables NC, SC, INT which are noise
%%%%% correlation, signal correlation, and interaction values respectively
%%%%% Each of them are cells like NC{L,C} where L referes to being either
%%%%% nonlabeles (L=1) or projection pairs (L=9) and C referes to the
%%%%% correctness (C=1 for correct trials and C=2 for incorrect trials)



%%%% nc for choice pools within and between 

    for label=[1 9]
        for co=1:2

            diff_Wnc(:,label,co)=bootstrp(1000,@nanmean,NC{label,co}(SC{label,1}>0))-mean(bootstrp(1000,@nanmean,NC{label,co}));
            diff_Bnc(:,label,co)=bootstrp(1000,@nanmean,NC{label,co}(SC{label,1}<0))-mean(bootstrp(1000,@nanmean,NC{label,co}));

            Wnc(:,label,co)=bootstrp(1000,@nanmean,NC{label,co}(SC{label,1}>0));
            Bnc(:,label,co)=bootstrp(1000,@nanmean,NC{label,co}(SC{label,1}<0));

        end
    end

    hk1=figure(1)
    for co=1:2
        subplot(2,2,co)
        bar(1,squeeze(mean(diff_Bnc(:,1,co))),'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(2,squeeze(mean(diff_Wnc(:,1,co))),'facecolor',[0 0 0],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        bar(3.5,squeeze(mean(diff_Bnc(:,9,co))),'facecolor',[1 0 1],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(4.5,squeeze(mean(diff_Wnc(:,9,co))),'facecolor',[1 0 1],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        errorbar([1 2 3.5 4.5],[squeeze(mean(diff_Bnc(:,1,co))) squeeze(mean(diff_Wnc(:,1,co))) squeeze(mean(diff_Bnc(:,9,co))) squeeze(mean(diff_Wnc(:,9,co)))],[squeeze(std(diff_Bnc(:,1,co)))  squeeze(std(diff_Wnc(:,1,co))) squeeze(std(diff_Bnc(:,9,co))) squeeze(std(diff_Wnc(:,9,co)))]/2,'.')
        xlim([0.5 5])
        ylim([-9 9]*1e-3)
        set(gca, 'XTick', []);
        ax = gca;
        tickPositions = [1 2 3.5 4.5];
        colors = {'blue', 'red', 'blue', 'red'};
        ticks={'Between', 'Within', 'Between', 'Within'};
        for i = 1:length(tickPositions)
            text(tickPositions(i), ax.YLim(1)-0.01/5, ticks{i}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', colors{i}, 'Rotation', 45); % Added 'Rotation' property for vertical text
        end
        box off
        if co==1
            ylabel('relative noise corr (choice pools)','fontsize',11)
            title('correct trials')
        else
            title('incorrect trials')
        end
        
    end
    set(hk1,'PaperPositionMode','Auto','Units','normalized','Position',[0 0 1 1]/3);

        
    hk1=figure(1)
    for co=1:2
        subplot(2,2,co+2)
        bar(1,squeeze(mean(Bnc(:,1,co))),'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(2,squeeze(mean(Wnc(:,1,co))),'facecolor',[0 0 0],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        bar(3.5,squeeze(mean(Bnc(:,9,co))),'facecolor',[1 0 1],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(4.5,squeeze(mean(Wnc(:,9,co))),'facecolor',[1 0 1],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        errorbar([1 2 3.5 4.5],[squeeze(mean(Bnc(:,1,co))) squeeze(mean(Wnc(:,1,co))) squeeze(mean(Bnc(:,9,co))) squeeze(mean(Wnc(:,9,co)))],[squeeze(std(diff_Bnc(:,1,co)))  squeeze(std(Wnc(:,1,co))) squeeze(std(Bnc(:,9,co))) squeeze(std(Wnc(:,9,co)))]/2,'.')
        xlim([0.5 5])

        set(gca, 'XTick', []);
        ax = gca;
        tickPositions = [1, 2, 3.5, 4.5];
        colors = {'blue', 'red', 'blue', 'red'};
        ticks = {'Between', 'Within', 'Between', 'Within'};

        % Loop to place each text at the specified position and color
        for i = 1:length(tickPositions)
            text(tickPositions(i), ax.YLim(1)-0.01/3, ticks{i}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', colors{i}, 'Rotation', 45); % Added 'Rotation' property for vertical text
        end
        box off
        if co==1
            ylabel('noise corr (choice pools)','fontsize',11)
            title('correct trials')
        else
            title('incorrect trials')
        end
                
        ylim([0 0.04])
    end
    set(hk1,'PaperPositionMode','Auto','Units','normalized','Position',[0 0 1/2 2]/3);


%%%% interaction for choice pools within and between 

    for label=[1 9]
        for co=1:2

            diff_Wnc(:,label,co)=bootstrp(1000,@nanmean,INT{label,co}(SC{label,1}>0))-mean(bootstrp(1000,@nanmean,INT{label,co}));
            diff_Bnc(:,label,co)=bootstrp(1000,@nanmean,INT{label,co}(SC{label,1}<0))-mean(bootstrp(1000,@nanmean,INT{label,co}));

            Wnc(:,label,co)=bootstrp(1000,@nanmean,INT{label,co}(SC{label,1}>0));
            Bnc(:,label,co)=bootstrp(1000,@nanmean,INT{label,co}(SC{label,1}<0));

        end
    end

    hk2=figure(2)
    for co=1:2
        subplot(2,2,co)
        bar(1,squeeze(mean(diff_Bnc(:,1,co))),'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(2,squeeze(mean(diff_Wnc(:,1,co))),'facecolor',[0 0 0],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        bar(3.5,squeeze(mean(diff_Bnc(:,9,co))),'facecolor',[1 0 1],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(4.5,squeeze(mean(diff_Wnc(:,9,co))),'facecolor',[1 0 1],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        errorbar([1 2 3.5 4.5],[squeeze(mean(diff_Bnc(:,1,co))) squeeze(mean(diff_Wnc(:,1,co))) squeeze(mean(diff_Bnc(:,9,co))) squeeze(mean(diff_Wnc(:,9,co)))],[squeeze(std(diff_Bnc(:,1,co)))  squeeze(std(diff_Wnc(:,1,co))) squeeze(std(diff_Bnc(:,9,co))) squeeze(std(diff_Wnc(:,9,co)))]/2,'.')
        xlim([0.5 5])
        ylim([-8 8]*1e-5)
        set(gca, 'XTick', []);
        ax = gca;
        tickPositions = [1 2 3.5 4.5];
        colors = {'blue', 'red', 'blue', 'red'};
        ticks={'Between', 'Within', 'Between', 'Within'};
        for i = 1:length(tickPositions)
            text(tickPositions(i), ax.YLim(1), ticks{i}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', colors{i}, 'Rotation', 45); % Added 'Rotation' property for vertical text
        end
        box off
        if co==1
            ylabel('relative interaction inf (choice pools)','fontsize',11)
            title('correct trials')
        else
            title('incorrect trials')
        end
        
    end

        
    hk2=figure(2)
    for co=1:2
        subplot(2,2,co+2)
        bar(1,squeeze(mean(Bnc(:,1,co))),'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(2,squeeze(mean(Wnc(:,1,co))),'facecolor',[0 0 0],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        bar(3.5,squeeze(mean(Bnc(:,9,co))),'facecolor',[1 0 1],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(4.5,squeeze(mean(Wnc(:,9,co))),'facecolor',[1 0 1],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        errorbar([1 2 3.5 4.5],[squeeze(mean(Bnc(:,1,co))) squeeze(mean(Wnc(:,1,co))) squeeze(mean(Bnc(:,9,co))) squeeze(mean(Wnc(:,9,co)))],[squeeze(std(diff_Bnc(:,1,co)))  squeeze(std(Wnc(:,1,co))) squeeze(std(Bnc(:,9,co))) squeeze(std(Wnc(:,9,co)))]/2,'.')
        xlim([0.5 5])

        set(gca, 'XTick', []);
        ax = gca;
        tickPositions = [1 2 3.5 4.5];
        colors = {'blue', 'red', 'blue', 'red'};
        ticks={'Between', 'Within', 'Between', 'Within'};
        for i = 1:length(tickPositions)
            text(tickPositions(i), ax.YLim(1), ticks{i}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', colors{i}, 'Rotation', 45); % Added 'Rotation' property for vertical text
        end
        box off
        if co==1
            ylabel('interaction inf (choice pools)','fontsize',11)
            title('correct trials')
        else
            title('incorrect trials')
        end
                
        ylim([-1e-5 3e-4])
    end
    set(hk2,'PaperPositionMode','Auto','Units','normalized','Position',[0 0 1/2 2]/3);


%%%% signal correlation for choice pools within and between 

    for label=[1 9]
        for co=1:2

            diff_Wnc(:,label,co)=bootstrp(1000,@nanmean,SC{label,co}(SC{label,1}>0))-mean(bootstrp(1000,@nanmean,SC{label,co}));
            diff_Bnc(:,label,co)=bootstrp(1000,@nanmean,SC{label,co}(SC{label,1}<0))-mean(bootstrp(1000,@nanmean,SC{label,co}));

            Wnc(:,label,co)=bootstrp(1000,@nanmean,SC{label,co}(SC{label,1}>0));
            Bnc(:,label,co)=bootstrp(1000,@nanmean,SC{label,co}(SC{label,1}<0));

        end
    end

    hk3=figure(3)
    for co=1:2
        subplot(2,2,co)
        bar(1,squeeze(mean(diff_Bnc(:,1,co))),'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(2,squeeze(mean(diff_Wnc(:,1,co))),'facecolor',[0 0 0],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        bar(3.5,squeeze(mean(diff_Bnc(:,9,co))),'facecolor',[1 0 1],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(4.5,squeeze(mean(diff_Wnc(:,9,co))),'facecolor',[1 0 1],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        errorbar([1 2 3.5 4.5],[squeeze(mean(diff_Bnc(:,1,co))) squeeze(mean(diff_Wnc(:,1,co))) squeeze(mean(diff_Bnc(:,9,co))) squeeze(mean(diff_Wnc(:,9,co)))],[squeeze(std(diff_Bnc(:,1,co)))  squeeze(std(diff_Wnc(:,1,co))) squeeze(std(diff_Bnc(:,9,co))) squeeze(std(diff_Wnc(:,9,co)))]/2,'.')
        xlim([0.5 5])
        ylim([-5 5]*1e-2)
        set(gca, 'XTick', []);
        ax = gca;
        tickPositions = [1 2 3.5 4.5];
        colors = {'blue', 'red', 'blue', 'red'};
        ticks={'Between', 'Within', 'Between', 'Within'};
        for i = 1:length(tickPositions)
            text(tickPositions(i), ax.YLim(1), ticks{i}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', colors{i}, 'Rotation', 45); % Added 'Rotation' property for vertical text
        end
        box off
        if co==1
            ylabel('relative signal corr (choice pools)','fontsize',11)
            title('correct trials')
        else
            title('incorrect trials')
        end
        
    end

        
    hk3=figure(3)
    for co=1:2
        subplot(2,2,co+2)
        bar(1,squeeze(mean(Bnc(:,1,co))),'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(2,squeeze(mean(Wnc(:,1,co))),'facecolor',[0 0 0],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        bar(3.5,squeeze(mean(Bnc(:,9,co))),'facecolor',[1 0 1],'facealpha',0.2,'edgecolor','b','linewidth',2)
        hold on
        bar(4.5,squeeze(mean(Wnc(:,9,co))),'facecolor',[1 0 1],'facealpha',1,'edgecolor','r','linewidth',2)
        hold on
        errorbar([1 2 3.5 4.5],[squeeze(mean(Bnc(:,1,co))) squeeze(mean(Wnc(:,1,co))) squeeze(mean(Bnc(:,9,co))) squeeze(mean(Wnc(:,9,co)))],[squeeze(std(diff_Bnc(:,1,co)))  squeeze(std(Wnc(:,1,co))) squeeze(std(Bnc(:,9,co))) squeeze(std(Wnc(:,9,co)))]/2,'.')
        xlim([0.5 5])

        set(gca, 'XTick', []);
        ax = gca;
        tickPositions = [1 2 3.5 4.5];
        colors = {'blue', 'red', 'blue', 'red'};
        ticks={'Between', 'Within', 'Between', 'Within'};
        for i = 1:length(tickPositions)
            text(tickPositions(i), ax.YLim(1), ticks{i}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', colors{i}, 'Rotation', 45); % Added 'Rotation' property for vertical text
        end
        box off
        if co==1
            ylabel('signal corr (choice pools)','fontsize',11)
            title('correct trials')
        else
            title('incorrect trials')
        end
                
        ylim([-0.05 0.05])
    end
    set(hk3,'PaperPositionMode','Auto','Units','normalized','Position',[0 0 1/2 2]/3);




