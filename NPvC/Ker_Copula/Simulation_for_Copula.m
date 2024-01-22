


switch simulation_num
    case 1
        clear x
        x(:,1)=rand(2000,1);%[0:0.01:1 0:0.01:1 0:0.01:1 0:0.01:1];
        x(1:1000,2)=x(1:1000,1)+rand(1000,1)/20;%+rand(100,1)/4;
        x(1001:2000,2)=1-x(1001:2000,1)+rand(1000,1)/20;%+rand(100,1)/4;
        
        % DDD=500*abs(randn(2000,1));
        % X(:,1) = DDD+randn(2000,1)*10;
        % X(1:1000,2) = 1*X(1:1000,1)+randn(1000,1)*10;
        % X(1001:2000,2) = 1*(X(1001:2000,1)+randn(1000,1)*10);
        % X=X+1;
        % x=X/500;
        
        f=x;
        x(:,1)=real(x(:,1)*1);
        
        X1 = x2fx(f(1:2000,[1 2])/1,'interaction');
        X1(:,1)=1;
        X2 = x2fx(x(1:2000,[1 2])/1,'interaction');
        X2(:,1)=1;
        
        mu(1,find(x(:,1)>median(x(:,1)))) = (X2(find(x(:,1)>median(x(:,1))),[1 2 3 4])*[1;1;0;0.0] + 1)+(X2(find(x(:,1)>median(x(:,1))),[1 2 3 4])*[1;0;2;0.0] + 1);
        mu(1,find(x(:,1)<=median(x(:,1)))) = (X2(find(x(:,1)<=median(x(:,1))),[1 2 3 4])*[1;-2;0;0.0] + 1)+ (X2(find(x(:,1)<=median(x(:,1))),[1 2 3 4])*[1;0;-2;0.0] + 1);
        y = poissrnd(mu);
        Xglm=X1(:,[1 2 3 4]);
        X=x;
        
        
        X=cat(2,y',X);1
        
    case 2
        
        
        for casee=2
            
            %% Construct 4D mixed copula vine
            disp('Constructing vine copula ...');
            d = 4; % Dimension
            vine.type = 'c-vine'; % Canonical vine type
            % Set margins
            
            switch casee
                case 1
                    vine.margins = cell(d,1);
                    % Standard normal margin
                    vine.margins{1}.dist = 'norm';
                    vine.margins{1}.theta = [1;2];
                    vine.margins{1}.iscont = true; % Continuous margin
                    % Gamma margin
                    vine.margins{2}.dist = 'gam';
                    vine.margins{2}.theta = [4;2];
                    vine.margins{2}.iscont = true; % Continuous margin
                    % Poisson margin
                    vine.margins{3}.dist = 'norm';
                    vine.margins{3}.theta = [1;3];
                    vine.margins{3}.iscont = true; % Discrete margin
                    % Binomial margin
                    vine.margins{4}.dist = 'norm';
                    vine.margins{4}.theta = [2;0.5];
                    vine.margins{4}.iscont = true; % Discrete margin
                    % Set copula families
                    vine.families = cell(d);
                    vine.theta = cell(d);
                    % Gaussian copula family
                    vine.families{1,2} = 'gaussian';
                    vine.theta{1,2} = 0.5;
                    % Student copula family
                    vine.families{1,3} = 'student';
                    vine.theta{1,3} = [0.7;1];
                    % Clayton copula family
                    vine.families{1,4} = 'clayton';
                    vine.theta{1,4} = 5;
                    % Clayton copula family rotated 90° clockwise
                    vine.families{2,3} = 'claytonrot090';
                    vine.theta{2,3} = 8;
                    % Clayton copula family survival transformed
                    vine.families{2,4} = 'student';
                    vine.theta{2,4} = [-1 3];
                    % Independence
                    vine.families{3,4} = 'claytonrot270';
                    vine.theta{3,4} = 2;
                    fprintf('\n');
                case 2
                    vine.margins = cell(d,1);
                    % Standard normal margin
                    vine.margins{1}.dist = 'norm';
                    vine.margins{1}.theta = [2;0.5];
                    vine.margins{1}.iscont = true; % Continuous margin
                    % Gamma margin
                    vine.margins{2}.dist = 'gam';
                    vine.margins{2}.theta = [2;1];
                    vine.margins{2}.iscont = true; % Continuous margin
                    % Poisson margin
                    vine.margins{3}.dist = 'norm';
                    vine.margins{3}.theta = [4;2];
                    vine.margins{3}.iscont = true; % Discrete margin
                    % Binomial margin
                    vine.margins{4}.dist = 'gam';
                    vine.margins{4}.theta = [1;2];
                    vine.margins{4}.iscont = true; % Discrete margin
                    % Set copula families
                    vine.families = cell(d);
                    vine.theta = cell(d);
                    % Gaussian copula family
                    vine.families{2,3} = 'gaussian';
                    vine.theta{2,3} = 0.1;
                    % Student copula family
                    vine.families{1,4} = 'student';
                    vine.theta{1,4} = [0.9;0.2];
                    % Clayton copula family
                    vine.families{2,4} = 'clayton';
                    vine.theta{2,4} = 5;
                    % Clayton copula family rotated 90° clockwise
                    vine.families{1,2} = 'claytonrot090';
                    vine.theta{1,2} = 3;
                    % Clayton copula family survival transformed
                    vine.families{3,4} = 'student';
                    vine.theta{3,4} = [-0.9;3];
                    % Independence
                    vine.families{1,3} = 'student';
                    vine.theta{1,3} = [0.3;1];
                    fprintf('\n');
            end
            
            %% Test probability density function
            disp('Calculating probability density function on a grid...');
            % Calculate probability density function on lattice
            x1gv = linspace(-3,3,50);
            x2gv = linspace(0.5,25,50);
            x3gv = linspace(-3,3,50);
            x4gv = linspace(-3,3,50);
            [x1,x2,x3,x4] = ndgrid(x1gv,x2gv,x3gv,x4gv);
            p = mixedvinepdf(vine,[x1(:),x2(:),x3(:),x4(:)]);
            p = reshape(p,size(x1));
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
            cases = 3000;
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
            
            
            
            if casee==2
                X=x;
            else
                X=cat(1,X,x);
            end
            if casee==2
                x=X;
                figure('Name','Scatter','Position',[0,0,1600,1000]);
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
            end
            
        end
        
        
        
        
        
        
        
    case 3
        %% Construct 4D mixed copula vine
        disp('Constructing vine copula ...');
        d = 2; % Dimension
        vine.type = 'c-vine'; % Canonical vine type
        % Set margins
        vine.margins = cell(d,1);
        
        vine.margins{1}.dist = 'norm';
        vine.margins{1}.theta = [1;0.5];
        vine.margins{1}.iscont = true; % Continuous margin
        
        vine.margins{2}.dist = 'poiss';
        vine.margins{2}.theta = 5%[-1;1];
        vine.margins{2}.iscont = true; % Continuous margin
        
        vine.margins{3}.dist = 'norm';
        vine.margins{3}.theta = [-1;0.5];
        vine.margins{3}.iscont = true; % Discrete margin
        % Set copula families
        vine.families = cell(d);
        vine.theta = cell(d);
        % Gaussian copula family
        vine.families{1,2} = 'gaussian';
        vine.theta{1,2} = 0.5;
        % Student copula family
        vine.families{1,3} = 'gaussian';
        vine.theta{1,3} = -0.9;
        % % Clayton copula family
        vine.families{2,3} = 'gaussian';
        vine.theta{2,3} = 0.5;
        
        %% Test probability density function
        disp('Calculating probability density function on a grid...');
        % Calculate probability density function on lattice
        x1gv = linspace(-3,3,100);
        x2gv = 0:15;%linspace(-7,7,100);
        x3gv = linspace(-3,3,100);
        % x4gv = 0:15;
        [x1,x2,x3] = ndgrid(x1gv,x2gv,x3gv);
        p = mixedvinepdf(vine,[x1(:),x2(:),x3(:)]);
        p = reshape(p,size(x1));
        % Plot 2D margins
        figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);
        subplot(2,3,1);
        margin12 = reshape(sum(p,3),[length(x1gv),length(x2gv)]);
        imagesc(x1gv,x2gv,margin12');
        colormap('hot');
        set(gca,'YDir','normal');
        xlabel('Margin 1');
        ylabel('Margin 2');
        subplot(2,3,2);
        margin13 = reshape(sum(p,2),[length(x1gv),length(x3gv)]);
        imagesc(x1gv,x3gv,margin13');
        colormap('hot');
        set(gca,'YDir','normal');
        xlabel('Margin 1');
        ylabel('Margin 3');
        subplot(2,3,3);
        margin23 = reshape(sum(p,1),[length(x2gv),length(x3gv)]);
        imagesc(x2gv,x3gv,margin23');
        colormap('hot');
        set(gca,'YDir','normal');
        xlabel('Margin 2');
        ylabel('Margin 3');
        %% Test sampling
        disp('Sampling from mixed copula vine...');
        % Draw samples
        cases = 10000;
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
        scatter(x(:,2),x(:,3),20,[0 0 0],'filled');
        xlabel('Margin 2');
        ylabel('Margin 3');
        
        % figure;imagesc(margin23/max(margin23(:))-margin13/max(margin13(:)));
        
        
        
        figure
        imagesc(squeeze(p(:,3,:)/max(max(p(:,3,:)))-p(:,13,:)/max(max(p(:,13,:)))))
        
        X=x;
        
        
        
    case 4
        
        pref1=12;
        pref2=14;
        cases=2000;
        
        clear X
        
        for neu=1:1
            c1=[];c2=[];b=[];
            for rr=1:0.5:20
                rrd=20*normpdf(rr,pref1,3)/normpdf(pref1,pref1,3);
                y=20-rr;
                y=20*normpdf(rr,pref2,3)/normpdf(pref2,pref2,3);
                MU(1)=rrd;
                SIG(:,1)=[2 1.3];
                MU(2)=y;
                SIG(:,2)=[1.3 2];
                r = abs(mvnrnd(MU,SIG,cases));
                c1=cat(1,c1,r(:,1));
                c2=cat(1,c2,r(:,2));
                b=cat(1,b,rr*ones(cases,1));
            end
            
            figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);
            subplot(2,3,1)
            scatter(c1,c2,20,[0 0 0],'filled')
            subplot(2,3,2)
            scatter(b,c1,20,[0 0 0],'filled')
            subplot(2,3,3)
            scatter(b,c2,20,[0 0 0],'filled')
            
            
            %         figure
            %         scatter(c1(randsample(b12,numel(b12))),c2(randsample(b12,numel(b12))),20,[1 0 0],'filled')
            %         hold on
            %         scatter(c1(randsample(b11,numel(b11))),c2(randsample(b11,numel(b11))),20,[0 0 1],'filled')
            
            %%%12 11
            for i=16%12%11%1:3:20
                b12=find(b==i);
                b11=find(b==i+1);
                X(:,(neu-1)*2+1)=[c1(b11);c1(b12)];
                X(:,(neu)*2)=[c2(b11);c2(b12)];
                
                figure
                scatter(c1(b==i),c2(b==i),20,[1 0 0],'filled')
                hold on
                scatter(c1(b==i+1),c2(b==i+1),20,[0 0 1],'filled')
                
                %                         X(:,1)=[c1(randsample(b11,numel(b11)));c1(randsample(b12,numel(b12)))];
                %                         X(:,2)=[c2(randsample(b11,numel(b11)));c2(randsample(b12,numel(b12)))];
                %                         X(:,3)=[ones(size(b11));2*ones(size(b12))];
                
                figure
                scatter(c1(randsample(b12,numel(b12))),c2(randsample(b12,numel(b12))),20,[1 0 0],'filled')
                hold on
                scatter(c1(randsample(b11,numel(b11))),c2(randsample(b11,numel(b11))),20,[0 0 1],'filled')
                
            end
            X(:,neu*2+1)=[ones(size(b11));2*ones(size(b12))];
            
            
        end
        
        
        
        
        
    case 5
        test_pop
        
        
    case 6
        
        %% Construct 4D mixed copula vine
        disp('Constructing vine copula ...');
        d = 2; % Dimension
        vine.type = 'c-vine'; % Canonical vine type
        % Set margins
        vine.margins = cell(d,1);
        
        for i=1:2
            vine.margins{1}.dist = 'norm';
            vine.margins{1}.theta = [1+i;1];
            vine.margins{1}.iscont = true; % Continuous margin
            
            vine.margins{2}.dist = 'norm';
            vine.margins{2}.theta = [-1;1];
            vine.margins{2}.iscont = true; % Continuous margin
            
            vine.margins{3}.dist = 'norm';
            vine.margins{3}.theta = [-1;0.5];
            vine.margins{3}.iscont = true; % Discrete margin
            % Set copula families
            vine.families = cell(d);
            vine.theta = cell(d);
            % Gaussian copula family
            vine.families{1,2} = 'gaussian';
            vine.theta{1,2} = 0.5;
            % Student copula family
            vine.families{1,3} = 'gaussian';
            vine.theta{1,3} = -0.9;
            % % Clayton copula family
            vine.families{2,3} = 'gaussian';
            vine.theta{2,3} = 0.5;
            
            %% Test probability density function
            disp('Calculating probability density function on a grid...');
            % Calculate probability density function on lattice
            x1gv = linspace(-3,3,100);
            x2gv = 0:15;%linspace(-7,7,100);
            x3gv = linspace(-3,3,100);
            % x4gv = 0:15;
            [x1,x2,x3] = ndgrid(x1gv,x2gv,x3gv);
            p = mixedvinepdf(vine,[x1(:),x2(:),x3(:)]);
            p = reshape(p,size(x1));
            % Plot 2D margins
            figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);
            subplot(2,3,1);
            margin12 = reshape(sum(p,3),[length(x1gv),length(x2gv)]);
            imagesc(x1gv,x2gv,margin12');
            colormap('hot');
            set(gca,'YDir','normal');
            xlabel('Margin 1');
            ylabel('Margin 2');
            subplot(2,3,2);
            margin13 = reshape(sum(p,2),[length(x1gv),length(x3gv)]);
            imagesc(x1gv,x3gv,margin13');
            colormap('hot');
            set(gca,'YDir','normal');
            xlabel('Margin 1');
            ylabel('Margin 3');
            subplot(2,3,3);
            margin23 = reshape(sum(p,1),[length(x2gv),length(x3gv)]);
            imagesc(x2gv,x3gv,margin23');
            colormap('hot');
            set(gca,'YDir','normal');
            xlabel('Margin 2');
            ylabel('Margin 3');
            %% Test sampling
            disp('Sampling from mixed copula vine...');
            % Draw samples
            cases = 2000;
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
            scatter(x(:,2),x(:,3),20,[0 0 0],'filled');
            xlabel('Margin 2');
            ylabel('Margin 3');
            
            % figure;imagesc(margin23/max(margin23(:))-margin13/max(margin13(:)));
            
            
            
            figure
            imagesc(squeeze(p(:,3,:)/max(max(p(:,3,:)))-p(:,13,:)/max(max(p(:,13,:)))))
            
            X{i}=x;
        end
        close all
        figure;plot3(X{1}(:,1),X{1}(:,2),X{1}(:,3),'.')
        hold on
        plot3(X{2}(:,1),X{2}(:,2),X{2}(:,3),'.r')
        
        
        
        
        
    case 7
        
        for i=1:2
            %% Construct 4D mixed copula vine
            disp('Constructing vine copula ...');
            d = 2; % Dimension
            vine.type = 'c-vine'; % Canonical vine type
            % Set margins
            vine.margins = cell(d,1);
            
            vine.margins{1}.dist = 'norm';
            vine.margins{1}.theta = [1+i*1;0.2];
            vine.margins{1}.iscont = true; % Continuous margin
            
            vine.margins{2}.dist = 'norm';
            vine.margins{2}.theta = [1+i*1;0.2];
            vine.margins{2}.iscont = true; % Continuous margin
            
            %                 vine.margins{3}.dist = 'norm';
            %                 vine.margins{3}.theta = [-1;0.5];
            %                 vine.margins{3}.iscont = true; % Discrete margin
            % Set copula families
            vine.families = cell(d);
            vine.theta = cell(d);
            % Gaussian copula family
            vine.families{1,2} = 'ind';
            vine.theta{1,2} = 0.99;
            % Student copula family
            %                 vine.families{1,3} = 'gaussian';
            %                 vine.theta{1,3} = -0.9;
            %                 % % Clayton copula family
            %                 vine.families{2,3} = 'gaussian';
            %                 vine.theta{2,3} = 0.5;
            
            %% Test probability density function
            disp('Calculating probability density function on a grid...');
            % Calculate probability density function on lattice
            x1gv = linspace(-3,3,100);
            x2gv = 0:15;%linspace(-7,7,100);
            x3gv = linspace(-3,3,100);
            % x4gv = 0:15;
            [x1,x2] = ndgrid(x1gv,x2gv);
            p = mixedvinepdf(vine,[x1(:),x2(:)]);
            p = reshape(p,size(x1));
            % Plot 2D margins
            figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);
            subplot(2,3,1);
            margin12 = reshape(sum(p,3),[length(x1gv),length(x2gv)]);
            imagesc(x1gv,x2gv,margin12');
            colormap('hot');
            set(gca,'YDir','normal');
            xlabel('Margin 1');
            ylabel('Margin 2');
            subplot(2,3,2);
            %                 margin13 = reshape(sum(p,2),[length(x1gv),length(x3gv)]);
            %                 imagesc(x1gv,x3gv,margin13');
            colormap('hot');
            set(gca,'YDir','normal');
            xlabel('Margin 1');
            ylabel('Margin 3');
            subplot(2,3,3);
            %                 margin23 = reshape(sum(p,1),[length(x2gv),length(x3gv)]);
            %                 imagesc(x2gv,x3gv,margin23');
            colormap('hot');
            set(gca,'YDir','normal');
            xlabel('Margin 2');
            ylabel('Margin 3');
            %% Test sampling
            disp('Sampling from mixed copula vine...');
            % Draw samples
            cases = 10000;
            x = mixedvinernd(vine,cases);
            % Plot samples in 2D
            figure('Name','Mixed vine copula samples in 2D','Position',[0,0,1600,1000]);
            subplot(2,3,1);
            scatter(x(:,1),x(:,2),20,[0 0 0],'filled');
            xlabel('Margin 1');
            ylabel('Margin 2');
            subplot(2,3,2);
            %                 scatter(x(:,1),x(:,3),20,[0 0 0],'filled');
            xlabel('Margin 1');
            ylabel('Margin 3');
            subplot(2,3,3);
            %                 scatter(x(:,2),x(:,3),20,[0 0 0],'filled');
            xlabel('Margin 2');
            ylabel('Margin 3');
            
            % figure;imagesc(margin23/max(margin23(:))-margin13/max(margin13(:)));
            
            
            
            figure
            imagesc(squeeze(p(:,3,:)/max(max(p(:,3,:)))-p(:,13,:)/max(max(p(:,13,:)))))
            
            X{i}=x;
        end
        
        
    case 8
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
        for i=1:2
            x(:,1)=unifrnd(1,6,10000,1);
            if i==1
                %            x(:,2)=unifrnd(1+i,4+2*i,10000,1);
                x(:,2)=unifrnd(2,7,10000,1);
                %            x(:,2)=unifrnd(2,5,10000,1);
            else
                %            x(:,2)=unifrnd(1+i,3+2*i,10000,1);
                x(:,2)=unifrnd(2,7,10000,1);
                %            x(:,2)=unifrnd(1+i,6,10000,1);
            end
            x(:,1)=x(:,1)-(-2+x(:,2)).^2+5*x(:,2);
            %            x(:,1)=x(:,1)+(-2+x(:,2)).^2;
            X{i}=x;
        end
        figure;plot(X{2}(:,2),X{2}(:,1),'.r');hold on;plot(X{1}(:,2),X{1}(:,1),'.')
        xlabel('B')
        ylabel('n')
        
        
    case 9
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
        for i=1:2
            
            clear x xx
            
            if i==1
                x(:,2)=normrnd(10,8,1000,1);
                %            x(1:2000,1)=x(1:2000,2)+normrnd(4,2,2000,1)+4;
                %            x(2001:10000,1)=x(2001:10000,2)+normrnd(-4,2,8000,1)+4;
                x(:,1)=x(:,2)+normrnd(0,5,1000,1)+4;
            else
                x(:,2)=normrnd(10,8,1000,1);
                x(1:500,1)=x(1:500,2)+normrnd(7,2,500,1)+4;
                x(501:1000,1)=x(501:1000,2)+normrnd(-7,2,500,1)+4;
            end
            X{i}=x;
        end
        
        figure;plot(X{2}(:,2),X{2}(:,1),'.r');hold on;plot(X{1}(:,2),X{1}(:,1),'.')
        xlabel('B')
        ylabel('n')
        axis equal
        
        
        clear Z
        nnn=1:1000;%randsample(1:10000,1000);
        Z{1}(:,1:2)=X{1}(nnn,:);
        Z{1}(:,3)=[ones(size(X{1}(nnn,:),1),1)];
        Z{2}(:,1:2)=X{2}(nnn,:);
        Z{2}(:,3)=[2*ones(size(X{2}(nnn,:),1),1)];
        figure
        plot3(Z{1}(:,1),Z{1}(:,2),Z{1}(:,3),'.b')
        hold on
        plot3(Z{2}(:,1),Z{2}(:,2),Z{2}(:,3),'.r')
        %%%% marginal f(B|C=1) and f(B|C=2))
        
        %         nnn=randsample(1:10000,5000);
        %         X{1}=X{1}(nnn,:);
        %         X{2}=X{2}(nnn,:);
        %
        
    case 10
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
        for i=1:2
            x(:,1)=unifrnd(1,6,1000,1);
            if i==1
                x(:,2)=unifrnd(2,8,1000,1);
            else
                x(:,2)=unifrnd(8,11,1000,1);
            end
            x(:,1)=x(:,1)-(-2+x(:,2)).^2+10*x(:,2);
            x(:,1)=x(:,1)/5;
            
            xx(:,1)=unifrnd(1,6,10000,1);
            if i==1
                xx(:,2)=unifrnd(8,11,10000,1);
            else
                xx(:,2)=unifrnd(2,8,10000,1);
            end
            xx(:,1)=xx(:,1)-(-2+xx(:,2)).^2+10*xx(:,2);
            xx(:,1)=xx(:,1)/5;
            
            %            x(:,1)=x(:,1)+(-2+x(:,2)).^2;
            X{i}=cat(1,x,xx);
        end
        figure;plot(X{2}(:,2),X{2}(:,1),'.r');hold on;plot(X{1}(:,2),X{1}(:,1),'.')
        xlabel('B')
        ylabel('n')
        %%%% marginal f(B|C=1) and f(B|C=2))
        
    case 11
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
        rng(1234)
        for i=1:2
            clear x xx
            x(:,1)=unifrnd(1,6,2000,1);
            
            noise=normrnd(0,10,2000,1);
            if i==1
                x(:,2)=unifrnd(0,12,2000,1);
            else
                x(:,2)=unifrnd(12,17,2000,1);
            end
            x(:,1)=normrnd(30-(-10+x(:,2)).^2,10,2000,1);
            x(:,1)=x(:,1)/5+20;
            
            xx(:,1)=unifrnd(1,6,8000,1);
            noise=normrnd(0,10,8000,1);
            if i==1
                xx(:,2)=unifrnd(12,17,8000,1);
            else
                xx(:,2)=unifrnd(0,12,8000,1);
            end
            xx(:,1)=normrnd(30-(-10+xx(:,2)).^2,10,8000,1);
            xx(:,1)=xx(:,1)/5+20;
            
            X{i}=cat(1,x,xx);
        end
        figure;plot(X{2}(:,2),X{2}(:,1),'.r');hold on;plot(X{1}(:,2),X{1}(:,1),'.')
        xlabel('B')
        ylabel('n')
        
        clear Z
        nnn=randsample(1:10000,1000);
        Z{1}(:,1:2)=X{1}(nnn,:);
        Z{1}(:,3)=[ones(size(X{1}(nnn,:),1),1)];
        Z{2}(:,1:2)=X{2}(nnn,:);
        Z{2}(:,3)=[2*ones(size(X{2}(nnn,:),1),1)];
        figure
        plot3(Z{1}(:,1),Z{1}(:,2),Z{1}(:,3),'.b')
        hold on
        plot3(Z{2}(:,1),Z{2}(:,2),Z{2}(:,3),'.r')
        %%%% marginal f(B|C=1) and f(B|C=2))
        
        nnn=randsample(1:10000,5000);
        X{1}=X{1}(nnn,:);
        X{2}=X{2}(nnn,:);
        
    case 12
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
        rng(1234)
        
        
        for i=1:2
            clear x
            
            if i==1
                mu = [5 10];
                SIGMA = [3 0; 0 3];
                x(:,[1 3]) = mvnrnd(mu,SIGMA,1000);
                %                                 x(:,2)=x(:,1)+(13-x(:,3)).^3/5+normrnd(0,2,2500,1);
                x(:,2)=normrnd(x(:,1)+(13-x(:,3)).^3/5,2);
                %                 x(:,2)=x(:,1)+x(:,1)+normrnd(0,2,2500,1);
            else
                mu = [5 10];
                SIGMA = [3 0; 0 3];
                x(:,[1 3]) = mvnrnd(mu,SIGMA,1000);
                %                                 x(:,2)=x(:,1)+(7-x(:,3)).^3/5+normrnd(0,2,2500,1);
                x(:,2)=normrnd(x(:,1)+(7-x(:,3)).^3/5,2);
                %                 x(:,2)=x(:,1)+x(:,1)+normrnd(0,2,2500,1);
            end
            
            X{i}=x;
        end
        
        figure;plot(X{2}(:,2),X{2}(:,1),'.r');hold on;plot(X{1}(:,2),X{1}(:,1),'.')
        xlabel('B')
        ylabel('n')
        axis equal
        %
        %         figure;plot(X{1}(:,1),X{1}(:,3),'.');hold on;plot(X{2}(:,1),X{2}(:,3),'.r')
        %         figure;plot(X{1}(:,1),X{1}(:,2),'.');hold on;plot(X{2}(:,1),X{2}(:,2),'.r')
        %         figure;plot(X{1}(:,2),X{1}(:,3),'.');hold on;plot(X{2}(:,2),X{2}(:,3),'.r')
        
        %         clear Z
        %         nnn=randsample(1:5000,1000);
        %         Z{1}(:,1:2)=X{1}(nnn,[1 2]);
        %         Z{1}(:,3)=[ones(size(X{1}(nnn,[1 2]),1),1)];
        %         Z{2}(:,1:2)=X{2}(nnn,[1 2]);
        %         Z{2}(:,3)=[2*ones(size(X{2}(nnn,[1 2]),1),1)];
        %         figure
        %         plot3(Z{1}(:,1),Z{1}(:,2),Z{1}(:,3),'.b')
        %         hold on
        %         plot3(Z{2}(:,1),Z{2}(:,2),Z{2}(:,3),'.r')
        %
        %         clear Z
        %         nnn=randsample(1:5000,1000);
        %         Z{1}(:,1:2)=X{1}(nnn,[2 3]);
        %         Z{1}(:,3)=[ones(size(X{2}(nnn,[2 3]),1),1)];
        %         Z{2}(:,1:2)=X{2}(nnn,[2 3]);
        %         Z{2}(:,3)=[2*ones(size(X{2}(nnn,[2 3]),1),1)];
        %         figure
        %         plot3(Z{1}(:,1),Z{1}(:,2),Z{1}(:,3),'.b')
        %         hold on
        %         plot3(Z{2}(:,1),Z{2}(:,2),Z{2}(:,3),'.r')
    case 13
        
        EI_network_simulation_test
        
        
        
    case 14    %%% to test the time-dependent information comp.
        
        
        rng(1234)
        
        for i=1:2
            clear x
            
            if i==1
                T=[sort(10*rand(1,1000)) sort(10*rand(1,2000))+10];
                mu_n=5*cos(2*pi*T/40-pi/2)+3;
                x(:,3)=normrnd(mu_n,2);
                x(:,1)=x(:,3)+normrnd(0*mu_n,0.5)';
                x(:,2)=T;
                
            else
                T=[sort(10*rand(1,2000)) sort(10*rand(1,1000))+10];
                mu_n=0*cos(2*pi*T/40-pi/2)+3;
                x(:,3)=normrnd(mu_n,2);
                x(:,1)=x(:,3)+normrnd(0*mu_n,0.5)';
                x(:,2)=T;
            end
            
            X{i}=x;
        end
        
        figure;plot(X{2}(:,2),X{2}(:,1),'.r');hold on;plot(X{1}(:,2),X{1}(:,1),'.')
        xlabel('B')
        ylabel('n')
        axis equal
        
        
        
        
    case 15
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
        rng(1234)
        
        
        for i=1:2
            clear x
            
            if i==1
                mu = [5 5 10];
                SIGMA = [5 1.5 1.5;1.5 3 2;1.5 2 3];
                x1 = mvnrnd(mu,SIGMA,1500);
                mu = [5 5 10];
                SIGMA = [5 -1.5 -1.5;-1.5 3 -2;-1.5 -2 3];
                x2 = mvnrnd(mu,SIGMA,500);
                x=cat(1,x1,x2);
                %                                 x(:,2)=x(:,1)+(13-x(:,3)).^3/5+normrnd(0,2,2500,1);
                %                 x(:,2)=normrnd(x(:,1)+(13-x(:,3)).^3/5,2);
                %                 x(:,2)=x(:,1)+x(:,1)+normrnd(0,2,2500,1);
            else
                mu = [5 5 10];
                SIGMA = [5 1.5 1.5;1.5 3 -2;1.5 -2 3];
                x1 = mvnrnd(mu,SIGMA,1500);
                mu = [5 5 10];
                SIGMA = [5 -1.5 -1.5;-1.5 3 2;-1.5 2 3];
                x2 = mvnrnd(mu,SIGMA,500);
                x=cat(1,x1,x2);
                
                %                                 x(:,2)=x(:,1)+(7-x(:,3)).^3/5+normrnd(0,2,2500,1);
                %                 x(:,2)=normrnd(x(:,1)+(7-x(:,3)).^3/5,2);
                %                 x(:,2)=x(:,1)+x(:,1)+normrnd(0,2,2500,1);
            end
            
            X{i}=x;
        end
        
        figure
        subplot(1,3,1)
        plot(X{2}(:,1),X{2}(:,2),'.r');hold on;plot(X{1}(:,1),X{1}(:,2),'.')
        xlabel('n')
        ylabel('B1')
        axis square
        subplot(1,3,2)
        plot(X{2}(:,1),X{2}(:,3),'.r');hold on;plot(X{1}(:,1),X{1}(:,3),'.')
        xlabel('n')
        ylabel('B2')
        axis square
        subplot(1,3,3)
        plot(X{2}(:,2),X{2}(:,3),'.r');hold on;plot(X{1}(:,2),X{1}(:,3),'.')
        xlabel('B1')
        ylabel('B2')
        axis square
        
        
    case 16
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
        rng(1234)
        
        
        for i=1:2
            clear x x1
            if i==1
                vineS.margins{1}.dist = 'norm';
                vineS.margins{1}.theta = [0;5];
                vineS.margins{1}.iscont = true; % Continuous margin
                vineS.margins{2}.dist = 'norm';
                vineS.margins{2}.theta = [0;5];
                vineS.margins{2}.iscont = true; % Continuous margin
                vineS.families{1,2} = 'student';
                vineS.theta{1,2} = [0.5;1];
                vineS.families{1,2} = 'claytonrot090';
                vineS.theta{1,2} = 50;  
%                 vineS.families{1,2} = 'gaussian';
%                 vineS.theta{1,2} = -0.8;  
                x1(:,[2 3]) = mixedvinernd(vineS,2000);
                
                
                
                x(:,[2 3]) = (x1(:,[2 3]))*[sin(pi/4) cos(pi/4);-cos(pi/4) sin(pi/4)];
                x(:,[3])=abs(x(:,3));
                x(:,[2 3]) = (x(:,[2 3]))*[sin(-pi/4) cos(-pi/4);-cos(-pi/4) sin(-pi/4)];
                
%                 x(:,[3]) = ((x1(:,[3])-0).^1);%*[sqrt(2) sqrt(2);-sqrt(2) sqrt(2)];
                
                x(:,1)=normrnd(x(:,2)+1*x(:,3),2);
            else
                vineS.margins{1}.dist = 'norm';
                vineS.margins{1}.theta = [0;5];
                vineS.margins{1}.iscont = true; % Continuous margin
                vineS.margins{2}.dist = 'norm';
                vineS.margins{2}.theta = [0;5];
                vineS.margins{2}.iscont = true; % Continuous margin
                vineS.families{1,2} = 'student';
                vineS.theta{1,2} = [-0.5;1];
                vineS.families{1,2} = 'claytonrot270';
                vineS.theta{1,2} = 10;
%                 vineS.families{1,2} = 'gaussian';
%                 vineS.theta{1,2} = -0.8;
                x1(:,[2 3]) = mixedvinernd(vineS,2000);
                x(:,[2 3]) = (x1(:,[2 3]))*[sin(-pi/4) cos(-pi/4);-cos(-pi/4) sin(-pi/4)];
                x(:,[2])=abs(x(:,2));
                x(:,[2 3]) = (x(:,[2 3]))*[sin(pi/4) cos(pi/4);-cos(pi/4) sin(pi/4)];

                
%                 x(:,[3]) = ((x1(:,[3])-0).^1);%*[sqrt(2) sqrt(2);-sqrt(2) sqrt(2)];
                x(:,1)=normrnd(x(:,2)+1*x(:,3),2)+0;
            end
            X{i}=x;
        end
        figure
        subplot(1,3,1)
        plot(X{2}(:,1),X{2}(:,2),'.r');hold on;plot(X{1}(:,1),X{1}(:,2),'.')
        xlabel('n')
        ylabel('B1')
        axis square
        subplot(1,3,2)
        plot(X{2}(:,1),X{2}(:,3),'.r');hold on;plot(X{1}(:,1),X{1}(:,3),'.')
        xlabel('n')
        ylabel('B2')
        axis square
        subplot(1,3,3)
        plot(X{2}(:,2),X{2}(:,3),'.r');hold on;plot(X{1}(:,2),X{1}(:,3),'.')
        xlabel('B1')
        ylabel('B2')
        axis square
        
            case 17
        
        %%%%%% the case with two stimulus which I(n;C|B) is zero but GLM
        %%%%%% shows nonzero. because the c(n,B) is nonlinear
        
   rng(1234)
        for i=1:2
            clear x1 x2 Z
            
            if i==1            
%             x1(:,2)=unifrnd(0,15,1000,1);
            x1(:,2)=normrnd(7.5,3,2000,1);
            x1(:,1)=normrnd(30-(-10+x1(:,2)).^2,8)/10;
            Z1(:,[2 3])=(x1-mean(x1)).^1;
%             x2(:,2)=unifrnd(10,15,1000,1);
            x2(:,2)=normrnd(12.5,3,2000,1);
            x2(:,1)=normrnd(30-(-10+x2(:,2)).^2,8)/10;
            Z2(:,[2 3])=(x2-mean(x1)).^1;
            Z=cat(1,Z1,Z2);
%             Z=Z1;%cat(1,Z1,Z2);
            else
%             x1(:,2)=unifrnd(0,15,1000,1);
            x1(:,2)=normrnd(7.5,3,2000,1);
            x1(:,1)=-1*normrnd(30-(-10+x1(:,2)).^2,8)/10;
            Z1(:,[2 3])=(x1-mean(x1)).^1;
%             x2(:,2)=unifrnd(10,15,1000,1);
            x2(:,2)=normrnd(12.5,3,2000,1);
            x2(:,1)=-1*normrnd(30-(-10+x2(:,2)).^2,8)/10;
            Z2(:,[2 3])=(x2-mean(x1)).^1;
            Z=cat(1,Z1,Z2);
            Z(:,[2 3])=(Z(:,[2 3])-[12 0]).^1;
            end
            
            
            Z(:,1)=normrnd(2*Z(:,2)-2*Z(:,3),10);            
            X{i}=Z;
        end
        
        X{1}(:,[2 3])=X{1}(:,[2 3])*[sin(pi/4) cos(pi/4);-cos(pi/4) sin(pi/4)];       
        X{2}(:,[2 3])=X{2}(:,[2 3])*[sin(pi/4) cos(pi/4);-cos(pi/4) sin(pi/4)];       
        X{1}(:,1)=normrnd(0.9*X{1}(:,2)+0.9*X{1}(:,3),2);            
        X{2}(:,1)=normrnd(0.9*X{2}(:,2)+0.9*X{2}(:,3),2);            
        
        figure
        subplot(1,3,1)
        plot(X{2}(:,1),X{2}(:,2),'.r');hold on;plot(X{1}(:,1),X{1}(:,2),'.')
        xlabel('n')
        ylabel('B1')
        axis square
        subplot(1,3,2)
        plot(X{2}(:,1),X{2}(:,3),'.r');hold on;plot(X{1}(:,1),X{1}(:,3),'.')
        xlabel('n')
        ylabel('B2')
        axis square
        subplot(1,3,3)
        plot(X{2}(:,2),X{2}(:,3),'.r');hold on;plot(X{1}(:,2),X{1}(:,3),'.')
        xlabel('B1')
        ylabel('B2')
        axis square
        
        
        
    case 18 %%%%%%%% time dependent signal

    om=2;
    tsamp1=[0:0.05:3];
    dt1=diff(tsamp1);
    tsamp1=tsamp1(1:end-1);
    tsamp2=[0:0.1:3];
    dt2=diff(tsamp2);
    tsamp2=tsamp2(1:end-1);
    
    dt2=dt1;
    tsamp2=tsamp1;
    
    clear X 
    L1=[];    
    T1=[];    
    L2=[];    
    T2=[];    
    for trial=1:100
    t1=tsamp1+rand(size(tsamp1)).*dt1;    
    L1=[L1 10*sin(2*pi/om*t1)+10];  
    T1=[T1 t1]; 
    t2=tsamp2+rand(size(tsamp2)).*dt2;    
    L2=[L2 10*sin(2*pi/om*t2)+10];  
    T2=[T2 t2]; 
    end
    
    for tt=1:numel(L1)
        vineS.margins{1}.dist = 'norm';
        vineS.margins{1}.theta = [0;3];
        vineS.margins{1}.iscont = true;
        vineS.margins{2}.dist = 'norm';
        vineS.margins{2}.theta = [L1(tt);3];
        vineS.margins{2}.iscont = true;
        vineS.families{1,2} = 'ind';
        vineS.theta{1,2} = 0.001;
        X{1}(tt,[1 3]) = mixedvinernd(vineS,1);
    end
 
    for tt=1:numel(L2)
        vineS.margins{1}.dist = 'norm';
        vineS.margins{1}.theta = [0;3];
        vineS.margins{1}.iscont = true;
        vineS.margins{2}.dist = 'norm';
        vineS.margins{2}.theta = [1.3*L2(tt);3];
        vineS.margins{2}.iscont = true;
        vineS.families{1,2} = 'ind';
        vineS.theta{1,2} = 0.001;
        X{2}(tt,[1 3]) = mixedvinernd(vineS,1);
    end
    X{1}(:,1)=X{1}(:,1)+X{1}(:,3);
    X{2}(:,1)=X{2}(:,1)+X{2}(:,3);
    

    X{1}(:,2)=T1;%+rand(size(T1))*mean(diff(unique(T1)));
    X{2}(:,2)=T2;%+rand(size(T2))*mean(diff(unique(T2)));
     figure
     subplot(1,2,1)
     plot(X{1}(:,2),X{1}(:,1),'.')
     hold on
     plot(X{2}(:,2),X{2}(:,1),'.r')
     axis tight
     subplot(1,2,2)
     plot(X{1}(:,2),X{1}(:,3),'.')
     hold on
     plot(X{2}(:,2),X{2}(:,3),'.r')
     axis tight
     
     
     
    case 19 %%%%%%%% time dependent signal NON_LINEARITY LNP
      
        clear B n X TI
        
        t=0:0.1:10;
        t=t(1:end-1);
        
        for trial=1:100
            T=t+rand(size(t))*0.1;
            X0=0*sin(2*pi/10*T)'*10;
            
            B{1}=normrnd(0,10,100,1);
            B{2}=normrnd(0,10,100,1);
            
            n{1}=B{1}+rand(size(B{1}))*0.05;
            n{1}(B{1}<0)=-n{1}(B{1}<0);
%             n{1}(B{1}<0 & T'>=5)=0;
            
            n{2}=B{2}+rand(size(B{2}))*0.05;
            n{2}(B{2}<0)=0;
            
            
            X{1}((trial-1)*100+1:trial*100,1)=n{1}+X0;
            X{1}((trial-1)*100+1:trial*100,2)=T;
            X{1}((trial-1)*100+1:trial*100,3)=B{1};
            
            X{2}((trial-1)*100+1:trial*100,1)=n{2}+X0;
            X{2}((trial-1)*100+1:trial*100,2)=T;
            X{2}((trial-1)*100+1:trial*100,3)=B{2};
            
            TI((trial-1)*100+1:trial*100)=T+(trial-1)*10;
            
        end
        
        figure
        subplot(2,2,1)
        plot(TI,X{1}(:,1))
        hold on
        plot(TI,X{2}(:,1),'r')
        axis tight
        xlabel('Time')
        ylabel('n')
        subplot(2,2,2)
        plot(TI,X{1}(:,3))
        hold on
        plot(TI,X{2}(:,3),'r')
        axis tight
        xlabel('Time')
        ylabel('B')
        
        subplot(2,2,3)
        plot(X{1}(:,3),X{1}(:,1),'.')
        axis tight
        xlabel('B')
        ylabel('n')
        title('C=1')
        subplot(2,2,4)
        plot(X{2}(:,3),X{2}(:,1),'.')
        axis tight
        xlabel('B')
        ylabel('n')
        title('C=2')

      case 20
        
          mu_1=[0 0 0 0 0];
          s_1=[1 1 10 1 10];
          co_1=[0.7 -0.2 0.9 0.7 -0.8 0.8 0.2 -0.9 0.1 0.1];
          mm=0;
          for i1=1:4
              sig_1(i1,i1)=1/s_1(i1)^2;
              for i2=i1+1:5
                  mm=mm+1;
                  sig_1(i1,i2)=co_1(mm)./(s_1(i1)*s_1(i2));
                  sig_1(i2,i1)=co_1(mm)./(s_1(i1)*s_1(i2));
              end
          end
          eig(sig_1)

         X = mvnrnd(mu_1,sig_1,5000);     

         figure
         mm=0;
         for i1=1:4
             for i2=i1+1:5
                 mm=mm+1;
                 subplot(5,5,(i1-1)*5+i2)
                 plot(X(:,i1),X(:,i2),'.')
             end
         end
     
end




return


%%% this is for paper figure



x1=normrnd(0,1,10000,1);
x2=normrnd(3,2,10000,1);

[y1 y2]=hist(x1,100);
[y3 y4]=hist(x2,100);


figure;
plot(y2,y1)
hold on
plot(y4,y3,'r')




z1=normrnd(0,3,10000,1);
z2=normrnd(4,2,10000,1);

[y1 y2]=hist(z1,100);
[y3 y4]=hist(z2,100);


figure;
plot(y2,y1)
hold on
plot(y4,y3,'r')


MU=[0 0];
SIGMA=[1 0.2;0.2 3]*2;
[xx1,xx2]=meshgrid(-15:0.1:15,-6:0.1:12);
zz(:,1)=xx1(:);
zz(:,2)=xx2(:);
y = mvnpdf(zz,MU,SIGMA);

figure
imagesc(xx1(1,:),xx2(:,1),reshape(y,size(xx1,1),size(xx1,2)))
set(gca,'YDir','normal')
set(gca,'XDir','normal')


MU=[3 4];
SIGMA=[2 1.5;1.5 2]*2;
[xx1,xx2]=meshgrid(-15:0.1:15,-6:0.1:12);
zz(:,1)=xx1(:);
zz(:,2)=xx2(:);
y = mvnpdf(zz,MU,SIGMA);

figure
imagesc(xx1(1,:),xx2(:,1),reshape(y,size(xx1,1),size(xx1,2)))
set(gca,'YDir','normal')
set(gca,'XDir','normal')









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 figure for paper/poster:


                vineS.margins{1}.dist = 'norm';
                vineS.margins{1}.theta = [0;2];
                vineS.margins{1}.iscont = true; % Continuous margin
                vineS.margins{2}.dist = 'norm';
                vineS.margins{2}.theta = [0;5];
                vineS.margins{2}.iscont = true; % Continuous margin
                vineS.families{1,2} = 'student';
                vineS.theta{1,2} = [0.5;1];
                
x1gv = linspace(-10,10,500);
x2gv = linspace(-10,10,500);
[x1,x2] = ndgrid(x1gv,x2gv);
p = mixedvinepdf(vineS,[x1(:),x2(:)]);
                
pp=reshape(p,500,500);
cases = 2000;
x = mixedvinernd(vineS,cases);


range(:,1)=min(x);
range(:,2)=max(x);

[vine_C]=prep_copula(x,{'kernel','kernel'},{'kercop' 'kercop';'kercop' 'kercop'},[0 0],'c-vine','rand',range([1 2],:))
vine_C.condition=0;
for i=1:size(range,1)
    for j=1:size(range,1)
        vine_C.METH{i,j}=[1 1];
    end
end
% [f_pC,f_C,copula_C,~] = Fit_vCopula(vine_C,x(1,:),TL,'LL1',1,0,'rand',vineest{stim},knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package
[f_pC,f_C,copula_C,~] = Fit_vCopula(vine_C,x(1,:),TL,'LL1',1,0,'rand',[],knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package
[f_pC,f_C,copula_C,~] = Fit_vCopula(vine_C,x(1,:),TL,'LL1',-1,copula_C,'rand',[],knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package


                              
Grid=GRID_Bands(0,knots,'LL1');gr=Grid.u(1:knots,1);

Y=Grid.u(1:200:end,2);
X=Grid.u(1:200,1)

Xx{1}=linspace(min(vine_C.margins{1}.ker),max(vine_C.margins{2}.ker+eps),500);
par.fit=0;par.s=copula_C{1,1}.MarginS.s;par.p=copula_C{1,1}.MarginS.p;
[pdfC1,Mar]=kernelpdf(Xx{1},2,Xx{1},par); 
Xx{2}=linspace(min(vine_C.margins{2}.ker),max(vine_C.margins{2}.ker+eps),500);
par.fit=0;par.s=copula_C{2,1}.MarginS.s;par.p=copula_C{2,1}.MarginS.p;
[pdfC2,Mar]=kernelpdf(Xx{2},2,Xx{2},par); 


figure(50)
subplot(1,4,4)
imagesc(reshape(pp,500,500))
axis square
xlabel('x_1')
ylabel('x_2')
set(gca,'xtick',[],'ytick',[],'fontsize',14)


subplot(1,4,1)
plot(Xx{1},pdfC1)
axis square
xlabel('X_1')
set(gca,'xtick',[],'ytick',[],'fontsize',14)

subplot(1,4,2)
plot(Xx{1},pdfC2)
axis square 
xlabel('X_2')
set(gca,'xtick',[],'ytick',[],'fontsize',14)

figure(60)
% subplot(1,4,3)
% figure
surf(X(25:end-25),Y(25:end-25),copula_C{1,2}.f_grid(25:end-25,25:end-25))
axis square
xlabel('u_1')
ylabel('u_2')
colormap(jet)
caxis([0 2]);
set(gca,'xtick',[],'ytick',[],'fontsize',14)
