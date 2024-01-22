function plot_copula(cop)

[trn strn]=size(cop);

n=0;
figure
for tr=1:trn
    for str=1:strn
        n=n+1;
        if ~isempty(cop{tr,str}) & isfield(cop{tr,str},'Grid')
                
        m=round(sqrt(size(cop{tr,str}.Grid.u,1)));
        subplot(trn,strn,n)
        surf(cop{tr,str}.Grid.u(1:m,1),cop{tr,str}.Grid.u(1:m:end,2),cop{tr,str}.f_grid,'FaceColor','interp','edgecolor','none')
        end
    end
end


n=0;
figure
for tr=1:trn
    for str=1:strn
        n=n+1;
        if ~isempty(cop{tr,str}) & isfield(cop{tr,str},'Grid')
                
        m=round(sqrt(size(cop{tr,str}.Grid.u,1)));
        subplot(trn,strn,n)
        imagesc(cop{tr,str}.Grid.u(1:m,1),cop{tr,str}.Grid.u(1:m:end,2),cop{tr,str}.f_grid)
        colormap('hot');

        end
    end
end



n=0;
figure
title('u space')
for tr=1:trn
    for str=1:strn
        n=n+1;
        if ~isempty(cop{tr,str}) & isfield(cop{tr,str},'Grid')
        subplot(trn,strn,n)
        plot(cop{tr,str}.data.u(:,1),cop{tr,str}.data.u(:,2),'.')
        end
    end
end


n=0;
figure
title('S space')
for tr=1:trn
    for str=1:strn
        n=n+1;
        if ~isempty(cop{tr,str}) & isfield(cop{tr,str},'Grid')
        subplot(trn,strn,n)
        plot(cop{tr,str}.data.X(:,1),cop{tr,str}.data.X(:,2),'.')
        end
    end
end


    


