function plot_empiricalCDF(vine)

[trn strn]=size(vine.margins);

if numel(vine.margins{1}.ker)>4000
sa=randsample(1:numel(vine.margins{1}.ker),4000);
else
sa=1:numel(vine.margins{1}.ker);    
end


for t=1:trn
   x=vine.margins{t}.ker(sa);
   x=x+1e-10*rand(size(x)); 
    
   C{t}=kernelcdf(x,1,1,x);
end





n=0;
figure
title('u space')
for tr=1:trn
    for str=1:trn
        n=n+1;
        if tr<str
        subplot(trn,trn,n)
        plot(C{tr},C{str},'.')
        xlim([0 1]);ylim([0 1])       
        end
    end
end

