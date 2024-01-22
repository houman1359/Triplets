function vine_sh=shuffle_vine(vine)


rng('shuffle')

N=numel(vine);
M=numel(vine{1}.margins{1}.ker);

vine_sh=vine;
for j=1:numel(vine{1}.margins)
    D{j}=[];
end

for n=1:N
    for j=1:numel(vine{n}.margins)
        D{j}=cat(1,D{j},vine{n}.margins{j}.ker);
        
    end
end


a(1)=0;
for n=1:N
a(n+1)=a(n)+numel(vine{n}.margins{1}.ker);
end

SH=randsample(1:numel(D{j}),numel(D{j}));


for n=1:N
    for j=1:numel(vine{n}.margins)
        vine_sh{n}.margins{j}.ker=D{j}(SH(a(n)+1:a(n+1)));
        vine{1}.theta{j,j}=D{j}(SH(a(n)+1:a(n+1)));
    end
end

% 
% figure
% subplot(1,2,1)
% plot(vine{1}.margins{2}.ker,vine{1}.margins{1}.ker,'.')
% hold on
% plot(vine{2}.margins{2}.ker,vine{2}.margins{1}.ker,'.r')
% subplot(1,2,2)
% plot(vine_sh{1}.margins{2}.ker,vine_sh{1}.margins{1}.ker,'.')
% hold on
% plot(vine_sh{2}.margins{2}.ker,vine_sh{2}.margins{1}.ker,'.r')
% 
