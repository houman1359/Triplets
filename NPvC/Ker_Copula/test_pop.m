

%         return



clear POP_st2 POP_st1 neu
% close all

% figure
mu = [0 0];
Sigma = [.25 .3; .3 1];
x1 = -3:.2:3; x2 = -3:.2:3;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));
% surf(x1,x2,F);

mu1 = [0 0]+5;
Sigma = [2 -1.6; -1.6 2];
r1=mvnrnd(mu1,Sigma,5000)';
mu2 = [1.5 1.5]+5+1;
% mu2 = [0.0 0.0]+5.1;
% mu2 = [0 0]+5;
Sigma = [2 -1.6; -1.6 2];
r2=mvnrnd(mu2,Sigma,5000)';

clear perf_in perf_single post_sing_pls post_sing_neg post_pop_pls post_pop_neg
for SH=1:3

for PP=1:20
[PP SH]
for po=1:PP
    
% if mod(po,2)==1
% if  PP==1 & po==1
%     mu1 = [0 0];
%     Sigma = [2 1.6; 1.6 2];
%     r1=mvnrnd(mu1,Sigma,20000)';
%     mu2 = [1.5 1.5]*2;
%     Sigma = [2 1.6; 1.6 2];
%     r2=mvnrnd(mu2,Sigma,20000)';
% end
% if SH<=5
% neu{po}(1,:)=r1(1,randsample(1:20000,20000));
% neu{po}(2,:)=r1(2,randsample(1:20000,20000));
% else
% neu{po}(1,:)=r2(1,:);
% neu{po}(2,:)=r2(2,:);
% end
% else
% if SH<=5
% neu{po}(1,:)=r1(1,randsample(1:20000,20000));
% neu{po}(2,:)=r1(2,randsample(1:20000,20000));
% else
% neu{po}(1,:)=r2(1,:);
% neu{po}(2,:)=r2(2,:);
% end
% end

if SH<-2
neu{po}(1,:)=r1(1,randsample(1:size(r1,2),size(r1,2)));
neu{po}(2,:)=r2(1,randsample(1:size(r1,2),size(r1,2)));
else
neu{po}(1,:)=r1(1,:);
neu{po}(2,:)=r2(1,:);
end

POP_st1{SH,PP}(po,:)=neu{po}(1,:)+0.5*rand(size(neu{po}(1,:)));
POP_st2{SH,PP}(po,:)=neu{po}(2,:)+0.5*rand(size(neu{po}(2,:)));
    
end


for po=1:PP
    
L{po}=linspace(min([neu{po}(1,:) neu{po}(2,:)]),max([neu{po}(1,:) neu{po}(2,:)]),50);

[yy xx]=hist([neu{po}(1,:) neu{po}(2,:)],L{po});
[YY{po} XX{po}]=histc([neu{po}(1,:) neu{po}(2,:)],L{po});
pn{po}=yy'/sum(yy);

[yy xx]=hist([neu{po}(1,:)],L{po});
[YY1{po} XX1{po}]=histc([neu{po}(1,:)],L{po});
likpn1{po}=yy'/sum(yy);

[yy xx]=hist([neu{po}(2,:)],L{po});
[YY2{po} XX2{po}]=histc([neu{po}(2,:)],L{po});
likpn2{po}=yy'/sum(yy);

pn{po}=likpn1{po}/2+likpn2{po}/2;

[class{po},err,post_sing{po},logpn{po}] = classify(L{po}',[neu{po}(1,:) neu{po}(2,:)]',[ones(size(neu{po}(1,:))) 2*ones(size(neu{po}(2,:)))]');

% likpn1{po}=2*post_sing{po}(:,1).*pn{po};
% likpn2{po}=2*post_sing{po}(:,2).*pn{po};
end

p1=1;p2=1;pnin=1;pnin1=1;pnin2=1;dd=1;
for po=1:PP
pnin=pnin.*pn{po};
pnin1=pnin1.*likpn1{po};
pnin2=pnin2.*likpn2{po};

pnin1=pnin1/sum(pnin1);
pnin2=pnin2/sum(pnin2);
% pnin=pnin/sum(pnin);
end
likpop1=pnin1/sum(pnin1);%2*post(:,1).*pnin;
likpop2=pnin2/sum(pnin2);%2*post(:,2).*pnin;
pnin=pnin/sum(pnin);


clear perf_sin post_sing_pl post_sing_ne post_singlee
for pop=1:PP    
    n1=log2(likpn1{pop}(XX1{po}));
    n2=log2(likpn2{pop}(XX2{po}));
    hx=log2(pn{pop}(XX{po}));
% perf_sin(pop)=0.5* nanmean (n1(n1~=Inf & n1~=-Inf))+0.5* nanmean (n2(n2~=Inf & n2~=-Inf))- nanmean(hx(hx~=Inf & hx~=-Inf));

n1=likpn1{pop}.*log2(likpn1{pop}./pn{pop});
n2=likpn2{pop}.*log2(likpn2{pop}./pn{pop});
perf_sin(pop)=0.5* nansum (n1(n1~=Inf & n1~=-Inf))+0.5* nansum (n2(n2~=Inf & n2~=-Inf));

post_sing_pl(pop)=mean(post_sing{pop}(1:5,1));
post_sing_ne(pop)=mean(post_sing{pop}(6:10,1));

post_singlee(pop,:)=0.5*likpn1{pop}./pn{pop};%./likpn2{pop};%0.5*likpn1{pop}./pn{pop};
post_singlee(pop,isnan(post_singlee(pop,:)) | isinf(post_singlee(pop,:)))=0;
end
% post_sing_pls(SH,PP)=mean(post_sing_pl);
% post_sing_neg(SH,PP)=mean(post_sing_ne);
perf_single(SH,PP)=mean(perf_sin);

if PP==1
    post_single(SH,PP,:)=squeeze(post_singlee);
else
    post_single(SH,PP,:)=squeeze(mean(post_singlee));
end

    n1=log2(likpop1(XX1{po}));
    n2=log2(likpop2(XX2{po}));
    hx=log2(pnin(XX{po}));
% perf_ind(SH,PP)=0.5* nanmean ( n1(n1~=Inf & n1~=-Inf))+0.5* nanmean ( n2(n2~=Inf & n2~=-Inf))- nanmean(hx(hx~=Inf & hx~=-Inf));

pnin=likpop1/2+likpop2/2;
n1=likpop1.*log2(likpop1./pnin);
n2=likpop2.*log2(likpop2./pnin);
perf_ind(SH,PP)=0.5* nansum ( n1(n1~=Inf & n1~=-Inf))+0.5* nansum ( n2(n2~=Inf & n2~=-Inf));

post_ind(SH,PP,:)=0.5*likpop1./pnin;%./likpop2;%0.5*likpop1./pnin;
post_ind(SH,PP,isnan(post_ind(SH,PP,:)) | isinf(post_ind(SH,PP,:)))=0;
% post_pop_pls(SH,PP)=mean(post(1:5,1));
% post_pop_neg(SH,PP)=mean(post(6:10,1));


end
end

figure
errorbar(mean(perf_single),std(perf_single)/10)
% hold on
% errorbar(mean(perf_ind(1,:),1),std(perf_ind(1,:),1)/10,'r')
% hold on
% errorbar(mean(perf_ind(2,:),1),std(perf_ind(2,:),1)/10,'k')
hold on
plot(mean(perf_ind(1,:),1),'r')
hold on
plot(mean(perf_ind(2,:),1),'k')

figure
plot(squeeze(mean(mean(post_ind,1),3)))
hold on
plot(squeeze(mean(mean(post_single,1),3)),'r')


figure;plot(squeeze((mean(post_ind,1))))


% figure
% plot(mean(perf_single),'.')
% hold on
% plot(mean(perf_ind(1:5,:)),'.r')
% hold on
% plot(mean(perf_ind(6:10,:)),'.k')
% figure
% errorbar(mean(post_sing_pls),std(post_sing_pls)/10,'--b')
% hold on
% errorbar(mean(post_sing_neg),std(post_sing_neg)/10,'--r')
% hold on
% errorbar(mean(post_pop_pls),std(post_pop_pls)/10,'b')
% hold on
% errorbar(mean(post_pop_neg),std(post_pop_neg)/10,'r')
% 
% 
% figure
% plot(diff(mean(post_pop_pls))>0,'.')

return
