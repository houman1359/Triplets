function I=CONF_INFO(REAL,PREDICT)

CM=confusionmat(REAL,PREDICT);

pssp=CM/sum(sum(CM));


for i=1:size(CM,1) 
    ps(i)=squeeze(sum(pssp(i,:)));
    psp(i)=squeeze(sum(pssp(:,i)));
end

I=0;

for i=1:size(CM,1)
    for j=1:size(CM,2)
    if pssp(i,j)~=0 & ps(i)~=0 & psp(j)~=0    
I=I+pssp(i,j)*log2(pssp(i,j)/(ps(i)*psp(j)));
    end
    end
end