function [fcon,f,y_vector1]=OneD_marginalize(copul,vine_stim,DIM,N,Tn,B,tvec)

%%%%% N  = the component of time =2 usually
%%%%% TN = the time point which we compute the 2D densities

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            samples=vinecopula_sample(vine_stim,copul,tvec,100,1000,2);    %%%% nn is time

            y1=vine_stim.margins{DIM(1)}.ker;
            y_vector1 =linspace(min(y1),max(y1)+eps,B(1));%unique(y);
          
            
            ne=0;
            TN=find(samples(:,N)>tvec(Tn-1) & samples(:,N)<tvec(Tn));
            poi=zeros(numel(y_vector1)*numel(y_vector2)*sum(samples(:,2)>tvec(Tn-1) & samples(:,2)<tvec(Tn)),numel(vine_stim.margins));
            for j=1:numel(TN)
                    for i1=1:numel(y_vector1)
                        ne=ne+1;
                        poi(ne,DIM)=[y_vector1(i1)];
                        poi(ne,setdiff(1:numel(vine_stim.margins),DIM))=samples(TN(j),setdiff(1:numel(vine_stim.margins),DIM));
                    end
            end
                        
            [f_points_train,~,~,~] = Fit_vCopula(vine_stim,poi,1000,[],-1,copul,'rand',[]);
            
            ff=reshape(f_points_train,B(1),numel(TN));
            f=squeeze(sum(ff,3));
            ffn=sum(f);
            ffn=repmat(ffn,size(f,1),1);
            
            fcon=f./ffn;