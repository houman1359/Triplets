function [copccdf coppdf copccdf_ker coppdf_ker co]=kernel_all(u,KER,dim,bw,METH)


% generate a random name for the temporary files that will be used to
% exchange data between MATLAB and R.



% if ~exist('/n/scratch2/TEMP','dir')
% mkdir('/n/scratch2/TEMP/tmp');
% end


% temp_filename = ['/n/scratch2/TEMP' tempname];
temp_filename = ['/home/hs258/Data_Folder/TEMP' tempname];


data_filename_input = [temp_filename, '.in.mat'];
data_filename_output = [temp_filename, '.out.mat'];

% prepare R options depending on the debug flag
r_opts = '--no-save --no-restore';
r_out_file = [temp_filename, '.log'];
r_err_file = [temp_filename, 'error.log'];

% save data and parameters for fit to temporary file and call R

if nargin<2
    dim=2;
end

copccdf=u(:,1)*0;
copcdf=u(:,1)*0;
coppcdf=u(:,1)*0;

n=1;
V=floor(size(u,1)/n);

% system('module load gcc/6.2.0')


if sum(METH.server(1:5)=='login')==5
setenv('LD_LIBRARY_PATH', '/n/app/gcc/6.2.0/lib64/');
system('module load gcc/6.2.0');
system('module load harfbuzz/1.3.4 cairo/1.14.6');
system('module load R/3.2.5');
else
% setenv('LD_LIBRARY_PATH', '/opt/openblas/0.2.14/lib/');
% setenv('LD_LIBRARY_PATH', '/opt/openblas/0.2.14/lib/:/opt/gcc/4.8.5/lib64/');
system('module load stats/R/3.3.1');
end
%         save(['fff00.mat'])

% HOST=char(getHostName(java.net.InetAddress.getLocalHost));
% setenv('LD_LIBRARY_PATH', '');
% HOST=getenv('HOSTNAME');


for ii=1:n
    
    if 1==1%rand<0.6
        if sum(METH.server(1:5)=='login')==5
%             command = sprintf('sbatch -p short -n 1 -t 0-4:00 --mem-per-cpu=7G --wrap="R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s"',...
%                 r_opts, [temp_filename num2str(ii)], r_out_file);
            command = sprintf('R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
        else
            if 1==1%rand<0.5
%             command = sprintf('bsub -q short -n 1 -W 20 R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
%                 r_opts, [temp_filename num2str(ii)], r_out_file);
            command = sprintf('R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
            else
            command = sprintf('bsub -q mcore -n 2 -W 20 R CMD BATCH %s''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
                r_opts, [temp_filename num2str(ii)], r_out_file);   
            end
        end
    else
        if sum(METH.server(1:5)=='login')==5
            command = sprintf('sbatch -p short -n 1 -t 0-5:00 --mem-per-cpu=10G --wrap="R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s"',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
        else
            command = sprintf('bsub -q mcore -n 2 -W 10:0 R CMD BATCH %s''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
        end
    end

    data_filename_output = [temp_filename, num2str(ii),'.out.mat'];
    data_filename_input = [temp_filename, num2str(ii),'.in.mat'];
if ii<n
    DD=(ii-1)*V+1:ii*V;
else
    DD=(n-1)*V+1:size(u,1);   
end
    e=u(DD,1:2);
    ke=KER(:,1:2);
    method=METH.method;%'TLL2';%'TLL1nn';    %%%%%%% this is used for Matthias
%     method='TLL1';  %%%% try for Selmaan
%     method='TTCV';  %%%% 
    
    save(data_filename_input, 'e','dim','ke','method','bw');
        
    system(command);
   
end

for ii=1:n
        
    data_filename_output = [temp_filename, num2str(ii),'.out.mat'];
    data_filename_input  = [temp_filename, num2str(ii),'.in.mat'];
    
    du=1;
    ERR=0;
    DTA={'ee','eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'};

    while (exist(data_filename_output,'file')==0 & du==1)
        du=1;
        if exist(r_out_file,'file')~=0
            DTA = importdata(r_out_file, ',');
            if numel(DTA)~=0
            if numel(cell2mat(DTA(numel(DTA))))==16
            if sum(cell2mat(DTA(numel(DTA)))=='Execution halted')==16
            du=0;
            ERR=1;
            end
            end
            end
        end
    end


    

if ii<n
    DD=(ii-1)*V+1:ii*V;
else
    DD=(n-1)*V+1:size(u,1);
end

if ERR==0
    kpdf = load(data_filename_output);

else
    
    data_filename_output = [temp_filename, num2str(ii),'.out.mat'];
    data_filename_input = [temp_filename, num2str(ii),'.in.mat'];
if ii<n
    DD=(ii-1)*V+1:ii*V;
else
    DD=(n-1)*V+1:size(u,1);   
end


    e=u(DD,1:2);
    ke=KER(:,1:2);
    method='T';  %%%% try for Selmaan
    
    if 1==1%rand<0.6
        if sum(METH.server(1:5)=='login')==5
%             command = sprintf('sbatch -p short -n 1 -t 0-3:00 --mem-per-cpu=7G --wrap="R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s"',...
%                 r_opts, [temp_filename num2str(ii)], r_out_file);
            command = sprintf('R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
        else
            command = sprintf('R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
        end
    else
        if sum(METH.server(1:5)=='login')==5
            command = sprintf('sbatch -p short -n 1 -t 0-1:00 --mem-per-cpu=10G --wrap="R CMD BATCH %s ''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s"',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
        else
            command = sprintf('bsub -q mcore -n 2 -W 10:0 R CMD BATCH %s''--args %s'' /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno/kern_all.R %s',...
                r_opts, [temp_filename num2str(ii)], r_out_file);
        end
    end
    
    
    save(data_filename_input, 'e','dim','ke','method','bw');

    system(command);
    
    data_filename_output = [temp_filename, num2str(ii),'.out.mat'];
    data_filename_input  = [temp_filename, num2str(ii),'.in.mat'];
    
   
    while (exist(data_filename_output,'file')==0)
        du=1;
    end
% kpdf = load(data_filename_output);
% else
% kpdf.pdf_ker=ones(size(KER,1),1)/size(KER,1);
% kpdf.pdf=ones(size(u,1),1)/size(KER,1);
% 
% [nn1 cc1]=sort(KER(:,1));
% [nn2 cc2]=sort(KER(:,2));
% 
% kpdf.ccdf=[];kpdf.ccdf_ker=[];
% for i=1:size(KER,1)
% kpdf.ccdf_ker(i,1)=cc1(i)*cc2(i)/(size(KER,1))^2;
% end
% 
% for i=1:size(u(DD,:),1)
%     [m1 n1]=min(abs(u(i,1)-KER(:,1)));
%     [m1 n2]=min(abs(u(i,1)-KER(:,1)));
% kpdf.ccdf(i,1)=cc1(n1)*cc2(n2)/(size(KER,1))^2;
% end


kpdf = load(data_filename_output);


end


    copccdf(DD)=kpdf.ccdf;
    coppdf(DD)=kpdf.pdf;
    if ii==n
    copccdf_ker=kpdf.ccdf_ker;
    coppdf_ker=kpdf.pdf_ker;
    end
    
display('one copula done')

%     % clean up temporary files
delete(data_filename_input);
delete(data_filename_output);
if ii==n
delete(r_out_file);
end

end

co=1;

