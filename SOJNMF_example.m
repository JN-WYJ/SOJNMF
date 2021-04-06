clear;clc;
cd('D:/SOJNMF')
%import data
GE = importdata('co_GE.txt');
ME = importdata('co_ME.txt');
DM = importdata('co_DM.txt');

genes = GE.textdata;
miRNAs = ME.textdata;
methylations = DM.textdata;

X1 = GE.data;
X2 = ME.data;
X3 = DM.data;

[n,m1]=size(X1);
[n,m2]=size(X2);
[n,m3]=size(X3);

nloop=100; % 100; 
verbose=1;
K = 238;alpha = 0.001;lambda = 0.1;
maxiter=1000; speak=1;

bestW=zeros(n,K);
bestH1=zeros(K,m1);
bestH2=zeros(K,m2);
bestH3=zeros(K,m3);

bestobj1=1000000000;
bestobj2=1000000000;
bestobj3=1000000000;

for iloop=1:nloop;
    if verbose 
        fprintf(1,' iteration %d\n',iloop); 
    end 
    
    [W,H1,H2,H3]=SOJNMF(X1, X2, X3, K, alpha, lambda, maxiter, speak);

    % compute residue
    newobj1 = sum(sum((X1-W*H1).^2));
    newobj2 = sum(sum((X2-W*H2).^2));
    newobj3 = sum(sum((X3-W*H3).^2));
    
    if (newobj1<bestobj1)||(newobj2<bestobj2)||(newobj3<bestobj3)          
        bestobj1 = newobj1;
        bestobj2 = newobj2;  
        bestobj3 = newobj3; 
        bestW=W;
        bestH1=H1;
        bestH2=H2;
        bestH3=H3;
    end
end

% compute the modules according to bestW, bestH1 and bestH2
W=bestW; H1=bestH1; H2=bestH2; H3=bestH3;

% Output rule
 tt0 = 2; tt1 = 2.5; tt2 = 2.5; tt3 = 2.5;
 
[ Co_module, Subpattern1, Suppattern2, Suppattern3]= SOJNMF_module(X1, X2, X3, W, H1, H2, H3, tt0, tt1, tt2, tt3);

module_file = 'Co_module_SOJNMF';
if ~isdir(module_file)
    mkdir(module_file);
end
cd('./Co_module_SOJNMF')
save SOJNMF_Comodule.mat

Co_genes = {};
for igenes=1:K;
    Co_genes(igenes,1:length(Co_module{igenes,2})) = genes(1,Co_module{igenes,2}+1);
end
xlswrite('Co_genes.xlsx', Co_genes);

Co_miRNA = {};
for imiRNAs=1:K;
    Co_miRNA(imiRNAs,1:length(Co_module{imiRNAs,3})) = miRNAs(1,Co_module{imiRNAs,3}+1);
end
xlswrite('Co_miRNA.xlsx', Co_miRNA);

Co_methy = {};
for imethy=1:K;
    Co_methy(imethy,1:length(Co_module{imethy,4})) = methylations(1,Co_module{imethy,4}+1);
end
xlswrite('Co_methy.xlsx', Co_methy);
