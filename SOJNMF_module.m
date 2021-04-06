function [ Co_module, Subpattern1, Subpattern2, Subpattern3]=SOJNMF_module(X1, X2, X3, W, H1, H2, H3, tt0, tt1, tt2, tt3) 
%
% Compute the mean and std in columns of W and rows in H1, H2 and H3 to determine
% the module member and output the Co-module based on W and H1, H2, H3 matrices.
% Co_module is the list of sample and SNPs, Genes in a Co-module.
% Subpattern1, Subpattern2 and Subpattern3 is the corresponding sub-matrices for each Co-module.
% 
m1 = size(H1,2);
m2 = size(H2,2);
m3 = size(H3,2);
n = size(W,1);
K = size(W,2);

MW =mean(W,1);      MH1 =mean(H1,2); MH2 =mean(H2,2); MH3 =mean(H3,2);

VW =std(W,0,1);     VH1 =std(H1,0,2);  VH2 =std(H2,0,2); VH3 =std(H3,0,2);

% Co-Module
for i=1:K
   c1=find(H1(i,:)> MH1(i) + tt1*VH1(i));
   module1{i,1}=c1'; 

   c2=find(H2(i,:)> MH2(i) + tt2*VH2(i));
   module2{i,1}=c2'; 
   
   c3=find(H3(i,:)> MH3(i) + tt3*VH3(i));
   module3{i,1}=c3'; 
   
   r=find(W(:,i)> MW(i) + tt0*VW(i));
   Co_module{i,1}=r'; Co_module{i,2}=c1'; Co_module{i,3}=c2'; Co_module{i,4}=c3';

   Subpattern1{i}=X1(r,c1');
   Subpattern2{i}=X2(r,c2');
   Subpattern3{i}=X3(r,c3');
end

