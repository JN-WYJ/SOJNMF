function [W,H1,H2,H3]=SOJNMF(X1, X2, X3, K, alpha, lambda, maxiter, speak)
%
% INPUT:
% X1 (N,M1): N (dimensionallity) x M1 (samples) non negative input matrix
% X2 (N,M2): N (dimensionallity) x M2 (samples) non negative input matrix
% X3 (N,M3): N (dimensionallity) x M3 (samples) non negative input matrix
% K        : Number of components
% lambda   : sparse parameter
% maxiter  : Maximum number of iterations to run
% speak    : prints iteration count and changes in connectivity matrix 
%            elements unless speak is 0
%
% OUTPUT:
% W        : N x K matrix
% H1       : K x M1 matrix
% H2       : K x M2 matrix
% H3       : K x M3 matrix
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_iter = 20; % iterations between print on screen and convergence test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X1 X2 and X3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (min(min(X1)) < 0) || (min(min(X2)) < 0) || (min(min(X3)) < 0)
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for same rows in X1 X2 and X3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,m1]=size(X1);
[n2,m2]=size(X2);
[n3,m3]=size(X3);

if (n1~=n2) || (n1~=n3)   
    error('Input matrices should have the same rows');
    return
end
n=n1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W and H1 ,H2, H3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W=rand(n,K);
H1=rand(K,m1);
H2=rand(K,m2);
H3=rand(K,m3);

% use W*H to test for convergence
Xr_old1=W*H1;
Xr_old2=W*H2;
Xr_old3=W*H3;

for iter=1:maxiter
    % Euclidean multiplicative method
    H1 = H1.*(W'*X1+alpha*H1)./((W'*W)*H1+2*alpha*H1*H1'*H1+lambda+eps);     
    H2 = H2.*(W'*X2+alpha*H2)./((W'*W)*H2+2*alpha*H2*H2'*H2+lambda+eps);   
    H3 = H3.*(W'*X3+alpha*H3)./((W'*W)*H3+2*alpha*H3*H3'*H3+lambda+eps);    
    W = W.*([H1 H2 H3]*[X1 X2 X3]')'./(W*([H1 H2 H3]*[H1 H2 H3]')+eps);   
  % iter
  % print to screen
    if (rem(iter,print_iter)==0) & speak,
        Xr1 = W*H1;            
        Xr2 = W*H2;
        Xr3 = W*H3;
        diff = sum(sum(abs(Xr_old1-Xr1)))+sum(sum(abs(Xr_old2-Xr2)))+sum(sum(abs(Xr_old3-Xr3)));
        Xr_old1 = Xr1;
        Xr_old2 = Xr2;
        Xr_old3 = Xr3;
        eucl_dist1 = nmf_euclidean_dist(X1,W*H1);
        eucl_dist2 = nmf_euclidean_dist(X2,W*H2);
        eucl_dist3 = nmf_euclidean_dist(X3,W*H3);
        eucl_dist = eucl_dist1 + eucl_dist2 + eucl_dist3;
        errorx1 = mean(mean(abs(X1-W*H1)))/mean(mean(X1));
        errorx2 = mean(mean(abs(X2-W*H2)))/mean(mean(X2));
        errorx3 = mean(mean(abs(X3-W*H3)))/mean(mean(X3));
        errorx = errorx1 + errorx2 + errorx3;
        disp(['Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end

function err = nmf_euclidean_dist(X,Y)

err = sum(sum((X-Y).^2));
