
% --- diffusionKernel function
% Written by R. Coifman & S. Lafon.
function [Y] = diffusionKernel (X,sigmaK,alpha,d)
D = L2_distance(X',X',1);
K = exp(-(D/sigmaK).^2);
p = sum(K);
p = p(:);
K1 = K./((p*p').^alpha);
v = sqrt(sum(K1));
v = v(:);
A = K1./(v*v');
if sigmaK >= 0.5
    thre = 1e-7;  
    M = max(max(A));
    A = sparse(A.*double(A>thre*M));
    [U,S,V] = svds(A,d+1);   %Sparse version.
    U = U./(U(:,1)*ones(1,d+1));
else
    [U,S,V] = svd(A,0);   %Full version.
    U = U./(U(:,1)*ones(1,size(U,1)));
end;
Y = U(:,2:d+1);


% --- mgs function: Modified Gram-Schmidt
% Used by HLLE function.
%%
function [Q, R] = mgs(A);
[m, n] = size(A);   % Assume m>=n.
V = A;
R = zeros(n,n);
for i=1:n
    R(i,i) = norm(V(:,i));
    V(:,i) = V(:,i)/R(i,i);
    if (i < n)
        for j = i+1:n
            R(i,j) = V(:,i)' * V(:,j);
            V(:,j) = V(:,j) - R(i,j) * V(:,i);
        end;
     end;
 end;
 Q = V;
% --- L2_distance function
% Written by Roland Bunschoten, University of Amsterdam, 1999
function d = L2_distance(a,b,df)
if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end
aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
d = real(d); 
if (df==1)
  d = d.*(1-eye(size(d)));
end
