function LabImage = h_SpectralClusteringSegmentation(BinaryImage,K,nFrame)
% global nFrame;
    warning('off','stats:kmeans:EmptyCluster');
    
    if( K == 1 )
        LabImage = double(BinaryImage>0);
        return;
    end

    LabImage = zeros(size(BinaryImage));
    for i=1:nFrame
       binImage = BinaryImage(:,:,i);
       
       [r,c] = find(binImage);
       if (length(r)<50);continue;end
       % compute distance matrix %
       X = [r c];
       W = pdist2(X,X);
       W = exp(-(W.^2)/4);
       n = length(r);
       W = sparse(W);
       [k l s] = find(W);
%        s = s.*(s>1e-5);
       s = s.*(s>0.6);
       W = sparse(k,l,s,n,n);
       % degree and regularization %
       d = sum(abs(W),2);
       D = diag(d);
       D = sparse(D); 
       L = D - W;
       % compute normalized laplacian %
       d(d==0) = eps;
       D = spdiags(1./d,0,size(D,1),size(D,2));
       L = D * L;

       % compute the k-smallest eigenvalues
       if(isinf(condest(L)));continue;end
       sigma = eps;
       [U,~] = eigs(L,K,sigma);
       U = real(U);
       % normalize rowise %
       U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
%        for jj = 1:n
%            temp = U(jj,:);
%            U(jj,:) = U(jj,:)/norm(temp);
%        end
       
       
       warning('off','stats:kmeans:FailedToConvergeRep')
       % perform kmeans on the eigenvectors %
       [label,~] = kmeans(U,K,'replicates',10,'emptyaction','singleton',...
                            'start','cluster');
       
       % relabel %
       for j=1:length(label)
           LabImage(r(j),c(j),i) = label(j);
       end

       
    end





end
