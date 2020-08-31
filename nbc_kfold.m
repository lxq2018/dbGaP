  function [yhat,post,C,Pre,Rec,Acc,Spe,Mcc,r] = nbc_kfold(X,y,kfold)
  
% default setting
  nc = size(X,1);
  r_y = max(y);
  if (nargin < 3)
      kfold = nc;
  end
  
% k-fold
  if (kfold == nc)
      K = 1:nc;
  else
      K = nan(nc,1);
      for i = 1:r_y
          I_i = find(y == i)';
          nc_i = length(I_i);
          K_i = crossvalind('kfold',nc_i,kfold);
          K(I_i) = K_i;
      end
  end
  
% nbc
  yhat = nan(nc,1);
  post = nan(nc,r_y);
  for i = 1:kfold
      I = (K == i);
      M = nbc_fit(X(~I,:),y(~I));
      [post(I,:),yhat(I)] = nbc_post(M,X(I,:));
      %M = NaiveBayes.fit(X(~I,:),y(~I),'Distribution','mvmn');
      %post(I,:) = posterior(M,X(I,:));
      %yhat(I) = predict(M,X(I,:));
%       C_i = confusionmat(yhat(I),y(I));
%       TP = C_i(2,2);              % number of true  positives
%       FP = C_i(2,1);              % number of false positives
%       FN = C_i(1,2);              % number of false negatives
%       TN = C_i(1,1);              % number of true  negatives
      TP = sum(y(I) == 2 & yhat(I) == 2);
      FP = sum(y(I) == 1 & yhat(I) == 2);
      FN = sum(y(I) == 2 & yhat(I) == 1);
      TN = sum(y(I) == 1 & yhat(I) == 1);
      r.Pre(i) = TP/(TP+FP);         % precision
      r.Rec(i) = TP/(TP+FN);         % recall (=sensitivity)
      r.Acc(i) = (TP+TN)/(TP+FP+FN+TN);%sum(C_i(:));  % accuracy
      r.Spe(i) = TN/(TN+FP);         % specification
      r.Mcc(i) = (TP*TN-FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FP)*(TN+FN));
  end
  if (nargout > 2)
      if (r_y ~= 2)
          error('Not 2-class data!')
      end
      C = confusionmat(yhat,y); % confusion matrix
      TP = C(2,2);              % number of true  positives
      FP = C(2,1);              % number of false positives
      FN = C(1,2);              % number of false negatives
      TN = C(1,1);              % number of true  negatives
      Pre = TP/(TP+FP);         % precision
      Rec = TP/(TP+FN);         % recall (=sensitivity)
      Acc = (TP+TN)/sum(C(:));  % accuracy
      Spe = TN/(TN+FP);         % specification
      C = confusionmat(3-yhat,3-y);
      Mcc = (TP*TN-FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FP)*(TN+FN));
  end
  
% --------
  function M = nbc_fit(X,y)
  %isnan_y = isnan(y);
  %X(isnan_y,:) = [];
  %y(isnan_y) = [];
  %if (isempty(y))
  %    error('!!!')
  %end
  r_y = max(y);
  nv = size(X,2);
  p_y = MIToolboxMex(1,y)';
  p_X = cell(r_y,nv);
  %E = ones(3,r_y);
  for j = 1:nv
      x_j = X(:,j);
      y_j = y;
      %isnan_x = isnan(x_j);
      %x_j(isnan_x) = [];
      %y_j(isnan_x) = [];
      nc_j = length(y_j);
      %if (nc_j == 0)
      %    error('!!!')
      %end
      min_x = min(x_j);
      max_x = max(x_j);
      r_x = max_x-min_x+1;
      q_x = [ones(min_x-1,r_y); ...
            reshape(MIToolboxMex(2,x_j,y_j),r_x,r_y)*nc_j+1; ...
            ones(3-max_x,r_y)];
      %q_x = [E(1:min_x-1,:); ...
      %       reshape(MIToolboxMex(2,x_j,y_j),r_x,r_y)*nc_j+1; ...
      %       E(1:3-max_x,:)];
      for i = 1:r_y
          q_x_i = q_x(:,i);
          p_X{i,j} = q_x_i/sum(q_x_i);
      end
  end
  M.NClasses = r_y;
  M.NDims = nv;
  M.Params = p_X;
  M.Prior = p_y;
  
% --------
  function [post,yhat] = nbc_post(M,x)
      
% For two very small positive numbers, x and y, compute z = x/(x+y).
% Letting a = x/max(x,y) and b = y/max(x,y), then z = a/(a+b).
% However, a and b cannot be computed directly. Instead, using:
% a = exp(log(x)/c), b = exp(log(y)/c), with c = max{log(x),log(y)}.
  
  [nc,nv] = size(x);
  %nc = size(x,1);
  logPrior = log(M.Prior);
  params = M.Params;
  logPost = nan(nc,M.NClasses);
  for k = 1:nc
      x_k = x(k,:);
      J_k = find(~isnan(x_k));
      logPost_k = logPrior;
      for i = 1:M.NClasses
          params_i = ones(1,nv);%
          for j = J_k
              params_i(j) = params{i,j}(x_k(j));%
              %logPost_k(i) = logPost_k(i)+log(params{i,j}(x_k(j)));
          end
          logPost_k(i) = logPost_k(i)+sum(log(params_i));%
      end
      logPost(k,:) = logPost_k;
  end
  RF = max(logPost,[],2); % RF: regularization factor
  p = exp(bsxfun(@minus,logPost,RF));
  post = bsxfun(@rdivide,p,sum(p,2));
  [~,yhat] = max(post,[],2);
  