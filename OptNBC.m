  function R = OptNBC(varargin)
  
% OptNBC: optimal NB classifier.
%
% CopyRight: Xu-Qing Liu (liuxuqing688@163.com) // Aug. 08, 2019
  
% begin
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('begin: %s\n',time)
  
% default and user-defined arguments
  folder = 'data_breast_v2';
  rsOnly = false;
  delMax = false;
  args = varargin;
  for i = 1:2:length(args)
      switch args{i}
          case 'folder', folder = args{i+1};
          case 'rsOnly', rsOnly = args{i+1};
          case 'delMax', delMax = args{i+1};
      end
  end
  if (delMax)
      delMax_flag = '_delMax';
  else
      delMax_flag = '_useMax';
  end
  
% paths
  datadir = [cd '\' folder];
  if (~rsOnly)
      To = [datadir '\D_OptNBC\a_snpAll\'];
      rsOnly_flag = '_snpAll';
  else
      To = [datadir '\D_OptNBC\b_rsOnly\'];
      rsOnly_flag = '_rsOnly';
  end
  
% data pre-processing
  load([datadir '\SNP'],'SNP')
  load([datadir '\X_IterMMPC.mat'],'f')
  snp = SNP(f);
  clearvars SNP
  nv = length(snp);
  if (rsOnly)
      idx_rs = true(1,nv);
      for i = 1:nv
          if (~strcmp(snp{i}(1:2),'rs'))
              idx_rs(i) = false;
          end
      end
      snp = snp(idx_rs);
      f = f(idx_rs);
      R.idx = find(idx_rs);
  else
      R.idx = 1:nv;
  end
  R.var_in_all = f;
  R.snp_in_all = snp;
  
% load X-data and y-data
  load([datadir '\X_IterMMPC.mat'],'X')
  load([datadir '\y.mat'],'y')
  R.y = y;
  nc = length(y);
  nv = length(snp);
  
% prior and posterior based on "Laplace Smoothing (LS)"
  [~,~,L,M] = LS(X,y);
  L_y = L;
  
% FW (forward search)
  I = false(1,nv);
  A = [];
  s_old = -inf;
  a = 0;
  while a < 0.999
      Ls = cell(1,nv);
      Ss = -inf(1,nv);
      for j = find(~I)
          L_j = zeros(nc,2);
          for i = 1:nc
              for k = 1:2
                  L_j(i,k) = L(i,k)*M(i,j,k);
              end
          end
          Ls{j} = L_j;
          Ss(j) = sum(log(L_j(:,1))-log(sum(L_j,2)));
      end
      if (delMax)
          [~,j] = max(Ss); Ss(j) = -inf;
      end
      [s_new,j] = max(Ss);
      if (s_new > s_old)
          L = Ls{j};
          I(j) = true;
          [~,~,C] = nbc_kfold(X(:,I),y,nc);
          %[~,~,C] = nbc_kfold(X(:,I),y,10);
          a = trace(C)/nc
          A = [A a];
          s_old = s_new;
      else
          break
      end
  end
  R.FW_var_in_iMMPC = R.idx(I);
  R.FW_X = X(:,R.FW_var_in_iMMPC);
  R.FW_var_in_all = f(I);
  R.FW_snp_in_all = snp(I);
  R.FW_accuracy = A;
  save(sprintf('%s%s%s%s%s',datadir,'\R',rsOnly_flag,delMax_flag),'R')
  
% BW (backward search)
  while a > 0.995
      Ss = -inf(1,nv);
      for j = find(I)
          L_j = zeros(nc,2);
          F_j = find(I);
          F_j(F_j == j) = [];
          for i = 1:nc
              for k = 1:2
                  L_j(i,k) = L_y(i,k)*prod(M(i,F_j,k));
              end
          end
          Ss(j) = sum(log(L_j(:,1))-log(sum(L_j,2)));
      end
      [s_new,j] = max(Ss);
      if (s_new >= s_old)
          I(j) = false;
          [~,~,C] = nbc_kfold(X(:,I),y,nc);
          %[~,~,C] = nbc_kfold(X(:,I),y,10);
          a = trace(C)/nc
          A = [A a];
          s_old = s_new;
      else
          break
      end
  end
  R.BW_var_in_iMMPC = R.idx(I);
  R.BW_X = X(:,R.BW_var_in_iMMPC);
  R.BW_var_in_all = f(I);
  R.BW_snp_in_all = snp(I);
  R.BW_accuracy = A;
  save(sprintf('%s%s%s%s%s',datadir,'\R',rsOnly_flag,delMax_flag),'R')
  
% end
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('e n d: %s\n',time)
  
% --------
  function [P_Xy,p_y,priorProb,postProb] = LS(X,y)
  [nc,nv] = size(X);
  
% joint distribution based on "Laplace Smoothing (LS)"
  P_Xy = cell(1,nv);
  for j = 1:nv
      P_j = reshape(MIToolboxMex(2,X(:,j),y),2,2);
      P_Xy{j} = (P_j*nc+1)/(nc+4);
  end
  
% prior probabilities of each y(i) based on "LS"
  %p_y = MIToolboxMex(1,y);
  p_y = (MIToolboxMex(1,y)*nc+1)/(nc+2);
  priorProb = nan(nc,2);
  for j = 1:2
      priorProb(y == j,1) = p_y(j);
      priorProb(y == j,2) = 1-p_y(j);
  end
  
% posterior probabilities of X(i,j) conditioned on y(i)
  postProb = nan(nc,nv,2);
  for j = 1:nv
      p_j = P_Xy{j};
      for i = 1:nc
          j1 = y(i);
          j2 = 3-y(i);
          postProb(i,j,1) = p_j(X(i,j),j1)/sum(p_j(:,j1));
          postProb(i,j,2) = p_j(X(i,j),j2)/sum(p_j(:,j2));
      end
  end
  