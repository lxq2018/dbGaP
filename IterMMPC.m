  function Fs = IterMMPC(varargin)
  
% IterMMPC: iterated version of MMPC algorithm used
%           for a problem with very high dimension.
%
% CopyRight: Xu-Qing Liu (liuxuqing688@163.com) // Aug. 08, 2019
  
% begin
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('begin: %s\n',time)
  
% default and user-defined arguments
  folder = 'data_breast_v2';
  alpha = 0.1;
  nPerP = 10;
  k_max = 2;
  args = varargin;
  for i = 1:2:length(args)
      switch args{i}
          case 'folder', folder = args{i+1};
          case 'alpha',  alpha  = args{i+1};
          case 'nPerP',  nPerP  = args{i+1};
          case 'k_max',  k_max  = args{i+1};
      end
  end
  opt.threshold = alpha;
  opt.epc = 5;
  opt.use_card_lim = 0;
  opt.max_card = 0;
  opt.maxK = k_max;
  
% paths
  datadir = [cd '\' folder];
  X_path = [datadir '\B_varVer\c_snp2binVer\'];
  y_path = datadir;
  To = [datadir '\C_IterMMPC\'];
  
% load y-data
  load([y_path '\y.mat'],'y')
  nc = length(y);
  c = 2*nc*log(2);
  cv = chi2inv(1-alpha,(max(y)-1)*(1:2))/c;
  
% preliminary-filtering: 0-order
  xNames = dir([X_path 'x_*.mat']);
  nv = length(xNames);
  F = true(1,nv);
  I = nan(1,nv);
  for j = 1:nv
      [j nv]
      load([X_path xNames(j).name],'x')
      I_xy = mi(x,y);
      I(j) = I_xy;
      if (I_xy <= cv(1))
          F(j) = false;
      end
  end
  Fs.pre_0 = find(F)
  save([datadir '\xNames.mat'],'xNames')
  save([datadir '\I_all.mat'],'I')
  save([datadir '\Fs.mat'],'Fs')
  
% preliminary-filtering: 1-order
  [~,j_cond] = max(I);
  load([X_path xNames(j_cond).name],'x')
  x_cond = x;
  for j = find(F)
      load([X_path xNames(j).name],'x')
      I_xyz = cmi(x,y,x_cond);
      if (I_xyz <= cv(2))
          F(j) = false;
      end
  end
  F(j_cond) = true;
  Fs.pre_1 = find(F)
  save([datadir '\Fs.mat'],'Fs')
  
% further-filtering: based on MMPC (randomly take one)
  s_rng = rng(123456);
  n_rem = ones(1,100);
  while mean(n_rem) <= 0.8*nPerP %sum(F) > 1e4
      f = find(F);
      f_tmp = f(randperm(length(f),nPerP));
      I_del = true(1,nPerP);
      X_tmp = nan(nc,nPerP);
      for j_tmp = 1:nPerP
          j = f_tmp(j_tmp);
          load([X_path xNames(j).name],'x')
          X_tmp(:,j_tmp) = x;
      end
      Xy = [X_tmp y]-1;
      tg = nPerP+1;
      ns = max(Xy)+1;
      pc = Causal_Explorer('MMPC',Xy,tg,ns,'MMPC',opt);
      n_rem = [n_rem(2:end) length(pc)];
      I_del(pc) = false;
      F(f_tmp(I_del)) = false;
      [length(pc) nPerP sum(F)]
  end
  rng(s_rng)
  Fs.fur = find(F)
  save([datadir '\Fs.mat'],'Fs')
  
% final-filtering: based on MMPC (randomly partition all)
  s_rng = rng(123456);
  while 1
      n_F = sum(F);
      f = find(F);
      k_par = floor(n_F/nPerP);
      K = crossvalind('kfold',n_F,k_par)';
      for k = 1:k_par
          f_k = f(K == k);
          n_k = length(f_k);
          x_k = nan(nc,n_k);
          I_k = true(1,n_k);
          for j_k = 1:n_k
              j = f_k(j_k);
              load([X_path xNames(j).name],'x')
              x_k(:,j_k) = x;
          end
          Xy = [x_k y]-1;
          tg = n_k+1;
          ns = max(Xy)+1;
          pc = Causal_Explorer('MMPC',Xy,tg,ns,'MMPC',opt);
          I_k(pc) = false;
          F(f_k(I_k)) = false;
          [length(pc) n_k]
      end
      %[sum(F) n_F]
      if (sum(F)/n_F > 0.99)
          break
      end
  end
  rng(s_rng)
  Fs.final = find(F)
  save([datadir '\Fs.mat'],'Fs')
  
% copy final-filtered data
  f = find(F);
  nv = length(f);
  X = nan(nc,nv);
  for j_tmp = 1:nv
      j = f(j_tmp);
      load([X_path xNames(j).name],'x')
      X(:,j_tmp) = x;
      copyfile([X_path xNames(j).name],[To xNames(j).name])
  end
  save([datadir '\X_IterMMPC.mat'],'X','f')
  
% end
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('e n d: %s\n',time)
  