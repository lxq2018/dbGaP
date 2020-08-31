  function obs2var_snp2bin(varargin)
  
% begin
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('begin: %s\n',time)
  
% default and user-defined arguments
  folder = 'data_breast_v2';
  args = varargin;
  for i = 1:2:length(args)
      switch args{i}
          case 'folder', folder = args{i+1};
      end
  end
  
% paths
  datadir = [cd '\' folder];
  AB = [datadir '\B_varVer\a_one2oneVer\'];
  Fr = [datadir '\B_varVer\b_grp2idxVer\'];
  To = [datadir '\B_varVer\c_snp2binVer\'];
  
% read grp2idx-based variable files
  varFiles = dir([Fr '*.mat']);
  
% read y-data
  load([datadir '\y.mat'],'y')
  
% transformation
  nv = length(varFiles);
  for j = 1:nv
      load([Fr varFiles(j).name],'x','g')
      coded_as_0 = true(1,length(g)); % length(g) == max(x)
      [x,coded_as_1] = Snp2Bin(x,y);
      coded_as_0(coded_as_1) = false;
      g_0 = g(coded_as_0); % alleles of taking 0
      g_1 = g(coded_as_1); % alleles of taking 1
      save([To sprintf('x_%07.0f.mat',j)],'x','g_0','g_1')
  end
  
% end
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('e n d: %s\n',time)
  