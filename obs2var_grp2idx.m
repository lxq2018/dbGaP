  function obs2var_grp2idx(varargin)
  
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
  Fr = [datadir '\B_varVer\a_one2oneVer\'];
  To = [datadir '\B_varVer\b_grp2idxVer\'];
  
% read one2one-based variable files
  varFiles = dir([Fr '*.mat']);
  
% transformation
  nv = length(varFiles);
  for j = 1:nv
      load([Fr varFiles(j).name],'x')
      [x,~,g] = grp2idx(x);
      save([To sprintf('x_%07.0f.mat',j)],'x','g')
  end
  
% end
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('e n d: %s\n',time)
  