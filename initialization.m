  function initialization(varargin)
  
% begin
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('begin: %s\n',time)
  
% default and user-defined arguments
  %folder = 'data_breast_v2';
  folder = 'data_nhs_subjects';
  args = varargin;
  for i = 1:2:length(args)
      switch args{i}
          case 'folder', folder = args{i+1};
      end
  end
  datadir = [cd '\' folder];
  
% generate "A_obsVer" folder and its subfolders
  obsVer = [datadir '\A_obsVer'];
  if (exist(obsVer,'dir') ~= 7)
      mkdir(obsVer)
  end
  f = [obsVer '\a_txtVer'];
  if (exist(f,'dir') ~= 7)
      mkdir(f)
  end
  f = [obsVer '\b_matVer'];
  if (exist(f,'dir') ~= 7)
      mkdir(f)
  end
  
% generate "B_varVer" folder and its subfolders
  varVer = [datadir '\B_varVer'];
  if (exist(varVer,'dir') ~= 7)
      mkdir(varVer)
  end
  f = [varVer '\a_one2oneVer'];
  if (exist(f,'dir') ~= 7)
      mkdir(f)
  end
  f = [varVer '\b_grp2idxVer'];
  if (exist(f,'dir') ~= 7)
      mkdir(f)
  end
  f = [varVer '\c_snp2binVer'];
  if (exist(f,'dir') ~= 7)
      mkdir(f)
  end
  
% generate "C_IterMMPC" folder
  f = [datadir '\C_IterMMPC'];
  if (exist(f,'dir') ~= 7)
      mkdir(f)
  end
  
% generate "D_OptNBC" folder
  f = [datadir '\D_OptNBC'];
  if (exist(f,'dir') ~= 7)
      mkdir(f)
  end
  
% end
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('e n d: %s\n',time)
  