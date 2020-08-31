  function obs2var_one2one(varargin)
  
% begin
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('begin: %s\n',time)
  
% default and user-defined arguments
  %folder = 'data_breast_v2';
  %n_vars = 561466;
  folder = 'data_nhs_subjects';
  n_vars = 546646;
  n_PerP = 100;
  args = varargin;
  for i = 1:2:length(args)
      switch args{i}
          case 'folder', folder = args{i+1};
          case 'n_vars', n_vars = args{i+1};
          case 'n_PerP', n_PerP = args{i+1};
      end
  end
  
% paths
  datadir = [cd '\' folder];
  Fr = 'C:\dataBC\';
  To = [datadir '\B_varVer\a_one2oneVer\'];
  
% list all SNPs and save them
  fid = fopen([Fr folder(6:end) '.bim'],'r');
  data = textscan(fid,'%d%s%d%d%s%s');
  fclose(fid);
  SNP = data{2};
  clearvars data
  save([datadir '\SNP'],'SNP')
  
% extract y-data
  fid = fopen([Fr folder(6:end) '.fam'],'r');
  data = textscan(fid,'%s%s%d%d%d%d');
  fclose(fid);
  sex = double(data{5}); % (1 = male, 2 = female, other = unknown)
  y = double(data{6});
  clearvars data
  save([datadir '\y'],'y','sex')
  nc = length(y);
  
% extract X-data (SNPs) using "plink"
  s_extract = ['plink --bfile ' Fr folder(6:end) ' --extract ' Fr ...
               'snp.txt' ' --make-bed --out ' Fr 'snp'];
  s_recode = ['plink --bfile ' Fr 'snp --recode --out '  Fr 'snp'];
  n_Par = ceil(n_vars/n_PerP);
  N_Par = [n_PerP*ones(1,n_Par-1) n_vars-n_PerP*(n_Par-1)];
  for k = 1:n_Par
      [k n_Par]
      
    % delete all snp-files
      delete([Fr 'snp.*'])
      snp_files = dir([Fr 'snp.*']);
      if (~isempty(snp_files))
          error('!!!')
      end
      
    % define formation of "textscan"
      s_format = '%s%s%d%d%d%d';
      for i = 1:N_Par(k)
          s_format = [s_format '%s%s'];
      end
      
    % create "snp.txt" and write j-th SNP into it
      fid = fopen([Fr 'snp.txt'],'w');
      j_0 = sum(N_Par(1:k-1));
      for j = 1:N_Par(k)
          fprintf(fid,'%s\n',SNP{j_0+j});
      end
      fclose(fid);
      
    % plink: extract SNP info about all samples and recode it
      dos(s_extract)
      dos(s_recode)
      
    % extract genotypes and transfer them into numeric values (1-to-1)
      fid = fopen([Fr 'snp.ped'],'r');
      data = textscan(fid,s_format);
      fclose(fid);
      for j = 1:N_Par(k)
          snp_1 = char(data{6+2*j-1});
          snp_2 = char(data{6+2*j});
          snp = cellstr([snp_1 snp_2]);
          x = nan(nc,1);
          for i = 1:nc
              x(i) = snp2num(snp{i});
          end
          save(sprintf('%sx_%07.0f.mat',To,j_0+j),'x')
      end
      clearvars data
  end
  
% end
  time = datestr(now,'yyyy-mm-dd HH:MM:SS');
  fprintf('e n d: %s\n',time)
  