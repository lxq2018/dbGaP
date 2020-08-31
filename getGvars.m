  function getGvars(studyID)
  if (nargin < 1)
      studyID = 'phs000147';
  end
  switch studyID
      case 'phs000147'
           f = [cd '\data_orig\GenotypeFiles\phg000419\genotype-calls-matrixfmt\breast_v2.bim'];
           %f = [cd '\data_orig\gdata\genotype-calls-matrixfmt\nhs_subjects.bim'];
           fid = fopen(f,'r');
           data = textscan(fid,'%d%s%d%d%s%s');
           fclose(fid);
           chr = data{1};
           snp = data{2};
           pos = data{4};
           clearvars data;
           load([cd '\data_nhs_subjects\SNP.mat'])
           gvars.snp = SNP;
           clearvars SNP
           gvars.chr = nan(length(gvars.snp),1);
           gvars.pos = nan(length(gvars.snp),1);
           [isMem,idx] = ismember(gvars.snp,snp);
           gvars.chr(isMem) = chr(idx(isMem));
           gvars.pos(isMem) = pos(idx(isMem));
           save([cd '\data_nhs_subjects\gvars.mat'],'gvars');
      case 'phs000517'
           f = [cd '\data_orig\gdata\genotype-calls-matrixfmt\BreastCancer_MultiethnicCohort_c2.map'];
           fid = fopen(f,'r');
           data = textscan(fid,'%d%s%d%d%s%s');
           fclose(fid);
           chr = data{1};
           snp = data{2};
           pos = data{4};
           clearvars data;
           load([cd '\data_phs000517\SNP.mat'])
           gvars.snp = SNP;
           clearvars SNP
           gvars.chr = nan(length(gvars.snp),1);
           gvars.pos = nan(length(gvars.snp),1);
           [isMem,idx] = ismember(gvars.snp,snp);
           gvars.chr(isMem) = chr(idx(isMem)); gvars.chr(gvars.chr == 0) = nan;
           gvars.pos(isMem) = pos(idx(isMem)); gvars.pos(gvars.pos == 0) = nan;
           save([cd '\data_phs000517\gvars.mat'],'gvars');
      case 'phs000634'
           f = [cd '\data_orig\gdata\genotype-calls-matrixfmt\nslc_c1.bim'];
           fid = fopen(f,'r');
           data = textscan(fid,'%d%s%d%d%s%s');
           fclose(fid);
           chr = data{1};
           snp = data{2};
           pos = data{4};
           clearvars data;
           load([cd '\data_nslc_c1\SNP.mat'])
           gvars.snp = SNP;
           clearvars SNP
           gvars.chr = nan(length(gvars.snp),1);
           gvars.pos = nan(length(gvars.snp),1);
           [isMem,idx] = ismember(gvars.snp,snp);
           gvars.chr(isMem) = chr(idx(isMem)); gvars.chr(gvars.chr == 0) = nan;
           gvars.pos(isMem) = pos(idx(isMem)); gvars.pos(gvars.pos == 0) = nan;
           save([cd '\data_nslc_c1\gvars.mat'],'gvars');
      case 'phs000753'
           f = [cd '\data_orig\gdata\genotype-calls-matrixfmt\lung_gwas_risk_to_dbGAP.bim'];
           fid = fopen(f,'r');
           data = textscan(fid,'%d%s%d%d%s%s');
           fclose(fid);
           chr = data{1};
           snp = data{2};
           pos = data{4};
           clearvars data;
           load([cd '\data_lung_gwas_risk\SNP.mat'])
           gvars.snp = SNP;
           clearvars SNP
           gvars.chr = nan(length(gvars.snp),1);
           gvars.pos = nan(length(gvars.snp),1);
           [isMem,idx] = ismember(gvars.snp,snp);
           gvars.chr(isMem) = chr(idx(isMem)); gvars.chr(gvars.chr == 0) = nan;
           gvars.pos(isMem) = pos(idx(isMem)); gvars.pos(gvars.pos == 0) = nan;
           save([cd '\data_lung_gwas_risk\gvars.mat'],'gvars');
  end