% plink --bfile C:\dataBC\breast_v2 --recode --out C:\dataBC\breast_v2 --noweb
% plink --bfile C:\dataBC\nhs_subjects --recode --out C:\dataBC\nhs_subjects --noweb
% plink --bfile C:\data_phs000147\phs000147 --recode --out C:\data_phs000147\phs000147 --noweb
  
  clear; clc; close all; format compact; format shortg
  %f = 'data_breast_v2';
  f = 'data_nhs_subjects';
  R = draw_classificationPerformance('folder',f,'rsOnly',false,'delMax',false);
  R = draw_classificationPerformance('folder',f,'rsOnly',false,'delMax',true);
  
  return
  
  %f = 'data_phs000147';
  initialization('folder',f)
  obs2var_one2one('folder',f,'n_vars',561466)
  %obs2var_one2one('folder',f,'n_vars',546646)
  obs2var_grp2idx('folder',f)
  obs2var_snp2bin('folder',f)
  Fs = IterMMPC('folder',f,'alpha',0.1,'nPerP',10,'k_max',2);
  R = OptNBC('folder',f,'rsOnly',false,'delMax',false);
  R = OptNBC('folder',f,'rsOnly',false,'delMax',true);
  