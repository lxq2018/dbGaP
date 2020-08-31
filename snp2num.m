  function num = snp2num(snp)
  switch snp
      case {'--','00'}, num = 0;
      case 'AA', num = 1;
      case 'AT', num = 2;
      case 'AG', num = 3;
      case 'AC', num = 4;
      case 'TA', num = 5;
      case 'TT', num = 6;
      case 'TG', num = 7;
      case 'TC', num = 8;
      case 'GA', num = 9;
      case 'GT', num = 10;
      case 'GG', num = 11;
      case 'GC', num = 12;
      case 'CA', num = 13;
      case 'CT', num = 14;
      case 'CG', num = 15;
      case 'CC', num = 16;
      case 'DD', num = 17;
      case 'DI', num = 18;
      case 'ID', num = 19;
      case 'II', num = 20;
      otherwise, num = 21;
  end
  