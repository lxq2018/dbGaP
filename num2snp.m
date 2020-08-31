  function snp = num2snp(num)
  switch num
      case { 0, '0'}, snp = '--';
      case { 1, '1'}, snp = 'AA';
      case { 2, '2'}, snp = 'AT';
      case { 3, '3'}, snp = 'AG';
      case { 4, '4'}, snp = 'AC';
      case { 5, '5'}, snp = 'TA';
      case { 6, '6'}, snp = 'TT';
      case { 7, '7'}, snp = 'TG';
      case { 8, '8'}, snp = 'TC';
      case { 9, '9'}, snp = 'GA';
      case {10,'10'}, snp = 'GT';
      case {11,'11'}, snp = 'GG';
      case {12,'12'}, snp = 'GC';
      case {13,'13'}, snp = 'CA';
      case {14,'14'}, snp = 'CT';
      case {15,'15'}, snp = 'CG';
      case {16,'16'}, snp = 'CC';
      case {17,'17'}, snp = 'DD';
      case {18,'18'}, snp = 'DI';
      case {19,'19'}, snp = 'ID';
      case {20,'20'}, snp = 'II';
      otherwise, error('!!!')
  end
  