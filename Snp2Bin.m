  function [z,J] = Snp2Bin(x,y)
  
% crosstab of x and y (note: x == grp2idx(X) and |y| == 2)
  n_x = max(x);
  n_y = 2;
  T = nan(n_x,n_y);
  for i = 1:n_x
      for j = 1:n_y
          T(i,j) = sum(x == i & y == j);
      end
  end
% crosstab of x and y: far faster than "T = crosstab(x,y)"
  
% transform x to a binary variable
  [~,I] = max(T,[],2);
  J = find(I == 2)';
  z = false(length(x),1);
  for j = J%(:)' % each x(j) is transformed to 2
      z(x == j) = true;
  end
  z = z+1;
  