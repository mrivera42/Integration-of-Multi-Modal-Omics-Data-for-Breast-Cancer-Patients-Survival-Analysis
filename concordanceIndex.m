function [ci] = concordanceIndex(targets, predicts)
  [~,i] = sort(targets,'descend');
  predicts = predicts(i);
  n = length(targets);
  total  = 0;
  norm_z = 0;
  for j=1:n
      for k=j:n
          if j ~= k
              h = stepFunction(predicts(j) - predicts(k));
              total = total + h;
              norm_z = norm_z + 1;
          end
      end
  end
  ci = total / norm_z;
end
