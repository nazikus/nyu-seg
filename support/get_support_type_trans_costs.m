% Returns the costs of transitioning between regions when a support
% relationship exists from below (TC_V) or from behind (TC_H).
%
% Args:
%   params - the parameters struct.
%
% Returns:
%   TC_V - transition costs (for support-from-below relationships)
%   TC_H - transition costs (for support-from-behind relationships).
function [TC_V, TC_H] = get_support_type_trans_costs(params)
  Consts;

  switch params.regionSrc
    case consts.REGION_SRC_LABELS

      TC_V = [0.78, 0.11, 0.11, 0.00;
              0.82, 0.08, 0.10, 0.00;
              0.93, 0.02, 0.02, 0.03;
              0.13, 0.19, 0.60, 0.08];

      TC_H = [0.00, 0.00, 0.00, 0.00;
              0.00, 0.99, 0.01, 0.00;
              0.01, 0.72, 0.26, 0.01;
              0.00, 0.89, 0.07, 0.04;];

    case consts.REGION_SRC_BOTTOM_UP
      
      TC_V = [0.93, 0.07, 0.00, 0.00;
              0.80, 0.09, 0.10, 0.01;
              0.90, 0.02, 0.05, 0.03;
              0.15, 0.18, 0.57, 0.10];

      TC_H = [0.00, 0.00, 0.00, 0.00;
              0.00, 1.00, 0.00, 0.00;
              0.00, 0.92, 0.06, 0.02;
              0.00, 0.92, 0.05, 0.03;];
      
  end
  
  TC_V = -log(TC_V + eps);
  TC_H = -log(TC_H + eps);
end
