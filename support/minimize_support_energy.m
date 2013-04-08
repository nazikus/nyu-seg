% Minimizes the energy of the IP formulation for support and
% structure-class inference.
%
% Args:
%   unaryCostSupport - the costs of assigning each region a support label
%                      from the following label set: {1..R, h, \emptyset},
%                      an R x (2*(R+1)+1) matrix.
%   unaryCostStructure - costs of assigning each region a structure class
%                        label from the label set {floor, structure,
%                        furniture, prop}, an Rx4 matrix.
%   TC_V - the transition costs from structure-class to structure-class
%          when a region is suported from below by another region: a 4x4
%          matrix.
%   TC_H - the transition costs from structure-class to structure-class
%          when a region is supported from behind by another region: a 4x4
%          matrix.
%   minHeights - the minimum heights (absolute) of each region.
%   maxHeights - the maximum heights (absolute) of each region.
%   minHorzDists - the minimum horizontal distances between each region.
%   params - the parameters struct. See Params.m
%
% Returns:
%   S - the support assignments, an R x (2*(R+1)+1) matrix.
%   M - the structure class assignments, an Rx4 matrix.
%   E - the energy of the solution.
function [S, M, E] = minimize_support_energy(unaryCostSupport, ...
    unaryCostStructure, TC_V, TC_H, minHeights, maxHeights, minHorzDists, ...
    params)

  Consts;
  
  TC_V = double(TC_V);
  TC_H = double(TC_H);
  minHeights = double(minHeights);
  maxHeights = double(maxHeights);
  minHorzDists = double(minHorzDists);
  
  %% Formalize the energy as costs and (in)equalities.
  R = size(unaryCostStructure,1);
  
  switch params.regionSrc
    case consts.REGION_SRC_BOTTOM_UP
      coeffs = params.support.coeffsSeg;
    case consts.REGION_SRC_LABELS
      coeffs = params.support.coeffsGt;
  end
  
  [f, A, b, Aeq, beq] = ...
      mex_setup_lp(unaryCostSupport', unaryCostStructure', ...
        TC_V', TC_H', minHeights, maxHeights, minHorzDists', coeffs);
      
  % For each region, consider only the possible vertical and horizontal
  % supports.
  validS = true(size(unaryCostSupport));
  validM = true(size(unaryCostStructure));
  validUnknownsWV = true(R, R * 16 + 4,1);
  validUnknownsWH = true(R, R * 16 + 4,1);
  
  validS = validS';
  validM = validM';
  
  validUnknownsWV = validUnknownsWV';
  validUnknownsWH = validUnknownsWH';
  
  A = A';
  Aeq = Aeq';

  %% Run solver.
  if consts.useGurobi
  
    model.A = [A; Aeq];

    model.obj = f(:)';
    model.rhs = [b; beq];

    model.sense = [repmat('<', [1 numel(b)]) repmat('=', [1, numel(beq)])];

    if params.support.infMethod == consts.SUP_INF_LP
      model.vtype = 'C';
    elseif params.support.infMethod == consts.SUP_INF_IP
      model.vtype = 'B';
    else
      error('Unknown support inference method');
    end

    model.modelsense = 'min';

    params2.LogToConsole = 0;
    params2.outputflag = 0;

    result = gurobi(model, params2);
    x = result.x;
    E = result.objval;
  else
    if params.support.infMethod == consts.SUP_INF_IP
      error('Not supported, switch to Gurobi');
    end
    
    options = optimset;
    options.LargeScale = 'off';
    options.Simplex = 'on';
    [x, E] = linprog(f, A', b, Aeq', beq, [], [], [], options);
  end
  
  [S, M] = recover_results_from_unknowns(x, unaryCostSupport, unaryCostStructure, validS, validM, ...
      validUnknownsWV, validUnknownsWH);
end

function [S, M] = recover_results_from_unknowns(x, D, T, validS, ...
    validM, validUnknownsWV, validUnknownsWH)
  
  numUnknownsS = nnz(validS);
  numUnknownsM = nnz(validM);
  numUnknownsWV = nnz(validUnknownsWV);
  numUnknownsWH = nnz(validUnknownsWH);
  totalUnknowns = numUnknownsS + numUnknownsM + numUnknownsWV + numUnknownsWH;
  
  if numel(x) ~= totalUnknowns
    error('Num Unknowns doesnt match');
  end
  
  offsetS = 0;
  offsetM = offsetS + numUnknownsS;
  offsetWV = offsetM + numUnknownsM;
  
  S = zeros(size(D'));
  S(validS) = x(offsetS+1 : offsetM);
  S = S';
  
  M = zeros(size(T'));
  M(validM) = x(offsetM+1 : offsetWV);
  M = M';
  
  % As an approximation, to obtain the final values, just take the max in
  % each row.
  R = size(M,1);
  for ii = 1 : R
    [~, ind] = max(S(ii,:));
    S(ii,:) = 0;
    S(ii,ind) = 1;
    
    [~, ind] = max(M(ii,:));
    M(ii,:) = 0;
    M(ii,ind) = 1;
  end
end
