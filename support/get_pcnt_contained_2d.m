% Determines whether or not the given prop object is contained within the
% given support object.
%
% Args:
%   propPoints - Nx3 point cloud.
%   supportPoints - Mx3 point cloud.
%   thresh - the threshold (0..1) for the percentage of prop points that
%            must fall inside the support object.
%
% Returns:
%   isContained - whether or not the prop object is contained by the
%                 support object.
function [pcntInside] = get_pcnt_contained_2d(propPoints, supportPoints)
  assert(size(propPoints, 2) == 2);
  assert(size(supportPoints, 2) == 2);
  propPoints = double(propPoints);
  supportPoints = double(supportPoints);
  
  % To speed up computation, take a sample of each.
  convHullInds = convhull(supportPoints);
  supportPoints = supportPoints(convHullInds(1:end-1), :);
  
  % Most points won't even overlap, so first check that ANY prop points are
  % found inside the axis aligned cuboid.
  if ~isInsideAxisAligned(propPoints, supportPoints)
    pcntInside = 0;
    return;
  end
  
  numProps = size(propPoints, 1);
  
%   d = delaunay(supportPoints(:,1), supportPoints(:,2));
  
  isInside = false(numProps, 1);
%   isInside2 = false(numProps, 1);
  
  for ii = 1 : numProps
    % To determine whether the point is inside or outside each triangle,
    % calculate the convex hull with the given prop point. Unless the point
    % is exactly inline with one of the trangular segments, we'll observe 4
    % points when the prop point is outside the triangle.
    
    K = convhulln([supportPoints; propPoints(ii,:)]);
    isInside(ii) = 1 - ismember(size(supportPoints,1)+1,K);

%     if isInside(ii)
%       sfigure(10); clf
% 
%       vis_point_cloud(supportPoints, 'r');
%       hold on;
%       vis_point_cloud(propPoints(ii,:), 'b');
%       hold off;
%       title(sprintf('isinside=%d', isInside(ii)));
%       keyboard;
%     end
    
    
%       sfigure(10); clf;
%       triplot(d, supportPoints(:,1), supportPoints(:,2));
%       hold on;
%       vis_point_cloud(triPoints);
%       vis_point_cloud(propPoints(ii,:), 'g');
%       hold off;
%       keyboard;
    
%     for jj = 1 : size(d,1)
%       triPoints = supportPoints(d(jj,:),:);
%       
% 
%       
%       ch = unique(convhull([triPoints; propPoints(ii,:)]));
%       
%       if ~any(ch == 4)
%         isInside2(ii) = 1;
%         continue;
%       end
%     end
%     
%     assert(all(isInside == isInside2));
  end
  
%   isInside = false(size(propPoints, 1), 1); 
%   
%   for ll = 1 : size(propPoints, 1)
%     propPoint = propPoints(ll, :);
%     
%     dirVecs = repmat(propPoint, [numSupport 1]) - supportPoints;
%     
%     signs = dirVecs >= 0;
% 
%     cnt = 0;
%     
%     % top-left.
%     if any(signs(:,1) & ~signs(:,2))
%       cnt = cnt + 1;
%     end
%     
%     % top-right.
%     if any(signs(:,1) & signs(:,2))
%       cnt = cnt + 1;
%     end
%     
%     % bottom-left.
%     if any(~signs(:,1) & ~signs(:,2))
%       cnt = cnt + 1;
%     end
%     
%     % bottom-right.
%     if any(~signs(:,1) & signs(:,2))
%       cnt = cnt + 1;
%     end
% 
%     if cnt >= 3
%       isInside(ll) = 1;
%     end
%   end
%   
%   % pcntBorder isnt enough. We need to make sure they overlap (birds
%   % eye view) as well. Otherwise, wall will be continuously supported
%   % by objects that occlude it.
  pcntInside = nnz(isInside) / numel(isInside);
end

function isInside = isInsideAxisAligned(propPoints, supportPoints)
  if max(propPoints(:,1)) < min(supportPoints(:,1)) || ...
      min(propPoints(:,1)) > max(supportPoints(:,1)) || ...
      max(propPoints(:,2)) < min(supportPoints(:,2)) || ...
      min(propPoints(:,2)) > max(supportPoints(:,2))
    
    isInside = 0;
  else
    isInside = 1;
  end
end
