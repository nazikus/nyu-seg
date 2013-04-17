function wseg = pb2wseg(pb, maxsp)
% pb - 4-directional gradient (4D-vector per each pixel)
% maxsp - threshold as MAXimum number of SuperPixels

nsp = Inf; % Number of SuperPixels
c = 1;     % median filter coefficient
pb = max(pb, [], 3);  % sqauasing 4D-gradeint vector to 1-D based on max-rule

% iteratively decrease number of watershed regions by enforcing median filter (c-coeff)
while nsp > maxsp
    wseg = watershed(medfilt2(pb, [c c]));
    nsp = max(wseg(:));  
    c = c + 2;
end
wseg = uint16(wseg);