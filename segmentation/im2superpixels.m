function [bndinfo, pball] = im2superpixels(im, initseg)

  if isa(im, 'uint8')
    im = im2double(im);
  end

  % NCS: running into numerical issues here.
  tmp = imfilter(im, fspecial('gaussian', 13, 2));
  tmp(tmp > 1) = 1;
  tmp(tmp < 0) = 0;
  
  Consts; Params; global ii_;
  DEBUG_ = params.debug;
  
  %MARKER pbCGTG_nonmax.m Compute probability of boundary using CG and TG.
  pball = pbCGTG_nonmax(tmp);
  %pball(pball<0.025) = 0;
  
  %DEBUG
  if DEBUG_
      h_ = figure('Visible','off');
      imshow(max(pball,[],3));
      if params.degub_fig
        saveas(h_, sprintf('%s%06d_a_pb.fig',consts.watershedDir,ii_), 'fig');
      end
      print(h_, '-dpng', sprintf('%s%06d_a_pb.png',consts.watershedDir,ii_));
      close(h_);
  end
    
  %MARKER compute watershed regions
  wseg = pb2wseg(pball, 4000);
  wseg = double(wseg);
  
  %DEBUG
  if DEBUG_
      h_ = figure('Visible','off');
      imshow(wseg);
      if params.degub_fig
        saveas(h_, sprintf('%s%06d_b_watershed.fig',consts.watershedDir,ii_), 'fig');
      end
      print(h_, '-dpng', sprintf('%s%06d_b_watershed.png',consts.watershedDir,ii_));
      close(h_);
  end
  
  imgEdges = wseg == 0;  % edges in watershed-map are marked with 0s
  %MARKER force consistency with initseg
  if exist('initseg', 'var') && ~isempty(initseg)
      wseg = wseg + max(wseg(:))*(initseg-min(initseg(:)));
      [us, ui, uj] = unique(wseg(:));
      wseg = reshape(uj-1, size(im, 1), size(im, 2));
      wseg(imgEdges) = 0;
  end
  
  %DEBUG
  if DEBUG_
      h_ = figure('Visible','off');
      imshow(wseg,[]);
      if params.degub_fig
        saveas(h_, sprintf('%s%06d_c_watershed_consist_gs.fig',consts.watershedDir,ii_), 'fig');
      end
      print(h_, '-dpng', sprintf('%s%06d_c_watershed_consist_gs.png',consts.watershedDir,ii_));
      close(h_);
      
      h_ = figure('Visible','off');
      %imagesc(wseg);
      imshow(wseg,[],'ColorMap',colormap('Jet'));
      if params.degub_fig
        saveas(h_, sprintf('%s%06d_c_watershed_consist.fig',consts.watershedDir,ii_), 'fig');
      end
      print(h_, '-dpng', sprintf('%s%06d_c_watershed_consist.png',consts.watershedDir,ii_));
      close(h_);
  end

  
  %MARKER seg2fragments.m Finds the boundary fragments for an oversegmentation of an image
  [edges, ~, neighbors, wseg] = seg2fragments(double(wseg), im, 25, 3, ii_);
  %MARKER processBoundaryInfo.m Adds several fields to the bndinfo structure based on bndinfo.wseg
  bndinfo = processBoundaryInfo(wseg, edges, neighbors);
  
  % Just doing this to be consistent with the rest of the codebase. Better
  % to fix it here than in processBoundaryInfo that is part of Derek's
  % iccv07 code.
  bndinfo.imgRegions = wseg;
  bndinfo = rmfield(bndinfo, 'wseg');
end
