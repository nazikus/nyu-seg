%SUBFIG   Pajdla: Opens figure in the specified part of the screen
%
%
%       function  f = subfig(rows,cols,ord,f)
%                           [left bottom height width]
%
%	rows	= number of row figures per screen
%	cols	= nimber of column figures per screen
%	ord	= order of figure (row oriented index), see SUBPLOT.
%	f	= figure handle
%
%       See also FIGURE, SUBPLOT.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%			08/03/94 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.1, (c) MathWorks
%       Last change  : $Id$
%       Status       : Ready ($Source$)  			 
%
function f = subfig(rows,cols,ord,f)
 
  if nargin < 4
    f	 = figure('Visible','off');
  end

 if nargin > 2

  screen = get(0, 'ScreenSize');
  sW     = screen(3);
  sH     = screen(4);
  fW     = sW/cols;
  fH     = sH/rows;
  i      = ceil(ord/cols);
  j      = rem(ord-1,cols);

  left   =      j * fW;
  bottom = sH - i * fH;

 else
  left   = rows(1);
  bottom = rows(2);
  fH     = rows(3);
  fW     = rows(4);
 end 

 % top pixel margin (so that window name is visible)
 margin = .04;
 marginpx = sH*margin;
 
 set(f,'Position',[left bottom-marginpx fW fH-21]);
 set(f,'Visible','on');

return
