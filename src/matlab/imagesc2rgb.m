% IMAGESC2RGB Scale data and return as RGB image
% [I,range] = imagesc2rgb(g,range,colormap,nancolor) applies color map 
% 'colormap' to scalar image g, clamping values below range(1) and above
% range(2) (just as imagesc) and returns an RGB image and the actual value
% range used. NaNs in g are represented by the color given by nancolor
% (default: [0,0,0]);
% 
% 'colormap' is set to 'jet()' by default
%
% if 'range' is empty, or contains 'nan', the minimum or maximum values of
% g are applied
%
% example:
% I = imagesc2rgb(g,[0,nan],jet(1000),[1,1,1]);
%

function [I,range] = imagesc2rgb(g,range,colormap,nanColor)
  if nargin < 2 || isempty(range)
    range = [nan,nan];
  end
  
  if nargin < 3
    colormap = jet();
  end
  
  if nargin < 4
    nanColor = [1,1,1];
  end  
  
  if isnan(range(1))
    range(1) = min(g(:));
  end
  if isnan(range(2))
    range(2) = max(g(:));
  end
  
  
  n = size(colormap,1);
  
  invalid = isnan(g);
  
  % scale range to [0,1]
  g = g - range(1);
  g = g / ( range(2) - range(1) );
  % scale range to [1,n]
  g = g * (n-1) + 1;
  % round to integer
  g = round(g);
  % clamp to [1,n]
  g = max(1, min(g,n));
  % convert to RGB using the colormap
  I = ind2rgb(g,colormap);
  
  for i = 1:3
    Ii = I(:,:,i);
    Ii(invalid) = nanColor(i);
    I(:,:,i) = Ii;
  end
end
  