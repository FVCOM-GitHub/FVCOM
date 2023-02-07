%
clear all
close all

% read data
[xseg1,xseg2,yseg1,yseg2,seg,type] = textread('check.dat','%f %f %f %f %d %d');
narrows = prod(size(seg));
seg = seg+1;
nsegs = max(seg);

% put x's on midpoints of open boundary edges
ii = 0;
for i=1:narrows
  if(type(i)==1);
    ii = ii + 1;
    xmid(ii) = .5*(xseg1(i)+xseg2(i));
    ymid(ii) = .5*(yseg1(i)+yseg2(i));
  end;
end;
nobc = prod(size(xmid));


% assign colors
cmap = colormap;
dims = size(cmap);
ncolors = dims(1);
scale = ncolors/nsegs;

icnt = zeros(nsegs,1);
for i=1:nsegs   
  k = ceil(i*scale);
  for j=1:narrows
    if(seg(j)==i)
      colors(j,1:3) = cmap(k,1:3);
    end;
  end;
end;

% plot arrows
scale = 1;
arrow([xseg1,yseg1],[xseg2,yseg2], ...
      'Length',scale)  ; hold on;


% plot x on open boundary segments
if(nobc > 0); 
  plot(xmid,ymid,'b+'); hold on;
end;
