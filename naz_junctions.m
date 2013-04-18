% first, put a breakpoint at the end of im2superpixels.m
juncs = bndinfo.junctions.position;
plot(juncs(:,1),juncs(:,2), 'ob', 'MarkerSize', 3);
