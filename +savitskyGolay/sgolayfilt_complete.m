function [filtered] = sgolayfilt_complete(signal,order,framelen)


%order = 3;
%framelen = 113;


lx = length(signal);
B = savitzkyGolay([ceil(-(framelen/2)):1:floor(framelen/2)],order,0,[ceil(-(framelen/2)):1:floor(framelen/2)],[],0);

m = (framelen -1)/2;
steady = conv(signal,B(m+1,:),'same');
ybeg = B(1:m,:)*signal(1:framelen);
yend = B(framelen-m+1:framelen,:)*signal(lx-framelen+1:lx);
cmplt = steady;
cmplt(1:m) = ybeg;
cmplt(lx-m+1:lx) = yend;

filtered  = cmplt;

end