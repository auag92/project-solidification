flist = dir('phi_*.dat');
N = size(flist,1);

for k = 1 : N-1

   fname =  sprintf('phi_%d0000.dat',k);
   A= load(fname);
   h=plot(A(:,1),A(:,2));
   xlabel('meshx');
   ylabel('phi');
   legend(fname);
   drawnow

end
