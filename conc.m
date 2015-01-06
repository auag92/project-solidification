flist = dir('Concentration_*.dat');
N = size(flist,1);

for k = 1 : N-1

   fname =  sprintf('Concentration_%d0.dat',k);
   A= load(fname);
   h=plot(A(:,1),A(:,2));
   xlabel('meshx');
   ylabel('Concentration');
   legend(fname);
   drawnow

end
