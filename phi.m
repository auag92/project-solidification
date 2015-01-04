flist = dir('phi_*.dat');

N = size(flist,1);

for k = 1:N

   A= load(flist(k).name);
   h=plot(A(:,1),A(:,2));
   xlabel('meshx');
   ylabel('phi');
   legend(flist(k).name);
   drawnow
   fname =  sprintf('phi_%d',k);
   if mod(k,10) == 0
	saveas(h,fname,'png');
   end
end

