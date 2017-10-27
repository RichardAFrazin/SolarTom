function status = sparwrite(mat,directory,fname_ext)
%function status = sparwrite(mat,directory,fname_ext)
%  this takes a matlab sparse matrix and writes it
%   to my sparse format in the directory (inclue final slash)
%   and with the filename_extension

[nr,nc] = size(mat);

[nmat,imat,vmat] = sparfix(mat);

whos nmat imat vmat

y = zeros(nr,1);

fidn = fopen([directory,'n',fname_ext],'wb')
fidi = fopen([directory,'i',fname_ext],'wb')
fidv = fopen([directory,'v',fname_ext],'wb')
fidy = fopen([directory,'y',fname_ext],'wb')
fidd = fopen([directory,'delta_',fname_ext],'wb')

fwrite(fidn,nmat,'int32')
fwrite(fidi,imat,'int32')
fwrite(fidv,vmat,'float32')
fwrite(fidy,y,'float32')
fwrite(fidd,y,'float32')

fclose(fidn)
fclose(fidi)
fclose(fidv)
fclose(fidy)
fclose(fidd)

status = 0;
