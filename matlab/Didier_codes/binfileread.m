function v = binfileread(prefix,base,type,num)
%base is the filename suffix, directory included if necessary
%prefix is the filname prefix, e.g., 'y'
%type is the data type.  usually 'int32' or 'float32'
%    this won't work on the 'info' file
%num 
%
%little endian (or at least consistency) is assumed
%D.Vibert & R.Frazin 7/9/2011


fid = fopen([prefix,base],'rb');
if (fid < 0)
    disp('Bad filename:');
    disp([prefix,base]);
    v = -1;
    return;
end

if (nargin == 3)
    v = fread(fid,inf,type);
else 
    v = fread(fid,num,type);
end

return;

