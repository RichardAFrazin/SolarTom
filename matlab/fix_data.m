% fix_data.m

files = {'x_nov2_15_mk_median';
         'x_nov9_22_mk_median';
         'x_oct12_25_mk_median';
         'x_oct19_nov1_mk_median';
         'x_oct26_nov8_mk_median';
         'x_oct5_18_mk_median'};

for i=1:length(files)
  file_i = char(files(i));
  file_new_i = ['/home/butala/src/srt/bindata/', file_i, '_fixed'];
  %file_new_i = ['/tmp/', file_i, '_fixed'];
  
  fix_r(file_i, file_new_i);
end
