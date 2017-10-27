% test.m

list = dir('~/src/srt/eit_data/fits_cp/*.mat');
file_list = {list.name};

cp_min = Inf;
cp_max = -Inf;
for i=1:length(file_list)
  fname_i = ['~/src/srt/eit_data/fits_cp/', char(file_list(i))];
  
  load(fname_i);
  
  cp = circ_proj(:);
  
  if (min(cp) < cp_min)
    cp_min = min(cp);
  end
  
  if (max(cp) > cp_max)
    cp_max = max(cp);
  end
end