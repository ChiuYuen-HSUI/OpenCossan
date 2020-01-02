Xds = cat(1,Tinput.eigenfrequencies);
Cnames = fieldnames(Toutput);
Cdata = num2cell(Xds.Mdata);
Toutput = cell2struct(Cdata,Cnames,2);