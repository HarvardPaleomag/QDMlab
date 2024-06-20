function Bfile2MagIC(Bdata,fID_out,meas_label,pixel_size,header)
arguments
    Bdata;
    fID_out;
    meas_label;
    
    pixel_size = 4.7*10^(-6); % standard 4x4 binning
    header = false;
end

% print header first if supplied
if header
    fprintf(header);
end

% convert t
Bn = 0;
for Bi = 1:size(Bdata,1)
    for Bj = 1:size(Bdata,2)
        measLine = ['\r\n',meas_label,'-',num2str(Bn),'\t',...
                num2str(Bdata(Bi,Bj),'%0.8e'),'\t',...
                num2str((Bi-1)*pixel_size,'%0.6e'),'\t',...
                num2str((Bj-1)*pixel_size,'%0.6e')];
        fprintf(fID_out, measLine);
        Bn = Bn + 1;
    end
end

fclose(fID_out);

end