function ferroorpara = ask_ferro_para()
%[ferroorpara] = ask_ferro_para()
ferroorpara=input('[ferro] or [para] map? [ferro]: ','s');

if isempty(ferroorpara)
    ferroorpara = 'ferro';
end

if ferroorpara(1:4) == 'para' || ferroorpara(1) == 'p'
    ferroorpara ='para ';
end