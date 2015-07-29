function y = surfid(x)
    %determine Gambit surface id from node numbering
    %called by gambit_hex2tet_with_BC
    if x=='123'
        y=1;
    elseif  x=='124'
        y=2;
    elseif  x=='134'
         y=4;
    else
        error('unkown surface %s',x)
    end
end
