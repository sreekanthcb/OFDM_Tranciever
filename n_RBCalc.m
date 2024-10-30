function nRB = n_RBCalc(scs,bw)
    % based on release 16 of 5GNR
    scs     = scs/1e3;
    bw      = bw/1e6;

    bw_vec  = [5 10 15 20 25 40 50 60 80 100];
    scs_vec = [15 30 60];

    if ~ismember(scs,scs_vec) || ~ismember(bw,bw_vec)
        error('I dont know how to handle this, Bad Sreekanth!!')
    end

    RB_vec  = [25 52 79 106 133 216 270 999 999 999;
               11 24 38 51 65 106 133 162 217 273;
               999 11 18 24 31 51 65 79 107 135];

    nRB     = RB_vec(scs == scs_vec,bw == bw_vec);

    if nRB == 999
        error('Wrong SCS and BW configuration')
    end
end
