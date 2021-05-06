%% Hypergeometric Analytic

% From Karch (2020), Collabra: Psychology

function r = hypergeom_analytic(c,z)

if(z == 0)
    r = 1;
elseif(z == 1)
    r = (c - 1) / (c - 2);
elseif(mod(c,1) == 0)
    shared = (z - 1) / z;
    the_sum = 0;
    for i = 2:c-1
        the_sum = the_sum + shared^i /(c - i);
    end
    r = ((c - 1)*z*(z - 1)^(-2))*(the_sum - shared^c * log(1 - z));
else
    cur_r = 1/(1 - z)*(1 + (sqrt(z)*asin(sqrt(z))))/sqrt(1 - z);
    for i = 1.5:c
        cur_c = i - 1;
        cur_r = (cur_c - cur_c*(1 - z)*cur_r)/(z*(cur_c - 1));
    end
    r = cur_r;
end

end

