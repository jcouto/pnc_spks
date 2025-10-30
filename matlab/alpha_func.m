function Y = alpha_function(amplitude, t_rise, t_decay, srate, N)

    t_rise = t_rise/1000.0;
    t_decay = t_decay/1000.0;
    
    fun_max  = (t_rise*t_decay/(t_decay-t_rise)) * log(t_rise-t_decay);
    normalization_factor = 1;%(exp(-fun_max/t_rise) - exp(-fun_max/t_decay))/(t_rise-t_decay);
    Y = zeros(1,N);
    for ii = 1:N
        %
        Y(ii) = amplitude*(1.0/(normalization_factor*(t_decay-t_rise))) * (exp(-((ii/srate)/t_decay)) - exp(-(ii/srate)/t_rise));
    end
