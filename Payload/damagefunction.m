% improvements drag
function [v_i, q_cr, P_f] = damage_analysis(C_D, k, ro_exp, ro_case, type, V_exp, V_case, E, B, D, d, ro_air, E_cr, R, A_T)
    % Calculating mass of case and explosive, kg
    m_exp = ro_exp*V_exp/1000
    m_case = ro_case*V_case/1000;
    
    % Initial fragment velocity for cylinder and sphere, m/s
    if type == 1
        v_i = E*((m_case/m_exp)+0.5)^(-0.5);
    elseif type == 2
        v_i = E*((m_case/m_exp)+0.6)^(-0.5);
    end
    
    % Mott distribution calculations, inches and lbs
    t = (D-d)/2;
    Mk = B*(t^(5/16))*(d^(1/3))*(1+t/d);
    N_t = (m_case*2.205)/2/Mk^2;
    %M_av = 2*Mk^2;

    % Final velocity of fragments, m and m/s
    L = (2*(k^(2/3))/C_D/ro_air)*100; %m
    v_f = v_i*exp(-R/L); %m/s

    
    %kg
    M_cr_low = 2*E_cr/(v_f.^2);
    M_cr_high = (2*E_cr/9.81/L)^(0.75);

    if M_cr_low < M_cr_high
        M_cr = M_cr_low;
    else
        M_cr = M_cr_high;
    end
    
    %pieces/m^2
    Q0 = N_t/4/pi;
    q_cr = Q0/(R^2)*exp(-(2*M_cr/m_case)^0.5); % Density of critical fragments, piece/m^2
    %Nm = m_case/2*(Mk^2)*exp(-(M_cr^0.5)/Mk);

    P_f = 1-exp(-q_cr*A_T);

end