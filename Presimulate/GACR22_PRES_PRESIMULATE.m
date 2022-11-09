function [t_presim] = GACR22_PRES_PRESIMULATE(vec)
        %GACR22_PRES_model([2000, 20000, 35, 23, 3, 800, 0.2]) - je reseni s f_cost = 0
        %vhodny rozsah parametru:
        %parametr 1 - 
        %parametr 2 - 
        %
        
%         fw = waitbar(0, 'Starting');
        T_PCM_initial = 20;
        lambda_PCM = vec(7);%0.2
        rho_PCM = vec(6);% 800
        Tpch_PCM_heat = vec(3);%35
        c0_PCM_heat = vec(1); %2000
        c1_PCM_heat = vec(2);% 20000
        sigma_PCM_heat_1 = vec(4);% 15 23
        sigma_PCM_heat_2 = vec(5);% 8 3
        %sigma_PCM_heat_2 = sigma_PCM_heat_1;
        Tpch_PCM_cool = Tpch_PCM_heat;
        c0_PCM_cool = c0_PCM_heat;
        c1_PCM_cool = c1_PCM_heat; %24632 a 20124
        sigma_PCM_cool_1 = sigma_PCM_heat_1; %50
        sigma_PCM_cool_2 = sigma_PCM_heat_2; %0.6071

        d_panel_thickness = 0.0039; %0.0041
        N_nodes_in_panel = 7;
        L_panel_length = 1.5;
        mdot = 0.078; %0.071086237 a taky 0.078, presim byl 0.04
        P_unit_perimeter = 2.12;
        M_panel_sections_along_flow_direction = 10;
        delka = L_panel_length/(M_panel_sections_along_flow_direction);
        %R_wall_resistance = 2;
        N_num_CSM_panels = 20;
        h_in_panel_wall_u = 6;
        h_in_panel_wall_l = 3.5;
        rho_vz = 1.178513127;
       
        cp_air = 1005;
        H_leading_edge = 0.45;
        DT = 60; % casovy krok "velke DT z TRNSYSu"

%         load GACR22_ASTFE_presim.mat
        %[T_amb, T_air_in, T_air_out_exp, size_min] = read_data_ENB();

        %nn = size_min;
        T_ambient_air_VEC = 25*ones(300);
        T_air_in_VEC(1:150) = 55;
        T_air_in_VEC(1:5) = 25;
        T_air_in_VEC(151:300) = 15;

        T_PCM_initial = 25; 

        x_var = 0.0293; 
        xxx = 2;
        base = 0.05;
        vaha(1) = base + sin(-5/5*pi/xxx)*x_var;       
        vaha(2) = base + sin(-4/5*pi/xxx)*x_var;
        vaha(3) = base + sin(-3/5*pi/xxx)*x_var;
        vaha(4) = base + sin(-2/5*pi/xxx)*x_var;
        vaha(5) = base + sin(-1/5*pi/xxx)*x_var;
        vaha(6) = base + sin(1/5*pi/xxx)*x_var;
        vaha(7) = base + sin(2/5*pi/xxx)*x_var;
        vaha(8) = base + sin(3/5*pi/xxx)*x_var;
        vaha(9) = base + sin(4/5*pi/xxx)*x_var;
        vaha(10) = base + sin(5/5*pi/xxx)*x_var;
        vaha(11) = vaha(10);
        vaha(12) = vaha(9);
        vaha(13) = vaha(8);
        vaha(14) = vaha(7);
        vaha(15) = vaha(6);
        vaha(16) = vaha(5);
        vaha(17) = vaha(4);
        vaha(18) = vaha(3);
        vaha(19) = vaha(2);
        vaha(20) = vaha(1);
% 
%         figure
%         plot(vaha)

        a = 0.45;
        b = 0.02;
        k_vz = 0.0286;

        D_h = 4*a*b/ (2*a + 2*b);
        Pr = 0.701;
        D = sqrt(0.009/pi);

        for k = 1:20
            m_air_in(k) = vaha(k)*mdot;
        end

        dt_max = delka/(m_air_in(10)/(rho_vz*0.009));
       
        for i = 1:20
            v(i) = m_air_in(i)/(rho_vz*0.009);
            Re(i) = v(i) *D_h /(1.79 * 10^-5);
            Nu(i) = 3.66 + (0.065 * Re(i) * Pr * D/D_h)/(1 + 0.04 * (Re(i) * Pr * D/D_h)^(2/3));
            Nu(i) = 8;
            h_VEC(i) = Nu(i) * k_vz / D_h;

            length(i) = v(i) * dt_max;
            w2(i) = length(i)/delka;
            w1(i) = 1 - w2(i);
        end

        T_air_out = zeros(numel(T_air_in_VEC), 1);	

%         Q_loss = zeros(numel(T_air_in_VEC), 1);
%         %        Q_al_1 = zeros(numel(T_air_in_VEC), 1);	
%         %        Q_al_2 = zeros(numel(T_air_in_VEC), 1);	
% 
%         Q_stored = zeros(numel(T_air_in_VEC), 1);	
% 
%         T_middle = zeros(numel(T_air_in_VEC), 1);

        weight_shift = zeros(2,1);

        S = 0.020106193;
        V_vz = mdot/rho_vz;
        w = V_vz/S;

        %dt_max = delka/w
        %        shift = w*dt_max
        %        delka
        %        weight_shift(2) = shift/delka
        %        weight_shift(1) = 1 - weight_shift(2)
       
       if (weight_shift(2)>1)
           print("Moc velky shift")
           weight_shift(2) = 1;
       end
       
       if (weight_shift(1)<0)
           print("Moc maly shift")
           weight_shift(1) = 0;
       end
       
       qR = 0;
        
       K_c = 5; %11.78
       K_d = 5; %12.26
       R_c = (P_unit_perimeter*L_panel_length)/K_c;
       R_d = (P_unit_perimeter*L_panel_length)/K_d;
       R = R_c;
       C_al = 100000/M_panel_sections_along_flow_direction;%J(K
       h_al = 1; %4
       %R = ((1 / h_in_panel_wall) + R_wall_resistance + (1 / h_out))
       
       %((1 / h_in_panel_wall) + R_wall_resistance + (1 / h_out))

        N_nodes_in_panel_INTEGER = N_nodes_in_panel;
        M_panel_sections_along_flow_direction_INTEGER = M_panel_sections_along_flow_direction;
        N_num_CSM_panels_INTEGER = N_num_CSM_panels;
        L_section = L_panel_length / M_panel_sections_along_flow_direction_INTEGER;
        dx = (d_panel_thickness) / (N_nodes_in_panel_INTEGER - 1);


        for i = 1:N_nodes_in_panel_INTEGER
            for j = 1:M_panel_sections_along_flow_direction_INTEGER
                for k=1:20
                    T(i, j, k) = T_PCM_initial;
                    TT(i, j, k) = 0;
                end
            end
        end


        for j = 1:M_panel_sections_along_flow_direction_INTEGER
            for k = 1:20
                T_air(j,k) = T_PCM_initial;
            end
        end



        ntsteps_dt = ceil(DT / dt_max);
        dt = DT / ntsteps_dt;
        ntsteps_DT = numel(T_air_in_VEC);
        charging = 1;

        f_cost = 0;

        for N = 1:ntsteps_DT
%             waitbar(N/ntsteps_DT, fw, sprintf('Progress: %d %%', floor(N/ntsteps_DT*100)));
            %HEATING/COOLING SWITCH IF EXISTS
            if (N==100000) %OFF
                Tpch_PCM_heat = Tpch_PCM_cool;
                c0_PCM_heat = c0_PCM_cool;
                c1_PCM_heat = c1_PCM_cool;
                sigma_PCM_heat_1 = sigma_PCM_cool_1;
                sigma_PCM_heat_2 = sigma_PCM_cool_2;
                R=R_d;
                charging = 0;
            end    

            %HEAT LOSSES NEGLECTED
            Q_PCM_total = 0;
            Q_loss_total = 0;
            Q_al_total = 0;

            T_ambient_air = T_ambient_air_VEC(N);

            T_air(1,:) = T_air_in_VEC(N);
            qu = zeros(M_panel_sections_along_flow_direction_INTEGER,20);
            ql = zeros(M_panel_sections_along_flow_direction_INTEGER,20);
            Q_PCM = zeros(M_panel_sections_along_flow_direction_INTEGER,20);
            Q_al = zeros(M_panel_sections_along_flow_direction_INTEGER,20);

            for n = 1:ntsteps_dt

                for j = 1:M_panel_sections_along_flow_direction_INTEGER

                    for k = 1:20

                        if (k>=2 && k<20)  
                            qu(j,k) = (h_VEC(k)) * (T_air(j,k) - T(1,j,k)); 
                            ql(j,k) = (h_VEC(k)) * (T_air(j,k+1) - T(N_nodes_in_panel_INTEGER,j,k));
                            Q_PCM(j,k) = L_section * H_leading_edge * dt * (qu(j,k) + ql(j,k-1));
                        end

                        if (k==1)
                            qu(j,k) = (h_VEC(k)) * (T_air(j,k) - T(1,j,k)); 
                            ql(j,k) = (h_VEC(k)) * (T_air(j,k+1) - T(N_nodes_in_panel_INTEGER,j,k));  
                            Q_PCM(j,k) = L_section * H_leading_edge * dt * (qu(j,k));
                        end

                        if (k==20)
                            qu(j,k) = (h_VEC(k)) * (T_air(j,k) - T(1,j,k)); 
                            ql(j,k) = 0; 
                            Q_PCM(j,k) = L_section * H_leading_edge * dt * (ql(j,k));
                        end

                        Q_loss(j,k) = L_section * P_unit_perimeter * dt * (T_air(j,k) - T_ambient_air) / R / 20; 
                        %ztraty delim na 20 casti. dalo by se to pres prumer teploty, da se jeste vyladit, nahore je toho vic a vnitrni maji jen vyrezy, vrch a spodek maji vic a konec a zacatek maji vic
                        Q_loss(j,k) = 0;
                        %Q_al(j,k) =  0.01 * h_al * dt * (T_air(j,k) - T_al(j,k));
                        Q_al(j,k) = 0;

                        %// pokud jsem v sekci 2 a dale, tak si nejdrive vypocitam teplotu vzduchu v mezere

                        if (j > 1)
                            T_air(j, k) = (T_air(j - 1, k) - ((Q_PCM(j - 1, k) + Q_loss(j - 1, k) + Q_al(j - 1, k)) / (m_air_in(k) * cp_air * dt)))*w2(k) + T_air(j-1)*w1(k);
                            %T_al(j, k) = T_al(j, k) + Q_al(j,k)/(C_al/20); % delim 20 i kapacitu 
                        end

                        TT(1, j, k) = T(1, j, k) + ...
                           (2 * dt / (dx * rho_PCM * ceff(T(1, j, k), Tpch_PCM_heat, c0_PCM_heat, c1_PCM_heat, sigma_PCM_heat_1,sigma_PCM_heat_2,charging))) * ( ...
                           qu(j,k) + lambda_PCM * ((T(2, j, k) - T(1, j, k)) / dx));

                        for i = 2:N_nodes_in_panel_INTEGER - 1
                            TT(i, j, k) = T(i, j, k) + ...
                                (dt / (dx * rho_PCM * ceff(T(i, j, k), Tpch_PCM_heat, c0_PCM_heat, c1_PCM_heat, sigma_PCM_heat_1,sigma_PCM_heat_2,charging))) * ( ...
                                 lambda_PCM * ((T(i - 1, j, k) - T(i, j, k)) / dx) + ...
                                 lambda_PCM * ((T(i + 1, j, k) - T(i, j, k)) / dx));
                        end

                        TT(N_nodes_in_panel_INTEGER, j, k) = T(N_nodes_in_panel_INTEGER, j, k) + ...
                           (2 * dt / (dx * rho_PCM * ceff(T(N_nodes_in_panel_INTEGER, j, k), Tpch_PCM_heat, c0_PCM_heat, c1_PCM_heat, sigma_PCM_heat_1,sigma_PCM_heat_2,charging))) * ( ...
                           lambda_PCM * (T(N_nodes_in_panel_INTEGER - 1, j, k) - T(N_nodes_in_panel_INTEGER, j, k)) / dx + ql(j,k));
                    end
                end


                for i = 1:N_nodes_in_panel_INTEGER
                    for j = 1:M_panel_sections_along_flow_direction_INTEGER
                        for k=1:20
                            T(i, j, k) = TT(i, j, k);
                        end
                    end
                end


            end %} %// end for n

            T_air_out(N) = sum(vaha .* T_air(M_panel_sections_along_flow_direction_INTEGER,:));

            T_air_m(N,:) = T_air(M_panel_sections_along_flow_direction_INTEGER,:);

%             f_cost = f_cost + (T_air_out(N)-t_presim(N))^2;

            t_presim(N) = T_air_out(N);


        end % end for N

%         figure()
%         plot(t_out)
%         hold on
%         plot(t_presim)
%         hold on
%         t_out
%         t_presim

%         close(fw)

save GACR22_PRES_presim.mat t_presim

end %// of main function


function [ r ] = ceff(T, Tpch, c0, c1, sigma1, sigma2, char)  

    if (T<Tpch)
        r =  c0 + c1 .* exp(- ((T - Tpch) .* (T - Tpch)) ./ sigma1); 
    end
    
    if (T>=Tpch)
        r =  c0 + c1 .* exp(- ((T - Tpch) .* (T - Tpch)) ./ sigma2); 
    end
    
end
