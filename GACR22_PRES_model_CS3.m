function [f_cost] = GACR22_PRES_model_3(vec)
        %GACR22_PRES_model_3([3000, 50000, 40, 5, 20, 800, 0.17]) - je reseni s f_cost = 0
        %vhodny rozsah parametru:
        %parametr 1 - c0 - [1500, 5000] - specific heat capacity solid/liquid
        %parametr 2 - c1 - [10000, 100000] - ceff peak hight
        %parametr 3 - Tpch - [10, 60] - peak phase change temperature
        %parametr 4 - sigma_1 - [0.5, 50] - solid ceff Gaussian curve sharpness
        %parametr 5 - sigma_2 - [0.5, 50] - liquid ceff Gaussian curve sharpness
        %parametr 6 - rho - [200, 2000] - density
        %parametr 7 - lambda - [0.1, 10] - thermal conductivity
        %nektere hustoty s tepelnou vodivosti (spolu s delta T modelu)
        %muzou porusit stabilitu modelu
        
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

        %load GACR22_PRES_presim.mat
        t_presim = [25	25	25	25	25	37.8550598729511	38.6263535941795	39.3665442419980	40.0768631836747	40.7575524221989	41.4003893442290	41.9836178307395	42.4939138277621	42.9299083030761	43.2957364667931	43.5973893892192	43.8424954895664	44.0402058751912	44.1999075470499	44.3303120215946	44.4389223501812	44.5316227803837	44.6126260511535	44.6848201264557	44.7501896144507	44.8101315482541	44.8656589003571	44.9175274546795	44.9663157036010	45.0124764229011	45.0563709191053	45.0982924081897	45.1384823882292	45.1771423739909	45.2144424839835	45.2505278425125	45.2855233003656	45.3195362280114	45.3526605954650	45.3849796394642	45.4165670696476	45.4474884028160	45.4778016730556	45.5075587673090	45.5368069531076	45.5655891266583	45.5939438986837	45.6219063085194	45.6495080281142	45.6767779913788	45.7037430544936	45.7304279822998	45.7568553969491	45.7830461197989	45.8090193476264	45.8347929020427	45.8603835054444	45.8858068243473	45.9110775014303	45.9362091509948	45.9612145793902	45.9861058276513	46.0108943175620	46.0355909520964	46.0602060715494	46.0847495416862	46.1092307344282	46.1336585480589	46.1580415747487	46.1823881208510	46.2067062327098	46.2310036689911	46.2552879521374	46.2795664900067	46.3038465432172	46.3281352031940	46.3524393942581	46.3767659064048	46.4011215447578	46.4255131447968	46.4499475648907	46.4744316455364	46.4989722173975	46.5235762148123	46.5482507313814	46.5730029985347	46.5978403753559	46.6227704031007	46.6478008766946	46.6729398746665	46.6981957662778	46.7235772568451	46.7490934636237	46.7747539606255	46.8005687895425	46.8265485188974	46.8527043526624	46.8790482214808	46.9055928657006	46.9323519121937	46.9593399874037	46.9865728584612	47.0140675900468	47.0418427162744	47.0699184397524	47.0983168683069	47.1270622923411	47.1561815070210	47.1857041858024	47.2156633122780	47.2460956759788	47.2770424337425	47.3085497296533	47.3406693497828	47.3734593569947	47.4069845963632	47.4413168694570	47.4765344291114	47.5127202345873	47.5499581526246	47.5883261103963	47.6278854115557	47.6686665611712	47.7106544825431	47.7537793422443	47.7979205461599	47.8429272660251	47.8886498237899	47.9349695961252	47.9818169842707	48.0291747282850	48.0770700945571	48.1255609766767	48.1747195825139	48.2246159770317	48.2753035456780	48.3268088663607	48.3791280188431	48.4322292976778	48.4860600993147	48.5405554894607	48.5956470769162	48.6512710620234	48.7073738196787	48.7639143462924	48.8208652562688	48.8782148129428	48.9359702139294	48.9941592383621	49.0528265850607	32.7235069915864	32.4375416773776	32.2410841105743	32.0958439635617	31.9783782568742	31.8775445536348	31.7877330732558	31.7058062745936	31.6298445619991	31.5585864926494	31.4911547190028	31.4269107383303	31.3653725989089	31.3061655618942	31.2489910554716	31.1936063104958	31.1398105159486	31.0874351117582	31.0363367983953	30.9863923848622	30.9374949140246	30.8895506963801	30.8424770033258	30.7962002480038	30.7506545325088	30.7057804743812	30.6615242487861	30.6178367992101	30.5746731812150	30.5319920122491	30.4897550067204	30.4479265801439	30.4064735089536	30.3653625665505	30.3245580636931	30.2840280425326	30.2437420716511	30.2036697910839	30.1637782889540	30.1240306883238	30.0843900526822	30.0448145431903	30.0052566237427	29.9656646503799	29.9259839485110	29.8861552755587	29.8461096484472	29.8057656036056	29.7650279508236	29.7237843946598	29.6818984880358	29.6391993215451	29.5954656387145	29.5503979043675	29.5035644606885	29.4542873741147	29.4013686975101	29.3423357111757	29.2714461546355	29.1821020490248	29.0844776926608	28.9862784007832	28.8867724087898	28.7830210860626	28.6696952824825	28.5444698773853	28.4082578992941	28.2622027323768	28.1075186790014	27.9462475271325	27.7779361031536	27.5960666919588	27.3983115593930	27.1875195791078	26.9663057992482	26.7327257045341	26.4823812481013	26.2165546137495	25.9380761396338	25.6450353660182	25.3356479194459	25.0113938356868	24.6717455457628	24.3177479753187	23.9539904243824	23.5833169838591	23.2073705811693	22.8300909809957	22.4541901346403	22.0824730679375	21.7184658685017	21.3651115980984	21.0260743481712	20.7036688877292	20.3978664698606	20.1079501465446	19.8331385302950	19.5726771875875	19.3258464408900	19.0919606029791	18.8703665636462	18.6604423807855	18.4615959197261	18.2732635424413	18.0949088456909	17.9260214469681	17.7661158170798	17.6147301581586	17.4714253258804	17.3357837946502	17.2074086645046	17.0859227084775	16.9709674591769	16.8622023333244	16.7593037930168	16.6619645424840	16.5698927591299	16.4828113576609	16.4004572861279	16.3225808527261	16.2489450822235	16.1793251009105	16.1135075489908	16.0512900193581	15.9924805217321	15.9368969711543	15.8843666998702	15.8347259916573	15.7878196376810	15.7435005129931	15.7016291728141	15.6620734677682	15.6247081772688	15.5894146602788	15.5560805226974	15.5245993006516	15.4948701589967	15.4667976043563	15.4402912120556	15.4152653663292	15.3916390132045	15.3693354254895	15.3482819793125	15.3284099416859	15.3096542685863	15.2919534130637	15.2752491429139	15.2594863674667	15.2446129730611	15.2305796667969];
        %[T_amb, T_air_in, T_air_out_exp, size_min] = read_data_ENB();

        %nn = size_min;
        T_ambient_air_VEC = 25*ones(300);
        T_air_in_VEC = zeros(1,300);
        T_air_in_VEC(1:150) = 55;
        T_air_in_VEC(1:5) = 25;
        T_air_in_VEC(151:300) = 15;

        T_PCM_initial = 25; 

        x_var = 0.0293; 
        xxx = 2;
        base = 0.05;
        vaha = zeros(1,20);
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
        m_air_in = zeros(1,20);
        for k = 1:20
            m_air_in(k) = vaha(k)*mdot;
        end

        dt_max = delka/(m_air_in(10)/(rho_vz*0.009));
        v = zeros(1,20); Re = zeros(1,20); Nu = zeros(1,20); h_VEC = zeros(1,20); length = zeros(1,20); w2 = zeros(1,20); w1 = zeros(1,20); 
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

        T = zeros(N_nodes_in_panel_INTEGER,M_panel_sections_along_flow_direction_INTEGER,20);
        TT =  zeros(N_nodes_in_panel_INTEGER,M_panel_sections_along_flow_direction_INTEGER,20);
        for i = 1:N_nodes_in_panel_INTEGER
            for j = 1:M_panel_sections_along_flow_direction_INTEGER
                for k=1:20
                    T(i, j, k) = T_PCM_initial;
                    TT(i, j, k) = 0;
                end
            end
        end

        T_air = zeros(M_panel_sections_along_flow_direction_INTEGER,20); 
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
        T_air_m = zeros(ntsteps_DT,20);
        t_out = zeros(1,ntsteps_DT);
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
            Q_loss = zeros(M_panel_sections_along_flow_direction_INTEGER,20);
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

                 T = TT;
%                 for i = 1:N_nodes_in_panel_INTEGER
%                     for j = 1:M_panel_sections_along_flow_direction_INTEGER
%                         for k=1:20
%                             T(i, j, k) = TT(i, j, k);
%                         end
%                     end
%                 end
                

            end %} %// end for n
            
            T_air_out(N) = sum(vaha .* T_air(M_panel_sections_along_flow_direction_INTEGER,:));

            T_air_m(N,:) = T_air(M_panel_sections_along_flow_direction_INTEGER,:);

            f_cost = f_cost + (T_air_out(N)-t_presim(N))^2;

            t_out(N) = T_air_out(N);


        end % end for N
        %         figure()
%         plot(t_out)
%         hold on
%         plot(t_presim)
%         hold on
%         t_out
%         t_presim

%         close(fw)

end %// of main function


function [ r ] = ceff(T, Tpch, c0, c1, sigma1, sigma2, char)  
    r = 1; %r = [];
    if (T<Tpch)
        r =  c0 + c1 .* exp(- ((T - Tpch) .* (T - Tpch)) ./ sigma1); 
    end
    
    if (T>=Tpch)
        r =  c0 + c1 .* exp(- ((T - Tpch) .* (T - Tpch)) ./ sigma2); 
    end
end
