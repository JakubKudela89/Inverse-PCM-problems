function [f_cost] = GACR22_PRES_model_5_v2(vec)
        %reseni nezname
        %vhodny rozsah parametru:
        %parametr 1 - c0_l - [1500, 5000] - specific heat capacity solid/liquid
        %parametr 2 - c1 - [10000, 100000] - ceff peak hight 
        %parametr 3 - Tpch - [10, 60] - peak phase change temperature
        %parametr 4 - sigma_1 - [0.5, 50] - solid ceff Gaussian curve sharpness
        %parametr 5 - sigma_2 - [0.5, 50] - liquid ceff Gaussian curve sharpness
        %parametr 6 - rho - [200, 2000] - density
        %parametr 7 - lambda - [0.1, 10] - thermal conductivity
        %parametr 8 - c0_u - [1500, 5000] - specific heat capacity solid/liquid

        %nektere hustoty s tepelnou vodivosti (spolu s delta T modelu)
        %muzou porusit stabilitu modelu
        
%         fw = waitbar(0, 'Starting');
        T_PCM_initial = 20;
        lambda_PCM = vec(7);%0.2
        rho_PCM = vec(6);% 800
        Tpch_PCM_heat = vec(3);%35
        c0_PCM_heat = vec(1); %2000
        c0_PCM_heat_up = vec(8);
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
        t_presim = [25	25	25	25	25	37.6311240749984	38.2151244452484	38.9029187882152	39.6313187844874	40.3311440268610	40.9758058592285	41.5591859818784	42.0810969649117	42.5444321399079	42.9537069609371	43.3143648117110	43.6319258078198	43.9118125017898	44.1591893702348	44.3787450239521	44.5746365052654	44.7505136334576	44.9095313558195	45.0542712711407	45.1869280199113	45.3092867256024	45.4228560166923	45.5288884105650	45.6284069108811	45.7222797015827	45.8112111916228	45.8958014052539	45.9765588854231	46.0539073547499	46.1282110035025	46.1997728831866	46.2688723967165	46.3357286497430	46.4005486150620	46.4634974205339	46.5247420398091	46.5844234746101	46.6426515520706	46.6995495570928	46.7552195144840	46.8097449445343	46.8632189140980	46.9157156237270	46.9673090611600	47.0180830616023	47.0680914174008	47.1174092260091	47.1661094495362	47.2142532162417	47.2619259480205	47.3091979894852	47.3561713684675	47.4029352185081	47.4496353097162	47.4964564639081	47.5437367845580	47.5921195013686	47.6435220843632	47.7071343142928	47.7849908132270	47.8611543812004	47.9341601413924	48.0046017806063	48.0752161011784	48.1523812783008	48.2350919719615	48.3239068078126	48.4152253180237	48.5073504710674	48.6039247598396	48.6969774786383	48.7871708303102	48.8795102358146	48.9827431853112	49.0929217635965	49.1997408898580	49.3054071603517	49.4139843735222	49.5191259769389	49.6221845235192	49.7308478963089	49.8521963366988	49.9726757115923	50.0877680435190	50.2024483154666	50.3190859510073	50.4336983096504	50.5494559196322	50.6745068947619	50.8037442638357	50.9286075295722	51.0492924898755	51.1691110144011	51.2883013222410	51.4098568631011	51.5435995622108	51.6845984783210	51.8209731316039	51.9533918737488	52.0846135277289	52.2187463398413	52.3520691899896	52.4855452995138	52.6212358353889	52.7574814116264	52.8869026523025	53.0093725194245	53.1308116091691	53.2520657094531	53.3808687742454	53.5040192376187	53.6170563422209	53.7239937718447	53.8296573280935	53.9295612291621	54.0286640382945	54.1283136848000	54.2200999552640	54.3025761235772	54.3761964423590	54.4420145803094	54.5010431411053	54.5538389666514	54.6010542982123	54.6432783609631	54.6810382323398	54.7148053899885	54.7450015705507	54.7720040117415	54.7961501425674	54.8177417797193	54.8370488820802	54.8543129098276	54.8697498297298	54.8835528038609	54.8958945950485	54.9069297188657	54.9167963688443	54.9256181387827	54.9335055635091	54.9405574972178	54.9468623464817	54.9524991732466	54.9575386815046	54.9620440998978	37.1388186243532	35.5037901614335	34.2090801321557	33.4049622115473	32.9773284680693	32.7543038429691	32.6037335115242	32.4767890525751	32.3618859816567	32.2543582868566	32.1516650503811	32.0522072857334	31.9548818499381	31.8588917780667	31.7636013507970	31.6685117269860	31.5731480399208	31.4771172625921	31.3800052305548	31.2813744342721	31.1807412551826	31.0775332227889	30.9710478538652	30.8603470687975	30.7441880715504	30.6208881494471	30.4884642989858	30.3448610869379	30.1885635957211	30.0191088840704	29.8365278442476	29.6405237099656	29.4312115480699	29.2094411103565	28.9765624279297	28.7363189854309	28.4907848502996	28.2358032931802	27.9702861143501	27.6970909126585	27.4181892309872	27.1335530052022	26.8425622289980	26.5439986727994	26.2401770428985	25.9332822288341	25.6223902237590	25.3069847176340	24.9867114104179	24.6644641817762	24.3418680972759	24.0180829229565	23.6937059844580	23.3696140477079	23.0486791178539	22.7318950509075	22.4186337642365	22.1096596795727	21.8075926856940	21.5154591455653	21.2327768697352	20.9590790122493	20.6949185546378	20.4412446304041	20.1997875318624	19.9700659159631	19.7509132072398	19.5428390726824	19.3469699729861	19.1621628175167	18.9869528992417	18.8209779263998	18.6640132114587	18.5162128122544	18.3762608328942	18.2425787260114	18.1144254022576	17.9916217017851	17.8742257316715	17.7621896187815	17.6555719409513	17.5543202295753	17.4582180826133	17.3670601978979	17.2806231589542	17.1986773975821	17.1209826560738	17.0472507781616	16.9772220207298	16.9106359244749	16.8472153301162	16.7867000997804	16.7288477366516	16.6734439256620	16.6203233724572	16.5693043433481	16.5202203041359	16.4729329914692	16.4273454527644	16.3833489314743	16.3408351135237	16.2997270108583	16.2599613174106	16.2214768048133	16.1842170634474	16.1481096712181	16.1130937698082	16.0791314382178	16.0461961808751	16.0142681781959	15.9832703649584	15.9531713702803	15.9239614759995	15.8956052491998	15.8680744909057	15.8413303305681	15.8153541912647	15.7901333162680	15.7656482002190	15.7418643672177	15.7187493810286	15.6962967344103	15.6744985813988	15.6533371485845	15.6327913847863	15.6128330152438	15.5934485252960	15.5746215024623	15.5563457851623	15.5386117954183	15.5213940028249	15.5046758543614	15.4884486351194	15.4726913388773	15.4573976705870	15.4425599832933	15.4281625238492	15.4141927502098	15.4006413269624	15.3874943624020	15.3747445452679	15.3623834575938	15.3503994705554	15.3387814924375	15.3275134142486	15.3165847686827	15.3059946584565	15.2957336630430	15.2857920681833	15.2761618853019];
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
                           (2 * dt / (dx * rho_PCM * ceff(T(1, j, k), Tpch_PCM_heat, c0_PCM_heat, c0_PCM_heat_up, c1_PCM_heat, sigma_PCM_heat_1,sigma_PCM_heat_2,charging))) * ( ...
                           qu(j,k) + lambda_PCM * ((T(2, j, k) - T(1, j, k)) / dx));

                        for i = 2:N_nodes_in_panel_INTEGER - 1
                            TT(i, j, k) = T(i, j, k) + ...
                                (dt / (dx * rho_PCM * ceff(T(i, j, k), Tpch_PCM_heat, c0_PCM_heat, c0_PCM_heat_up, c1_PCM_heat, sigma_PCM_heat_1,sigma_PCM_heat_2,charging))) * ( ...
                                 lambda_PCM * ((T(i - 1, j, k) - T(i, j, k)) / dx) + ...
                                 lambda_PCM * ((T(i + 1, j, k) - T(i, j, k)) / dx));
                        end

                        TT(N_nodes_in_panel_INTEGER, j, k) = T(N_nodes_in_panel_INTEGER, j, k) + ...
                           (2 * dt / (dx * rho_PCM * ceff(T(N_nodes_in_panel_INTEGER, j, k), Tpch_PCM_heat, c0_PCM_heat, c0_PCM_heat_up, c1_PCM_heat, sigma_PCM_heat_1,sigma_PCM_heat_2,charging))) * ( ...
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


function [ r ] = ceff(T, Tpch, c0_l, c0_u, c1, sigma1, sigma2,char)  
    r = 1; %r = [];
    if (T<Tpch)
        r =  c0_l + (c1-c0_l) .* exp(- ((T - Tpch) .* (T - Tpch)) ./ sigma1); 
    end
    
    if (T>=Tpch)
        r =  c0_u + (c1-c0_u) .* exp(- ((T - Tpch) .* (T - Tpch)) ./ sigma2); 
    end
end
