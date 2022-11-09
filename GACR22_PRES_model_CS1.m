function [f_cost] = GACR22_PRES_model(vec)
        %GACR22_PRES_model([2000, 20000, 35, 23, 3, 800, 0.2]) - je reseni s f_cost = 0
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
        t_presim = [25	25	25	25	25	37.9484606004569	38.6089553518044	39.0875722880408	39.4559968677376	39.7544964049070	40.0058962222669	40.2238674816984	40.4170808629460	40.5913206554563	40.7506145783029	40.8978702287715	41.0352515870879	41.1644116485374	41.2866417477793	41.4029706604002	41.5142323279909	41.6211133539752	41.7241870894636	41.8239419996102	41.9208357053755	42.0152785638036	42.1076916936061	42.1985355026757	42.2883528177020	42.3778376532085	42.4680080058366	42.5606286919788	42.6596981835691	42.7774402285259	42.9216473566247	43.0702080270549	43.2280213102168	43.4055533686179	43.5994574420895	43.8038677720004	44.0134535135343	44.2389699269329	44.4792791109021	44.7257616956456	44.9760383866484	45.2380337908507	45.5142129408828	45.7932777686337	46.0764523195707	46.3669592998020	46.6706486082577	46.9760430291288	47.2829300688482	47.5993354943632	47.9328442589558	48.2690865824682	48.6059770567248	48.9433112291771	49.2843779698420	49.6277637156807	49.9605087233715	50.2894225303817	50.6187626119733	50.9386453471951	51.2451715324896	51.5410116111587	51.8275654388359	52.1009796898059	52.3524182294892	52.5822681826267	52.7923353193146	52.9843059044327	53.1597242259238	53.3200041679201	53.4664400923771	53.6002168257507	53.7224188221382	53.8340385697837	53.9359843025098	54.0290870727125	54.1141072380185	54.1917404095216	54.2626229056626	54.3273367522700	54.3864142660083	54.4403422554744	54.4895658714083	54.5344921349394	54.5754931704358	54.6129091673703	54.6470510936244	54.6782031808259	54.7066252006341	54.7325545493369	54.7562081567058	54.7777842337441	54.7974638727623	54.8154125121105	54.8317812768808	54.8467082059608	54.8603193749608	54.8727299237500	54.8840449966152	54.8943606023876	54.9037644012767	54.9123364245889	54.9201497329923	54.9272710185212	54.9337611550753	54.9396757017770	54.9450653631797	54.9499764099892	54.9544510636521	54.9585278478821	54.9622419099401	54.9656253142428	54.9687073106623	54.9715145796751	54.9740714563406	54.9764001349193	54.9785208557886	54.9804520761743	54.9822106260850	54.9838118507211	54.9852697405211	54.9865970499082	54.9878054057111	54.9889054061486	54.9899067111921	54.9908181250512	54.9916476714626	54.9924026624052	54.9930897608103	54.9937150377881	54.9942840248451	54.9948017615284	54.9952728388925	54.9957014391546	54.9960913718672	54.9964461069133	54.9967688046009	54.9970623431087	54.9973293435161	54.9975721926272	54.9977930637825	54.9979939358336	54.9981766104436	54.9983427278568	54.9984937812761	54.9986311299666	37.3853237747908	35.8843719500663	34.4733482136390	33.1479959794781	31.9147948498699	30.8406194186350	29.9606649494407	29.2651355775671	28.7332032373169	28.3421896336469	28.0601543822673	27.8554452584911	27.7003544779019	27.5734978641209	27.4628937391544	27.3624952252085	27.2689330223601	27.1801325958251	27.0947118279056	27.0116960279587	26.9303657122514	26.8501661267300	26.7706580502089	26.6914776031119	26.6123123800650	26.5328871801524	26.4529436904643	26.3722313751296	26.2904932750132	26.2074546716662	26.1228083953464	26.0361990522003	25.9471996972685	25.8552806190629	25.7597665033635	25.6597793614118	25.5541783554720	25.4415597735150	25.3205111383071	25.1903620852168	25.0518686476175	24.9063629272546	24.7542629331229	24.5951538136499	24.4288406417710	24.2558123583442	24.0768657329932	23.8925610610723	23.7029702030367	23.5079683193805	23.3078582989010	23.1033593293659	22.8949881105129	22.6827953214307	22.4666458361438	22.2467421636169	22.0237040380135	21.7980117366084	21.5697265897567	21.3387697376288	21.1053211887367	20.8698796967753	20.6327998230234	20.3939975609490	20.1533246253063	19.9111781634321	19.6686108144965	19.4267675692757	19.1864663071324	18.9483615931963	18.7132450952370	18.4821137046030	18.2559139225370	18.0352822626019	17.8207548470107	17.6131571999546	17.4135166923583	17.2226338523675	17.0410202598784	16.8692078706450	16.7079450971462	16.5580126031316	16.4198006786058	16.2931311484149	16.1774149070917	16.0718759914185	15.9756952186742	15.8880774394644	15.8082773492339	15.7356068553520	15.6694349588763	15.6091848135583	15.5543298629039	15.5043898156206	15.4589267563981	15.4175415003522	15.3798702222176	15.3455813597704	15.3143727783755	15.2859691790169	15.2601197310834	15.2365959115604	15.2151895332684	15.1957109459842	15.1779873955342	15.1618615271534	15.1471900205571	15.1338423452298	15.1216996254299	15.1106536053059	15.1006057053578	15.0914661622330	15.0831532445455	15.0755925380425	15.0687162940235	15.0624628354512	15.0567760156763	15.0516047251443	15.0469024418556	15.0426268217227	15.0387393253038	15.0352048777017	15.0319915586991	15.0290703204571	15.0264147303399	15.0240007366421	15.0218064551904	15.0198119749706	15.0179991810932	15.0163515935598	15.0148542204279	15.0134934240948	15.0122567995359	15.0111330634342	15.0101119532314	15.0091841352190	15.0083411208624	15.0075751906261	15.0068793246303	15.0062471395296	15.0056728310585	15.0051511217372	15.0046772132782	15.0042467432703	15.0038557457616	15.0035006153877	15.0031780747325	15.0028851446280	15.0026191171319	15.0023775309422];
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
