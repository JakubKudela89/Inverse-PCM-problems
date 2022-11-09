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
        t_presim = [25	25	25	25	25	37.5469102171437	38.2231011650285	38.9143687521620	39.5624533882718	40.1586174277111	40.7017632492257	41.1930322332004	41.6351041902665	42.0317929456120	42.3874182634012	42.7062681109355	42.9926432326891	43.2504849479272	43.4833727363839	43.6946087184328	43.8871007788371	44.0633551711254	44.2255055218431	44.3754099355575	44.5146113981291	44.6445157673529	44.7661730753793	44.8805991110035	44.9886023921442	45.0909312866545	45.1881766432177	45.2808989001486	45.3695265308459	45.4544568879513	45.5360626469607	45.6146415419269	45.6904350590849	45.7637271472064	45.8346996629193	45.9035430229240	45.9704176291231	46.0354692685231	46.0988531418136	46.1606682982439	46.2210182082266	46.2800062902573	46.3377236654590	46.3942423762675	46.4496379645133	46.5039682289276	46.5572877433832	46.6096511727543	46.6611050016851	46.7116861868638	46.7614155251146	46.8103265525823	46.8584779012456	46.9062223412590	46.9543525478292	47.0013278465554	47.0502393155495	47.1035188957076	47.1795578857009	47.2538344935247	47.3248228619110	47.3936671092068	47.4606873223449	47.5330473559425	47.6124669677311	47.6977927275291	47.7875380476046	47.8763215291414	47.9720428264665	48.0637850203103	48.1517393531438	48.2410592473441	48.3411955319471	48.4504550161452	48.5561856489004	48.6604675984283	48.7705801650694	48.8752206527434	48.9765597549179	49.0818316071614	49.2029467870942	49.3240466764322	49.4395842944186	49.5544550206645	49.6734915774426	49.7875166343861	49.9015798886598	50.0238297279498	50.1559615338049	50.2826008696858	50.4055730271654	50.5270065244390	50.6482867399169	50.7672898535915	50.8974892878283	51.0415563466287	51.1818680083298	51.3172627560857	51.4509687372569	51.5885819251481	51.7256868103920	51.8596074972027	51.9989192366859	52.1436526214570	52.2840632185398	52.4152997636662	52.5444999419546	52.6714385198436	52.8062077058347	52.9435752626325	53.0720949911689	53.1925179436211	53.3118753768918	53.4278466450714	53.5369985551932	53.6503780219565	53.7629597164810	53.8671269918103	53.9625895519124	54.0501183590148	54.1304049355687	54.2043672038247	54.2721934511483	54.3348929977294	54.3926325938966	54.4453926542450	54.4935791202483	54.5375876242455	54.5777796484821	54.6144854664148	54.6480068312752	54.6786194335311	54.7065751470643	54.7321040821933	54.7554164621098	54.7767043378876	54.7961431559182	54.8138931904441	54.8301008527740	54.8448998877734	54.8584124673123	54.8707501895274	54.8820149919910	54.8922999861896	54.9016902200775	54.9102633748903	37.3748296504520	36.0021386008609	34.7570971166170	33.7462101749014	33.0172582369116	32.5274487129711	32.2375005582540	32.0717350139022	31.9682058524977	31.8792530262479	31.7942161147602	31.7100516686187	31.6256417103299	31.5405066401504	31.4545376614444	31.3675819791711	31.2796792577793	31.1906389440702	31.1004173450384	31.0088362908941	30.9156532001236	30.8207238802638	30.7236649026318	30.6241430412716	30.5217567681779	30.4159240947067	30.3060657474003	30.1913974779921	30.0710438917304	29.9440764508434	29.8095089509647	29.6664836048452	29.5143686312647	29.3529521113029	29.1822541342971	29.0020174931455	28.8124093573262	28.6137875104559	28.4065930954150	28.1915401077304	27.9697169495632	27.7464457444611	27.5198994571573	27.2852056080838	27.0418083654242	26.7934373127384	26.5417797169416	26.2871562418224	26.0286828337499	25.7652800241816	25.4956892741757	25.2235761515305	24.9514946531772	24.6763128721861	24.3985201974413	24.1171238260939	23.8334899200527	23.5524178747837	23.2719301824701	22.9908729127007	22.7104502162197	22.4308070167492	22.1555129168775	21.8855210948344	21.6186664112514	21.3555570001651	21.0967618527376	20.8460618366931	20.6057037608783	20.3724977073565	20.1460401945951	19.9274408762452	19.7166150734259	19.5150092512290	19.3230414840112	19.1381534946386	18.9602972657066	18.7906411406351	18.6310354685568	18.4786303090771	18.3320735794834	18.1921592883004	18.0582717317322	17.9318714795079	17.8113795053596	17.6950228692904	17.5823943651811	17.4735063389245	17.3686476534058	17.2680640098824	17.1718886287294	17.0802948912291	16.9932548026779	16.9106035636187	16.8322283072080	16.7579276164734	16.6875700915161	16.6210226816746	16.5581272570134	16.4987259185982	16.4425419764788	16.3894115714076	16.3391846670201	16.2916540981415	16.2467658566574	16.2042990053260	16.1640808624460	16.1259936921231	16.0898563275418	16.0555633501235	16.0229960636658	15.9920218744187	15.9625276782087	15.9343919821638	15.9075131918839	15.8817985246061	15.8571541509221	15.8334909511596	15.8107277366722	15.7887933753367	15.7676030784377	15.7471115354580	15.7272813280386	15.7080593834099	15.6894076660803	15.6712953431044	15.6536867958882	15.6365531451033	15.6198721281760	15.6036261319908	15.5877968511729	15.5723686595725	15.5573251213450	15.5426537607214	15.5283460298071	15.5143873677287	15.5007648005470	15.4874689874602	15.4744942512525	15.4618317257636	15.4494700455270	15.4374019734203	15.4256238540522	15.4141304211295	15.4029123732175	15.3919596474466	15.3812689996280	15.3708378164043	15.3606596391853	15.3507262722062];
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
