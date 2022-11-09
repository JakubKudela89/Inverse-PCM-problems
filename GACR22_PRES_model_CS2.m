function [f_cost] = GACR22_PRES_model_2(vec)
        %GACR22_PRES_model_2([2500, 31000, 30, 15, 15, 800, 0.22]) - je reseni s f_cost = 0
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
        t_presim = [25	25	25	25	25	37.2475855403868	37.4649177870157	37.6484700616800	37.8092731030520	37.9532189758091	38.0841634038152	38.2047872601970	38.3170345867636	38.4223625497692	38.5218934209866	38.6165115277775	38.7069275509229	38.7937225779518	38.8773791388311	38.9583036092668	39.0368427371472	39.1132960770365	39.1879255205415	39.2609627327126	39.3326150600325	39.4030703140140	39.4725007262505	39.5410662978797	39.6089177177578	39.6761989924399	39.7430499134552	39.8096084816266	39.8760134144918	39.9424068835635	40.0089376685704	40.0757649861756	40.1430633698821	40.2110291811564	40.2798896867279	40.3499162774012	40.4214446055808	40.4949067725248	40.5708854948558	40.6502101727711	40.7341344369373	40.8246573468439	40.9249238093141	41.0384357654522	41.1643429198995	41.2971059225161	41.4348306170911	41.5797918476420	41.7346343691706	41.8995115372502	42.0728087394444	42.2526743174159	42.4376993213480	42.6279580831168	42.8257112095648	43.0318927421702	43.2442687024687	43.4607365434919	43.6807440583849	43.9047509137596	44.1347871759117	44.3716641712249	44.6129936461920	44.8567830766182	45.1031413117975	45.3527406844390	45.6069626872553	45.8663892811893	46.1290005943060	46.3928399105964	46.6583765184758	46.9275490443876	47.2024888360934	47.4831859366542	47.7664460443921	48.0496190965409	48.3323452257624	48.6147023161076	48.8968222731693	49.1788334836617	49.4591890059304	49.7350611409124	50.0057384371265	50.2726423947584	50.5365190079062	50.7959064670161	51.0480228817434	51.2916813757059	51.5273443099106	51.7552087952279	51.9752938062347	52.1866051888814	52.3866602607164	52.5737874422370	52.7479058925627	52.9096707693730	53.0598992782863	53.1993979571104	53.3289261332208	53.4491915042068	53.5608522182156	53.6645200430126	53.7607635565550	53.8501111632909	53.9330539101240	54.0100481082928	54.0815177730191	54.1478568931906	54.2094315427567	54.2665818447748	54.3196237983009	54.3688509776279	54.4145361127182	54.4569325590725	54.4962756647074	54.5327840413890	54.5666607467748	54.5980943836600	54.6272601220959	54.6543206497503	54.6794270555084	54.7027196509695	54.7243287341698	54.7443752995648	54.7629716980233	54.7802222503273	54.7962238174284	54.8110663304851	54.8248332834985	54.8376021911649	54.8494450143833	54.8604285556870	54.8706148267089	54.8800613896449	54.8888216745411	54.8969452741053	54.9044782176234	54.9114632254504	54.9179399454455	54.9239451726214	54.9295130531929	54.9346752741246	54.9394612392008	54.9438982325698	54.9480115706475	54.9518247432033	37.6750905970732	36.4531233076513	35.2897372876459	34.1827035757263	33.1301501047819	32.1340097608809	31.2101215625103	30.3791450359428	29.6500141157331	29.0218151353596	28.4887884149769	28.0428075107986	27.6736416721439	27.3694636666372	27.1181170394650	26.9082276770122	26.7300151508060	26.5757494749457	26.4397069873887	26.3177758713747	26.2070047111016	26.1052429649767	26.0108904696695	25.9227302921347	25.8398176551332	25.7614052044114	25.6868917756934	25.6157865438877	25.5476834227861	25.4822424310859	25.4191758864491	25.3582380088856	25.2992169741878	25.2419287568532	25.1862122997029	25.1319256806891	25.0789430387205	25.0271520839327	24.9764520627682	24.9267520804257	24.8779697066035	24.8300298076167	24.7828635607080	24.7364076159240	24.6906033781594	24.6453963874885	24.6007357801436	24.5565738157735	24.5128654591586	24.4695680065304	24.4266407481692	24.3840446601297	24.3417421188204	24.2996966328034	24.2578725866015	24.2162349915267	24.1747492385742	24.1333808482640	24.0920952119215	24.0508573182463	24.0096314580588	23.9683808987551	23.9270675181266	23.8856513846296	23.8440902676950	23.8023390569043	23.7603490623502	23.7180671596071	23.6754347305714	23.6323863348576	23.5888480241439	23.5447351827301	23.4999497419263	23.4543765790198	23.4078788930477	23.3602924072079	23.3114185358672	23.2610175554814	23.2088051864280	23.1544613231400	23.0976688053985	23.0382059725185	22.9760888933118	22.9116661528818	22.8455138730542	22.7781579206868	22.7098604914129	22.6405948967136	22.5701359286058	22.4981814630478	22.4244755324241	22.3489041720919	22.2715195230021	22.1924952515367	22.1120632350972	22.0304502209926	21.9478323170809	21.8643297470371	21.7799992423052	21.6947912601846	21.6085315086693	21.5209903249310	21.4320209191555	21.3416779669813	21.2502165098267	21.1579558898104	21.0651217474206	20.9717913514377	20.8779381031356	20.7834764653097	20.6882631704290	20.5921132584037	20.4948872092902	20.3966099982742	20.2975105280488	20.1979116811394	20.0980505382737	19.9979841510778	19.8976290561260	19.7968518583756	19.6955220673229	19.5935277719064	19.4908120613542	19.3874380048327	19.2836116773701	19.1796008407601	19.0755927651547	18.9716036748330	18.8674889978940	18.7630149579567	18.6579387516663	18.5520845515072	18.4454276008988	18.3381575009674	18.2306465155556	18.1233073274288	18.0164392658449	17.9101725039245	17.8045184382308	17.6994623160850	17.5950359389904	17.4913406986235	17.3885372593293	17.2868350812045	17.1864773778736	17.0876987426619	16.9906626151016	16.8954157616895	16.8019022143632	16.7100414786966];
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
