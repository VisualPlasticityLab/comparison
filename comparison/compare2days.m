function compare2days(Day1,Day2)
%for Lab computer
        Day1=load('D:\Box\Enhancement-Andrea\1847\032917\000\data\1847_329_000_1\picked171_sp+ca\peakSI.mat');
        Day2=load('D:\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_1\picked455_sp+ca_3\peakSI.mat');
        
        p1=load('D:\Box\Enhancement-Andrea\1847\329&404\1847_329_000_1&1847_329_003_1&&1847_404_000_1&1847_404_001_1_plane1_selected142.mat');
%with anatomical only
%         Day1=load('D:\Box\Enhancement-Andrea\1847\032917\000&003\1847_329_000_1\picked49_sp+ca\peakSI.mat');
%         Day2=load('D:\Box\Enhancement-Andrea\1847\040417\000&001\1847_404_000_1\picked164_sp+ca\peakSI.mat');
%         
%         p1=load('D:\Box\Enhancement-Andrea\1847\329&404\1847_329_000&1847_404_000_plane_1_selected35.mat');


%        
%         %for home computer (ucl)
%         Day1=load('C:\Users\UCL053324\Box\Enhancement-Andrea\1847\032917\000\data\1847_329_000_1\picked171_sp+ca\peakSI.mat');
%         Day2=load('C:\Users\UCL053324\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_1\picked455_sp+ca_3\peakSI.mat');
% 
%         p1=load('C:\Users\UCL053324\Box\Enhancement-Andrea\1847\329&404\1847_329_000_1&1847_329_003_1&&1847_404_000_1&1847_404_001_1_plane1_selected142.mat');


%         %for home computer (surface)
%         Day1=load('C:\Users\Surface\Box\Enhancement-Andrea\1847\032917\000\data\1847_329_000_1\picked171_sp+ca\peakSI.mat');
%         Day2=load('C:\Users\Surface\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_1\picked455_sp+ca_3\peakSI.mat');
% 
%         p1=load('C:\Users\Surface\Box\Enhancement-Andrea\1847\329&404\1847_329_000_1&1847_329_003_1&&1847_404_000_1&1847_404_001_1_plane1_selected142.mat');




        pair1_final = p1.pair1; % 329 000&003 pair 1
        pair2_final=p1.pair2; %404 000&001 pair2

         [SI1.baselineR,SI1.baselineSDR,SI1.peakR,SI1.errorR,...
            SI1.baselineS,SI1.baselineSDS,SI1.peakS,SI1.errorS,cal_ER_option]= sigFcmp(Day1.sigF,Day1.win_sig,Day1.run.matrix);  % sigR:seg,1,Var,ncell
        [SI1.OSI_R,SI1.gOSI_R,SI1.DSI_R, SI1.gDSI_R,SI1.pref_ori_R,SI1.pref_dir_R]=calOSI(SI1.peakR);
        [SI1.OSI_S,SI1.gOSI_S,SI1.DSI_S, SI1.gDSI_S,SI1.pref_ori_S,SI1.pref_dir_S]=calOSI(SI1.peakS);
        
        [SI2.baselineR,SI2.baselineSDR,SI2.peakR,SI2.errorR,...
            SI2.baselineS,SI2.baselineSDS,SI2.peakS,SI2.errorS,cal_ER_option]= sigFcmp(Day2.sigF,Day2.win_sig,Day2.run.matrix);  % sigR:seg,1,Var,ncell
        [SI2.OSI_R,SI2.gOSI_R,SI2.DSI_R, SI2.gDSI_R,SI2.pref_ori_R,SI2.pref_dir_R]=calOSI(SI2.peakR);
        [SI2.OSI_S,SI2.gOSI_S,SI2.DSI_S, SI2.gDSI_S,SI2.pref_ori_S,SI2.pref_dir_S]=calOSI(SI2.peakS);
       
     

        OSI_change = SI2.OSI_R(pair2_final)- SI1.OSI_R(pair1_final); %OSI difference day5-day1
        pref_ori_change = SI2.pref_ori_R(pair2_final)- SI1.pref_ori_R(pair1_final); %prefered orinetation change from day5 to day1
        pref_dir_change = SI2.pref_dir_R(pair2_final)- SI1.pref_dir_R(pair1_final); %prefered direction change from day5 to day1
        SI1_mag= SI1.peakR(:,:,pair1_final); %peak responses for all orientaion day 1
        SI2_mag= SI2.peakR(:,:,pair2_final); %peak responses for all orientaion day 5
        all_mag= [SI1_mag; SI2_mag]; 
        SI1_max= max(SI1_mag(:,:,:)); % chose the max response to the prefered orinetaion
        SI2_max= max(SI2_mag(:,:,:));
        max_mag_change = log2(SI2_max./SI1_max); % magnitude change from day5 to day1
        max_mag_change = reshape(max_mag_change,1,142);
        mag_positive= find(max_mag_change>0); %find the positive change in magnitude for increased magnitude
        
      

        %select cells with positive or negative dOSI
        OSI_change_positive = find(OSI_change>0); 
        OSI_change_negative = find(OSI_change<0);
        


        %select cells with positive or negative dprefered
        %orientaion/direction
        pref_ori_change_positive =  find(pref_ori_change>0);
        pref_ori_nochange= find(pref_ori_change==0);%prefered orintation remains the same
        noorichangepoz=intersect(pref_ori_nochange,mag_positive);% cells with increased magnitude and the same prefered orination
        d0_pref_dir =  find(pref_dir_change==0); %prefered direction stays the same
        nodirchangepoz = intersect(mag_positive, d0_pref_dir); %cells with increased magnitude and the same prefered direction
        dorth_pref_dir = find(pref_dir_change==2);%prefered direction changes to orthogonal direction
        
        
        dpref_dir_0 = find(pref_dir_change == 0);%find cells with not chnage in the prefered direction
        dpref_dir_45_1 = find(pref_dir_change == 1);
        dpref_dir_45_3 = find(pref_dir_change == 3);
        dpref_dir_45_m3 = find(pref_dir_change == -3);
        dpref_dir_45_m1 = find(pref_dir_change == -1);
        dpref_dir_45= [dpref_dir_45_m1, dpref_dir_45_m3, dpref_dir_45_3, dpref_dir_45_1]; % cells with a 45 degree direction change
        dpref_dir_90_p = find(pref_dir_change == -2);
        dpref_dir_90_n = find(pref_dir_change == 2);
        dpref_dir_90= [dpref_dir_90_p dpref_dir_90_n];%cells with 90 degrees pref direction change
        
      
        % could not make  aviolin plot for all values as the array have
        % diffrent dimensions
        OSI_dir_change =[OSI_change(dpref_dir_0(1:42));OSI_change(dpref_dir_45(1:42));OSI_change(dpref_dir_90(1:42))];
        OSI_dir_change= OSI_dir_change';

        mag_dir_change = [max_mag_change(dpref_dir_0(1:42));max_mag_change(dpref_dir_45(1:42));max_mag_change(dpref_dir_90(1:42))];
        mag_dir_change = mag_dir_change';
        
        
        pref_ori_change_negative =  find(pref_ori_change<0);
        pref_dir_change_negative =  find(pref_dir_change<0);

        [OSI_mag_positive]=intersect(OSI_change_positive,mag_positive); %cells with increased magnitude and OSI

        [pref_dir0_mag_positive]=intersect(d0_pref_dir,mag_positive); %cells with increased magnitude and pref dir
        
        [OSI_pref_ori_positive]=intersect(OSI_change_positive,pref_ori_change_positive);
        [OSI_pref_dir_positive]=intersect(OSI_change_positive,d0_pref_dir);

        [OSI_mag_negative]=intersect(OSI_change_negative,mag_change_negative);%cells with decreased magnitude and OSI
        
        [OSI_pref_ori_negative]=intersect(OSI_change_negative,pref_ori_change_negative);
        [OSI_pref_dir_negative]=intersect(OSI_change_negative,pref_dir_change_negative);

        [OSI_positive_mag_negative]=intersect(OSI_change_positive,mag_change_negative);%cells with increased OSI and decreased magnitude 

        [OSI_positive_pref_ori_negative]=intersect(OSI_change_positive,pref_ori_change_negative);
        [OSI_positive_pref_dir_negative]=intersect(OSI_change_positive,pref_dir_change_negative);

        [OSI_negative_mag_positive]=intersect(OSI_change_negative,mag_positive);%cells with decreased OSI and increased magnitude

        [OSI_negative_pref_ori_positive]=intersect(OSI_change_negative,pref_ori_change_positive);
        [OSI_negative_pref_dir_positive]=intersect(OSI_change_negative,d0_pref_dir);

       
        Cor= [1:142];
        Cor1=[1:171];
        Cor2=[1:455];

        day1_picked=[134 28 31 40 112 38 113 76 92 68 146 33 109 54 86 49 98 153 161 71 22];
        day2_picked = [1 7 33 36 38 52 53 57 61 64 67 71 88 100 137 201 204 251 252 346 422 ];

        
  %% plot graphs      
  
        hsigF1=sigFplt(Day1.sigF,Day1.run.matrix,Day1.win_sig,Cor1(day1_picked),folder1);
        hsigF2=sigFplt(Day2.sigF,Day2.run.matrix,Day2.win_sig,Cor2(day2_picked),folder2);  % sigF:seg,rep,Var,ncell
  
        
         % violin plot
         figure ('Name','OSI vs prefered direction violin plot');
         title('OSI and prefered direction change violin plot');
         violinplot(OSI_dir_change);
         xlabel('prefered direction change');
         ylabel('OSI change');
         ax = gca;
         %ax.XAxisLocation = 'origin';
         %ax.YAxisLocation = 'origin';
         labelname={'0 degrees','45 degrees','90 degrees'};
         set(gca, 'XTick', 1:3, 'XTickLabels', labelname);
        
        
        figure('Name','OSI and magnitude change');
        title('magnitude and prefered direction change violin plot');
        violinplot(mag_dir_change);
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        ylabel('magnitude change');xlabel('prefered direction change');
        labelname={'0 degrees','45 degrees','90 degrees'};
         set(gca, 'XTick', 1:3, 'XTickLabels', labelname);
        
        figure('Name','magnitude and prefered direction change');
        subplot;hold on;
        title('magnitude and prefered direction Change')
        scatter(pref_dir_change,max_mag_change,'blue','*','jitter','on', 'jitterAmount', 0.15);hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        xlabel('magnitude change');ylabel('prefered direction change');
        direction = floor(size(sigF,3)/2);

        figure('Name','magnitude and OSI change');
        subplot;hold on;
        title('magnitude and OSI Change')
        scatter(OSI_change,max_mag_change,'blue','*','jitter','on', 'jitterAmount', 0.15);hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        xlabel('magnitude change');ylabel('OSI change');
        lsline
        R = corr2(OSI_change, max_mag_change);
   
    
        
        direction = floor(size(sigF,3)/2);
        
        figure('Name','OSI comparison');
        title('OSI_R change')
        histplt2(SI1.OSI_R(pair1_final),SI2.OSI_R(pair2_final))
        xlabel('OSI_R');ylabel('#cell');
        legend({'OSI 03/29','OSI 04/04'},'Location', 'best');


        
        
        figure('Name','OSI two dates comparison');
        subplot;hold on;
        title('OSI_R change')
        scatter(SI1.OSI_R(pair1_final),SI2.OSI_R(pair1_final),'b','jitter','on', 'jitterAmount', 0.15);
        plot([0 1],[0 1],'--');lsline;
        xlabel('OSI_R 03/29');ylabel('OSI_R 04/04');

        orientation = floor(size(sigF,3));

        figure('Name','Prefered orientation change');
        title('Prefered orientaion change')
        histplt2(SI1.pref_ori_R(pair1_final),SI2.pref_ori_R(pair2_final),.5:1:orientation+1.5)
        legend({'Pref Ori 03/29','Pref Ori 04/04'}, 'Location', 'bestoutside');
        legend('boxoff');
        xlabel('Pref Ori');ylabel('cell count');
        
        
        figure('Name','Prefered direction change');
        title('Prefered direction change')
        histplt2(SI1.pref_dir_R(pair1_final),SI2.pref_dir_R(pair2_final),.5:1:direction+1.5)
        legend({'Pref Dir 03/29','Pref Dir 04/04'}, 'Location', 'bestoutside');
        legend('boxoff');
        xlabel('Pref Dir');ylabel('cell count');

        figure('Name','Magnitude change');
        title('Magnitude change')
        histplt2(SI1_max,SI2_max)
        legend({'Magnitude 03/29','Magnitude 04/04'}, 'Location', 'bestoutside');
        legend('boxoff');
        xlabel('Magnitude');ylabel('cell count');


        picked = [4 20 19 42 50 49 56 70 64 96];
        polarplt(all_mag,Cor(nodirchangepoz));
        Name='Orientation tuning (polar)';
        legend ('03/29','04/04');   
%         
%         hpk(5)=figure('Name','Preferred Orientation & Direction');
%         subplot(1,2,1);histplt2(SI1.pref_ori_S,SI1.pref_ori_R,.5:1:direction*2+1.5);hold on;
%         subplot(1,2,1);histplt2(SI2.pref_ori_S,SI2.pref_ori_R,.5:1:direction*2+1.5);
%         legend({'still','run'},'orientation','horizontal')
%         xlabel('Orientation');ylabel('cell count');
%         subplot(1,2,2);histplt2(SI1.pref_dir_S,SI1.pref_dir_R,.5:1:direction+1.5);hold on;
%         subplot(1,2,2);histplt2(SI2.pref_dir_S,SI2.pref_dir_R,.5:1:direction+1.5);
%         legend({'still','run'},'orientation','horizontal')
%         xlabel('Direction');ylabel('cell count');
%         
%         [SI1.baseline,SI1.baselineSD,SI1.peak,SI1.error,cal_ER_option]=cal_ER(sigF,win_sig);  %peak:1*Var*ncell
%         [SI1.OSI,SI1.gOSI,SI1.DSI,SI1.gDSI,SI1.pref_ori,SI1.pref_dir]=calOSI(SI1.peak);
%         direction = floor(size(sigF,3)/2);
%         
%         hpk(6)=figure('Name','Overall Preferred Orientation&Direction');
%         subplot(1,2,1);h1 = histogram(SI1.pref_ori,.5:1:direction*2+1.5);
%         h1.FaceColor=[ 0.5 0.5 0.5];
%         legend('overall','Location','north','Orientation','horizontal');
%         xlabel('Orientation');ylabel('cell count');
%         subplot(1,2,2);h1 = histogram(SI1.pref_dir,.5:1:direction+1.5);
%         h1.FaceColor=[ 0.5 0.5 0.5];
%         legend('overall','Location','north','Orientation','horizontal');
%         xlabel('Direction');ylabel('cell count');
%         
% end