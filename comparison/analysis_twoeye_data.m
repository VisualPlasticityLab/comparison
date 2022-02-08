
clear all;close all;clc
% func_save_data2('ipsi')
% func_save_data2('contra')

[figname,figpath]=uigetfile('*.fig','find fig with the labeled cells');
open(fullfile(figpath,figname));

% botheye=uigetfile('*_traces.mat','open trace for botheye');
% both=load(botheye)
%%

[contraeye,path]=uigetfile('*signals.mat','select trace for contraeye');
contra=load(fullfile(path,contraeye))
contra.sigF= permute(contra.sigF,[4 3 2 1]);
[contraeye2,path]=uigetfile('peakSI.mat','select peakSI matfile for contraeye');
contra2=load(fullfile(path,contraeye2));
contra.finalvalue = permute(contra2.SI.peakS,[3 2 1]);
contra.stdofeachvalue = permute(contra2.SI.errorS,[3 2 1]);
contra.finalvalueR = permute(contra2.SI.peakR,[3 2 1]);
contra.stdofeachvalueR = permute(contra2.SI.errorR,[3 2 1]);


[ipsieye,path]=uigetfile('*signals.mat','open trace for ipsieye');
ipsi=load(fullfile(path,ipsieye))
ipsi.sigF= permute(ipsi.sigF,[4 3 2 1]);
[ipsieye2,path]=uigetfile('peakSI.mat','select peakSI matfile for ipsieye');
ipsi2=load(fullfile(path,ipsieye2));
ipsi.finalvalue = permute(ipsi2.SI.peakS,[3 2 1]); 
ipsi.stdofeachvalue = permute(ipsi2.SI.errorS,[3 2 1]);
ipsi.finalvalueR = permute(ipsi2.SI.peakR,[3 2 1]);
ipsi.stdofeachvalueR = permute(ipsi2.SI.errorR,[3 2 1]);
%% calculate ODI using preferred orientation
ncells = size(contra.sigF,1);
nSteps=size(contra.sigF,2);
framestocapture=size(contra.sigF,4);
nSteps1=nSteps-1;


finalvalue=cat(3,contra.finalvalue,ipsi.finalvalue); %finalvalue size: ncells,nSteps,1,1
blank=finalvalue(:,end,:);
response=finalvalue(:,1:end-1,:);
[~,pref_contra]=max(response(:,:,1),[],2);
[~,pref_ipsi]=max(response(:,:,2),[],2);
response_comb=(response(:,:,1)+response(:,:,2))/2;
[response_pref,pref]=max(response_comb,[],2);

pos=sub2ind(size(contra.finalvalue),1:ncells,pref');
% ODI_pref=(contra.finalvalue(pos)-ipsi.finalvalue(pos))./(ipsi.finalvalue(pos)+contra.finalvalue(pos));
ODI_prefR=(contra.finalvalueR(pos)-ipsi.finalvalueR(pos))./(ipsi.finalvalueR(pos)+contra.finalvalueR(pos));

%%
cellnum=1:ncells-1;
defaultgd=cellstr(num2str(cellnum));
selected=inputdlg('input good cell#','includes:',[4 80],defaultgd);
goodcell(str2num(selected{1}))=1;

figure;hold on;
plot(abs(mod(pref_contra(logical(goodcell))-pref(logical(goodcell))+1,4)-1),'ro')
plot(abs(mod(pref_ipsi(logical(goodcell))-pref(logical(goodcell))+1,4)-1),'b+')
ylim([-.5 2.5]);
g=gca;
g.XLabel.String='Neurons';
g.YTick=[0 1 2];
g.YTickLabel={'0','\pi/4','\pi/2'};
g.YLabel.String='Difference of preferred orientation';
legend('contra','ipsi')
saveas(gca,'preferred orientation.fig')


h=figure;
hold on;
seg=.3;
histogram(ODI_pref(logical(goodcell)),[-1.2:.3:1.2]);
hold on;
histogram(ODI_prefR(logical(goodcell)),[-1.2:.3:1.2]);

title(sprintf('total cell %d, using preferred orientation', sum(goodcell)));
ODIfigname= [figname(1:end-4) 'good' num2str(sum(goodcell)) '_ODI_prefori'];
saveas(h,[ODIfigname '.fig']);
saveas(h,[ODIfigname '.png']);
save([ODIfigname '.mat'] ,'ODI_pref','goodcell');

% sigF:ncells,nSteps, numberofcycles,framestocapture
%% hand pick cells or directl7 load from ,mat
[filename,path]=uigetfile('Pick the preselected good cell?','*.mat');
if filename
    temp=load(fullfile(path,filename));
    goodcell = temp.goodcell;
else
    
    goodcell=zeros(1,ncells);
    
    for j=1:ncells-1
        h1=figure('Position',[100 200 1000 600],'Name',['alltraces_cell#' num2str(j)])
        
        subplot(3,1,1);    hold on;
        errorbar(1:nSteps,contra.finalvalue(j,:),contra.stdofeachvalue(j,:),'r*--');
        errorbar(1:nSteps,ipsi.finalvalue(j,:),ipsi.stdofeachvalue(j,:),'b+--');
        errorbar(1:nSteps,contra.finalvalueR(j,:),contra.stdofeachvalueR(j,:),'r*-');
        errorbar(1:nSteps,ipsi.finalvalueR(j,:),ipsi.stdofeachvalueR(j,:),'b+-');
        legend('contra','ipsi','contra run','ipsi run')
        %     plot([1 nSteps-1],[contra.finalvalue_z(j,end) contra.finalvalue_z(j,end)],'r-')
        %     plot([1 nSteps-1],[ipsi.finalvalue_z(j,end) ipsi.finalvalue_z(j,end)],'b-')
        xlim([1 nSteps])
       % title(['stillODI: ' num2str(ODI_pref(j)) 'runODI: ' num2str(ODI_pref2(j))]);
        
        temp_contra=squeeze(contra.sigF(j,:,:,:));
        temp_ipsi=squeeze(ipsi.sigF(j,:,:,:));
        ymax=max([temp_contra(:);temp_ipsi(:)]);
        ymin=min([temp_contra(:);temp_ipsi(:)]);
        buffer=(ymax-ymin)*.1;
        for ori=1:nSteps
            %     a=subplot(3,nSteps,ori); hold on;
            %     plot(squeeze(temp_both(ori,:,:))');
            %     plot(squeeze(mean(temp_both(ori,:,:))),'r','linewidth',3);
            %     axis([1 framestocapture ymin-buffer ymax+buffer])
            %     title(num2str(ori));
            
            b=subplot(3,nSteps,nSteps+ori); hold on;
            plot(squeeze(temp_contra(ori,:,:))');
            plot(squeeze(mean(temp_contra(ori,:,:))),'r','linewidth',3);
            axis([1 framestocapture ymin-buffer ymax+buffer])
            title(num2str(contra.finalvalue(j,ori)));
            
            c=subplot(3,nSteps,2*nSteps+ori); hold on;
            plot(squeeze(temp_ipsi(ori,:,:))');
            plot(squeeze(mean(temp_ipsi(ori,:,:))),'r','linewidth',3);
            axis([1 framestocapture ymin-buffer ymax+buffer])
            title(num2str(ipsi.finalvalue(j,ori)))
            if ori==1
                %         a.YLabel.String='both eye';
                b.YLabel.String='contra eye';
                c.YLabel.String='ipsi eye';
            else
                b.YTickLabel=''
                c.YTickLabel=''
            end
        end
        drawnow;
        if menu('Good cell?','Yes','No') ==1
            goodcell(j)=1;
            savefig(h1,['alltraces_cell#' num2str(j)])
        end
    end
end
%%
figure;hold on;
plot(abs(mod(pref_contra(logical(goodcell))-pref(logical(goodcell))+1,4)-1),'ro')
plot(abs(mod(pref_ipsi(logical(goodcell))-pref(logical(goodcell))+1,4)-1),'b+')
ylim([-.5 2.5]);
g=gca;
g.XLabel.String='Neurons';
g.YTick=[0 1 2];
g.YTickLabel={'0','\pi/4','\pi/2'};
g.YLabel.String='Difference of preferred orientation';
legend('contra','ipsi')
saveas(gca,'preferred orientation.fig')

% 
% h=figure;
% hold on;
% seg=.3;
% histogram(ODI_pref(logical(goodcell)),[-1.2:.3:1.2]);
% hold on;
% histogram(ODI_pref2(logical(goodcell)),[-1.2:.3:1.2]);
% 
% title(sprintf('total cell %d, using preferred orientation', sum(goodcell)));
% ODIfigname= [figname(1:end-4) 'good' num2str(sum(goodcell)) '_ODI_prefori'];
% saveas(h,[ODIfigname '.fig']);
% saveas(h,[ODIfigname '.png']);
% save([ODIfigname '.mat'] ,'ODI_pref','ODI_pref2','goodcell');

%% calculating the difference
% oppo=mod(pref+nSteps1/2-1,nSteps1)+1
% response_oppo=response_comb(:,oppo);
%
% orth1=mod(pref+nSteps1/4-1,nSteps1)+1;
% orth2=mod(pref-nSteps1/4-1,nSteps1)+1;
% response_orth=(response_comb(:,orth1)+response_comb(:,orth2))/2;
%
% figure;
% subplot(1,2,1)
% hold on
% scatter(pref,pref_ipsi,'jitter','on','jitterAmount',0.02)
% scatter(pref,pref_contra,'r+','jitter','on','jitterAmount',0.02)
% legend({'pref_ ipsi','pref_ contra'})
% subplot(1,2,2)
% hold on
% scatter(pref_contra,pref_ipsi,'jitter','on','jitterAmount',0.02)
% P = polyfit(pref_contra,pref_ipsi,1);
% yfit = P(1)*pref_contra+P(2);
% plot(pref_contra,yfit,'r-.');
% title('pref_ contra vs pref_ ipsi')

