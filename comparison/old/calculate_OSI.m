clear all;close all;clc
[figname,path]=uigetfile('*.png','find fig with the labeled cells');
imshow(fullfile(path,figname));

%%
[contraeye,path]=uigetfile('*signals.mat','select trace for 1st');
eye1=load(fullfile(path,contraeye))
if isfield(eye1,'window')
    win_sig1 = [eye1.window(1) eye1.window(end)];
else
    win_sig1 = [eye1.win_sig(1) eye1.win_sig(end)];
end
if ~isfield(eye1,'matrix')
    eye1.matrix=eye1.run.matrix;
end
[~,eye1.SI.peakR,eye1.SI.errorR,~,~,~]= sigAcmp(eye1.sigF,win_sig1,eye1.matrix);  % sigR:seg,1,Var,ncell
% eye1.sigF= permute(eye1.sigF,[4 3 2 1]);
% eye1.finalvalueR = permute(eye1.SI.peakR,[3 2 1]);
% eye1.stdofeachvalueR = permute(eye1.SI.errorR,[3 2 1]);
% eye1.finalvalue = permute(eye1.SI.peakS,[3 2 1]);
% eye1.stdofeachvalue = permute(eye1.SI.errorS,[3 2 1]);


[ipsieye,path]=uigetfile('*signals.mat','select trace for 2nd');
eye2=load(fullfile(path,ipsieye))
if isfield(eye2,'window')
    win_sig2 = [eye2.window(1) eye2.window(end)];
else
    win_sig2 = [eye2.win_sig(1) eye2.win_sig(end)];
end
if ~isfield(eye2,'matrix')
    eye2.matrix=eye2.run.matrix;
end
[~,eye2.SI.peakR,eye2.SI.errorR,~,~,~]= sigAcmp(eye2.sigF,win_sig2,eye2.matrix);  % sigR:seg,1,Var,ncell
% eye2.sigF= permute(eye2.sigF,[4 3 2 1]);
% eye2.finalvalueR = permute(eye2.SI.peakR,[3 2 1]);
% eye2.stdofeachvalueR = permute(eye2.SI.errorR,[3 2 1]);
% eye2.finalvalue = permute(eye2.SI.peakS,[3 2 1]);
% eye2.stdofeachvalue = permute(eye2.SI.errorS,[3 2 1]);
% [~,pref_ori1]=nanmax(eye1.finalvalueR(:,1:end-1,:),[],2);
% [~,pref_ori2]=nanmax(eye2.finalvalueR(:,1:end-1,:),[],2);

[f,p]=uigetfile('*.mat','selected neurons');
load(fullfile(p,f));
ncell=numel(pair1);

% finalvalueR=(eye1.finalvalueR(pair1,:)+eye2.finalvalueR(pair2,:))/2;
% [~,pref_ori0]= nanmax(finalvalueR(:,1:end-1),[],2);
[OSI1,pref_ori1]=calOSI(eye1.SI.peakR);
[OSI2,pref_ori2]=calOSI(eye2.SI.peakR);

gOSI1=calgOSI(eye1.SI.peakR);
gOSI2=calgOSI(eye2.SI.peakR);
%%

pref_ori = nan(1,ncell);
ORI = nan(1,ncell);
ORI0= nan(1,ncell);
ORI1= nan(1,ncell);
ODI = nan(1,ncell);
ODI1 = nan(1,ncell);
ODI0= nan(1,ncell);
ODI_bi = nan(1,ncell);
ODI1_bi = nan(1,ncell);
ODI0_bi= nan(1,ncell);



%%
nSteps = size(eye1.finalvalueR,2);
for j=1:ncell
    nth1= pair1(j);
    nth2= pair2(j);
    h1=figure('Position',[100 200 1500 500],'Name',['alltraces_cell#' num2str(j)])
    subplot(3,3,1);    hold on;
    %     errorbar(1:nSteps,eye1.finalvalue(j,:),eye1.stdofeachvalue(j,:),'r*--');
    %     errorbar(1:nSteps,eye2.finalvalue(j,:),eye2.stdofeachvalue(j,:),'b+--');
    errorbar(1:nSteps,eye1.finalvalueR(nth1,:),eye1.stdofeachvalueR(nth1,:),'r*-');
    errorbar(1:nSteps,eye2.finalvalueR(nth2,:),eye2.stdofeachvalueR(nth2,:),'b+-');

    
    %[odi,odi_bi] = ODIcalculation(pref_ori,eye1.finalvalueR(nth1,:),eye2.finalvalueR(nth2,:));
    %text(0,-a1/10,sprintf('(%.2f)(%.2f)',odi,odi_bi))
    
    subplot(3,3,3);    hold on;
    polarplt(cat(1,eye1.finalvalueR(nth1,:),eye2.finalvalueR(nth2,:)));
    %     legend('contra','ipsi','contra run','ipsi run')
    %     plot([1 nSteps-1],[contra.finalvalue_z(j,end) contra.finalvalue_z(j,end)],'r-')
    %     plot([1 nSteps-1],[ipsi.finalvalue_z(j,end) ipsi.finalvalue_z(j,end)],'b-')
    % title(['stillODI: ' num2str(ODI_pref(j)) 'runODI: ' num2str(ODI_pref2(j))]);
    
    temp_1st=squeeze(eye1.sigF(nth1,:,:,:));
    temp_2nd=squeeze(eye2.sigF(nth2,:,:,:));
    ymax=max([temp_1st(:);temp_2nd(:)]);
    ymin=min([temp_1st(:);temp_2nd(:)]);
    buffer=(ymax-ymin)*.1;
    for ori=1:nSteps
        %     a=subplot(3,nSteps,ori); hold on;
        %     plot(squeeze(temp_both(ori,:,:))');
        %     plot(squeeze(mean(temp_both(ori,:,:))),'r','linewidth',3);
        %     axis([1 framestocapture ymin-buffer ymax+buffer])
        %     title(num2str(ori));
        
        b=subplot(3,nSteps,nSteps+ori); hold on;
        plot(squeeze(temp_1st(ori,eye1.matrix(:,ori),:))','b','linewidth',1);
        plot(squeeze(mean(temp_1st(ori,eye1.matrix(:,ori),:))),'r','linewidth',3);
        axis([win_sig1(1)-5 win_sig1(end)+5 ymin-buffer ymax+buffer])
        title(sprintf('%.3f',eye1.finalvalueR(nth1,ori)));
        b.XLabel.String=sprintf('%d',ori);
        
        c=subplot(3,nSteps,2*nSteps+ori); hold on;
        plot(squeeze(temp_2nd(ori,eye2.matrix(:,ori),:))','b','linewidth',1);
        plot(squeeze(mean(temp_2nd(ori,eye2.matrix(:,ori),:))),'r','linewidth',3);
        axis([win_sig2(1)-5 win_sig2(end)+5 ymin-buffer ymax+buffer])
        title(sprintf('%.3f',eye2.finalvalueR(nth2,ori)));
        if ori==1
            %         a.YLabel.String='both eye';
            b.YLabel.String=sprintf('cell# %d',nth1);
            c.YLabel.String=sprintf('cell# %d',nth2);
        else
            b.YTickLabel=''
            b.XTickLabel=''
            c.YTickLabel=''
            c.XTickLabel=''
        end
    end
    
    m=inputdlg('Good cell& best ori=','Selected cell?',1,{num2str(pref_ori)})
    if ~isempty(m)
        ORI(j) =str2num(m{1});
        ORI1(j) = pref_ori;
        ORI0(j) = ori0;
        savefig(h1,['alltraces_cell#' num2str(j)])
    end
    close(h1)
    
end
%%
for j=1:ncell
    if ~isnan(ORI1(j))
        nth1= pair1(j);
        nth2= pair2(j);
        [ODI(j),ODI_bi(j)]= ODIcalculation(ORI(j),eye2.finalvalueR(nth2,:),eye1.finalvalueR(nth1,:));
        [ODI1(j),ODI1_bi(j)]= ODIcalculation(ORI1(j),eye2.finalvalueR(nth2,:),eye1.finalvalueR(nth1,:));
        [ODI0(j),ODI0_bi(j)]= ODIcalculation(ORI0(j),eye2.finalvalueR(nth2,:),eye1.finalvalueR(nth1,:));
    end
end
%%
h=figure('Name',sprintf('ODI of total cell %d', numel(ODI(ODI>-10))));
seg=.3;

subplot(3,2,1);hold on;title('Handpicked')
histogram(ODI,[-1.2:seg:1.2]);
subplot(3,2,2);hold on;title('Handpicked-bidirectional')
histogram(ODI_bi,[-1.2:seg:1.2]);

subplot(3,2,3);hold on;title('Separate-fitted')
histogram(ODI1,[-1.2:seg:1.2]);
subplot(3,2,4);hold on;title('Separate-fitted-bidirectional')
histogram(ODI1_bi,[-1.2:seg:1.2]);

subplot(3,2,5);hold on;title('Combined-fitted')
histogram(ODI0,[-1.2:seg:1.2]);
subplot(3,2,6);hold on;title('Combined-fitted-bidirectional')
histogram(ODI0_bi,[-1.2:seg:1.2]);

ODIfigname= fullfile(p,[f(1:end-4) '_ODI']);
saveas(h,[ODIfigname '.fig']);
saveas(h,[ODIfigname '.png']);

%%
save(fullfile(p,f),'ORI','ORI0','ORI1','-append');

