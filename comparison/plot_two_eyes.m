%% Load file
[file,path]=uigetfile('*signals.mat','select sigF file');
signals1=load(fullfile(path,file))
matrix1 = signals1.matrix;
if isfield(signals1,'window')
    win_sig1 = [signals1.window(1) signals1.window(end)];
else
    win_sig1 = [signals1.win_sig(1) signals1.win_sig(end)];
end
sigF1 = signals1.sigF;
% [~,SI.peakR1,SI.errorR1,~,SI.peakS1,SI.errorS1]= sigFcmp(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
% peak1=fixpeak(SI.peakR1);

[file2,path2]=uigetfile('*signals.mat','select sigF file');
signals2=load(fullfile(path2,file2))
matrix2 = signals2.matrix;
if isfield(signals2,'window')
    win_sig2 = [signals2.window(1) signals2.window(end)];
else
    win_sig2 = [signals2.win_sig(1) signals2.win_sig(end)];
end
sigF2 = signals2.sigF;

[f,p]=uigetfile('*.mat','selected neurons');
load(fullfile(p,f));

%% 
% [~,SI.peakR1,SI.errorR1,~,SI.peakS1,SI.errorS1]= sigFcmp(sigF1,win_sig1,matrix1);  % sigR:seg,1,Var,ncell
% peak1=fixpeak(SI.peakR1);
% [pref1,pref_ori1]=nanmax(peak1(:,1:end-1,:),[],2);

% [~,SI.peakR2,SI.errorR2,~,SI.peakS2,SI.errorS2]= sigFcmp(sigF2,win_sig2,matrix2);  % sigR:seg,1,Var,ncell
% peak2=fixpeak(SI.peakR2);
% [pref2,pref_ori2]=nanmax(peak2(:,1:end-1,:),[],2);
%%
[peakR1,xR1,yR1]= sigFcmp(sigF1(:,:,1:end-1,pair1),win_sig1,matrix1);
[peakR2,xR2,yR2]= sigFcmp(sigF2(:,:,1:end-1,pair2),win_sig2,matrix2);
[pref1,pref_ori1]= nanmax(peakR1,[],2);
[pref2,pref_ori2]= nanmax(peakR2,[],2);
pref_ori = nan(size(pref_ori1));
ODI = nan(size(pref1));
    
%%
figure
ncell=numel(pair1);
nrow=ceil(sqrt(ncell));
ncol=nrow*2;
for i=1:ncell
    subplot(nrow,ncol,i*2-1);
    nth1=pair1(i);
    temp1= Fit2Gaussian(pref_ori1(1,:,nth1),xR1,yR1(:,nth1));
    
    subplot(nrow,ncol,i*2);
    nth2=pair2(i);
    temp2 = Fit2Gaussian(pref_ori2(1,:,nth2),xR2,yR2(:,nth2));
    
    
    
end

