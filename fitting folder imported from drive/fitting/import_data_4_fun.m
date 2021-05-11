function [freq, data] = import_data_4_fun(folder)


%New in version 4: 
%1. read parameters from paramsFile
%2. colormap according to power

% close all
% clear all
% hold all
% folder = 'Z:\Measurements\2020.dir\13_07_20 cooldown\Israel\Transmission and Reflections Initial\Nockit_port_10-2';
saveFigures = false;
%Circulators=-2; HEMT=+40; BAZ=+26; PwrSpl_4K=-9; Cable=-12; LNF=0;
attenuation = -20-10-12; %-20-20-10+Cable;
amplification= 0 ;%HEMT+BAZ+2*Circulators;%+PwrSpl_4K;

paramsFileName = 'paramsFile.txt';

firstFigOff = 1;

saveFolder = fullfile(folder,'figs');

 
fid = fopen(fullfile(folder,paramsFileName));

%read resonances and power levels
line = fgetl(fid);
while ischar(line)
    if ~isempty(strfind(line,'Resonances'))
        s = regexp(line,'\[');
        e = regexp(line,'\]');
        resonances = regexp(line(s+1:e-1),',','split');
        %remove ' and space characters
        for ri=1:length(resonances)
            resonances{ri}(strfind(resonances{ri},'''')) = [];
            resonances{ri}(strfind(resonances{ri},' ')) = [];
        end
    elseif ~isempty(strfind(line,'Power levels'))
            while isempty(strfind(line,'] dBm')) 
                line = [line fgetl(fid)];
            end
            s = regexp(line,'\[');
            e = regexp(line,'\]');
            pwrsStr = regexp(line(s+1:e-1),',','split');
            pwrs = str2num(cell2mat(pwrsStr));
            break;
    end
    line = fgetl(fid);
end
fclose(fid);

if length(pwrs)>1
%     ssetRainbowColors(pwrs(1:end-1));
end

if saveFigures
    mkdir(saveFolder);
end

for rn=1:1:length(resonances)
    resonance = resonances{rn};
    freq = importdata(fullfile(folder,[resonance,'_freq.out']));
    data_raw = importdata(fullfile(folder,[resonance,'_data.out']));
    data = data_raw(:,1:2:end-1)+1i*data_raw(:,2:2:end);
    magM = abs(data);
    %figures
%     h1(rn) = figure(2*rn+firstFigOff); hold all;
%     plot(freq/1e6,20*log10(magM)-attenuation-amplification);
%     plot(freq/1e6,20*log10(magM));
%     title(strcat(resonance,' GHZ'))  
%     xlabel('Frequency [MHz]'); ylabel('Transmission [dB]');
%     h2(rn) = figure(2*rn+1+firstFigOff); hold all;
%     for p=1:size(data,1)
%         plot(real(data(p,:)),imag(data(p,:)),'.');
%     end
%     title(strcat(resonance,' GHZ '))
%     xlabel('Re(S_{21})'); ylabel('Im(S_{21})');
%     h3(rn) = figure(1000+rn+firstFigOff); imagesc(freq/1e6,pwrs+attenuation,20*log10(magM)); 
%     axis xy; colorbar; 
%     xlabel('Frequency [MHz]'); ylabel('Power [dBm]'); title([resonance]); 
end


pwrsLgnd = {};
for c=1:1:length(pwrs)-1
    pwrsLgnd{end+1} = num2str(pwrs(c)-attenuation);
end

for rn=1:1:length(resonances)
%     figure(2*rn+firstFigOff)
%     legend(pwrsLgnd);
%     colormap jet; cb =colorbar; caxis([min(pwrs+attenuation) max(pwrs+attenuation)]);
%     cb.Label.String = 'Power [dBm]';
%     figure(2*rn+1+firstFigOff)
%     legend(pwrsLgnd);
    if saveFigures
        saveas(h1(rn),fullfile(saveFolder,['mag_' resonances{rn} '.fig']));
        saveas(h1(rn),fullfile(saveFolder,['mag_' resonances{rn} '.jpg']));    
        saveas(h2(rn),fullfile(saveFolder,['pol_' resonances{rn} '.fig']));
        saveas(h2(rn),fullfile(saveFolder,['pol_' resonances{rn} '.jpg']));
        saveas(h3(rn),fullfile(saveFolder,['2d_' resonances{rn} '.fig']));
        saveas(h3(rn),fullfile(saveFolder,['2d_' resonances{rn} '.jpg']));

    
    end
end

% set(0,'DefaultAxesColorOrder',get(0,'FactoryAxesColorOrder'));
end