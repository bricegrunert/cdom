function [dat]=process_acdom(file,pathlength,dat,lam_min,lam_max)

% process_acdom
%
% [dat]=process_acdom(file,pathlength,dat,lam_min,lam_max)
%
% Function uses measured CDOM absorbance, DI scans and pathlength to
% calculate CDOM absorption. The function is specific to the naming
% convention used in CUNY's Land-Ocean-Atmosphere Interctions Lab (PI:
% Prof. Maria Tzortziou) for sorting through absorbance headers (see
% below). Once absorbance scans are classified to a specific station, based
% on location, date/time and depth, an outlier analysis is used to ensure
% all scans used to calculate absorption are quality scans. Outliers are
% identified using Matlab's isoutlier function with default settings,
% meaning an outlier is defined as more than three scaled median absolute
% deviations away from the median. Individual absorbance scans are compared
% to the group median at 290-310, 390-410, and 490-510 nm, to ensure no
% unusual spectral deviations impact a scan. Alternatively, outliers can be
% identified as the integral of each curve; however, while this approach is
% inclusive of all wavelengths, it is less sensitive to deviations as the
% comparison product is a sum of a larger population.

% Sample naming convention
% cruise_date_station_depth_sample-type_replicate
%
% Example: LS_200309_1_5_ag_R1
%
% Inputs:
%
%     file               = user date file (xls format); alteratively, user
%                          can alter line XXX to use preferred file read 
%                          format/function
%     pathlength         = pathlength of cuvette used in cm (e.g., 5)
%     dat                = preexisting data structure to append new data
%                          to, if one exists
%     lam_min            = minimum wavelength to calculate absorption over 
%     lam_max            = maximum wavelength to calculate absorption over
%
% Returns:
% Data structure called dat, which contains the following fields
%
%     abs_raw            = mean of raw absorbance scans that passed outlier 
%                          test
%     di1                = nearest DI scan made before CDOM absorbance
%     di2                = nearest DI scan made after CDOM absorbance
%     di_avg             = mean DI scan used for correction, calculated as
%                          average of di1 and di2 (note: no QAQC done on DI
%                          scans, assumed this is done prior to reading in)
%     abs_corr           = raw absorbance with DI absorbance removed
%                          (DI corrected absorbance)
%     ag                 = CDOM absorption calculated from abs_corr scans,
%                          i.e., ag = 2.303*mean(abs_corr)/pathlength
%
% Recommended reading
% Green, Sarah A., and Neil V. Blough. "Optical absorption and fluorescence 
% properties of chromophoric dissolved organic matter in natural waters." 
% Limnology and Oceanography 39.8 (1994): 1903-1916.
%
%
% Function is compatible from Matlab 2017a
% copyright (c) 2020 Brice K. Grunert
% email: bricegrunert@gmail.com
%
% Last modified on 23 December 2020 by BG
%
% pending updates: change xlsread to readtable
%
%%%%%%%


if exist('dat','var')==0
    sampNum=0;
else
    sampNum=length(dat);
end

if exist('lam_min','var')==0
    lam_min=250;
end

if exist('lam_max','var')==0
    lam_max=750;
end

[NUM,TXT,RAW]=xlsread(file);

fishing=find(~isnan(NUM(:,1))==1);

NUM=NUM(fishing,:);

clear newhdr di_ind

for ii=1:size(TXT,2)
    if contains(TXT{1,ii},'_R')
        newhdr{ii}=extractBefore(TXT{1,ii},'_R');
    else
        newhdr{ii}=TXT{1,ii};
    end
end

cnt=0;
for ii=1:size(TXT,2)
    if contains(TXT{1,ii},'DI')
        cnt=cnt+1;
        di_ind(cnt)=ii;
    end
end

di_ind=di_ind+1;



[uu,aa,bb]=unique(newhdr);


for ii=1:length(uu)
    if length(uu{ii}) > 0 & contains(uu{ii},'line')==0 & contains(uu{ii},'DI')==0
        sampNum=sampNum+1;
        if contains(uu{ii},'_ag')
            dat(sampNum).station=extractBefore(uu{ii},'_ag');
        else
            dat(sampNum).station=uu{ii};
        end
        ind=find(bb==ii);
        lam_ind=find(NUM(:,ind(1)) >= lam_min & NUM(:,ind(1) <= lam_max));
        dat(sampNum).wavelength=NUM(lam_ind,ind(1));
        ind=ind+1;
        A=NUM(lam_ind,ind);
        
        if length(ind) > 2
            %look for outlier scans
            lind=find(dat(sampNum).wavelength >=490 & dat(sampNum).wavelength <= 510);
            igood=isoutlier(nanmean(A(lind,:)));
            
            lind=find(dat(sampNum).wavelength >=390 & dat(sampNum).wavelength <= 410);
            igood2=isoutlier(nanmean(A(lind,:)));
            
            lind=find(dat(sampNum).wavelength >=290 & dat(sampNum).wavelength <= 310);
            igood3=isoutlier(nanmean(A(lind,:)));
            
            igood=igood+igood2+igood3;
            igood=find(igood==0);
            
            A=nanmean(A(:,igood),2);
            
        else
            
            A=nanmean(A,2);
            
        end
        
        dat(sampNum).abs_raw=A;
        
        ind=ind-1;
        dat(sampNum).abs_good_scans=extractAfter(TXT(1,ind),'ag_');
        
        di1=find(di_ind < min(ind));
        di1=di1(end);
        
        di2=find(di_ind > max(ind));
        di2=di2(1);
        
        dat(sampNum).di1=NUM(lam_ind,di_ind(di1));
        dat(sampNum).di2=NUM(lam_ind,di_ind(di2));
        
        dat(sampNum).di_avg=(dat(sampNum).di1+dat(sampNum).di2)/2;
        
        dat(sampNum).abs_corr=dat(sampNum).abs_raw-dat(sampNum).di_avg;
        if dat(sampNum).wavelength(1) > dat(sampNum).wavelength(end)
            dat(sampNum).wavelength=flipud(dat(sampNum).wavelength);
            dat(sampNum).abs_raw=flipud(dat(sampNum).abs_raw);
            dat(sampNum).di1=flipud(dat(sampNum).di1);
            dat(sampNum).di2=flipud(dat(sampNum).di2);
            dat(sampNum).di_avg=flipud(dat(sampNum).di_avg);
            dat(sampNum).abs_corr=flipud(dat(sampNum).abs_corr);
        end
        dat(sampNum).ag=2.303.*dat(sampNum).abs_corr/pathlength;
    end
end
