# cdom
Various packages for processing and analyzing colored dissolved organic matter (CDOM) absorption

## **Primary Functions**

#### **process_acdom**
This function uses CDOM absorbance, deionized water (DI) scans and pathlength to calculate CDOM absorption. The function follows a specific station naming convention noted below. Once absorbance scans are classified to a specific station, an outlier analysis is used to ensure all scans used to calculate absorption are quality scans. Outliers are identified using Matlab's isoutlier function (Matlab 2017a) with default settings (outlier is defined as more than three scaled median absolute deviations away from median of all available scans).

Sample naming convention

cruise_date_station_depth_sample-type_replicate

Example: LS_200309_1_5_ag_R1

#### **cdom_spectral_slope**
This function fits an exponential model to measured CDOM absorption and calculates the spectral slope, or shape of the exponential curve, over a user-specified wavelength range. This function requires Matlab's curve fitting toolbox and the fit function, along with sub-functions cdom_model and cdom_model_noK. Model output contains coefficients in dot notation (e.g., model.s is the spectral slope).

### **Dependent Functions**

**cdom_model**

**cdom_model_noK**

#### Example Code


```
file = '~/data_folder/cdom_examples/raw_cdom/cdom_absorbance_examples.xlsx';

% note - change data_folder in file path to location of file

dat = process_acdom(file,0.05,[],200,750);

for ii = 1:length(dat)
    model = cdom_spectral_slope(dat(ii).wavelength,dat(ii).ag,275,295,350);
    S275295(ii) = model.s;
end

S275295 = round(S275295,4);
```
