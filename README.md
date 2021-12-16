# Wind speed forecasting using the stationary wavelet transform and quaternion adaptive-gradient methods
This repo contains the supported mfiles to reproduce the [paper](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9535164) results using the stationary wavelet transform and quaternion adaptive-gradient methods 

## ABSTRACT
Accurate wind speed forecasting is a fundamental requirement for advanced and economically  viable large-scale wind power integration. The hybridization of the quaternion-valued neural networks and stationary wavelet transform has not been proposed before. In this paper, we propose a novel wind-speed forecasting model that combines the stationary wavelet transform with quaternion-valued neural networks.  The proposed model represents wavelet subbands in quaternion vectors, which avoid separating the naturally correlated subbands. The model consists of three main steps. First, the wind speed signal is decomposed using the stationary wavelet transform into sublevels. Second, a quaternion-valued neural network is used to  forecast wind speed components in the stationary wavelet domain. Finally, the inverse stationary wavelet transform is applied to estimate the predicted wind speed. In addition, a softplus quaternion variant of the RMSProp learning algorithm is developed and used to improve the performance and convergence speed of the proposed model. The proposed model is tested on wind speed data collected from different sites in China and the United States, and the results demonstrate that it consistently outperforms similar models. In the meteorological terminal aviation routine (METAR) dataset experiment, the proposed wind speed forecasting model reduces the mean absolute error, and root mean squared error of predicted wind speed values by 26.5% and 33%, respectively, in comparison to several existing approaches.

 ![image](https://user-images.githubusercontent.com/78357759/146317399-47fe976c-d832-4c72-bc85-14935aa5d007.png)

FIGURE 1. The proposed architecture. A quaternion representation of stationary wavelet subbands is used to train a quaternion-valued neural network. The predicted wind speed is computed by using the inverse stationary wavelet transform

## Getting started
1. Download all files and put them in the same folder. 
2. Add the folder to you MATLAB's path. 
3. Run the mfile [main.m](https://github.com/LyesSaadSaoud/Wind_forecast/blob/main/main.m) 

Please cite it as: L. S. Saoud, H. Al-Marzouqi and M. Deriche, "Wind Speed Forecasting Using the Stationary Wavelet Transform and Quaternion Adaptive-Gradient Methods," in IEEE Access, vol. 9, pp. 127356-127367, 2021, doi: 10.1109/ACCESS.2021.3111667.
