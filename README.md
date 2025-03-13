# Adaptive-Beamforming-for-Directional-Signal-Enhancement
# Project Overview
This project focuses on the development of an advanced adaptive beamforming system aimed at enhancing signal reception in complex multi-source environments. The primary goals are to improve the accuracy of direction-of-arrival (DOA) estimation and effectively suppress interference. The system integrates AIC-SVD for DOA denoising, alongside the RAIS-LCMV algorithm for spatial spectrum estimation. By continuously adapting the beamforming weights in response to environmental changes, the system optimizes the quality of the received signal, ensuring robust performance even in the presence of dynamic interference and fluctuating source directions.
# Background
**Beamforming** is a crucial technology in modern wireless communication systems, enabling signals to be focused in specific directions while reducing interference. Accurate estimation of the directions-of-arrival (DOA) of signals and interference is essential for effective beamforming. This project addresses DOA estimation errors by using AIC-SVD to enhance both the stability and accuracy of the estimates.
# Approach
* Applied **AIC-SVD** to refine real-time DOA estimates for both desired signals and interference sources.
* Used **RAIS-LCMV** algorithm to generate spatial spectra , identify signal  and interference directions.
* Designed an **adaptive beamformer** that dynamically adjusts antenna array weights based on the filtered DOA information.
# Key Features
* Real-time DOA tracking using AIC-SVD.
* High-resolution spatial spectrum estimation via RAIS-LCMV algorithm.
* "Adaptation of beamforming weights in real-time to respond to dynamic changes in the environment.
* Improved Signal-to-Interference Ratio (SIR) and robustness in multi-source scenarios.
# Tools & Technologies
**MATLAB** (System simulation & algorithm implementation)
**AIC-SVD** (DOA signal smoothing & prediction)
**RAIS-LCMV** (Spatial spectrum estimation)
**Adaptive Beamforming** (Dynamic weight calculation)
# Results
* The adaptive beamforming system outperformed conventional LCMV beamforming in terms of **interference rejection** and **target signal clarity**, particularly in dynamic scenarios with time-varying directions-of-arrival (DOAs).
