# Reproduction code for: Quantifying Co-benefits of Water Quality Policies: An Integrated Assessment Model of Land and Nitrogen Management

## Description
This repository contains the code, scripts, and data required to reproduce the results presented in the paper.

## Computational Requirements
Before running the code, please ensure that you have the following software/libraries installed:
* GAMS
* R
* R studio

## Instructions for Data Preparation and Analysis
We used GAMS to generate calibration and simulation results and used R to produce related graphs. 

We processed land use data from Cropland Data Layer for the study site. Processed land use data can be found here: https://github.com/econweng/AJAEcobenefit/tree/main/Data/land%20use%20data

Part I: Calibration and Simulation using GAMS
* Calibration based on the observed baseline: https://github.com/econweng/AJAEcobenefit/blob/main/Code/GAMS/mendota_cobenefit_baseline.gms
* Policy simulation: https://github.com/econweng/AJAEcobenefit/blob/main/Code/GAMS/mendota_cobenefit_policysimuation.gms

Part II: Visualisation using R
*Related data can be found here: https://github.com/econweng/AJAEcobenefit/tree/main/Data/Simulated%20Results
*Related code: https://github.com/econweng/AJAEcobenefit/blob/main/Code/R/Co-benefits.Rmd

## License

This project is licensed under the MIT License.

Copyright (c) 2023 Weizhe Weng

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


