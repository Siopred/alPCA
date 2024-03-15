GPL-3.0-or-later

# alPCA
An automatic software for the selection and combination of forecasts in monthly series

In this alPCA.r version, the user can obtain forecasts for a 12-month planning horizon for any time series. To select the methods, only 3 percentiles are taken into account. The 5th percentile if we want few of them, the 95th percentile if we want to consider almost all of them, and the 50th percentile if we want an intermediate number of them. In this case only the weights obtained from the scores of each selected configuration are used. The associated function file in this case is f_alPCA.monthly.r.

The average monthly data of CO2 mole fraction (co2-mm-mlo), which represents CO2 emissions, have been taken as an illustrative example. The series is obtained from daily averages and covers from 01/01/2000 to 31/12/2019 (20 years and 240 observations). The data are from the U.S. Government’s Earth System Research Laboratory, Global Monitoring Division, and are available from the Trends in Atmospheric Carbon Dioxide website https://datahub.io/core/co2-ppm#data. When the algorithm is finished, it generates a 12x3 dimension prediction matrix under the name result alPCA, the first column contains the predictions obtained with the 5th percentile, the second column those obtained with the 50th percentile and the last column those obtained with the 95th percentile.

The previous results that led us to this version can be found in the paper García-Aroca C., Martínez-Mayoral MA., Morales-Socuéllamos J. and Segura Heras JV. (2024). An algorithm for automatic selection and combination of forecast models. Expert Systems with Applications. 237(121636). DOI: https://doi.org/10.1016/j.eswa.2023.121636. Source Code https://doi.org/10.24433/CO.4777598.v1
