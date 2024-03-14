GPL-3.0-or-later

# alPCA
An automatic software for the selection and combination of forecasts in monthly series

In this alPCA.r version, the user can obtain forecasts for a 12-month planning horizon for any time series. To select the methods, only 3 percentiles are taken into account. The 5th percentile if we want few of them, the 95th percentile if we want to consider almost all of them, and the 50th percentile if we want an intermediate number of them. In this case only the weights obtained from the scores of each selected configuration are used.

The previous results that led us to this version can be found in the paper García-Aroca C., Martínez-Mayoral MA., Morales-Socuéllamos J. and Segura Heras JV. (2024). An algorithm for automatic selection and combination of forecast models. Expert Systems with Applications. 237(121636). DOI: https://doi.org/10.1016/j.eswa.2023.121636. Source Code https://doi.org/10.24433/CO.4777598.v1

The data of the monthly series used as an example and the 6 functions defined for the calculations have been included in the code.
