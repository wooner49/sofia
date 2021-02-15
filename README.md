# Robust Factorization of Real-world Tensor Streams with Patterns, Missing Values, and Outliers (ICDE'21)
This repository contains the source code for the paper [Robust Factorization of Real-world Tensor Streams with Patterns, Missing Values, and Outliers](https://google.co.kr), by [Dongjin Lee](https://github.com/wooner49) and [Kijung Shin](https://kijungs.github.io/), presented at [ICDE 2021](https://icde2021.gr/).

In this work, we propose **SOFIA**, an online algorithm for factorizing real-world tensors that evolve over time with missing entries and outliers. By smoothly and tightly combining tensor factorization, outlier detection, and temporal-pattern detection, SOFIA achieves the following strengths over state-of-the-art competitors:
* **Robust and accurate**: SOFIA yields up to 76% and 71% lower imputation and forecasting error than its best competitors.
* **Fast**: Compared to the second-most accurate method, using SOFIA makes imputation up to 935X faster.
* **Scalable**: SOFIA incrementally processes new entries in a time-evolving tensor, and it scales linearly with the number of new entries per time step.

## Datasets
| Name          | Description                                       | Size                   | Granularity in Time | Processed Dataset | Original Source   |
| ------------- |:-----------------------------------------------:| :---------------------:| :----------:| :---: |:---:|
| Intel Lab Sensor   | locations x sensor x time         | 54 x 4 x 1152     | every 10 minutes      | [Dataset]() | [Link](http://db.csail.mit.edu/labdata/labdata.html) |
| Network Traffic | sources x destinations x time       | 23 x 23 x 2000       | hourly     | [Dataset]() | [Link](https://www.cs.utexas.edu/~yzhang/research/AbileneTM/) |
| Chicago Taxi | sources x destinations x time       | 77 x 77 x 2016    | hourly   | [Dataset]()  | [Link](https://data.cityofchicago.org/Transportation/Taxi-Trips/wrvz-psew) |
| NYC Taxi   | sources x destinations x time | 265 x 265 x 904 | daily    |  [Dataset]() | [Link](https://www1.nyc.gov/site/tlc/about/tlc-trip-record-data.page) |
