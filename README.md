## Synopsis

This software is used to analyze a simple flux model. The original Excel data file is used as input. Each experiment should be in one sheet. The results of relative flux value in each experiment are exported to one .csv file. In each experiment, the average and standard deviation of each condition are calculated separately.

## Requirements

This software is developed and tested on Python 3.6. It also relies on pandas, numpy and scipy. It has been tested on pandas 0.22, numpy 1.14 and scipy 1.0.

## Usages

The raw data are in `13C-Glucose_tracing_Mike.xlsx` file. To solve the MFA, just run `main.py`. Results will be returned by `.csv` files in the same path. Each `.csv` file correspond one sheet in Excel file.

## Contributors

**Shiyu Liu**

+ [http://github.com/liushiyu1994](http://github.com/liushiyu1994)

## License

This software is released under the [MIT License][opensource].