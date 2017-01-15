This repository contains full scripts to automate the process of constructing fertility density curves for US data. The worklow is as follows:

  * `download` data in archived files and `unzip` those; the unarchived yearly data files are quite large and can go up to `5.8 Gb`
  * optionally, download `documentation` files to understand data layout
  * `extract` only the columns of interest from yearly data files for further processing
  * construct an `SQLite` database with data from the extracted columns
  * do statistical analysis on the data by querying the `SQL` database

First three steps use only `bash` scripting. They are separated since each of the steps may be time consuming. The last two steps use `Python` and `R` respectively. Construction of an `SQL` database may take a few hours.

To reconstruct the data, one has to run the following commands.

```
chmod +x {download,unzip,extract}
./download && ./unzip && ./extract
python preprocessed/data.py
Rscript preprocessed/data.R
```
Be aware that some extra dependencies are necessary, notably `sqlite` module for `Python` as well as `RSQLite` library for `R`. Some additional libraries for `R` are useful to treat raw data and obtain plots.

Plurality data for 1969 and 1970, i.e. number of children at birth, are not available. `SQLite` database is therefore constructed from 1971 only.
