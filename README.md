This repository contains full scripts to automate the process of constructing fertility density curves for US data. The worklow is as follows:

  * `download` data in archived files and `unzip` those; the unarchived yearly data files are quite large and can go up to `5.8 Gb`
  * optionally, download `documentation` files to understand data layout
  * `extract` only the columns of interest from yearly data files for further processing
  * construct an `SQLite` database with data from the extracted columns
  * do statistical analysis on the data by querying the `SQL` database

First three steps use only `bash` scripting. The last two steps use `Python` and `R` respectively.
