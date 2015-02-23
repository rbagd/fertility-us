This repository contains full scripts to automate the process of constructing fertility density curves for US data. The worklow is as follows:

  * download data in archived files and unzip those; the unarchived yearly data files are quite large and can go up to `5.8 Gb`
  * extract only the columns of interest for further processing
  * construct an `SQLite` database with data from the extracted columns
  * do statistical analysis on the data by querying the `SQL` database

First two steps use only `bash` scripting. The last two steps use `Python` and `R` respectively.
