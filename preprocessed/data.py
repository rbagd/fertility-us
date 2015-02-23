import os
import sqlite3
import csv

# This is to be used with pre-processed raw data, i.e. extract the necessary
# columns with bash into files of type 'age_1976'. The rest of the code is 
# then a hack to put it all into a SQLite database.

os.chdir('/home/rytis/Downloads/downloads_birth/preprocessed')
conn = sqlite3.connect('data.db')
c = conn.cursor()
c.execute("CREATE TABLE data (id integer primary key autoincrement, age integer, children integer, year integer)")
conn.commit()

names = ['age_', 'child_']

years = range(1971, 2014)

for year in years:
    print(year) # Just to see where the loop is
    filenames = [ i + str(year) for i in names ]

    fa = open(filenames[0], 'r')
    fb = open(filenames[1], 'r')

    fa_reader = csv.reader(fa, delimiter="\n")
    fb_reader = csv.reader(fb, delimiter="\n")
    
    age = []
    child = []
    for line in fa_reader:
        age.append(int(line[0]))
    for line in fb_reader:
        child.append(int(line[0]))
    year_list = [year]*len(age)

    c.executemany("INSERT INTO data (age, children, year) VALUES (?,?,?)", zip(age, child, year_list))
    conn.commit()
    if year % 10 == 0:
        c.execute("VACUUM")

# Age coding for 2003 is different. It suffices to increment it by 13.
c.execute("UPDATE data SET age = age + 13 WHERE year = 2003")

conn.commit()
conn.close()
