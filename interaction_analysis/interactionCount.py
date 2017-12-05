import sqlite3
connection = sqlite3.connect("../db/test.sqlite")

cursor = connection.cursor()

count = 0;
gene_dict = {"ADGRG1":0,"MYH9":0,"RASA1":0,"PML":0,"CDK6":0,"OPHN1":0,"NF1":0,
"ACHE":0,"NRP1":0,"GLI1":0,"THY1":0,"AMOT":0,"HEY1":0,"PLG":0,
"SCG2":0,"TLE3":0,"CRMP1":0,"VEGFA":0,"ETS2":0,"LDB1":0,"SLIT1":0,
"TLE1":0,"RTN1":0,"NRCAM":0,"NRP2":0,"PTCH1":0,"SHH":0,"NKX6-1":0,
"VLDLR":0,"CDK5R1":0,"DPYSL2":0,"CELSR1":0,"HEY2":0,"L1CAM":0,
"CNTFR":0,"UNC5C":0}

for key in gene_dict.keys():
    cursor.execute('SELECT DISTINCT UniProtKBGeneNameID FROM UniProtKBmap WHERE Genename =\'{k}\';'.format(k=key));
    result = cursor.fetchall()
    for i in range(len(result)):
      cursor.execute('SELECT COUNT(*) FROM( SELECT DISTINCT interactorA, interactorB FROM Interactions WHERE interactorA = \'uniprotkb:{r}\' UNION SELECT DISTINCT interactorA, interactorB FROM Interactions WHERE interactorB = \'uniprotkb:{r}\') x'.format(r=result[i][0]));
      temp = cursor.fetchall()
      count = count + temp[0][0]
    gene_dict[key] = count
    count = 0  

      
print(gene_dict)    
    

connection.commit()

connection.close()
