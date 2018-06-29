# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 15:17:29 2017

@author: Karim
"""
import mysql.connector
cnx = mysql.connector.connect(host = '***', user = '***', password = '***', database = 'world')
cur = cnx.cursor()

genes_proteins = 'CREATE TABLE Genes_Proteins (Id INTEGER NOT NULL AUTO_INCREMENT,\
                                                    Gene TEXT,\
                                                    Protein TEXT,\
                                                    Organism TEXT,\
                                                    PRIMARY KEY (Id))'

cur.execute(genes_proteins)

sequences = 'CREATE TABLE Sequences (Id INTEGER NOT NULL AUTO_INCREMENT,\
                                        Sequence TEXT,\
                                        PRIMARY KEY (Id))'


cur.execute(sequences)

function = 'CREATE TABLE Functions (Id INT NOT NULL AUTO_INCREMENT,\
                                    Molecular_function TEXT,\
                                    Biological_process TEXT,\
                                    PRIMARY KEY (Id))'


cur.execute(function)


cur.close()
cnx.close()
