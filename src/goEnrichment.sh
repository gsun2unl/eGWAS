#!/bin/bash
find_enrichment.py HbTempLow.list HbPop.list All_ass.txt --pval=0.01 --method=fdr_bh --pval_field=fdr_bh --outfile=HBTemplowGO.xlsx

