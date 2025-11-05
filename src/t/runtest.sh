profit -f fit.prf

../profitcore -o1 1yqv.core -o2 8fab.core zones.txt pdb1yqv_0P.mar 8fab.fit
cat 1yqv.core 8fab.core | pdbchain >both.pdb
