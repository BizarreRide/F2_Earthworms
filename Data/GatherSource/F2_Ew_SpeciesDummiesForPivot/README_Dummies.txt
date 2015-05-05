Create a Dummy, simulating that each species was found in each hole at each field. 
Creation of dummies is performed with Dummies.txt in R, using the melt() function from reshape2 package.
Write out .csv
The .csv data has to be processed and added to "Alle Aufnahmen" in F2_EW_Data.xlsx.
Weight of the Dummy must stay empty in the Excel sheet, hence NAs must be removed.
Furthermore X.Leer has to be changes to empty cells.

Then Excel Pivot table can create weight sums taking zeros for empty fields, and reveals counts for weights as proxy for 
abundance data.
To gain real abundance data, all real non-dummy, non-weight Individulas must get weight value of zero.