# Structure 

The subfolders contain the FST files of the full model and all model simplifications (each row is a simulated parameter set).
Structure of the files:
| Column | Content | Variable name |
|:------------|:-------------:|-------------:|
| 1:30        | parameters          |          |
| 31:33       | FST Ini1, averaged over all alleles on RNAP,RNA,Protein level | avg_fsw1        |
| 34:36       | FST Ini2, averaged over all alleles on RNAP,RNA,Protein level | avg_fsw2        |
|37	nr_par+7 |	Mean of RNAP A distribution over last 50h and all alleles - Ini1 |	mean_pa1 |
|38	nr_par+8 |	Mean of RNAP A distribution over last 50h and all alleles - Ini2 |	mean_pa2 |
|39	nr_par+9 |	Mean of RNAP B distribution over last 50h and all alleles - Ini1 |	mean_pb1 |
|40	nr_par+10 |	Mean of RNAP B distribution over last 50h and all alleles - Ini2 |	mean_pb2 |
|41	nr_par+11 |	Mean of RNA A distribution over last 50h and all alleles - Ini1 |	mean_ra1 |
|42	nr_par+12 |	Mean of RNA A distribution over last 50h and all alleles - Ini2 |	mean_ra2 |
|43	nr_par+13 |	Mean of RNA B distribution over last 50h and all alleles - Ini1 |	mean_rb1 |
|44	nr_par+14 |	Mean of RNA B distribution over last 50h and all alleles - Ini2 |	mean_rb2 |
|45	nr_par+15 |	Mean of Protein A distribution over last 50h and all alleles - Ini1 |	mean_pra1 |
|46	nr_par+16 |	Mean of Protein A distribution over last 50h and all alleles - Ini2 |	mean_pra2 |
|47	nr_par+17 |	Mean of Protein B distribution over last 50h and all alleles - Ini1 |	mean_prb1 |
|48	nr_par+18 |	Mean of Protein B distribution over last 50h and all alleles - Ini2 |	mean_prb2 |
