efetch -db taxonomy -id 2559597 -format xml | \
xtract -pattern Taxon -first TaxId -element Taxon -block "*/Taxon"  \
-unless Rank -equals "no rank" -tab "," -sep "_" -element Rank,ScientificName
