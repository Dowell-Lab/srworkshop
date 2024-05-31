# Library QC homework answers
Author: Lynn Sanford, 2021

<h3>Library 1</h3>
This library has some duplication issues, but generally passes QC. The base content and adapter content  are fine, and while only about half of the library is unique, either due to enrichment or PCR duplication,  it’s probably still usable (depending on overall read depth). 

<h3>Library 2</h3>
This library has no adapter contamination, but most of the reads are poly-A sections, based on base composition. This could be due to adapter ligation problems (some library preps involve adding poly-A  sections to inserts) or PCR duplication of poly-A stretches.

<h3>Library 3</h3>
This library has low complexity, based on very high duplication levels. Depending on the application, this  may be expected, but it’s probably a PCR duplication problem that will lead to an unusable library.

<h3>Library 4</h3>
This library looks pretty good across all three metrics.

<h3>Library 5</h3>
This library has major adapter contamination based on the base composition and adapter content  graphs. More cleanup needs to be done on the library before sequencing. It may still be usable after  trimming adapters (more on that later in the week), but the read depth would be compromised.
