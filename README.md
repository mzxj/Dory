# Dory: a computational method for differential chromatin tracing analysis
- [Installation](#Installation)
- [Usage](#Usage)
- [Output](#Output)
- [Example](#Example)
- [Contact us](#Contact-us)
## Installation
1. Download Dory from GitHub and navigate into the directory:
```
git clone https://github.com/mzxj/Dory.git && cd Dory
```
2. Make the bash script executable:
```
chmod +x Dory
```
3. Install required R packages manually, one by one, for the most reliable setup, including `ggplot2`, `optparse`, `reshape2`, `dplyr`, `data.table`, `future`, `parallel`, and `furrr`. Alternatively, you can install all of them automatically using the provided script:
```
Rscript ./R/install_packages.R
```
or inside R: ```source("./R/install_packages.R")```. Note: If installation fails for some packages in all of the above methods, please update your R installation to a newer version (R ≥ 4.4 recommended).

4. To make Dory globally runnable, add its path to your environment's PATH variable. Append `export PATH=/path/to/Dory/:$PATH` to your `~/.bashrc` or `~/.bash_profile`. Replace `/path/to/Dory` with the actual path to the Dory directory. Then reload your shell:
```
source ~/.bashrc  # or source ~/.bash_profile
```
### Verify installation
Test the Dory is available globally:
```
Dory -h
```
You should see the help message printed in your terminal.




## Usage
### Compare Two Conditions, Cell Types, or Cell States
```
Dory -m $FileFormat -a $Input1 -b $Input2 -o $OutputPath -c $ChrNum
```
- `-m $FileFormat`: (*Required*) Format of the input data ["4DN" or "csv"]. 
    - "4DN": If the input data is directly downloaded from the [4DN data portal](https://data.4dnucleome.org/resources/data-collections/chromatin-tracing-datasets). 
    - "csv": If the input data is standard .csv file. The 3D_coordinates.csv file requires these columns:
        - **Trace_ID**: The ID for each chromatin trace.
        - **X**, **Y**, **Z**: 3D spatial coordinates of the genomic region.
        - **Chrom**, **Chrom_Start**, **Chrom_End** or **Region_ID**: Indicate genomic location either by chromosome name, start and end coordinate, or by the genomic region ID (e.g., 1,2, ...). If "Chrom, Chrom_Start Chrom_End" columns are not provided, then the "Region_ID" column is required. 
- `-a $Input1`: (*Required*). File containing 3D coordinates of genomic regions for foreground condition. 
- `-b $Input2`: (*Required*). File containing 3D coordinates of genomic regions for background condition.
- `-o $OutputPath`: (*Optional*). Output directory path. Defaults to the current directory "./" if not specified.
- `-c $ChrNum`: (*Optional*). Number of traced chromosomes ["one" or "more"]. Required if the 3D coordinate file uses "Region_ID" instead of "Chrom, Chrom_Start Chrom_End". Optional otherwise. 


### Compare Multiple Cell Types (More Than Two)
```
Dory -m $FileFormat -i $InputCoordFile -l $LabelFile -n $LabelColame -o $OutputPath -c $ChrNum
```
- `-m $FileFormat`: (*Required*) Format of the input data ["4DN" or "csv"]. 
    - "4DN": If input files are directly downloaded from the [4DN data portal](https://data.4dnucleome.org/resources/data-collections/chromatin-tracing-datasets). 
    - "csv": If the input data is standard .csv file. 
        - The 3D_coordinates.csv file requires these columns:
            - **Trace_ID**: The ID for each chromatin trace.
            - **X**, **Y**, **Z**: 3D spatial coordinates of the genomic region.
            - **Chrom**, **Chrom_Start**, **Chrom_End** or **Region_ID**: Indicate genomic location either by chromosome name, start and end coordinate, or by the genomic region ID (e.g., 1,2, ...). If "Chrom, Chrom_Start Chrom_End" columns are not provided,  then the "Region_ID" column is required. 
            - **Cell_ID**: The ID for the cell.
        - The cell_label.csv file requires these columns:
            - **Cell_ID**: The ID for the cell.
            - **Cell_Type**: The cell-type label or classification assigned to the corresponding cell.

- `-i $InputCoordFile`: (*Required*). File containing 3D coordinates of genomic regions for all cells. 
- `-l $LabelFile`: (*Optional*) File with cell-type labels or classifications for cells. Required if cell-type labels or classifications are not included in the input 3D coordinate file. 
- `-n $LabelColname`: (*Optional*). The column name for the cell-type label or classification. Defaults to "Cell_Type" if not specified.
- `-o $OutputPath`: (*Optional*). Output directory path. Defaults to the current directory "./" if not specified.
- `-c $ChrNum`: (*Optional*). Number of traced chromosomes ["one" or "more"]. Required if using "Region_ID" instead of "Chrom, Chrom_Start Chrom_End" in the input 3D coordinate file. Optional otherwise. 


## Output
If the traced genomic regions are located on one chromosome, the output includes:
- S0_DataInfo
    - Region.tsv: Contains traced genomic regions with columns for Chrom, Chrom_Start, Chrom_End and Region_ID.
    - TraceCount.pdf: Displays the number of traces for each condition.
- S1_Distance
    - RegionPairByTrace_CT_x.tsv: The Euclidean distances between region pairs (rows) across all chromatin traces (columns) for each condition (CT_x).
- S2_DiffScore
    - DiffScoreMatrix_CT_xVSCT_y.tsv: The DiffScore matrix generated by Dory for comparing condition CT_x (foreground) with CT_y (background).
    - DiffScoreHeatmap_CT_xVSCT_y.pdf: The heammap visualization of the DiffScore matrix.
### Special Case: Multiple Chromosomes
If the traced genomic regions span multiple chromosomes, the output directory structure changes slightly:
- Because Dory performs intra-chromosomal analysis, the results are separated by chromosome (**Chrom**).
- Each main output folder (S0_DataInfo, S1_Distance and S2_DiffScore) contains subdirectories named "chr_**Chrom**" (e.g., chr_chr1, chr_chr2, etc.), corresponding to each chromosome.
- Each "chr_**Chrom**" subdirectory contains the same types of output files as described above, but limited to the genomic regions on corresponding chromosome.

## Example
Using the two datasets [4DNFIHSXQZIV.csv](https://data.4dnucleome.org/files-processed/4DNFIHSXQZIV/) and [4DNFITZBMT6Q.csv](https://data.4dnucleome.org/files-processed/4DNFITZBMT6Q/) downloaded from the 4DN data portal as an example. Both files are located in the folder [tests](https://github.com/mzxj/Dory/tree/main/tests).
```
Dory -m 4DN -a /path/to/4DNFIHSXQZIV.csv -b /path/to/4DNFITZBMT6Q.csv -o output/to/path 
```


## Contact us
For any questions regarding Dory, please contact Zhaoxia Ma (prm2zt@virginia.edu).
