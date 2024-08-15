# Daniel Stratis - TIGeR Bioinformatician Interview Assignment

## Part 1: VCF filtering
### Script
**part1-VCF_filtering.sh**

This script requires GATK to run. I used the official [GATK Docker image](https://hub.docker.com/r/broadinstitute/gatk) from the Broad Insititute.

### Instructions to run the GATK Docker container and script
1. Install [Docker](https://docs.docker.com/engine/install/).

2. Pull the image.
```docker pull broadinstitute/gatk:4.1.3.0```

3. Start the container.
```docker run -v .:/gatk/my_data -it broadinstitute/gatk:4.1.3.0```

4. Change directory and run the script
```cd my_data```
```./part1-VCF_filtering.sh```

### Output
1. Updated VCF files in the ```filtered_vcfs/``` folder

2. Logs stored in the ```logs/``` folder

## Part 2: VCF summary and visualization
### Script
**part2-VCF_analysis.Rmd**

This Rmarkdown script can be re-run by knitting in Rstudio. The script will install all required packages.

### Output
1. part2-VCF_analysis.html

2. variant_count_distribution.pdf
