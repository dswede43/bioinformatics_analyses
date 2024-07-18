# Daniel Stratis - TIGeR Bioinformatician Interview Assignment

## Part 1: VCF filtering
This part of the assignment uses a script that requires GATK to run. I used the official [GATK Docker image](https://hub.docker.com/r/broadinstitute/gatk) from the Broad Insititute.

### Instructions to run the GATK Docker container and script
1. Install [Docker](https://docs.docker.com/engine/install/).

2. Pull the image.

```docker pull broadinstitute/gatk:4.1.3.0```

3. Start the container.

```docker run -v .:/gatk/my_data -it broadinstitute/gatk:4.1.3.0```

4. Change directory and run the script

```cd my_data```

```./part1-VCF_filtering.sh```

## Part 2: VCF summary and visualization
This part of the assignment uses an Rmarkdown script. Knitting this script in Rstudio will re-run it and install all necessary packages required. The output is an HTML file (```part2-VCF_analysis.html```) and PDF (```variant_count_distribution.pdf```).
