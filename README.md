# Generic qPCR Limit of Detection (LOD) / Limit of Quantification (LOQ) calculator

Christopher M. Merkes, Katy E. Klymus, Michael J. Allison, Caren Goldberg, Caren C. Helbing, Margaret E. Hunter, Craig A. Jackson, Richard F. Lance, Anna M. Mangan, Emy M. Monroe, Antoinette J. Piaggio, Joel P. Stokdyk, Chris C. Wilson, Catherine Richter

This script is designed to analyze qPCR data for many replicates of known concentration DNA standards and determine the limit of detection (LoD) and limit of quantification (LoQ) for use in environmental DNA applications. It is written in the hopes that users with even limited knowledge of R will be able to successfully use the code to analyze their own data in the same way as other eDNA researchers to get similar results and automatically generate plots to visualize the data. The code has 5 lines for user input and requires a few simple restrictions on the input data, but after that the user should be able to run the code without any further coding or interventions.

No R programming ability is required to run this script. However, the code does include many descriptive comments for those with moderate coding ability to understand what the commands are doing, and this file also includes tips for savvy users to make minor adjustments and gain additional functionality for refining their own analyses.

## Code files

This repository contains the following files:
- `README.md`: This file
- `LICENSE`: The standard USGS software license
- `LoD-calculator.R`: The script for automated generic LOD / LOQ analysis
- `Data.csv`: Example data that can be run for new users to do test runs

## Contact for code

Primary code developer: Chris Merkes (cmerkes@usgs.gov)

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www2.usgs.gov/visual-id/credit_usgs.html#copyright/).


This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

This software is provided "AS IS".

## Data Requirements

The code is written to expect a comma separated values (\*.csv) file with at least 3 specific columns: Target, Cq, and SQ. Additional columns can be included if they are helpful for the user to keep track of the data, but they will be ignored by the script.

If there is no column called "Target" and spelled exactly that way, the code will automatically search for a single column with the word "target" in it and rename it to "Target". If it is unable to find the column or finds multiple columns, it will print an error message in the Analysis Log.txt output file.

If there is no column called "Cq" and spelled exactly that way, the code will automatically search for a single column with "cq", "ct", or "cycle" and rename it to "Cq". If it is unable to find the column or finds multiple columns, it will print an error message in the Analysis Log.txt output file.

If there is no column called "SQ" and spelled exactly that way, the code will automatically search for a single column with "sq", "copies", "starting", or "quantity" and rename it to "SQ". If it is unable to find the column or finds multiple columns, it will print an error message in the Analysis Log.txt output file.

Negative reactions can be identified in a variety of ways. The criteria used are any value in the Cq column that is not read in as a numerical value are considered negative reactions. Blank space, NA, N/A, Undetermined, are all acceptable examples that will be detected as negative reactions. Nearly any text other than numbers in the Cq column are suitable and considered as negative reactions. If the data is manipulated to populate negative reaction Cq values with 0 or some other identifiable number, the code will analyze as if those are positive detects (although they will be flagged as potential outliers in most cases).

SQ column should be populated with numerical values for expected standard concentration. If estimated copy numbers per some calibration curve included on the plate are provided instead of the expected concentration in the SQ column, they will be assumed to all be different unique standards tested rather than replicates of the same standard, and that will disrupt all of the models and calculations performed accordingly. *Coming up with a generic way to detect this and automatically flag it is on my to-do list, and will likely be included in a future update.*

