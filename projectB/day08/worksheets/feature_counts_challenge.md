# Feature Counts Advanced Challenge
This is a extra challenge problem for the people who have completed the Day 8 worksheet. Because this is a challenge, we have intentially kept the instructions minimal and encourage you to read documentation and troubleshoot.

## Goal
Calculate FRiP scores from your day 8 MACS output (peaks) and BAM files using feature counts.

## Instructions
1. Run feature counts on your MACS data to generate counts for your narrow peaks.
    - You will need to modify your annotation file to work with feature counts.
    - Here is the Feature Counts [documentation](https://subread.sourceforge.net/featureCounts.html) and [User Manual](https://subread.sourceforge.net/SubreadUsersGuide.pdf) for reference.
    - Also refer to day07 scripts for how you ran feature counts previously.

2. Use the counts to calculate a FRiP score per peak.
    - See the [ChIP Video](https://www.youtube.com/watch?v=PBaC2qG4f4k) for information about what a FRiP score is.
    - Below is the equation for calculating a FRiP score:

$$FRiP = \frac{\text{No. of Reads in Peaks}}{\text{Total No. of Reads}}$$

