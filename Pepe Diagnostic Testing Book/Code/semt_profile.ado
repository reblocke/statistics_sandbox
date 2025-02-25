*! 1.0.0 03 Apr 2003
/* semt_profile.ado

   called by example do-files for examples from:

   "The Statistical Evaluation of Medical Tests for Classification and Prediction"

   place this file in your Stata personal ado directory
     (type -sysdir- in Stata to determine where this directory is for your system.)

   Edit the file to  point to the appropriate data directory and the desired log file target directory.

   global macros are defined which specify

    1) the data directory path where example datasets reside
        (these can be downloaded from  http://www.fhcrc.org/labs/pepe/book)

    2) the directory path specifying where you would like exercise log files to go.
       set this to "" if you wish log files to go into the current working directory.
*/
pro def semt_profile
    version 7.0
    global semt_data  "..\..\data\"    /* data directory: edit this line accordingly */
    global semt_log   "h:\tmp\"        /* directory for log files & result files: edit this line  */
end
