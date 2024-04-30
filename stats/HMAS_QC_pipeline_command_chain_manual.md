
## Generate confusion matrix after running hmas2 QC pipeline 
<br>


>  `python3 hmas2_confusion_matrix.py`  
>  `-i hmas2 QC pipeline output folder` (which contains subfolders for each sample)    
>  `-o output confusion_matrix file path`  
>  `-r common reference file for all those samples`  
>  `-m the metasheet file for all those samples`  (this is usually generated while extracting amplicon sequences)    
>  `-p mapping file` (mapping between sample and isolates. A sample might has multiple isolates in it)  
>  `-s the path for parse_count_table_confusion_matrix.py script`  

***note***  
1. Need to have blast loaded: ml ncbi-blast+/LATEST   
2. Possibly need to load a conda env which has the python pandas module installed  
3. the mapping file is a csv file, with header: `Sample	isolate_1	isolate_2	isolate_3`. If a sample has more than 3 isolates in it, you can add more columns to it. If a sample has only one isolate, you can leave the other 2 isolates column blank.  
4. There are some intermediary files generated in the working directory of this script, which you might want to clean up after the run.  

<br>



## Run stand-alone `parse_count_table_confusion_matrix.py` script
### <br>


**note**:    
1. Need to have blast loaded: ml ncbi-blast+/LATEST   
2. **this script has been updated so that it will take one single argument of confusion_matrix.ini file, and all required arguments are set in that config.ini file, as:** 
>  `python3 parse_count_table_confusion_matrix.py`  
>  `-c confusion_matrix.ini`   

<br>

---