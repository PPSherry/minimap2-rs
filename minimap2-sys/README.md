minimap2 public API:

example workflow presented by example.c
1. Initialization: 
    set mm_idxopt_t and mm_mapopt_t
2. Index Loading and Processing: 
    index get seperated into several parts, and load index part
3. Mapping:
    map each query to a single index part and then loop through all the index parts
4. Output: mm_reg1_t
5. Cleanup
