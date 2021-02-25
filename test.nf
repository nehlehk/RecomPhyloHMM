
params.out = '/home/nehleh/work/results/'
values = Channel.of( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

process tupleExample {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${x}_$filename" }
     maxForks 1
    
    
    input:
    tuple x, 'latin.txt' from values
    
    output:
        path "copy_${x}.txt" , emit: copy


    """
    cat - latin.txt > copy_${x}.txt
    """
}
