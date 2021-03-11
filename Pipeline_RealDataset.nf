nextflow.enable.dsl = 2



params.wholegenome = "${PWD}/wholegenome.fasta"
params.xml_file = "${PWD}/bin/GTR_template.xml"
params.out =  "${PWD}/PipelineResults" 


frequencies = Channel.value(' 0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
nu_hmm = Channel.of(0.03)
mix_prob = Channel.of(0.9)






process Gubbins {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "$filename" }
    conda  "bioconda::gubbins=2.4.1"
    errorStrategy 'retry'
    maxRetries 3
    label 'short'
        
    
     output:
        path "gubbins.final_tree.tre" , emit: gubbinstree
        
    
     """   
        run_gubbins.py  -r GTRGAMMA -p gubbins   ${params.wholegenome} 
          
     """
}


process get_raxml_tree {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "$filename" }
    maxForks 1
    conda "bioconda::raxml=8.2.12"
    errorStrategy 'retry'
    maxRetries 3
    label 'short'

   
    output:
        path 'RAxML_bestTree.tree', emit: myRaxML
    
    
    """
     raxmlHPC -m GTRGAMMA   -p 12345 -s ${params.wholegenome} -N 10 -n tree 
    """
}


process make_xml_seq {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "$filename" }
    conda "numpy bioconda::dendropy"
    errorStrategy 'retry'
    maxRetries 3
    label 'short'
    

    input:
        path myRaxML

     output:
        path 'originalSeq.xml' , emit: original_XML 


     """
       python3.8  $PWD/bin/MakeXMLSeq.py -t ${myRaxML} -a ${params.wholegenome}  -x ${params.xml_file}        
        
     """
}

process Beast_alignment {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "$filename" }
     conda "conda-forge::python=3.7  bioconda::beast2"
     errorStrategy 'retry'
     maxRetries 3
     label 'short'

     input:
         path original_XML
         
         
     output:    

         path 'wholegenome.trees' , emit:  beastSeqTree
     
     """
       beast  ${original_XML}        
     """
}


process treeannotator_alignment {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "$filename" } 
     conda "conda-forge::python=3.7  bioconda::beast2"
     errorStrategy 'retry'
     maxRetries 3
     label 'short'

     input:
     
         path beastSeqTree 
 
         
     output:     
         path 'beastSeqTree.nexus' , emit: SeqTree
     
     """
       treeannotator -b 10  ${beastSeqTree}  beastSeqTree.nexus      
     """
}


process convertor_SeqTree {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "$filename" }
     conda "conda-forge::python=3.8  bioconda::dendropy" 
     errorStrategy 'retry'
     maxRetries 3  
     label 'short'  

     
     input:
        path SeqTree
    
     output:
         path 'beastSeqTree.newick' , emit : beastTree
         
     """
       python3.8 $PWD/bin/NexusToNewick.py -t ${SeqTree} -o 'beastSeqTree.newick'
        
     """
}


process CFML {

   publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "$filename" }
   conda "bioconda::clonalframeml=1.12"
   errorStrategy 'retry'
   maxRetries 3
   label 'short'
    
    input:
        path myRaxML 
        
    
    output:
        path "CFML.labelled_tree.newick" , emit: CFMLtree
        path "CFML.importation_status.txt" , emit: CFML_recom
    
    """ 
    
     ClonalFrameML ${myRaxML} ${params.wholegenome} CFML >  CFML.result.txt  
      
    """
}



process phyloHMM {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas scikit-learn conda-forge::hmmlearn"
     errorStrategy 'retry'
     maxRetries 3
     label 'short'    

     input:
        path myRaxML 
        each nu_hmm
        each mix_prob


     output:
        path 'RecomPartial.xml' , emit: partial_XML 
        path 'PhyloHMM_Recombination.jpeg' , emit: phyloHMMFig
        path 'Recom_phyloHMM_Log.txt' , emit: phyloHMMLog
        val nu_hmm , emit: my_nu
        val mix_prob , emit: prob
        
        
     """
       python3.8  $PWD/bin/phyloHMM.py -t ${myRaxML} -a ${params.wholegenome}  -x ${params.xml_file} -nu ${nu_hmm} -p ${mix_prob} -st 0 
     
     """
}




process Beast_partial {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "nu_${nu_hmm}_prob_${mix_prob}_$filename" } 
     conda "conda-forge::python=3.7  bioconda::beast2"
     errorStrategy 'retry'
     maxRetries 3
     label 'short'



     input:
        path partial_XML
        each nu_hmm
        each mix_prob
         
         
     output:    
         path 'wholegenome.trees' , emit: beastPartialTree
     
     """
       beast  ${partial_XML}        
     """
}




process treeannotator_partial {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     conda "conda-forge::python=3.7  bioconda::beast2"
     errorStrategy 'retry'
     maxRetries 3
     label 'short'

     input:   
        path beastPartialTree
        each nu_hmm
        each mix_prob
         
         
     output:
     
         path 'beastOurTree.nexus' , emit: beastOurTree

     
     """
       treeannotator -b 10  ${beastPartialTree}  beastOurTree.nexus      
     """
}



process convertor_ourTree {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     conda "conda-forge::python=3.8 bioconda::dendropy=4.5.2"
     errorStrategy 'retry'
     maxRetries 3
     label 'short'
     
     input: 
        path beastOurTree
        each nu_hmm
        each mix_prob

    
     output:
         path 'beastHMMTree.newick' , emit:  beastHMMTree
         
     """
       python3.8 $PWD/bin/NexusToNewick.py -t ${beastOurTree} -o 'beastHMMTree.newick'
        
     """
}


process mergeTreeFiles {

    publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "nu_${nu_hmm}_prob_${mix_prob}_$filename" } 
    maxForks 1
    errorStrategy 'retry'
    maxRetries 3
    label 'short'

    input:
         path beastHMMTree
         path myRaxML 
         path beastTree
         path gubbinstree
         path CFMLtree
         each nu_hmm
         each mix_prob
         

    output:
         path 'allOtherTrees.newick' , emit: allOtherTrees
     
     """
       python3.8 $PWD/bin/mergeFiles.py ${beastHMMTree}  ${myRaxML}  ${beastTree}  ${CFMLtree}  ${gubbinstree} > allOtherTrees.newick
       
     """


}


process TreeCmp {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "nu_${nu_hmm}_prob_${mix_prob}_$filename" }    
     maxForks 1
     errorStrategy 'retry'
     maxRetries 3
     label 'short'

     input:
     
         path allOtherTrees
         each nu_hmm
         each mix_prob
         
         
     output:
         path 'TreeCmpResult.result' , emit: Comparison     
    
     """
       java -jar $PWD/bin/TreeCmp_v2.0-b76/bin/treeCmp.jar  -r ${clonaltree}  -i ${allOtherTrees} -d qt pd rf ms um rfw gdu -o TreeCmpResult.result -W
     """
}


process phyloHMM_beast {

     publishDir "${params.out}" , mode: 'copy' , saveAs:{ filename -> "num_${repeat_range}/num_${repeat_range}_nu_${nu_hmm}_prob_${mix_prob}_$filename" }
     maxForks 1  
     conda "numpy bioconda::dendropy=4.5.2 conda-forge::matplotlib=3.3.4 pandas scikit-learn conda-forge::hmmlearn"
     errorStrategy 'retry'
     maxRetries 3
     label 'short'

     input:
        path wholegenome
        path beastHMMTree 
        tuple val(repeat_range), path('recomlog') 
        path CFML_recom
        path CFMLtree
        each nu_hmm
        val repeat_range

//     output:
//        path 'rmse.rmse' , emit : rmse_partialBeast 
//        path 'PhyloHMM_Recombination.jpeg' , emit: phyloHMMFig_partialBeast

     """
       python3.8 $PWD/bin/phyloHMM.py -t ${beastHMMTree} -a ${wholegenome}  -l ${recomlog}  -nu ${nu_hmm} -ct ${CFMLtree} -c ${CFML_recom} -x ${params.xml_file}
       
        
     """
}



workflow {
    
    Gubbins()

    get_raxml_tree()
    
    make_xml_seq(get_raxml_tree.out.myRaxML)
    
    Beast_alignment(make_xml_seq.out.original_XML)

    treeannotator_alignment(Beast_alignment.out.beastSeqTree)

    convertor_SeqTree(treeannotator_alignment.out.SeqTree)  
   
    CFML(get_raxml_tree.out.myRaxML)
   
    phyloHMM(get_raxml_tree.out.myRaxML,nu_hmm,mix_prob)
      
    Beast_partial(phyloHMM.out.partial_XML,nu_hmm,mix_prob)
    
    treeannotator_partial(Beast_partial.out.beastPartialTree,nu_hmm,mix_prob)
    
    convertor_ourTree(treeannotator_partial.out.beastOurTree,nu_hmm,mix_prob)
    
//    mergeTreeFiles(convertor_ourTree.out.beastHMMTree,get_raxml_tree.out.myRaxML,convertor_SeqTree.out.beastTree,Gubbins.out.gubbinstree,CFML.out.CFMLtree,nu_hmm,mix_prob)

//    TreeCmp(BaciSim.out.clonaltree,mergeTreeFiles.out.allOtherTrees,nu_hmm)
    
//    phyloHMM_beast(seq_gen.out.wholegenome,convertor_ourTree.out.beastHMMTree,BaciSim.out.recomlog,CFML.out.CFML_recom,CFML.out.CFMLtree,nu_hmm,seq_gen.out.range)

    
       
}